# === 包加载 ===
library(tidyverse)
library(vegan)
library(ggpubr)
dir()
# === 数据读取与匹配 ===
kegg_css <- read.csv("CSS_relative_KEGG_bases_insitu.csv", row.names = 1)
metadata  <- read.csv("metadata.csv")
kegg_t    <- as.data.frame(t(kegg_css))
metadata  <- metadata %>% slice(match(rownames(kegg_t), Sample))

# === 1. 计算 Shannon 多样性 ===
shannon_df <- kegg_t %>%
  mutate(Sample = rownames(.)) %>%
  rowwise() %>%
  mutate(Shannon = diversity(c_across(-Sample), index = "shannon")) %>%
  ungroup() %>%
  left_join(metadata, by = "Sample")

# === 新增：计算统计量 ===
stats_summary <- shannon_df %>%
  group_by(Group) %>%  # 按组（IT/ST）分组计算
  summarise(
    Median = median(Shannon, na.rm = TRUE),
    Mean   = mean(Shannon, na.rm = TRUE),
    SD     = sd(Shannon, na.rm = TRUE),
    IQR    = IQR(Shannon, na.rm = TRUE),
    Min    = min(Shannon, na.rm = TRUE),
    Max    = max(Shannon, na.rm = TRUE)
  )

# 打印统计结果
print(stats_summary)
# === 2. 正态性 & 方差齐性 检验 ===
sw_IT    <- shapiro.test(shannon_df$Shannon[shannon_df$Group=="IT"])
sw_ST    <- shapiro.test(shannon_df$Shannon[shannon_df$Group=="ST"])
levene   <- car::leveneTest(Shannon ~ Group, data = shannon_df)  # 需要先 install.packages("car")

# === 3. 选择检验并提取结果 ===
if (sw_IT$p.value > 0.05 && sw_ST$p.value > 0.05 && levene[["Pr(>F)"]][1] > 0.05) {
  base_test   <- t.test(Shannon ~ Group, data = shannon_df, var.equal = TRUE)
  stat_method <- "t.test"
} else {
  base_test   <- wilcox.test(Shannon ~ Group, data = shannon_df)
  stat_method <- "wilcox.test"
}
p_val <- base_test$p.value

# === 4. 绘图 ===
p1 <- ggplot(shannon_df, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, size = 2, alpha = 0.7) +
  stat_compare_means(
    method   = stat_method,   # "t.test" 或 "wilcox.test"
    label    = "p.format",    # 预设之一，表示 “p = 0.001”
    label.x  = 1.5,
    label.y  = max(shannon_df$Shannon) * 1.05
  ) +
  scale_fill_manual(values = c(IT = "#E64B35", ST = "#4DBBD5")) +
  scale_color_manual(values = c(IT = "#E64B35", ST = "#4DBBD5")) +
  labs(
    y        = "Shannon Functional Diversity",
    subtitle = paste0(stat_method)  # P 值会自动在图上注
  ) +
  theme_classic(base_size = 14)


# === 5. 保存并展示 ===
ggsave("Shannon_diversity_CSS.png", p1, width = 6, height = 5, dpi = 300)
print(p1)

library(rstatix)
# === 2. 标准化 B 值计算 ===
# 过滤：至少 20% 样本出现
prev <- rowMeans(kegg_css > 0)
abun_filt <- kegg_css[prev >= 0.2, ]

compute_B <- function(mat, samples) {
  n <- length(samples)
  bs <- apply(mat[, samples, drop=FALSE], 1, function(x) {
    if (all(x==0)) return(NA)
    p_raw <- x / sum(x)
    B_raw <- 1 / sum(p_raw^2)
    (B_raw - 1) / (n - 1)
  })
  tibble(genus = names(bs), B = bs)
}

IT_samps <- metadata %>% filter(Group=="IT") %>% pull(Sample)
ST_samps <- metadata %>% filter(Group=="ST") %>% pull(Sample)

B_IT <- compute_B(abun_filt, IT_samps)  %>% mutate(Group="IT")
B_ST <- compute_B(abun_filt, ST_samps) %>% mutate(Group="ST")
B_data <- bind_rows(B_IT, B_ST) %>% drop_na()

write.csv(B_data, "Processed_Niche_Breadth_Data-genus.csv", row.names = FALSE)

cat("\nB值统计摘要（优化后）:\n")
print(B_data %>% group_by(Group) %>% 
        summarise(Median = median(B),
                  Mean = mean(B),
                  SD = sd(B),
                  IQR = IQR(B),
                  Min = min(B),
                  Max = max(B)))
# 组间比较
test_B <- wilcox_test(B ~ Group, data = B_data) %>% add_significance()

# β 多样性 (Bray–Curtis)
dist_mat <- vegdist(t(abun_filt), method="bray")
perm <- adonis2(dist_mat ~ Group, data = metadata, permutations = 9999)
disp <- betadisper(dist_mat, metadata$Group)
disp_test <- permutest(disp, pairwise=TRUE, permutations=9999)

# 可视化 B 值
p2 <- ggplot(B_data, aes(Group, B, fill = Group)) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  stat_compare_means(
    method   = "wilcox.test",
    label    = "p.format",
    label.x  = 1.5,
    label.y  = max(B_data$B) * 1.05
  ) +
  scale_fill_manual(values = c(IT = "#E64B35", ST = "#4DBBD5")) +
  labs(y = "Standardized Niche Breadth (B)") +
  theme_minimal(base_size = 14)

# 效应量
effsize <- B_data %>% wilcox_effsize(B ~ Group) %>% pull(effsize)

p2 <- p2 + annotate(
  "text",
  x     = 1.5,
  y     = max(B_data$B) * 1.1,
  label = paste0("r = ", signif(effsize, 2)),
  size  = 4
)

print(p2)

# === 输出 ===
print(p1)
print(p2)
print(perm)
print(disp_test)



















#分组计算功能生态位宽度（Shannon index，功能多样性）
# 加载必要包
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# 读取数据
kegg_css <- read.csv("CSS_relative_KEGG_bases_lab.csv", row.names = 1) # KO为行，样本为列
metadata <- read.csv("metadata.csv") # 包含Sample和Group列

# 转置数据（vegan要求样本为行）
kegg_t <- t(kegg_css)

# 检查样本顺序是否一致
stopifnot(rownames(kegg_t) == metadata$Sample)

# 计算Shannon多样性（使用CSS标准化后的相对丰度）
shannon <- diversity(kegg_t, index = "shannon")

# 合并结果
div_df <- data.frame(
  Sample = names(shannon),
  Shannon = shannon
) %>%
  left_join(metadata, by = "Sample")

# 检查组间样本量
table(div_df$Group)

# 正态性检验
shapiro.test(div_df$Shannon[div_df$Group == "IT"])
shapiro.test(div_df$Shannon[div_df$Group == "ST"])

# 方差齐性检验
bartlett.test(Shannon ~ Group, data = div_df)

# 选择适当检验（参数/非参数）
if(all(
  shapiro.test(div_df$Shannon[div_df$Group == "IT"])$p.value > 0.05,
  shapiro.test(div_df$Shannon[div_df$Group == "ST"])$p.value > 0.05,
  bartlett.test(Shannon ~ Group, data = div_df)$p.value > 0.05
)) {
  test_result <- t.test(Shannon ~ Group, data = div_df, var.equal = TRUE)
} else {
  test_result <- wilcox.test(Shannon ~ Group, data = div_df)
}

# 保存检验结果
test_method <- ifelse(exists("var.equal"), "Student's t-test", "Wilcoxon test")
p_value <- format.pval(test_result$p.value, digits = 3)


# 颜色设置（符合IT=红，ST=蓝的惯例）
group_colors <- c("IT" = "#E64B35", "ST" = "#4DBBD5")

# 箱线图+散点图
ggplot(div_df, aes(x = Group, y = Shannon, fill = Group)) +
  # 半透明箱线图
  geom_boxplot(
    width = 0.6,
    alpha = 0.3,
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  # 散点（添加抖动避免重叠）
  geom_jitter(
    width = 0.1,
    size = 2.5,
    alpha = 0.7,
    aes(color = Group),
    show.legend = FALSE
  ) +
  # 颜色映射
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  # 统计标注
  stat_compare_means(
    method = ifelse(grepl("t-test", test_method), "t.test", "wilcox.test"),
    label = "p.format",
    label.x = 1.5,
    label.y = max(div_df$Shannon) * 1.05,
    size = 3.5
  ) +
  # 坐标轴标签
  labs(
    x = "",
    y = "Shannon Functional Diversity Index",
    title = "Functional Niche Width Comparison",
    subtitle = paste0("CSS-normalized data (", test_method, ": p = ", p_value, ")")
  ) +
  # 主题美化
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  ) +
  # 扩展y轴空间
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

# 保存图形
ggsave("Shannon_diversity_CSS.pdf", width = 6, height = 5, dpi = 300)
ggsave("Shannon_diversity_CSS.png", width = 6, height = 5, dpi = 300)


cat(sprintf(
  "\nGroup Comparison (IT vs ST):\n- Method: %s\n- p-value: %s\n- IT median: %.2f\n- ST median: %.2f",
  test_method,
  p_value,
  median(div_df$Shannon[div_df$Group == "IT"]),
  median(div_df$Shannon[div_df$Group == "ST"])
))

# PERMANOVA检验（功能β多样性）
adonis2(kegg_t ~ Group, data = metadata, method = "bray")










####################################
# 加载包
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix)  # 用于更丰富的统计检验
library(ggsci)    # 用于科学期刊配色
dir()
# 数据过滤（更严格）
abundance_base <- read.csv("CSS_relative_KEGG_bases_lab.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv") %>% 
  filter(Sample %in% colnames(abundance_base))
# 过滤低丰度菌（0.1%）
abundance_df <- abundance_base[rowSums(abundance_base) > 0.001, ]

##########################################################################
# 检查过滤后保留的物种总丰度占比
total_abundance_before <- sum(abundance_base)
total_abundance_after <- sum(abundance_df)
cat("保留物种的总丰度占比:", round(100 * total_abundance_after / total_abundance_before, 2), "%\n")

# 过滤前后零值比例对比
cat("过滤前零值比例:", mean(abundance_base == 0), "\n")
cat("过滤后零值比例:", mean(abundance_df == 0), "\n")

# 数据验证
cat("\n数据维度验证:\n")
cat("物种表 - 样本数:", ncol(abundance_df), "物种数:", nrow(abundance_df), "\n")
cat("元数据 - 样本数:", nrow(metadata), "\n")

# 样本匹配
metadata <- metadata %>% filter(Sample %in% colnames(abundance_df))
cat("匹配后的有效样本数:", nrow(metadata), "\n")

# 3. 计算标准化B值（优化版）
# ----------------------
cat("\n=== STEP 3: 计算标准化B值（优化版） ===\n")

calculate_standardized_B <- function(abundance, group_samples, group_name) {
  B_values <- apply(abundance[, group_samples, drop = FALSE], 1, function(x) {
    # 关键改进1：包含零丰度样本（反映真实分布）
    p <- x / sum(x)  # 不过滤x > 0，保留零值
    
    # 计算原始B值
    B_raw <- 1 / sum(p^2)
    
    # 关键改进2：Hurlbert标准化（基于生境数n）
    n <- length(group_samples)
    B_standardized <- (B_raw - 1) / (n - 1)  # 标准化到[0,1]
    
    return(B_standardized)
  })
  
  result <- data.frame(
    genus = names(B_values),
    B_standardized = unname(B_values),
    Group = group_name,
    stringsAsFactors = FALSE
  ) %>% 
    filter(!is.na(B_standardized)) %>% 
    arrange(B_standardized)
  
  cat(group_name, "组计算完成，有效物种数:", nrow(result), "\n")
  return(result)
}

# 分组样本（保持不变）
IT_samples <- metadata %>% filter(Group == "IT") %>% pull(Sample)
ST_samples <- metadata %>% filter(Group == "ST") %>% pull(Sample)

# 计算B值
B_IT <- calculate_standardized_B(abundance_df, IT_samples, "IT")
B_ST <- calculate_standardized_B(abundance_df, ST_samples, "ST")

# 合并数据
B_data <- bind_rows(B_IT, B_ST) %>% 
  mutate(Group = factor(Group, levels = c("IT", "ST")))

# 保存处理后的数据
write.csv(B_data, "Processed_Niche_Breadth_Data-genus.csv", row.names = FALSE)

cat("\nB值统计摘要（优化后）:\n")
print(B_data %>% group_by(Group) %>% 
        summarise(Median = median(B_standardized),
                  Mean = mean(B_standardized),
                  SD = sd(B_standardized),
                  IQR = IQR(B_standardized),
                  Min = min(B_standardized),
                  Max = max(B_standardized)))

# ----------------------
# 4. 组间比较（更全面的统计检验）
# ----------------------
cat("\n=== STEP 4: 组间比较（优化版） ===\n")

# 正态性检验
cat("\n正态性检验（Shapiro-Wilk test）:\n")
print(B_data %>% group_by(Group) %>% shapiro_test(B_standardized))

# 方差齐性检验
cat("\n方差齐性检验（Levene's test）:\n")
print(B_data %>% levene_test(B_standardized ~ Group))

# 非参数检验（Wilcoxon rank-sum test）
cat("\nWilcoxon rank-sum test结果:\n")
wilcox_result <- B_data %>% 
  wilcox_test(B_standardized ~ Group) %>% 
  add_significance()
print(wilcox_result)

# 效应量计算
cat("\n效应量（r）计算:\n")
effect_size <- B_data %>% 
  wilcox_effsize(B_standardized ~ Group)
print(effect_size)

# ----------------------
# 5. 差异分析（更全面的β多样性分析）
# ----------------------
cat("\n=== STEP 5: β多样性分析 ===\n")

# 计算Bray-Curtis距离
dist_matrix <- vegdist(t(abundance_df), "bray")

# PERMANOVA分析（增加置换检验次数）
set.seed(123)
adonis_result <- adonis2(dist_matrix ~ Group, data = metadata, permutations = 9999)
cat("\nPERMANOVA结果:\n")
print(adonis_result)

# 组间离散度检验（β离散度）
dispersion <- betadisper(dist_matrix, metadata$Group)
cat("\n组间离散度检验结果:\n")
print(permutest(dispersion, pairwise = TRUE, permutations = 9999))

# ----------------------
# 6. 科学可视化（符合SCI顶刊格式）
# ----------------------
# 自定义颜色
group_colors <- c("IT" = "#E64B35", "ST" = "#4DBBD5")
#重置图形系统
try(dev.off(), silent = TRUE)
graphics.off()
# 箱线图+小提琴图
# 科学可视化（优化点颜色与组别匹配）
p_main <- ggplot(B_data, aes(x = Group, y = B_standardized, fill = Group, color = Group)) +
  # 小提琴图（填充色半透明，边框色深）
  geom_violin(alpha = 0.7, width = 0.8, trim = TRUE) +  # 移除小提琴边框
  # 箱线图（白色填充，细边框）
  geom_boxplot(
    width = 0.15, 
    alpha = 0.8, 
    outlier.shape = NA,
    color = "black"  # 箱线图边框色
       # 箱线图填充色
  ) +
  # 散点（使用半透明且与组别匹配的颜色）
  geom_jitter(
    width = 0.1, 
    size = 1.5, 
    alpha = 0.4,      # 适当降低透明度
    aes(color =Group ) # 颜色与组别匹配
  ) +
  # 统计检验标注
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    label.y = max(B_data$B_standardized) * 1.01,
    size = 4.5,
    vjust = 0.5,
    label.x = 1.5      # 将p值居中显示
  ) +
  # 颜色设置（使用您指定的颜色）
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +  # 点颜色与填充色一致
  labs(
    title = "Niche Width Comparison",
    y = "Standardized Niche Breadth (B)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text = element_text(color = "black"),
    axis.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5),
    plot.margin = unit(c(10, 10, 10, 10), "points")  # 增加边距
  )
print(p_main)
# 添加效应量标注（优化位置和颜色）
p_main <- p_main + 
  annotate("text", 
           x = 1.5, 
           y = max(B_data$B_standardized) * 1.1,
           label = paste0("Effect size (r) = ", round(effect_size$effsize, 3)),
           size = 4.5,
           color = "black",  # 使用黑色保证可读性
           fontface = "italic")

# 显示优化后的图形
print(p_main)

###################################################################
##################################################################
#相关检验
##################################################################
##################################################################
# 计算Shannon生态位宽度（验证趋势是否一致）
shannon_B <- apply(abundance_df, 1, function(x) {
  p <- x[x > 0] / sum(x)
  -sum(p * log(p))
})

# 比较与Levins' B的相关性
cor(B_data$B_standardized, shannon_B[B_data$genus], use = "complete.obs")


# 1. 绘制散点图观察关系
df_cor <- data.frame(
  Levins_B = B_data$B_standardized,
  Shannon = shannon_B[B_data$genus]
) %>% na.omit()

ggplot(df_cor, aes(x = Levins_B, y = Shannon)) +
  geom_point(alpha = 0.6, color = "#0072B2") +
  geom_smooth(method = "lm", color = "#D55E00", se = FALSE) +
  labs(
    x = "Standardized Levins' B",
    y = "Shannon Niche Width",
    title = paste("Correlation: r =", round(cor(df_cor$Levins_B, df_cor$Shannon), 3))
  ) +
  theme_minimal()

# 2. 检查是否存在非线性关系（如Spearman相关）
cor.test(B_data$B_standardized, shannon_B[B_data$genus], 
         method = "spearman", use = "complete.obs")
##############################################################
library(ggpubr)
library(ggrepel)

# 带物种标签的散点图
ggplot(df_cor, aes(x = Levins_B, y = Shannon)) +
  geom_point(aes(color = ifelse(Levins_B > 0.8 | Shannon > 2, "Highlight", "Normal")), 
             size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", linetype = "dashed") +
  geom_text_repel(
    data = subset(df_cor, Levins_B > 0.8 | Shannon > 2),
    aes(label = rownames(subset(df_cor, Levins_B > 0.8 | Shannon > 2))),
    size = 3, box.padding = 0.5
  ) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(
    x = "Standardized Levins' B", 
    y = "Shannon Niche Width",
    caption = paste("Pearson r =", round(cor(df_cor$Levins_B, df_cor$Shannon), 3))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")



###################################################################
#查看零值比例
zero_ratio <- apply(abundance_df[, IT_samples], 1, function(x) sum(x == 0) / length(x))
hist(zero_ratio, main = "Proportion of IT group bacteria with zero values")


zero_ratio <- apply(abundance_df[, ST_samples], 1, function(x) sum(x == 0) / length(x))
hist(zero_ratio, main = "Proportion of ST group bacteria with zero values")

#可视化代码
# 1. B值分布箱线图
ggplot(B_data, aes(x = Group, y = B_standardized, fill = Group)) +
  geom_boxplot() +
  labs(title = "Niche Breadth (B) by Habitat", y = "Standardized B")

#可视化验证
#菌属B值分布直方图
ggplot(B_data, aes(x = B_standardized, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  geom_vline(xintercept = c(0.1, 0.5), linetype = "dashed", color = "red") +
  labs(x = "Standardized B", y = "Count", title = "Distribution of Niche Breadth") +
  facet_wrap(~Group, scales = "free_y")
# 计算每个菌属在各组的零值比例
zero_ratio_data <- metadata %>%
  split(.$Group) %>%  # 按组拆分
  map_df(~ {
    group_samples <- .x$Sample
    apply(abundance_df[, group_samples], 1, function(x) {
      sum(x == 0) / length(x)  # 零值比例
    }) %>%
      data.frame(
        Species = names(.),
        Zero_Ratio = .,
        Group = unique(.x$Group),
        stringsAsFactors = FALSE
      )
  }, .id = "Group")

#可视化零值比例分布
library(ggplot2)
# 绘制直方图
ggplot(zero_ratio_data, aes(x = Zero_Ratio, fill = Group)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  facet_wrap(~Group) +
  labs(title = "Zero Abundance Ratio Distribution by Group",
       x = "Proportion of Samples with Zero Abundance",
       y = "Number of Taxa")

###########################################################################
ggplot(df_cor, aes(x = Levins_B, y = Shannon)) +
  geom_point(aes(color = ifelse(Levins_B > quantile(Levins_B, 0.9) | 
                                  Shannon > quantile(Shannon, 0.9), 
                                "Top 10%", "Normal")), 
             size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "grey40", se = TRUE) +  # 添加置信区间
  geom_text_repel(
    data = subset(df_cor, 
                  Levins_B > quantile(Levins_B, 0.9) | 
                    Shannon > quantile(Shannon, 0.9)),
    aes(label = rownames(subset(df_cor, 
                                Levins_B > quantile(Levins_B, 0.9) | 
                                  Shannon > quantile(Shannon, 0.9)))),
    size = 3, 
    box.padding = 0.5,
    max.overlaps = 50  # 增加可显示标签数
  ) +
  scale_color_manual(values = c("Top 10%" = "#E69F00", "Normal" = "#56B4E9")) +
  labs(
    x = "Standardized Levins' B", 
    y = "Shannon Niche Width",
    caption = paste("Pearson r =", round(cor(df_cor$Levins_B, df_cor$Shannon), 3),
                    "\nSpearman ρ = 0.753 (p < 0.001)")
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none")
#############################################################

##############################################################

         
         
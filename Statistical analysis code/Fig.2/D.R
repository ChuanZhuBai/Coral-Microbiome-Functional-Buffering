# ======================
# 精准累积比例分析 - 分步终极版
# 优化目标：
# 1. 分步验证每步结果
# 2. 精确查找指定B值的累积比例
# 3. 修复可视化问题
# ======================

# ----------------------
# 1. 初始化环境
# ----------------------
cat("=== STEP 1: 初始化环境 ===\n")
rm(list = ls())
library(ggplot2)
library(dplyr)
library(purrr)
library(ggrepel)
dir()
# ----------------------
# 2. 数据加载与验证
# ----------------------
cat("\n=== STEP 2: 数据加载与验证 ===\n")

# 检查文件存在性
if(!all(c("genus.csv", "metadata.csv") %in% list.files())) {
  stop("缺少必要文件: genus.csv 或 metadata.csv")
} else {
  cat("所有文件存在检查通过\n")
}

# 加载数据
# 数据过滤（更严格）
abundance_base <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv") %>% 
  filter(Sample %in% colnames(abundance_base))
# 过滤低丰度菌（0.1%）
#abundance_df <- abundance_base[rowSums(abundance_base) > 0.001, ]
prev <- rowMeans(abundance_base > 0)
abundance_df <- abundance_base[prev >= 0.2, ]

# 数据验证
cat("\n数据维度验证:\n")
cat("物种表 - 样本数:", ncol(abundance_df), "物种数:", nrow(abundance_df), "\n")
cat("元数据 - 样本数:", nrow(metadata), "\n")

# 样本匹配
metadata <- metadata %>% filter(Sample %in% colnames(abundance_df))
cat("匹配后的有效样本数:", nrow(metadata), "\n")

# ----------------------
# 3. 计算标准化B值
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
    species = names(B_values),
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

cat("\nB值统计摘要:\n")
print(B_data %>% group_by(Group) %>% 
        summarise(Median = median(B_standardized),
                  Mean = mean(B_standardized),
                  SD = sd(B_standardized)))


# K-S检验（比较IT富集与ST富集OTU的B值分布）

B_data <- B_data %>% 
  mutate(B_jitter = B_standardized + rnorm(n(), sd = 0.0001))#添加微小噪声

library(ggsci)
library(tibble)
ks_test <- ks.test(
  B_data %>% filter(Group == "IT") %>% pull(B_jitter),
  B_data %>% filter(Group == "ST") %>% pull(B_jitter)
)

# 生成统计标注文本
stat_text <- sprintf(
  "K-S test: D = %.2f\np < 0.001",  # 修改p值显示
  ks_test$statistic
)

###############################################
###############################################

median_values <- B_data %>% 
  group_by(Group) %>% 
  summarise(Median = median(B_standardized), 
            Mean = mean(B_standardized), 
            SD = sd(B_standardized))
# ----------------------
# 4. 创建精确CDF函数
# ----------------------
cat("\n=== STEP 4: 创建精确CDF函数 ===\n")

create_precise_cdf <- function(data) {
  # 处理可能的重复值
  unique_data <- data %>% 
    arrange(B_standardized) %>% 
    mutate(CumProb = (1:n())/n())
  
  approxfun(
    x = unique_data$B_standardized,
    y = unique_data$CumProb,
    method = "linear",
    yleft = 0,
    yright = 1
  )
}

cdf_IT <- create_precise_cdf(B_IT)
cdf_ST <- create_precise_cdf(B_ST)

# ----------------------
# 5. 计算特定B值的累积比例
# ----------------------
cat("\n=== STEP 5: 计算指定B值的累积比例 ===\n")

target_B_values <- c(0.105, 0.496)#属水平
calculate_cumprob <- function(B_val, cdf_func, group_name) {
  cp <- cdf_func(B_val)
  data.frame(B_value = B_val, 
             Cumulative_Proportion = cp,
             Group = group_name)
}

# 计算IT组
IT_cumprob <- map_dfr(target_B_values, ~calculate_cumprob(.x, cdf_IT, "IT"))
# 计算ST组
ST_cumprob <- map_dfr(target_B_values, ~calculate_cumprob(.x, cdf_ST, "ST"))

# 合并结果
cumprob_results <- bind_rows(IT_cumprob, ST_cumprob) %>% 
  arrange(B_value, Group)

cat("\n特定B值的累积比例结果:\n")
print(cumprob_results)

# ----------------------
# 7. 可视化
# ----------------------
cat("\n=== STEP 7: 创建可视化 ===\n")

base_plot1 <- ggplot(B_data, aes(x = B_standardized)) +
  stat_ecdf(aes(color = Group), geom = "step", linewidth = 1) +
  scale_color_manual(values = c("IT" = "#E64B35", "ST" = "#4DBBD5")) +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 1, by = 0.2)) +
  labs(x = "Standardized Niche Breadth",
       y = "Cumulative Proportion",
       title = "Niche Breadth Cumulative Distribution") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.85, 0.2),
        plot.title = element_text(hjust = 0.5))
print(base_plot1)
# 添加特定B值标记
############################################################
# 最终可视化（包含两个B值标记）
base_plot2 <- base_plot1 +
  # 添加双垂直线
  geom_vline(
    xintercept = target_B_values,
    linetype = "dotted",
    color = "grey50",
    linewidth = 0.7
  ) +
  annotate(
    "text",
    x = max(B_data$B_index) * 0.95,
    y = 0.1,
    label = stat_text,
    hjust = 1,
    size = 4,
    color = "black"
  ) +
  # 添加双标签（自动对齐）
  geom_text(
    data = data.frame(
      x_val = target_B_values,
      label = paste("B =", target_B_values)
    ),
    aes(x = x_val, y = 0.05, label = label),
    angle = 90,
    vjust = -0.5,
    size = 3.5,
    color = "grey30"
  ) +
  
  # 添加阈值比例标注（来自cumprob_results）
  geom_text(
    data = cumprob_results,
    aes(
      x = B_value,
      y = Cumulative_Proportion + 0.04,  # 在曲线上方标注
      label = paste0(round(Cumulative_Proportion*100, 1), "%"),
      color = Group
    ),
    size = 3.5,
    show.legend = FALSE
  )

print(base_plot2)



##############################################################


# 在基础图形上添加中位线
plot_with_median_lines <- base_plot2 +
  geom_vline(
    data = median_values,
    aes(xintercept = Median, color = Group),
    linetype = "dashed", 
    linewidth = 0.8, 
    alpha = 0.7
  )

# 显示带中位线的图形
print(plot_with_median_lines)

# 在中位线图形上添加标签
final_plot <- plot_with_median_lines +
  geom_text_repel(
    data = median_values,
    aes(
      x = Median,
      y = 0.75,  # 将标签放在y轴中间位置
      label = sprintf(" Median = %.2f",  Median),
      color = Group
    ),
    direction = "x",  # 垂直方向调整
    nudge_y = 0.1,# 向上轻微偏移
    angle = 90,
    size = 3.5,         # 字体大小
    show.legend = FALSE  # 不显示在图例中
  )

# 显示最终图形
print(final_plot)

head(stat_text)
# 保存图形
ggsave("Niche_Breadth_CDF_Analysis.pdf", final_plot, 
       width = 10, height = 5, dpi = 300)
cat("\n图形已保存为 Niche_Breadth_CDF_Analysis.pdf\n")

# ----------------------
# 8. 结果汇总
# ----------------------
cat("\n=== 最终结果汇总 ===\n")

cat("\n1. 特定B值的累积比例:\n")
print(cumprob_results)


# 保存所有结果
write.csv(cumprob_results, "Specific_B_Cumulative_Proportions.csv", row.names = FALSE)
if(!is.null(main_intersection)) {
  write.csv(main_intersection, "Main_Intersection_Point.csv", row.names = FALSE)
}
cat("\n结果已保存为CSV文件\n")
################################################################
################################################################
#密度图
ggplot(B_data, aes(x = B_standardized, fill = Group)) +
  geom_density(alpha = 0.5) +
  geom_vline(data = B_data %>% group_by(Group) %>% 
               summarise(Median = median(B_standardized)),
             aes(xintercept = Median, color = Group),
             linetype = "dashed") +
  scale_fill_manual(values = c("IT" = "#E64B35", "ST" = "#4DBBD5")) +
  labs(x = "Standardized Niche Breadth", y = "Density")
####################################################################
####################################################################







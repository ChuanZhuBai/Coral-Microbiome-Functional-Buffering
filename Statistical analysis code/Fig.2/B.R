###R代码：生态位宽度区间比例小提琴图
### R代码：分面生态位宽度比较图（顶刊标准版）
# 0. 初始化设置 ------------------------------------------------------------
rm(list = ls())
options(stringsAsFactors = FALSE)

# 1. 数据准备 -------------------------------------------------------------
cat("\n=== 1. 数据准备 ===\n")

# 1.1 加载必要包
required_packages <- c("ggplot2", "ggpubr", "rstatix", "dplyr", "tidyr")
invisible(lapply(required_packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}))
dir()

# 1.2 数据读取与处理
metadata <- read.csv("metadata.csv", header = TRUE)
b_data <- read.csv("Processed_Niche_Breadth_Data-insitu.csv", header = TRUE)

combined_data <- b_data %>%
  left_join(metadata, by = c("Group" = "Sample")) %>%
  select(genus,  B_standardized, Group = Group) %>%
  filter(!is.na(B_standardized))

# 1.3 生态位宽度分类
breaks_quantile <- quantile(combined_data$B_standardized, probs = c(0.25, 0.75))
combined_data$B_category <- cut(combined_data$B_standardized,
                                breaks = c(-Inf, breaks_quantile, Inf),
                                labels = c("Low", "Medium", "High"))

# 2. 分步绘制图形 --------------------------------------------------------
cat("\n=== 2. 分步绘制图形 ===\n")

# 2.1 定义通用绘图参数
my_colors <- c("IT" = "#E64B35", "ST" = "#4DBBD5") # Nature常用配色
base_size <- 12 # 基础字体大小

# 2.2 Low组绘图
cat("\n2.2 绘制Low组图形...\n")
low_data <- combined_data %>% filter(B_category == "Low")

library(rstatix)
# 确保使用rstatix包的函数
low_stats <- rstatix::wilcox_test(
  data = low_data, 
  formula = B_standardized ~ Group
) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Group")

# 效应量
low_effsize <- wilcox_effsize(low_data, B_standardized ~ Group)

# 基础图形
p1 <- ggplot(low_data, aes(x = Group, y = B_standardized)) +
  geom_violin(aes(fill = Group), alpha = 0.6, trim = TRUE, width = 0.7) +
  geom_boxplot(aes(fill = Group), width = 0.15, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(title = "Low Niche Breadth", 
       y = "Standardized Niche Breadth (B)",
       x = "") +
  theme_classic(base_size = base_size) +
  theme(plot.title = element_text( hjust = 0.5),
        legend.position = "none")

print(p1)
# 添加统计标注
p1 <- p1 + stat_pvalue_manual(low_stats, 
                              label = "p = {p} ({p.signif})",
                              tip.length = 0.01,
                              bracket.size = 0.3,
                              size = 3.5,
                              y.position = max(low_data$B_standardized) * 1.05)

# 添加效应量
p1 <- p1 + annotate("text", x = 1.5, 
                    y = max(low_data$B_standardized) * 1.12,
                    label = sprintf("r = %.2f", low_effsize$effsize),
                    size = 3.5)

# 添加样本量
sample_sizes <- table(low_data$Group)
p1 <- p1 + annotate("text", x = 1:2,
                    y = min(low_data$B_standardized),
                    label = paste0("n=", sample_sizes),
                    vjust = 2, size = 3.5)

print(p1)

# 2.3 Medium组绘图
cat("\n2.3 绘制Medium组图形...\n")
medium_data <- combined_data %>% filter(B_category == "Medium")

# 统计检验
medium_stats <- rstatix::wilcox_test(medium_data, B_standardized ~ Group) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Group")

# 效应量
medium_effsize <- wilcox_effsize(medium_data, B_standardized ~ Group)

# 基础图形
p2 <- ggplot(medium_data, aes(x = Group, y = B_standardized)) +
  geom_violin(aes(fill = Group), alpha = 0.6, trim = TRUE, width = 0.7) +
  geom_boxplot(aes(fill = Group), width = 0.15, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(title = "Medium Niche Breadth", 
       y = "Standardized Niche Breadth (B)",
       x = "") +
  theme_classic(base_size = base_size) +
  theme(plot.title = element_text( hjust = 0.5),
        legend.position = "none")

print(p2)
# 添加统计标注
p2 <- p2 + stat_pvalue_manual(medium_stats, 
                              label = "p = {p} ({p.signif})",
                              tip.length = 0.01,
                              bracket.size = 0.3,
                              size = 3.5,
                              y.position = max(medium_data$B_standardized) * 1.05)

# 添加效应量
p2 <- p2 + annotate("text", x = 1.5, 
                    y = max(medium_data$B_standardized) * 1.12,
                    label = sprintf("r = %.2f", medium_effsize$effsize),
                    size = 3.5)

# 添加样本量
sample_sizes <- table(medium_data$Group)
p2 <- p2 + annotate("text", x = 1:2,
                    y = min(medium_data$B_standardized),
                    label = paste0("n=", sample_sizes),
                    vjust = 2, size = 3.5)

print(p2)

# 2.4 High组绘图
cat("\n2.4 绘制High组图形...\n")
high_data <- combined_data %>% filter(B_category == "High")

# 统计检验
high_stats <- rstatix::wilcox_test(high_data, B_standardized ~ Group) %>%
  rstatix::add_significance() %>%
  rstatix::add_xy_position(x = "Group")

# 效应量
high_effsize <- wilcox_effsize(high_data, B_standardized ~ Group)

# 基础图形
p3 <- ggplot(high_data, aes(x = Group, y = B_standardized)) +
  geom_violin(aes(fill = Group), alpha = 0.6, trim = TRUE, width = 0.7) +
  geom_boxplot(aes(fill = Group), width = 0.15, outlier.shape = NA) +
  geom_jitter(aes(color = Group), width = 0.1, alpha = 0.5, size = 1.5) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  labs(title = "High Niche Breadth", 
       y = "Standardized Niche Breadth (B)",
       x = "") +
  theme_classic(base_size = base_size) +
  theme(plot.title = element_text( hjust = 0.5),
        legend.position = "none")

print(p3)
# 添加统计标注
p3 <- p3 + stat_pvalue_manual(high_stats, 
                              label = "p = {p} ({p.signif})",
                              tip.length = 0.01,
                              bracket.size = 0.3,
                              size = 3.5,
                              y.position = max(high_data$B_standardized) * 1.05)

# 添加效应量
p3 <- p3 + annotate("text", x = 1.5, 
                    y = max(high_data$B_standardized) * 1.12,
                    label = sprintf("r = %.2f", high_effsize$effsize),
                    size = 3.5)

# 添加样本量
sample_sizes <- table(high_data$Group)
p3 <- p3 + annotate("text", x = 1:2,
                    y = min(high_data$B_standardized),
                    label = paste0("n=", sample_sizes),
                    vjust = 2, size = 3.5)

print(p3)

# 3. 图形合并与导出 ------------------------------------------------------
cat("\n=== 3. 图形合并与导出 ===\n")

# 3.1 合并图形
final_plot <- ggarrange(p1, p2, p3,
                        ncol = 3,
                        common.legend = TRUE,
                        legend = "bottom",
                        labels = "AUTO")

# 3.2 添加全局标题和注释
final_plot <- annotate_figure(final_plot,
                              top = text_grob("Comparison of Niche Breadth Between Groups", 
                                               size = 13),
                              bottom = text_grob(paste("Classification thresholds:", 
                                                       round(breaks_quantile[1], 3), "(Low/Medium),",
                                                       round(breaks_quantile[2], 3), "(Medium/High)"),
                                                 color = "gray40", size = 10))

print(final_plot)

# 3.3 导出图形
ggsave("Niche_Breadth_Comparison_Separate.tiff",
       plot = final_plot,
       device = "tiff",
       dpi = 600,
       width = 12,
       height = 8,
       compression = "lzw")

cat("\n=== 分析完成 ===\n")
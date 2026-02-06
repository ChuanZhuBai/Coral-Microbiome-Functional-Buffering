library(tidyverse)
library(ggpubr)
library(patchwork)

# ==========================================
# 1. 读入预计算的数据
# ==========================================
# 1.1 读入全局功能冗余数据 (由之前代码生成)
fri_df <- read.csv("Global_FRI_per_sample.csv")

# 1.2 读入不同阈值下的广适者丰富度数据
# 假设列名为: Sample, Generalist_Richness
rich_10 <- read.csv("Sample_Generalist_Richness_10%.csv")
rich_20 <- read.csv("Sample_Generalist_Richness_20%.csv")

# ==========================================
# 2. 数据合并与准备
# ==========================================
# 合并 10% 数据
df_10 <- fri_df %>%
  inner_join(rich_10 %>% select(Sample, Generalist_Richness), by = "Sample") %>%
  mutate(Threshold = "Top 10% Threshold")

# 合并 20% 数据
df_20 <- fri_df %>%
  inner_join(rich_20 %>% select(Sample, Generalist_Richness), by = "Sample") %>%
  mutate(Threshold = "Top 20% Threshold")

# ==========================================
# 3. 定义统一的绘图函数 (顶级期刊风格)
# ==========================================
plot_robustness <- function(data, x_label, title_lab) {
  ggplot(data, aes(x = Generalist_Richness, y = Global_FRI)) +
    # 拟合线与置信区间
    geom_smooth(method = "lm", color = "black", fill = "gray90", alpha = 0.6, size = 1) +
    # 散点：IT红色，ST蓝色
    geom_point(aes(fill = Group), shape = 21, size = 4, stroke = 0.3, alpha = 0.8) +
    scale_fill_manual(values = c("IT" = "#E64B35FF", "ST" = "#4DBBD5FF")) +
    # 统计标注：R2 和 P值
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             label.x.npc = "left", label.y.npc = "top", size = 4.5) +
    labs(
      title = title_lab,
      x = x_label,
      y = "Mean Global FRI"
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.position = "none"
    )
}

# ==========================================
# 4. 绘制两个分面的图
# ==========================================
p1 <- plot_robustness(df_10, "Richness of Generalist Genera (Top 10%)", "A. Strict Classification")
p2 <- plot_robustness(df_20, "Richness of Generalist Genera (Top 20%)", "B. Moderate Classification")

# ==========================================
# 5. 组合图形并导出
# ==========================================
# 使用 patchwork 组合，并添加总标题
final_s6 <- p1 + p2 + 
  plot_annotation(
    title = "Sensitivity Analysis: Generalist Richness vs. Mean Functional Redundancy",
    subtitle = "Validation across different niche breadth thresholds",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, face = "italic")
    )
  )

print(final_s6)

# 导出符合出版要求的 PDF
ggsave("Figure_S6_Robustness_Analysis.pdf", final_s6, width = 10, height = 5)

# 导出合并后的汇总表供 SI 使用
summary_robustness <- bind_rows(df_10, df_20)
write.csv(summary_robustness, "Table_S_Robustness_Data_Combined.csv", row.names = FALSE)

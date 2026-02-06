library(tidyverse)
library(ggpubr)

# ==========================================
# 1. 数据读取与基础计算
# ==========================================
niche_data <- read.csv("Processed_Niche_Breadth_Data-insitu.csv")
genus_abundance <- read.csv("genus.csv")
fri_data <- read.csv("FRI_per_sample.csv")
metadata <- read.csv("metadata.csv")
paired_data <- read.csv("paired_data.csv")

# 识别 Generalists (B > 0.496)
generalists_list <- niche_data %>%
  filter(B_standardized > 0.496) %>%
  pull(genus)

# 计算每个样本中“广适者”的种类数 (Richness)
gen_richness_summary <- genus_abundance %>%
  pivot_longer(cols = -genus, names_to = "Sample", values_to = "Abundance") %>%
  filter(genus %in% generalists_list, Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(Generalist_Richness = n_distinct(genus))

# 计算每个样本的平均功能冗余度 (FRI)
sample_fri_avg <- fri_data %>%
  group_by(Sample) %>%
  summarise(Mean_FRI = mean(FRI_eff))

# 计算 Fv/Fm Retention Rate
fvfm_summary <- paired_data %>%
  mutate(Retention_Rate = 100 + ChangeRate) %>%
  select(Sample, Retention_Rate)

# 合并所有数据到主表
master_df <- gen_richness_summary %>%
  left_join(sample_fri_avg, by = "Sample") %>%
  left_join(fvfm_summary, by = "Sample") %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Group = factor(Group, levels = c("IT", "ST"))) %>%
  na.omit()

# ==========================================
# 2. 绘图函数定义 (统一 SA 风格)
# ==========================================
plot_regression <- function(data, x_var, y_var, x_label, y_label, title) {
  ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_smooth(method = "lm", color = "black", fill = "gray90", alpha = 0.6, size = 1) +
    geom_point(aes(fill = Group), shape = 21, size = 5, stroke = 0.5, alpha = 0.8) +
    scale_fill_manual(values = c("IT" = "#E64B35FF", "ST" = "#4DBBD5FF")) +
    # 同时标注 R2 和 P
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
             label.x.npc = "left", label.y.npc = "top", size = 5) +
    labs(x = x_label, y = y_label, title = title) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
}

# ==========================================
# 3. 生成三个核心图表
# ==========================================

# Figure 2G: 解释冗余来源 (Diversity -> Redundancy)
p2g <- plot_regression(master_df, "Generalist_Richness", "Mean_FRI",
                       "Richness of Generalist Genera", "Microbial Functional Redundancy (FRI)",
                       "A. Diversity Supports Redundancy")

# Figure 3E: 解释保护机制 (Redundancy -> Resilience)
p3e <- plot_regression(master_df, "Mean_FRI", "Retention_Rate",
                       "Microbial Functional Redundancy (FRI)", "Coral Fv/Fm Retention Rate (%)",
                       "B. Redundancy Underpins Resilience")

# Figure 3F: 终极预测结果 (Diversity -> Resilience)
p3f <- plot_regression(master_df, "Generalist_Richness", "Retention_Rate",
                       "Richness of Generalist Genera", "Coral Fv/Fm Retention Rate (%)",
                       "C. Diversity Predicts Thermotolerance")

# 展示并保存
print(p2g)
print(p3e)
print(p3f)

ggsave("Fig2G_Richness_vs_FRI.pdf", p2g, width = 6, height = 5)
ggsave("Fig3E_FRI_vs_FvFm.pdf", p3e, width = 6, height = 5)
ggsave("Fig3F_Richness_vs_FvFm.pdf", p3f, width = 6, height = 5)

# 导出 master_df 供分析使用
write.csv(master_df, "Master_Correlation_Data.csv", row.names = FALSE)

# 加载必要包
library(tidyverse)
library(vegan)
library(ggpubr)
library(ggsci)
library(rstatix) 

# 1. 读取数据
genus_data <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv")

# 2. 计算 Shannon 指数
# vegan::diversity 要求样本为行，物种为列，所以需要转置 t()
shannon_values <- diversity(t(genus_data), index = "shannon")

# 3. 整合数据表
alpha_df <- data.frame(
  Sample = names(shannon_values),
  Shannon = as.numeric(shannon_values)
) %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Habitat = factor(Habitat, levels = c("Intertidal", "Subtidal")))

# 4. 统计检验 (Wilcoxon test)
stat_test <- alpha_df %>%
  wilcox_test(Shannon ~ Habitat) %>%
  add_significance() %>%
  add_y_position()

# 5. 绘图 (修改后：仅 Boxplot + Jitter)
p_shannon <- ggplot(alpha_df, aes(x = Habitat, y = Shannon, fill = Habitat)) +
  # --- 已删除 geom_violin ---
  # 箱线图：展现分位数 (宽度调大至 0.5)
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.8, color = "black") +
  # 抖动点：展现真实样本分布
  geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
  # 颜色方案
  scale_fill_manual(values = c("Intertidal" = "#E64B35FF", "Subtidal" = "#4DBBD5FF")) +
  # 显著性标注
  stat_pvalue_manual(stat_test, label = "p.signif", tip.length = 0.01, size = 6) +
  labs(
    x = NULL,
    y = "Shannon Index (Genus level)",
    title = "Bacterial Alpha Diversity"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 11, color = "black"),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 显示结果
print(p_shannon)

# 6. 保存图片
ggsave("Figure_Genus_Shannon_NoViolin.pdf", p_shannon, width = 4, height = 4)

# 7. 导出数据表
write.csv(alpha_df, "Alpha_Diversity_Genus_Level.csv", row.names = FALSE)

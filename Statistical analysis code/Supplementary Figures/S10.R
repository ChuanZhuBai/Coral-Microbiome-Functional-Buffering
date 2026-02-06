library(tidyverse)
library(ggpubr)
library(rstatix)

# ==========================================
# 1. 数据准备
# ==========================================
ko_data <- read.csv("CSS_relative_KEGG_bases_lab.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv")

# 1.1 定义 Hub 基因及其严格的显示顺序
hub_genes_info <- data.frame(
  KO = c("K00836", "K02566", "K02703", "K05373", "K09835", "K07305", 
         "K05501", "K02315", "K01903", "K10235", "K14127", "K03594", "K00992"),
  Gene = c("ectB", "nagD", "psbA", "ybtX", "crtISO", "msrB", 
           "slmA", "dnaC", "sucC", "aglK", "mvhD", "bfr", "murU")
)

# 1.2 提取并转换长表
hub_data <- ko_data[rownames(ko_data) %in% hub_genes_info$KO, ]
hub_long <- hub_data %>%
  rownames_to_column("KO") %>%
  pivot_longer(-KO, names_to = "Sample", values_to = "Abundance") %>%
  left_join(metadata, by = "Sample") %>%
  left_join(hub_genes_info, by = "KO")

# ==========================================
# 2. 【核心修改】设置因子水平（Factor Levels）
# ==========================================
# 2.1 设置 X 轴顺序：IT对照 -> IT热 -> ST对照 -> ST热
group_order <- c("IT-C", "IT-H", "ST-C", "ST-H")
hub_long$Group <- factor(hub_long$Group, levels = group_order)

# 2.2 设置分面顺序：按照 hub_genes_info 定义的顺序
hub_long$Gene <- factor(hub_long$Gene, levels = hub_genes_info$Gene)

# ==========================================
# 3. 统计分析 (解决 0 值报错的稳健循环)
# ==========================================
comparisons_list <- list(c("IT-C", "IT-H"), c("ST-C", "ST-H"), c("IT-H", "ST-H"), c("IT-C", "ST-C"))
all_stats <- data.frame()

for (g in levels(hub_long$Gene)) {
  sub_df <- hub_long %>% filter(Gene == g)
  for (comp in comparisons_list) {
    test_data <- sub_df %>% filter(Group %in% comp)
    try({
      # 如果数据全为 0 且无变化，设 p=1，否则做 Wilcoxon
      if(length(unique(test_data$Abundance)) == 1) {
        p_val <- 1
      } else {
        res <- wilcox.test(Abundance ~ Group, data = test_data, exact = FALSE)
        p_val <- res$p.value
      }
      all_stats <- rbind(all_stats, data.frame(Gene = g, group1 = comp[1], group2 = comp[2], p = p_val))
    }, silent = TRUE)
  }
}

# 多重校正 (FDR)
all_stats <- all_stats %>%
  group_by(group1, group2) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  add_significance("p.adj") %>%
  # 也要给统计表的因子排序，确保与图对应
  mutate(Gene = factor(Gene, levels = hub_genes_info$Gene))

# ==========================================
# 4. 绘图：Figure S10
# ==========================================
# 自动计算显著性标注的位置
stat_for_plot <- all_stats %>%
  add_y_position(data = hub_long, formula = Abundance ~ Group, step.increase = 0.12)

# 定义顶刊配色
color_pal <- c("IT-C" = "#FC8D62", "IT-H" = "#E41A1C", "ST-C" = "#4DBBD5", "ST-H" = "#377EB8")

p_hubs <- ggplot(hub_long, aes(x = Group, y = Abundance, fill = Group)) +
  facet_wrap(~Gene, scales = "free_y", ncol = 4) +
  geom_boxplot(width = 0.6, alpha = 0.8, outlier.shape = NA, color = "black", size = 0.4) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1.2, color = "grey20") +
  scale_fill_manual(values = color_pal) +
  # 标注显著性 (隐藏非显著 ns)
  stat_pvalue_manual(stat_for_plot, label = "p.adj.signif", hide.ns = TRUE, tip.length = 0.01) +
  labs(
    title = "Relative Abundance of Identified Hub Genes",
    y = "CSS Normalized Relative Abundance",
    x = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95", color = "black"),
    strip.text = element_text(face = "bold.italic", size = 11), # 基因名加粗斜体
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.title.y = element_text(face = "bold", size = 12),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "none"
  )

print(p_hubs)

# 导出高质量图和表
ggsave("Figure_S10_Hubs_Ordered.pdf", p_hubs, width = 12, height = 10)
write.csv(all_stats, "Table_S6_Hub_Genes_Stats_Final.csv", row.names = FALSE)

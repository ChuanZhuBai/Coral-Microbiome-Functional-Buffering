# 加载必要包
library(vegan)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(rstatix)
library(patchwork)

# ==========================================
# 1. 定义稳定性计算与绘图函数
# ==========================================
analyze_stability <- function(data_path, metadata, title_prefix, y_label) {
  
  # 读取数据
  df <- read.csv(data_path, header = TRUE, row.names = 1, check.names = FALSE)
  
  # 确保样本顺序与 metadata 一致
  df <- df[, metadata$Sample]
  
  # 计算 Bray-Curtis 距离矩阵
  dist_mat <- as.matrix(vegdist(t(df), method = "bray"))
  
  # 将矩阵转为长表
  dist_long <- as.data.frame(dist_mat) %>%
    rownames_to_column("Sample1") %>%
    pivot_longer(-Sample1, names_to = "Sample2", values_to = "Distance") %>%
    left_join(metadata, by = c("Sample1" = "Sample")) %>%
    left_join(metadata, by = c("Sample2" = "Sample"), suffix = c(".1", ".2")) %>%
    filter(Sample1 != Sample2) # 去除自身对比
  
  # 筛选同一生境内的 Control vs Heat 对比
  # 逻辑：S1来自Control组，S2来自Heat组，且两者同属一个Habitat
  stability_data <- dist_long %>%
    filter(Habitat.1 == Habitat.2) %>%
    filter(Treatment.1 == "Control" & Treatment.2 == "Heat") %>%
    mutate(Habitat = factor(Habitat.1, levels = c("Intertidal", "Subtidal")))
  
  # 统计检验
  stat_test <- stability_data %>% 
    wilcox_test(Distance ~ Habitat) %>%
    add_significance() %>%
    add_y_position()
  
  # 效应量计算
  eff_size <- stability_data %>% wilcox_effsize(Distance ~ Habitat)
  
  # 绘图
  p <- ggplot(stability_data, aes(x = Habitat, y = Distance)) +
    geom_violin(aes(fill = Habitat), alpha = 0.2, color = NA) +
    geom_boxplot(aes(fill = Habitat), width = 0.2, outlier.shape = NA, alpha = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")) +
    stat_pvalue_manual(stat_test, label = "p = {p}{p.signif}", tip.length = 0.01) +
    labs(
      title = paste(title_prefix, "Stability"),
      subtitle = paste0("Effect size (r) = ", round(eff_size$effsize, 3)),
      x = "",
      y = y_label
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "none"
    )
  
  return(list(plot = p, data = stability_data, stats = stat_test))
}

# ==========================================
# 2. 准备元数据 (添加 Treatment 和 Habitat 区分)
# ==========================================
metadata <- read.csv("metadata.csv")
# 假设你的 Group 是 IT-C, IT-H 等，我们需要拆分它们以便函数识别
metadata <- metadata %>%
  mutate(
    Habitat = ifelse(grepl("IT", Group), "Intertidal", "Subtidal"),
    Treatment = ifelse(grepl("-C", Group), "Control", "Heat")
  )

# ==========================================
# 3. 运行分析
# ==========================================

# A. 属水平分类稳定性 (Genus)
res_genus <- analyze_stability(
  "Experimental_Genus_RelAbund_Final.csv", 
  metadata, 
  "Taxonomic", 
  "Bray-Curtis Dissimilarity\n(Control vs. Heat)"
)

# B. KO水平功能稳定性 (KO)
res_ko <- analyze_stability(
  "Experimental_KO_CSS_Normalized_Final.csv", 
  metadata, 
  "Functional", 
  "Bray-Curtis Dissimilarity\n(Control vs. Heat)"
)

# ==========================================
# 4. 拼图展示并保存
# ==========================================
combined_plot <- res_genus$plot + res_ko$plot + 
  plot_annotation(tag_levels = 'A') # 自动添加 A, B 标签

print(combined_plot)

ggsave("Figure_3_Stability_Comparison.pdf", combined_plot, width = 10, height = 5)

# 导出原始计算结果供 SI Table 使用
write.csv(res_genus$data, "Stability_Values_Genus.csv", row.names = FALSE)
write.csv(res_ko$data, "Stability_Values_KO.csv", row.names = FALSE)


library(tidyverse)
library(ggpubr)
library(rstatix)

# ==========================================
# 1. 参数设置与数据读入
# ==========================================
min_rel_abund <- 1e-4  # 过滤极低丰度噪声

genus_abun <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
taxon_kegg <- read.csv("extracted_taxon_kegg_cleaned.csv") 
kegg_paths <- read.csv("KEGG_paths.csv")
metadata   <- read.csv("metadata.csv")
# 读入之前计算好的包含 Generalist_Richness 的 master 表
# 如果你手头没有这个csv，代码后面有重新计算的逻辑
master_data <- read.csv("Master_Correlation_Data.csv") 

# ==========================================
# 2. KEGG 通路清洗 (全路径)
# ==========================================
kegg_paths_all <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>% 
  mutate(paths = str_trim(paths)) %>%
  separate(paths, into = c("level1","level2","level3"), sep = "; ", fill = "right") %>%
  filter(!is.na(level3)) %>% # 剔除没有三级分类的
  select(KO, path = level3) %>%
  distinct()

# 构建 “属-通路” 全局映射
taxon_path_all <- taxon_kegg %>%
  inner_join(kegg_paths_all, by = "KO", relationship = "many-to-many") %>%
  select(genus, path) %>%
  distinct()

# ==========================================
# 3. 计算 Global FRI
# ==========================================
# A. 计算样本内相对丰度
rel_abun <- genus_abun %>%
  as.matrix() %>% { sweep(., 2, colSums(.), "/") } %>%
  as.data.frame()

# B. 长表化并滤噪
abun_long <- rel_abun %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "RelAbund") %>%
  filter(RelAbund >= min_rel_abund)

# C. 定义加权冗余计算函数 (Hill q=2)
hill_q2 <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(0)
  p <- x / sum(x)
  1 / sum(p^2)
}

# D. 计算所有样本在所有通路上的 FRI
fri_all_paths <- abun_long %>%
  inner_join(taxon_path_all, by = "genus", relationship = "many-to-many") %>%
  group_by(Sample, path) %>%
  summarise(eff_contrib = hill_q2(RelAbund), .groups = "drop")

# E. 计算每个样本的全局平均 FRI (Global FRI)
sample_global_fri <- fri_all_paths %>%
  group_by(Sample) %>%
  summarise(Global_FRI = mean(eff_contrib), .groups = "drop") %>%
  left_join(metadata, by = "Sample")

# ==========================================
# 4. 分析一：IT vs ST 全局冗余差异 (Boxplot)
# ==========================================
stat_test <- sample_global_fri %>%
  wilcox_test(Global_FRI ~ Group) %>%
  add_significance()

p_box <- ggplot(sample_global_fri, aes(x = Group, y = Global_FRI, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  scale_fill_manual(values = c("IT" = "#E64B35FF", "ST" = "#4DBBD5FF")) +
  stat_pvalue_manual(stat_test, label = "p", y.position = max(sample_global_fri$Global_FRI)*1.05) +
  labs(title = "Global Functional Redundancy", x = "Habitat", y = "Mean Global FRI") +
  theme_bw() + theme(panel.grid = element_blank(), legend.position = "none")

print(p_box)
# ==========================================
# 5. 分析二：Richness vs. Global FRI 回归 (Figure 2G)
# ==========================================
# 合并数据
regression_df <- sample_global_fri %>%
  left_join(master_data %>% select(Sample, Generalist_Richness), by = "Sample") %>%
  na.omit()

p_reg <- ggplot(regression_df, aes(x = Generalist_Richness, y = Global_FRI)) +
  geom_smooth(method = "lm", color = "black", fill = "gray90", alpha = 0.6) +
  geom_point(aes(fill = Group), shape = 21, size = 5, stroke = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("IT" = "#E64B35FF", "ST" = "#4DBBD5FF")) +
  # 标注 R2 和 P
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), 
           label.x.npc = "left", label.y.npc = "top", size = 5) +
  labs(x = "Richness of Generalist Genera", 
       y = "Global Functional Redundancy Index",
       title = "Diversity of Generalists Drives Global Redundancy") +
  theme_bw() + theme(panel.grid = element_blank())

print(p_reg)
# ==========================================
# 6. 结果导出
# ==========================================
# 组合图展示
library(patchwork)
p_final <- p_box + p_reg + plot_layout(widths = c(1, 2.5))
print(p_final)

# 导出表格
write.csv(sample_global_fri, "Global_FRI_per_sample.csv", row.names = FALSE)
write.csv(regression_df, "Table_Fig2G_Global_Regression.csv", row.names = FALSE)

ggsave("Figure_2G_Global_Logic.pdf", p_final, width = 10, height = 5)
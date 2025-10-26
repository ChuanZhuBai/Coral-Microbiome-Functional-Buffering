#最终的版本
# ---------------------------
# FRI for thermal-tolerance pathways (robust version)
# ---------------------------
library(tidyverse)
library(rstatix)

# ========== 0) 参数 ==========
# 在样本内的最小相对丰度阈值（过滤极低丰度噪声）；可设为 0 以完全保留
min_rel_abund <- 1e-4
group_levels  <- c("IT-H","ST-H")   # 只比较热胁迫组；若要 IT vs ST，把这里改回 c("IT","ST")

# 明确与耐热相关的 KEGG 三级通路（以 Pathway 名称为准）
key_paths_strict <- c(
  # 蛋白稳态/去折叠蛋白清除
  "Chaperones and folding catalysts",  # 若你的 level3 不含此 BRITE 名，可删去或改用 Protein export/Proteasome
  "Proteasome",
  #"Protein export",
  
  # DNA 损伤应答/修复
  "Base excision repair",
  "Mismatch repair",
  #"Nucleotide excision repair",
  "Homologous recombination",
  #"DNA replication",
  
  # 能量代谢与氧化应激耦合
  "Oxidative phosphorylation",
  #"Glycolysis / Gluconeogenesis",
  "Pentose phosphate pathway",   # 供给 NADPH 的重要来源
  
  # 抗氧化与红氧稳态
  "Glutathione metabolism",
  # 如你的注释有 Peroxisome，可加上：
  # "Peroxisome",
  
  # 膜脂重塑（热敏感表型常见）
  #"Fatty acid biosynthesis",
  "Biosynthesis of unsaturated fatty acids",
  #"Glycerophospholipid metabolism",
  
  # 类异戊二烯衍生抗氧化/电子载体
  "Terpenoid backbone biosynthesis",
  "Carotenoid biosynthesis"
  #"Ubiquinone and other terpenoid-quinone biosynthesis"
  
  # 微生物环境感知/转运（可选；如需更严格可移除）
  # ,"Two-component system", "ABC transporters"
)

# ========== 1) 读入数据 ==========
genus_abun <- read.csv("genus.csv", row.names = 1, check.names = FALSE)   # 行=属，列=样本
taxon_kegg <- read.csv("extracted_taxon_kegg_cleaned.csv")                # 列示例: genus, KO
kegg_paths <- read.csv("KEGG_paths.csv")                                   # 列示例: KO, paths (可能是 "level1; level2; level3")
metadata   <- read.csv("metadata.csv")                                     # 列示例: Sample, Group

# 将 IT / ST 改为 IT-H / ST-H（如果你只想看热胁迫组）
metadata <- metadata %>%
  mutate(Group = case_when(
    Group == "IT" ~ "IT-H",
    Group == "ST" ~ "ST-H",
    TRUE          ~ Group
  ))

# 仅保留我们要比较的组
metadata <- metadata %>% filter(Group %in% group_levels)
genus_abun <- genus_abun[, colnames(genus_abun) %in% metadata$Sample, drop = FALSE]

# ========== 2) 规范化 KEGG 通路层级，并筛选热相关通路 ==========
kegg_paths_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%       # 如果一个 KO 映射多条通路，以 " | " 分隔
  mutate(paths = str_trim(paths)) %>%
  separate(paths, into = c("level1","level2","level3"), sep = "; ", fill = "right")

# 检查 key_paths 是否都能在文件中找到
unmatched <- setdiff(key_paths_strict, unique(kegg_paths_long$level3))
if (length(unmatched) > 0) {
  message("未在 KEGG_paths.csv$level3 中找到以下通路（请核对命名或以 KO ID 方式筛选）: ",
          paste(unmatched, collapse = "; "))
}

# 只保留热相关通路
kegg_target <- kegg_paths_long %>%
  filter(level3 %in% key_paths_strict) %>%
  select(KO, level3) %>%
  distinct() %>%
  rename(path = level3)

# 如果你的 KEGG_paths.csv 存的是 koID（例如 ko00190）而非 level3 名字：
# 可改为： kegg_target <- kegg_paths_long %>% filter(ko %in% c("ko00190", ...)) %>% select(KO=ko, path=ko)

# ========== 3) 构建 “属-通路” 二分映射 ==========
taxon_path <- taxon_kegg %>%
  inner_join(kegg_target, by = "KO", relationship = "many-to-many") %>%
  select(genus, path) %>%
  distinct()

# ========== 4) 相对丰度与滤噪 ==========
# 样本内相对丰度
rel_abun <- genus_abun %>%
  as.matrix() %>% { sweep(., 2, colSums(.), "/") } %>%
  replace(is.na(.), 0) %>%
  as.data.frame()

# 长表
abun_long <- rel_abun %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "RelAbund") %>%
  left_join(metadata %>% select(Sample, Group), by = "Sample") %>%
  filter(Group %in% group_levels)

# 过滤极低丰度（减少偶然性“出现”）
abun_long_filt <- abun_long %>%
  filter(RelAbund >= min_rel_abund)

# ========== 5) 计算三种 FRI ==========
# 定义一个安全的 Hill 数 q=2（有效贡献者数）：1/sum(p_i^2)，p_i 为通路内各属的相对贡献（在通路内归一化）
hill_q2 <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) == 0) return(0)
  p <- x / sum(x)
  1 / sum(p^2)
}

# 生成 Sample × path 的 FRI
fri_long <- abun_long_filt %>%
  inner_join(taxon_path, by = "genus", relationship = "many-to-many") %>%
  group_by(Sample, Group, path) %>%
  summarise(
    # 通路内“贡献属”的集合与计数
    n_genera     = n_distinct(genus),
    # 通路内各属的丰度（用于加权 FRI）
    eff_contrib  = hill_q2(RelAbund),             # 有效贡献者数（加权）
    .groups = "drop"
  )

# 样本的属数（用于归一化）
sample_richness <- abun_long_filt %>%
  group_by(Sample) %>%
  summarise(sample_genus_rich = n_distinct(genus), .groups = "drop")

fri_long <- fri_long %>%
  left_join(sample_richness, by = "Sample") %>%
  mutate(
    FRI_count = n_genera,
    FRI_norm  = ifelse(sample_genus_rich > 0, n_genera / sample_genus_rich, NA_real_),
    FRI_eff   = eff_contrib
  ) %>%
  select(Sample, Group, path, FRI_count, FRI_norm, FRI_eff)

# ========== 6) 统计检验（IT-H vs ST-H），多指标 + 多重校正 ==========
# 长表汇总
fri_stats_input <- fri_long %>%
  pivot_longer(cols = starts_with("FRI_"), names_to = "Metric", values_to = "Value")

# 按通路与指标做 Wilcoxon 与效应量（Cliff's delta），BH 调整
stat_res <- fri_stats_input %>%
  group_by(path, Metric) %>%
  wilcox_test(Value ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj") %>%
  left_join(
    fri_stats_input %>%
      group_by(path, Metric) %>%
      wilcox_effsize(Value ~ Group, alternative = "two.sided") %>%
      select(path, Metric, effsize, magnitude),
    by = c("path","Metric")
  ) %>%
  ungroup()

# 组内均值±SE（用于绘图或表格）
summary_tbl <- fri_long %>%
  pivot_longer(cols = starts_with("FRI_"), names_to = "Metric", values_to = "Value") %>%
  group_by(path, Metric, Group) %>%
  summarise(
    mean = mean(Value, na.rm = TRUE),
    se   = sd(Value,  na.rm = TRUE) / sqrt(sum(is.finite(Value))),
    n    = sum(is.finite(Value)),
    .groups = "drop"
  )

# ========== 7) 导出结果 ==========
write.csv(fri_long,     "FRI_per_sample.csv", row.names = FALSE)
write.csv(summary_tbl,  "FRI_summary_by_group.csv", row.names = FALSE)
write.csv(stat_res,     "FRI_stats_wilcox_BH.csv", row.names = FALSE)

# ========== 8) 可选绘图（示例：加权 FRI 的柱图）==========
# —— 选择要画的指标： "FRI_eff" / "FRI_count" / "FRI_norm"
metric_to_plot <- "FRI_eff"

plot_df <- summary_tbl %>% filter(Metric == metric_to_plot)

# 为当前指标计算 y 位置
y_max <- max(plot_df$mean + plot_df$se, na.rm = TRUE) * 1.08

ggplot(plot_df, aes(x = path, y = mean, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.8), width = 0.2) +
  # 显著性标注要与当前 metric 对齐
  geom_text(
    data = stat_res %>% filter(Metric == metric_to_plot),
    aes(x = path, y = y_max, label = p.adj.signif),
    inherit.aes = FALSE, size = 4.5, fontface = "bold"
  ) +
  labs(x = NULL,
       y = case_when(
         metric_to_plot == "FRI_eff"  ~ "FRI (effective contributors, q=2)",
         metric_to_plot == "FRI_norm" ~ "FRI (normalized count)",
         TRUE                         ~ "FRI (count)"
       )) +
  scale_fill_manual(values = c("IT-H" = "#E64B35", "ST-H" = "#4DBBD5")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))











#######################################
######################################
library(tidyverse)
library(ggpubr)

# 1. 数据加载
genus_abun <- read.csv("genus.csv", row.names = 1)  # 菌属丰度表
taxon_kegg <- read.csv("extracted_taxon_kegg_cleaned.csv")  # 菌属-KO关联
kegg_paths <- read.csv("KEGG_paths.csv")  # KO-通路映射
metadata <- read.csv("metadata.csv")  # 样本分组信息

# 如果 metadata$Group 只有 "IT" 和 "ST"，改为 "IT-H" / "ST-H"
metadata <- metadata %>%
  mutate(Group = case_when(
    Group == "IT" ~ "IT-H",
    Group == "ST" ~ "ST-H",
    TRUE          ~ Group
  ))
# 2. 定义关键抗热通路
key_paths <- c(
  # 蛋白稳态
  "Chaperones and folding catalysts", "Proteasome", "Protein export",
  # DNA修复（具体）
  "Base excision repair", "Mismatch repair", "Nucleotide excision repair", "Homologous recombination", "DNA replication",
  # 能量与抗氧化
  "Oxidative phosphorylation", "Glycolysis / Gluconeogenesis", "Glutathione metabolism",
  # 膜脂重塑
  "Fatty acid biosynthesis", "Biosynthesis of unsaturated fatty acids", "Glycerophospholipid metabolism",
  # 异戊二烯/类胡萝卜素/醌
  "Terpenoid backbone biosynthesis", "Carotenoid biosynthesis", "Ubiquinone and other terpenoid-quinone biosynthesis",
  # 微生物感知与转运（可选）
  "Two-component system", "ABC transporters"
)


# 3. 处理KEGG通路数据（处理多级分类）
kegg_paths_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(paths = str_trim(paths)) %>%
  separate(paths, into = c("level1", "level2", "level3"), sep = "; ", fill = "right") %>%
  filter(level3 %in% key_paths) %>%
  select(KO, level3) %>%
  distinct()

# 4. 构建菌属-通路对应表（明确处理多对多关系）
taxon_path <- taxon_kegg %>%
  inner_join(kegg_paths_long, by = "KO", relationship = "many-to-many") %>%
  select(genus, paths = level3) %>%
  distinct()

# 5. 计算每个样本每条通路的FRI（菌属数量）
fri_raw <- genus_abun %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%  # 只计入丰度>0的菌属
  inner_join(taxon_path, by = "genus", relationship = "many-to-many") %>%
  group_by(Sample, paths) %>%
  summarise(FRI = n_distinct(genus), .groups = "drop")

# 6. 合并分组信息
fri_df <- fri_raw %>%
  left_join(metadata %>% select(Sample, Group), by = "Sample") %>%
  filter(Group %in% c("IT-H", "ST-H"))  # 只保留热胁迫组

# 7. 计算均值和标准误
plot_df <- fri_df %>%
  group_by(paths, Group) %>%
  summarise(
    mean_FRI = mean(FRI),
    se_FRI = sd(FRI)/sqrt(n()),
    .groups = "drop"
  )


write.csv(plot_df, "FRI.csv")
library(rstatix)  # 提供add_xy_position()函数
# 8. 统计检验（处理警告）
# 8. 统计检验（修正版）
# 统计检验
stat_test <- fri_df %>%
  group_by(paths) %>%
  wilcox_test(FRI ~ Group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

# 计算标记位置
max_height <- max(plot_df$mean_FRI + plot_df$se_FRI)
y_position <- max_height * 1.1

# 绘图
ggplot(plot_df, aes(x = paths, y = mean_FRI, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_FRI - se_FRI, ymax = mean_FRI + se_FRI),
    position = position_dodge(width = 0.8), width = 0.2
  ) +
  geom_text(
    data = stat_test,
    aes(x = paths, y = y_position, label = p.adj.signif),
    inherit.aes = FALSE,
    vjust = -0.5,
    size = 5,
    fontface = "bold"  
  ) +
  scale_fill_manual(values = c("IT-H" = "#E64B35", "ST-H" = "#4DBBD5")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +  # 为标记留出空间
  labs(x = NULL, y = "FRI") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))









# ───────────────────────────────────────────────────────────────────────────────
# Panel B: 抗热通路功能冗余指数（FRI）— IT-H vs ST-H 分组柱状图
# ───────────────────────────────────────────────────────────────────────────────

# 0. 加载必要的包
library(tidyverse)   # 数据整形 + ggplot2
library(ggpubr)      # stat_compare_means
dir()
# 1. 读取数据
genus_abun   <- read.csv("genus.csv",    row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
taxon_kegg   <- read.csv("extracted_taxon_kegg_cleaned.csv", stringsAsFactors = FALSE)
kegg_paths   <- read.csv("KEGG_paths.csv",        stringsAsFactors = FALSE)
metadata     <- read.csv("metadata.csv",          stringsAsFactors = FALSE)



# 如果 metadata$Group 只有 "IT" 和 "ST"，改为 "IT-H" / "ST-H"
metadata <- metadata %>%
  mutate(Group = case_when(
    Group == "IT" ~ "IT-H",
    Group == "ST" ~ "ST-H",
    TRUE          ~ Group
  ))

# 2. 定义关键抗热通路
key_paths <- c(
  "Chaperones and folding catalysts",
  "Glutathione metabolism",
  "Base excision repair",
  "DNA repair and recombination proteins"
)


# 3. 拆分 paths → 提取 level3 → 筛选关键抗热通路
library(tidyr)  # for separate_rows()

kegg_paths_long <- kegg_paths %>%
  # 先按 " | " 拆成多行
  separate_rows(paths, sep = " \\| ") %>%
  # 去掉前后空白
  mutate(paths = str_trim(paths)) %>%
  # 按 "; " 拆成 3 列：level1、level2、level3
  separate(paths, into = c("level1", "level2", "level3"), sep = "; ") %>%
  # 只保留我们关心的第三级通路
  filter(level3 %in% key_paths) %>%
  select(KO, level3)

# 4. 构建 “菌属 → 通路(level3)” 对应表
taxon_path <- taxon_kegg %>%
  inner_join(kegg_paths_long, by = "KO") %>%
  select(genus, paths = level3) %>%
  distinct()


# 5. 计算每个样本每条通路的 FRI
fri_raw <- genus_abun %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%                            # 只计入丰度>0的菌属
  inner_join(taxon_path, by = "genus") %>%            # 关联到通路
  group_by(Sample, paths) %>%
  summarise(FRI = n_distinct(genus), .groups = "drop")# FRI = 菌属数

# 6. 合并分组信息
fri_df <- fri_raw %>%
  left_join(metadata %>% select(Sample, Group), by = "Sample")

# 7. 统计每条通路 IT-H vs ST-H 的 Wilcoxon 检验
stat_df <- fri_df %>%
  group_by(paths) %>%
  summarise(
    p.value = wilcox.test(FRI ~ Group)$p.value,
    .groups = "drop"
  )

# 8. 计算均值与标准误，用于柱状图
plot_df <- fri_df %>%
  group_by(paths, Group) %>%
  summarise(
    mean_FRI = mean(FRI),
    se_FRI   = sd(FRI) / sqrt(n()),
    .groups = "drop"
  )
head(plot_df)
# 9. 绘图
ggplot(plot_df, aes(x = paths, y = mean_FRI, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_FRI - se_FRI, ymax = mean_FRI + se_FRI),
                position = position_dodge(width = 0.8), width = 0.2) +
  # 添加 Wilcoxon 显著性标注
  stat_compare_means(
    data        = fri_df,
    aes(x = paths, y = FRI),
    method      = "wilcox.test",
    comparisons = list(c("IT-H", "ST-H")),
    label       = "p.signif",
    label.y     = max(plot_df$mean_FRI + plot_df$se_FRI) * 1.1
  ) +
  scale_x_discrete(limits = key_paths) +
  labs(
    x    = NULL,
    y    = "Functional Redundancy Index (FRI)",
    fill = "Group"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    text        = element_text(size = 12)
  )

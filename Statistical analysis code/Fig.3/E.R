# —— 0. 载入包 —— 
library(tidyverse)
library(readr)

# —— 1. 定义 Stress Response 通路与提取对应 KO —— 
heat_paths <- c(
 "Signal transduction",
  "Replication and repair","Metabolism of terpenoids and polyketides","Carbohydrate metabolism", "Energy metabolism",
 "Lipid metabolism",       "Nucleotide metabolism",
 "Amino acid metabolism", "Metabolism of other amino acids",
 "Glycan biosynthesis and metabolism",
 "Metabolism of cofactors and vitamins",
 "Biosynthesis of other secondary metabolites",
 "Xenobiotics biodegradation and metabolism"
 
  
  
)

kegg_paths <- read_csv("KEGG_paths.csv", show_col_types = FALSE) %>%
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  mutate(Level2 = str_split(paths, ";\\s*", simplify = TRUE)[,2])

stress_kos <- kegg_paths %>%
  filter(Level2 %in% heat_paths) %>%
  distinct(KO) %>% pull(KO)

# —— 2. 读取丰度矩阵（行=KO，列=所有样本） —— 
abun <- read_csv("CSS_relative_KEGG_bases_lab.csv", show_col_types = FALSE) %>%
  column_to_rownames("KO")

# —— 3. 计算每个 KO 的“专属伪计数” —— 
#     取该 KO 在所有样本非零值的最小值的一半
pc <- abun[stress_kos, ] %>%
  apply(1, function(x){
    nz <- x[x > 0]
    if(length(nz)>0) min(nz)/2 else 1e-5
  }) %>% set_names(stress_kos)

# —— 4. 自动配对：加热样本 vs 对照样本 —— 
all_samps <- colnames(abun)
heated <- all_samps[str_detect(all_samps, "-H\\d+$")]
pairs <- tibble(heat = heated) %>%
  mutate(ctrl = sub("-H", "-C", heat)) %>%
  filter(ctrl %in% all_samps)

# —— 5. 向量化计算每对样本的 FRR —— 
fr_table <- pairs %>%
  mutate(
    FRR_val = map2_dbl(heat, ctrl, ~{
      h <- abun[stress_kos, .x] %>% unlist()
      c <- abun[stress_kos, .y] %>% unlist()
      log2fc <- log2((h + pc) / (c + pc))
      mean(abs(log2fc) <= 1, na.rm = TRUE) * 100
    }),
    FRR = sprintf("%.2f%%", FRR_val)
  )

# 用 base R 提取并重命名
fr_table <- data.frame(
  Sample = fr_table$heat,
  FRR    = fr_table$FRR,
  stringsAsFactors = FALSE
)
print(fr_table)

# —— 6. 分别提取 IT-H 和 ST-H 系列，并保存 —— 
fr_it_h <- fr_table %>% filter(str_starts(Sample, "IT-H"))
fr_st_h <- fr_table %>% filter(str_starts(Sample, "ST-H"))

# 输出到控制台
print(fr_it_h)
print(fr_st_h)

# 保存为 CSV
write_csv(fr_it_h, "IT-H_per_pair_StressResponse_FRR.csv")
write_csv(fr_st_h, "ST-H_per_pair_StressResponse_FRR.csv")


######################################################
#表型相关
# —— 加载必要的包 —— 
library(tidyverse)
library(ggpubr)
library(ggsci)
library(scales)

# —— 第一步：读入 FRR 和表型数据 —— 

# 1. Stress‐Response 通路保留率（FRR）
it_frr <- read_csv("IT-H_per_pair_StressResponse_FRR.csv", show_col_types = FALSE) %>%
  rename(FRR = FRR) %>%
  mutate(FRR_val = parse_number(FRR))

st_frr <- read_csv("ST-H_per_pair_StressResponse_FRR.csv", show_col_types = FALSE) %>%
  rename(FRR = FRR) %>%
  mutate(FRR_val = parse_number(FRR))

frr_data <- bind_rows(it_frr, st_frr)

# 2. 珊瑚表型数据
paired <- read_csv("paired_data.csv", show_col_types = FALSE) 

# 3. 合并 FRR 与表型
combined_data <- paired %>%
  left_join(frr_data, by = "Sample") %>%
  # 重新标记分面用的 Group
  mutate(
    Group = case_when(
      Group == "IT-H" ~ "Intertidal (IT)",
      Group == "ST-H" ~ "Subtidal (ST)"
    ),
    Group = factor(Group, levels = c("Intertidal (IT)", "Subtidal (ST)"))
  )

# —— 第二步：统计分析 —— 

group_stats <- combined_data %>%
  group_by(Group) %>%
  summarise(
    # 线性回归：ChangeRate ~ FRR_val
    lm_r2 = summary(lm(ChangeRate ~ FRR_val, data = .))$r.squared,
    lm_p  = summary(lm(ChangeRate ~ FRR_val, data = .))$coefficients[2,4],
    # Spearman 相关
    rho   = cor.test(FRR_val, ChangeRate, method = "spearman")$estimate,
    rho_p = cor.test(FRR_val, ChangeRate, method = "spearman")$p.value,
    .groups = "drop"
  )

# —— 第三步：绘图 —— 

base_plot <- ggplot(combined_data,
                    aes(x = FRR_val, y = ChangeRate,
                        color = Group, fill = Group)) +
  # 散点
  geom_point(shape = 21, size = 3.5, stroke = 0.8, alpha = 0.9, show.legend = FALSE) +
  # 回归线
  geom_smooth(method = "lm", formula = y ~ x,
              se = TRUE, level = 0.95,
              linewidth = 1.2, alpha = 0.2) +
  # 颜色
  scale_color_manual(values = c("#E64B35", "#4DBBD5")) +
  scale_fill_manual(values  = c("#E64B35", "#4DBBD5")) +
  # 坐标轴
  scale_x_continuous(name = "Stress Response FRR (%)",
                     limits = c(0, 100),
                     breaks = seq(0, 100, by = 20)) +
  scale_y_continuous(name = "ΔFv/Fm Change Rate (%)",
                     breaks = pretty_breaks(n = 6)) +
  # 主题
  theme_pubr(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    axis.title      = element_text(face = "bold"),
    legend.position = "top",
    legend.title    = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.background = element_rect(fill = "white")
  )

final_plot <- base_plot +
  # 添加 Pearson R² 和 p
  stat_cor(
    aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),
    method = "pearson",
    label.x = 20, label.y = c(-10, -20),
    size = 3.5, show.legend = FALSE
  ) +
  # 添加 Spearman ρ 和 p
  geom_text(
    data = group_stats,
    aes(x = 20, y = c(-30, -40),
        label = sprintf("Spearman's ρ = %.2f\np = %.3f", rho, rho_p)),
    size = 3.5, hjust = 0, color = "black"
  ) +
  # 分面
  facet_wrap(~Group, ncol = 2) +
  labs(
    title    = "Stress Response FRR Predicts Coral Thermal Tolerance",
    subtitle = "Intertidal vs Subtidal Habitats"
  )

# 打印
print(final_plot)

# —— 第四步：保存 —— 
ggsave("FRR_vs_ChangeRate.tiff", plot = final_plot,
       device = "tiff",
       width = 18, height = 10, units = "cm",
       dpi = 600, compression = "lzw")

ggsave("FRR_vs_ChangeRate.pdf", plot = final_plot,
       device = "pdf",
       width = 18, height = 10, units = "cm")


##################################################













#最终代码-生成代谢通路保留率
# ===================================================
# 核心/热响应 Level-2 通路纯保留率分析 & 可视化
# ===================================================

# 0. 预设：自己关注的 Level-2 通路列表（基础+热响应）
basic_L2_keep <- c(
  "Carbohydrate metabolism","Energy metabolism","Lipid metabolism",
  "Nucleotide metabolism","Amino acid metabolism","Metabolism of other amino acids",
  "Glycan biosynthesis and metabolism","Metabolism of cofactors and vitamins",
  "Metabolism of terpenoids and polyketides","Xenobiotics biodegradation and metabolism"
)

stress_L2_keep <- c(
  "Replication and repair","Folding, sorting and degradation","Signal transduction"
)

target_L2 <- c(basic_L2_keep, stress_L2_keep)


# 1. 加载所需 R 包
library(tidyverse)
library(readr)
library(ggpubr)
library(lme4)
library(broom.mixed)
library(ggrepel)

# 2. 读取并合并差异结果
it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C.csv") %>% mutate(Group="IT-H")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C.csv") %>% mutate(Group="ST-H")
kegg_paths <- read_csv("KEGG_paths.csv") %>% 
  mutate(KO = str_extract(KO, "K\\d{5}"))

# 3. 拆行提取 Level2，打上分类标签，并筛选感兴趣通路
combined <- bind_rows(it_h, st_h) %>%
  left_join(kegg_paths, by = "KO") %>%
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  mutate(Level2 = str_split(paths, ";\\s*", simplify = TRUE)[,2]) %>%
  filter(!is.na(Level2), Level2 %in% target_L2) %>%
  mutate(
    PathType = if_else(Level2 %in% basic_L2_keep, "Core Metabolic", "Stress Response"),
    Unchanged = as.integer(is.na(PValue) | is.na(logFC) | !(PValue < 0.01 & abs(logFC) > 1)),
    Increased = as.integer(PValue < 0.01 & logFC >  1),
    Decreased = as.integer(PValue < 0.01 & logFC < -1)
  )


combined %>% 
  filter(Unchanged == 1) %>%
  count(Level2, KO) %>%
  filter(n > 1)

# 4. 计算纯保留率等指标
pathway_stats <- combined %>%
  group_by(Group, PathType, Level2) %>%
  summarise(
    ## 总 KO 数：对 KO 去重
    Total_KO     = n_distinct(KO),
    ## Unchanged_KO：对 Unchanged==1 的 KO 去重
    Unchanged_KO = n_distinct(KO[Unchanged == 1]),
    Increased_KO = n_distinct(KO[Increased == 1]),
    Decreased_KO = n_distinct(KO[Decreased == 1]),
    .groups      = "drop"
  ) %>%
  filter(Total_KO >= 5) %>%
  mutate(
    FRR      = 100 * Unchanged_KO / Total_KO,
    ActRate  = 100 * Increased_KO / Total_KO,
    LossRate = 100 * Decreased_KO / Total_KO
  )

glimpse(pathway_stats)

# 5. 配对比较 IT-H vs ST-H 并排序
retention_comp <- pathway_stats %>%
  dplyr::select(Group, Level2, FRR) %>%         # 明确调用 dplyr::select
  tidyr::pivot_wider(
    names_from  = Group,
    values_from = FRR
  ) %>%
  mutate(
    FRR_diff    = `IT-H` - `ST-H`,
    Delta_Label = ifelse(
      FRR_diff > 0,
      paste0("+", round(FRR_diff, 1), "%"),
      paste0(round(FRR_diff, 1), "%")
    )
  ) %>%
  filter(!is.na(FRR_diff))

ordered_paths <- retention_comp %>% arrange(FRR_diff) %>% pull(Level2)

# 6. 统计检验
wilcox_res <- wilcox.test(retention_comp$`IT-H`, retention_comp$`ST-H`,
                          paired = TRUE, exact = FALSE)
message("Wilcoxon p-value: ", signif(wilcox_res$p.value,3))

model_data <- retention_comp %>%
  pivot_longer(cols=c(`IT-H`,`ST-H`), names_to="Group", values_to="FRR")

lmm <- lmer(FRR ~ Group + (1|Level2), data = model_data)
lmm_tidy <- tidy(lmm)
print(lmm_tidy)

# 7. 可视化
# 7.1 ΔFRR 柱状图 + 按 PathType 分色
p_main <- retention_comp %>%
  left_join(
    dplyr::select(pathway_stats, Level2, PathType) %>% distinct(),
    by = "Level2"
  ) %>%
  mutate(Level2 = factor(Level2, levels = ordered_paths)) %>%
  ggplot(aes(x = Level2, y = FRR_diff, fill = PathType)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = Delta_Label),
            hjust = ifelse(retention_comp$FRR_diff > 0, -0.1, 1.1),
            size = 3) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Core Metabolic"    = "#4DBBD5",
      "Stress Response"   = "#E64B35",
      "Other"             = "grey80"
    )
  ) +
  labs(
    x = "KEGG Level-2 Pathway", y = "ΔFRR (IT-H − ST-H, %)",
    fill = "Pathway Type"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.major.y = element_blank())

print(p_main)


# 7.2 热图 & 箱线图 同前，可视化仅关注 target_paths 的结果
## 7.2 热图：IT-H vs ST-H FRR
p_heatmap <- pathway_stats %>%
  filter(Group %in% c("IT-H","ST-H")) %>%
  mutate(Level2=factor(Level2, levels=ordered_paths)) %>%
  ggplot(aes(x=Group, y=Level2, fill=FRR)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=round(FRR,1)), size=3) +
  scale_fill_gradient2(low="#3288BD", mid="#FFFFBF", high="#D53E4F",
                       midpoint=50, limits=c(0,100)) +
  labs(x=NULL,y=NULL,fill="Pure retention rate (%)",
       title="IT-H vs ST-H Pure Retention Rate Heatmap") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45,hjust=1))

## 7.3 激活率 & 丢失率箱线图
p_rates <- pathway_stats %>%
  pivot_longer(cols=c(ActRate, LossRate),
               names_to="Metric", values_to="Value") %>%
  mutate(Metric = recode(Metric, ActRate="Activation", LossRate="Loss")) %>%
  ggplot(aes(x=Metric, y=Value, fill=Group)) +
  geom_boxplot() +
  labs(x=NULL,y="%", title="Activation rate vs Loss rate") +
  theme_minimal()
print(p_rates)
## 合并
final_plot <- (p_main / (p_heatmap | p_rates)) + 
  plot_layout(heights=c(2,2))

# 8. 保存
write_csv(pathway_stats, "core_pathway_pure_FRR_rates-1.csv")
write_csv(retention_comp, "core_pathway_FRR_diff-1.csv")
ggsave("core_pathway_retention_overview-1.pdf", final_plot,
       width=12, height=10, dpi=300)

print(final_plot)

# ===================================================
# 核心/热响应 Level-2 & 指定 Level-3 通路功能保留率分析 & 可视化（终极优化）
# ===================================================

# 0. 通路列表 ---------------------------------------------------------------
basic_paths <- c(
  "Carbohydrate metabolism", "Energy metabolism",
  "Lipid metabolism", "Nucleotide metabolism",
  "Amino acid metabolism", "Metabolism of other amino acids",
  "Glycan biosynthesis and metabolism",
  "Metabolism of cofactors and vitamins",
  "Metabolism of terpenoids and polyketides",
  "Biosynthesis of other secondary metabolites",
  "Xenobiotics biodegradation and metabolism",
  "Chemical structure transformation maps"
)
heat_paths2 <- c(
  "Glutathione metabolism", "Oxidative phosphorylation",
  "Carotenoid biosynthesis", "Protein processing in endoplasmic reticulum",
  "Proteasome", "Ubiquitin mediated proteolysis",
  "Base excision repair", "Nucleotide excision repair",
  "MAPK signaling pathway", "Two-component system",
  "Glutathione metabolism"
)
key_l3 <- c(
  "Base excision repair",
  "Chaperones and folding catalysts",
  "DNA repair and recombination proteins",
  "Glutathione metabolism"
)
heat_paths <- unique(c(heat_paths2, key_l3))
target_paths <- unique(c(basic_paths, heat_paths))

# 1. 加载包 ----------------------------------------------------------------
library(tidyverse)
library(readr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(patchwork)
library(ggrepel)

# 2. 数据读取 --------------------------------------------------------------
it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C.csv") %>% mutate(Group = "IT-H")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C.csv") %>% mutate(Group = "ST-H")
kegg_paths <- read_csv("KEGG_paths.csv")

# 3. KO 格式化与拆分 --------------------------------------------------------
extract_ko <- function(x) str_extract(x, "K\\d{5}")
kegg_paths$KO <- extract_ko(kegg_paths$KO)
it_h$KO        <- extract_ko(it_h$KO)
st_h$KO        <- extract_ko(st_h$KO)

combined <- bind_rows(it_h, st_h) %>%
  left_join(kegg_paths, by = "KO") %>%
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  mutate(
    parts  = str_split(paths, ";\\s*", simplify = TRUE),
    Level2 = trimws(parts[,2]),
    Level3 = trimws(parts[,3]),
    Path   = if_else(Level3 %in% key_l3, Level3, Level2)
  ) %>%
  filter(!is.na(Path) & Path %in% target_paths) %>%
  mutate(
    PathType  = case_when(
      Path %in% basic_paths ~ "Core Metabolic",
      Path %in% heat_paths   ~ "Stress Response",
      TRUE                   ~ "Other"
    ),
    Unchanged = as.integer(
      is.na(PValue) | is.na(logFC)
      | !(PValue < 0.01 & abs(logFC) > 1)
    ),
    Increased = as.integer(PValue < 0.01 & logFC >  1),
    Decreased = as.integer(PValue < 0.01 & logFC < -1)
  )

# 4. FRR 等指标 ------------------------------------------------------------
pathway_stats <- combined %>%
  group_by(Group, PathType, Path) %>%
  summarise(
    Total_KO     = n_distinct(KO),
    Unchanged_KO = n_distinct(KO[Unchanged == 1]),
    Increased_KO = n_distinct(KO[Increased == 1]),
    Decreased_KO = n_distinct(KO[Decreased == 1]),
    .groups       = "drop"
  ) %>%
  filter(Total_KO >= 5) %>%
  mutate(
    FRR      = 100 * Unchanged_KO / Total_KO,
    ActRate  = 100 * Increased_KO / Total_KO,
    LossRate = 100 * Decreased_KO / Total_KO
  )

# 5. 排序 & 全局检验 --------------------------------------------------------
retention_comp <- pathway_stats %>%
  dplyr::select(Group, Path, FRR) %>%
  pivot_wider(names_from = Group, values_from = FRR) %>%
  mutate(
    FRR_diff    = `IT-H` - `ST-H`,
    Delta_Label = sprintf("%+.1f%%", FRR_diff)
  ) %>%
  filter(!is.na(FRR_diff))
ordered_paths <- retention_comp %>% arrange(FRR_diff) %>% pull(Path)

wilcox_res <- wilcox.test(retention_comp$`IT-H`, retention_comp$`ST-H`,
                          paired = TRUE, exact = FALSE)
message("Global Wilcoxon p-value: ", signif(wilcox_res$p.value,3))
model_data <- retention_comp %>%
  pivot_longer(cols = c(`IT-H`,`ST-H`), names_to = "Group", values_to = "FRR")

# 6. 单通路检验 & FDR 校正 --------------------------------------------------
compare_by_pathway <- model_data %>%
  group_by(Path) %>%
  summarise(
    p_raw = wilcox.test(FRR[Group=="IT-H"], FRR[Group=="ST-H"], exact=FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj     = p.adjust(p_raw, method="fdr"),
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE          ~ ""
    )
  )

# 7. 合并星号 & 计算位置 ----------------------------------------------------
retention_plot_data2 <- model_data %>%
  left_join(compare_by_pathway %>% dplyr::select(Path, sig_label), by = "Path") %>%
  group_by(Path) %>%
  mutate(
    max_rate   = max(FRR, na.rm = TRUE),
    y_position = max_rate * 1.1
  ) %>%
  ungroup()

# 8. 可视化 ---------------------------------------------------------------

## 8.1 ΔFRR 柱状图 + 星号
p_main <- retention_plot_data2 %>%
  mutate(Path = factor(Path, levels = ordered_paths)) %>%
  ggplot() +
  # 主柱状
  geom_col(aes(x = Path, y = FRR, fill = Group),
           position = position_dodge(0.8), width = 0.7) +
  # ΔFRR 文本
  geom_text(data = retention_comp %>% mutate(Path = factor(Path, levels = ordered_paths)),
            aes(x = Path, y = pmax(`IT-H`, `ST-H`) * 1.05, label = Delta_Label),
            inherit.aes = FALSE, size = 4) +
  # 显著性星号
  geom_text(data = retention_plot_data2 %>% filter(Group == "ST-H"),
            aes(x = Path, y = y_position, label = sig_label),
            inherit.aes = FALSE,
            position = position_dodge(0.8), vjust = 0, size = 6, fontface = "bold") +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x        = NULL,
    y        = "Retention Rate (%)",
    title    = "Functional Retention Comparison per Pathway",
    subtitle = paste0("Global Wilcoxon p = ", signif(wilcox_res$p.value,3))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y         = element_text(size = 12),
    legend.position     = "top"
  )


print(p_main)
## 8.2 Heatmap
p_heatmap <- pathway_stats %>%
  filter(Group %in% c("IT-H","ST-H")) %>%
  mutate(Path = factor(Path, levels = ordered_paths)) %>%
  ggplot(aes(x = Group, y = Path, fill = FRR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(FRR,1)), size = 3) +
  scale_fill_gradient2(low="#3288BD",mid="#FFFFBF",high="#D53E4F",
                       midpoint = 50, limits = c(0,100)) +
  labs(x = NULL, y = NULL, fill = "FRR (%)") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 9.3 Activation vs Loss Rate
p_rates <- pathway_stats %>%
  pivot_longer(cols=c(ActRate, LossRate), names_to="Metric", values_to="Val") %>%
  ggplot(aes(x = Metric, y = Val, fill = Group)) +
  geom_boxplot(width=0.6, outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(0.2), size=1, alpha=0.6) +
  scale_fill_brewer(palette="Set1") +
  labs(x=NULL, y="%", title="Activation vs Loss Rate") +
  theme_minimal(base_size=14)

# 10. 合并与导出 ----------------------------------------------------------
final_plot <- (p_main / (p_heatmap | p_rates)) + plot_layout(heights=c(3,2))

write_csv(pathway_stats, "pathway_stats.csv")
write_csv(compare_by_pathway, "pathway_pvalues.csv")

ggsave("retention_overview.pdf", final_plot, width=14, height=12, dpi=300)
print(final_plot)




# ===================================================
# 核心/热响应 Level-2 & 指定 Level-3 通路功能保留率分析 & 可视化
# ===================================================

# 0. 预设关注的通路列表
basic_paths <- c(
  "Carbohydrate metabolism", "Energy metabolism",
  "Lipid metabolism",       "Nucleotide metabolism",
  "Amino acid metabolism", "Metabolism of other amino acids",
  "Glycan biosynthesis and metabolism",
  "Metabolism of cofactors and vitamins",
  "Metabolism of terpenoids and polyketides",
  "Biosynthesis of other secondary metabolites",
  "Xenobiotics biodegradation and metabolism",
  "Chemical structure transformation maps"
)
heat_paths2 <- c(
    "Oxidative phosphorylation",
  "Carotenoid biosynthesis",
  "Protein processing in endoplasmic reticulum",
  "Proteasome",
  "Ubiquitin mediated proteolysis",
  "Base excision repair",
  "Nucleotide excision repair",
  "MAPK signaling pathway",
  "Two-component system"
)
# 指定关注的 Level-3 通路
key_l3 <- c(
  "Base excision repair",
  "Chaperones and folding catalysts",
  "DNA repair and recombination proteins",
  "Glutathione metabolism"
)
# 合并为最终热响应列表
heat_paths <- unique(c(heat_paths2, key_l3))
# 全部目标通路
target_paths <- unique(c(basic_paths, heat_paths))

# 1. 加载所需 R 包
library(tidyverse)
library(readr)
library(ggpubr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(patchwork)
library(ggrepel)

# 2. 读取数据并做质量检查
it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C.csv") %>% mutate(Group = "IT-H")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C.csv") %>% mutate(Group = "ST-H")
kegg_paths <- read_csv("KEGG_paths.csv")

# KO 格式化
extract_ko <- function(x) str_extract(x, "K\\d{5}")
kegg_paths$KO <- extract_ko(kegg_paths$KO)
it_h$KO        <- extract_ko(it_h$KO)
st_h$KO        <- extract_ko(st_h$KO)

# 3. 拆行提取 Level2/Level3，并生成 Path 字段
combined <- bind_rows(it_h, st_h) %>%
  left_join(kegg_paths, by = "KO") %>%
  separate_rows(paths, sep = "\\s*\\|\\s*") %>%
  mutate(
    parts   = str_split(paths, ";\\s*", simplify = TRUE),
    Level2  = trimws(parts[,2]),
    Level3  = trimws(parts[,3]),
    Path    = if_else(Level3 %in% key_l3, Level3, Level2)
  ) %>%
  filter(!is.na(Path) & Path %in% target_paths) %>%
  mutate(
    PathType  = case_when(
      Path %in% basic_paths ~ "Core Metabolic",
      Path %in% heat_paths   ~ "Stress Response",
      TRUE                   ~ "Other"
    ),
    Unchanged = as.integer(
      is.na(PValue) | is.na(logFC)
      | !(PValue < 0.05 & abs(logFC) > 1)
    ),
    Increased = as.integer(PValue < 0.05 & logFC >  1),
    Decreased = as.integer(PValue < 0.05 & logFC < -1)
  )

# 4. 统计每条 Path 的 KO 数及纯保留率等指标
pathway_stats <- combined %>%
  group_by(Group, PathType, Path) %>%
  summarise(
    Total_KO     = n_distinct(KO),
    Unchanged_KO = n_distinct(KO[Unchanged == 1]),
    Increased_KO = n_distinct(KO[Increased == 1]),
    Decreased_KO = n_distinct(KO[Decreased == 1]),
    .groups = "drop"
  ) %>%
  filter(Total_KO >= 5) %>%
  mutate(
    FRR      = 100 * Unchanged_KO / Total_KO,
    ActRate  = 100 * Increased_KO / Total_KO,
    LossRate = 100 * Decreased_KO / Total_KO
  )

# 5. 配对比较 IT-H vs ST-H 并排序
retention_comp <- pathway_stats %>%
  dplyr::select(Group, Path, FRR) %>%
  pivot_wider(names_from = Group, values_from = FRR) %>%
  mutate(
    FRR_diff    = `IT-H` - `ST-H`,
    Delta_Label = ifelse(FRR_diff > 0,
                         paste0("+",round(FRR_diff,1),"%"),
                         paste0(round(FRR_diff,1),"%"))
  ) %>%
  filter(!is.na(FRR_diff))
ordered_paths <- retention_comp %>% arrange(FRR_diff) %>% pull(Path)

# 6. 全局 Wilcoxon & LMM 检验
wilcox_res <- wilcox.test(retention_comp$`IT-H`, retention_comp$`ST-H`,
                          paired = TRUE, exact = FALSE)
message("Global Wilcoxon p-value: ", signif(wilcox_res$p.value,3))

model_data <- retention_comp %>%
  pivot_longer(cols = c(`IT-H`,`ST-H`), names_to = "Group", values_to = "FRR")
lmm <- lmer(FRR ~ Group + (1|Path), data = model_data)
print(tidy(lmm))

# 7. 每 Path 的组间检验与 FDR 校正
compare_by_pathway <- model_data %>%
  group_by(Path) %>%
  summarise(
    p_raw = wilcox.test(FRR[Group=="IT-H"], FRR[Group=="ST-H"], exact=FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj     = p.adjust(p_raw, method = "fdr"),
    sig_label = case_when(
      p_adj<0.001 ~ "***",
      p_adj<0.01  ~ "**",
      p_adj<0.05  ~ "*",
      TRUE         ~ ""
    )
  )

# 8. 合并星号、计算标注位置
retention_plot_data2 <- model_data %>%
  left_join(compare_by_pathway %>% dplyr::select(Path, sig_label), by = "Path") %>%
  group_by(Path) %>%
  mutate(
    max_rate   = max(FRR, na.rm = TRUE),
    y_position = max_rate + 5
  ) %>%
  ungroup()

# 9. 可视化: ΔFRR 柱状图+星号
p_main <- retention_plot_data2 %>%
  mutate(Path = factor(Path, levels = ordered_paths)) %>%
  ggplot(aes(x = Path, y = FRR, fill = Group)) +
  geom_col(position = position_dodge(width=0.8), width=0.7, color="black") +
  geom_text(
    data = retention_plot_data2 %>% filter(Group=="ST-H"),
    aes(y=y_position, label=sig_label),
    position = position_dodge(width=0.8), vjust=0, size=5
  ) +
  coord_flip() +
  scale_fill_manual(values=c("IT-H"="#E41A1C","ST-H"="#377EB8")) +
  labs(x="Pathway", y="Retention Rate (%)",
       title="Per-Pathway Retention (IT-H vs ST-H)") +
  theme_minimal(base_size=13) + theme(legend.position="top")

# 10. 热图 & 箱线图
p_heatmap <- pathway_stats %>%
  filter(Group %in% c("IT-H","ST-H")) %>%
  mutate(Path=factor(Path, levels=ordered_paths)) %>%
  ggplot(aes(x=Group,y=Path,fill=FRR)) +
  geom_tile(color="white",linewidth=0.5) +
  geom_text(aes(label=round(FRR,1)),size=3) +
  scale_fill_gradient2(low="#3288BD",mid="#FFFFBF",high="#D53E4F",
                       midpoint=50,limits=c(0,100)) +
  labs(fill="FRR (%)",title="FRR Heatmap") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,hjust=1))

p_rates <- pathway_stats %>%
  pivot_longer(cols=c(ActRate,LossRate),names_to="Metric",values_to="Val") %>%
  mutate(Metric=recode(Metric,ActRate="Activation",LossRate="Loss")) %>%
  ggplot(aes(x=Metric,y=Val,fill=Group)) + geom_boxplot() +
  labs(y="%",title="Activation vs Loss Rate") + theme_minimal()

final_plot <- (p_main / (p_heatmap | p_rates)) + plot_layout(heights=c(2,2))

# 导出
write_csv(pathway_stats, "pathway_stats.csv")
write_csv(compare_by_pathway, "pathway_pvalues.csv")
ggsave("retention_overview.pdf", final_plot, width=12, height=10, dpi=300)
print(final_plot)


















###########################################################
#6.0
# 1. 数据准备与Level2功能归类 - 增强质量检查 --------------------------------
library(tidyverse)
library(readr)
library(lme4)       # 用于更稳健的混合效应模型
library(multcomp)   # 用于多重比较校正
dir()
# 读取数据并检查完整性
kegg_paths <- read_csv("KEGG_paths.csv")
it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C.csv")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C.csv")

# 数据质量检查
message(paste("KEGG_paths contains", nrow(kegg_paths), "rows with", sum(!is.na(kegg_paths$KO)), "non-missing KOs"))
message(paste("IT-H data contains", nrow(it_h), "rows with", sum(!is.na(it_h$KO)), "non-missing KOs"))
message(paste("ST-H data contains", nrow(st_h), "rows with", sum(!is.na(st_h$KO)), "non-missing KOs"))

# 标准化K编号格式并检查匹配率
extract_ko <- function(x) str_extract(x, "K\\d{5}")  # 更严格的K编号匹配

kegg_paths$KO <- extract_ko(kegg_paths$KO)
it_h$KO <- extract_ko(it_h$KO)
st_h$KO <- extract_ko(st_h$KO)

# 报告K编号匹配情况
shared_kos <- intersect(kegg_paths$KO, intersect(it_h$KO, st_h$KO))
message(paste(length(shared_kos), "KOs are shared across all datasets"))
message(paste("Coverage:", round(length(shared_kos)/nrow(kegg_paths)*100, 1), "% of KEGG_paths KOs are analyzed"))

# 2. 改进的核心代谢通路选择 - 基于生物学重要性 ----------------------------
# 定义真正的核心代谢通路 (根据珊瑚-微生物组代谢互作文献)
strict_core_pathways <- c(
  "Carbohydrate metabolism", "Energy metabolism",
  "Lipid metabolism",       "Nucleotide metabolism",
  "Amino acid metabolism", "Metabolism of other amino acids",
  "Glycan biosynthesis and metabolism",
  "Metabolism of cofactors and vitamins",
    "Biosynthesis of other secondary metabolites",
  "Xenobiotics biodegradation and metabolism",
  "Chemical structure transformation maps",
   
   "Glutathione metabolism",     # 抗氧化
  "Oxidative phosphorylation",
  "Protein processing in endoplasmic reticulum",
  "Proteasome",                # 受损蛋白降解
  "Ubiquitin mediated proteolysis",
  "Carotenoid biosynthesis",
  "Two-component system",
  "MAPK signaling pathway",
  "Base excision repair",       # DNA 修复
  "Nucleotide excision repair",
   
   "Signal transduction",
  "Metabolism of terpenoids and polyketides"
  
)

# 定义应激响应通路 (用于辅助分析)
stress_pathways <- c(
  "Metabolism of other amino acids",
  "Signal transduction",
  "Replication and repair","Metabolism of terpenoids and polyketides",
  "Metabolism of cofactors and vitamins","Glutathione metabolism",     # 抗氧化
  "Oxidative phosphorylation",
  "Carotenoid biosynthesis",
  "Protein processing in endoplasmic reticulum",
  "Proteasome",                # 受损蛋白降解
  "Ubiquitin mediated proteolysis",
  "Base excision repair",       # DNA 修复
  "Nucleotide excision repair",
  "MAPK signaling pathway",
  "Two-component system"
)

# 拆分Level2并分类
kegg_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(
    Level2 = str_split(paths, "; ", simplify = TRUE)[,2],
    PathwayType = case_when(
      Level2 %in% strict_core_pathways ~ "Core Metabolic",
      Level2 %in% stress_pathways ~ "Stress Response",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(Level2))

# 3. 改进的功能状态分类系统 ----------------------------
classify_functional_status <- function(df) {
  df %>%
    mutate(
      Status = case_when(
        # 条件1：显著上调（P<0.01且logFC>1）
        PValue < 0.01 & logFC > 1 ~ "Significantly Increased",
        
        # 条件2：显著下调（P<0.01且logFC<-1）
        PValue < 0.01 & logFC < -1 ~ "Significantly Decreased",
        
        # 其他所有情况（包括PValue不显著、logFC介于-1和1之间、或NA值）均归类为"Unchanged"
        TRUE ~ "Unchanged"
      ),
      # 固定因子顺序（可选，便于后续绘图）
      Status = factor(Status, levels = c("Significantly Increased", "Unchanged", "Significantly Decreased"))
    )
}


classify_functional_status <- function(df) {
  df %>%
    mutate(
      Status = case_when(
        is.na(PValue) | is.na(logFC) ~ "Unchanged",
        PValue < 0.01 & logFC < -1 ~ "Significantly Decreased",
        PValue < 0.01 & logFC > 1 ~ "Significantly Increased",
        TRUE ~ "Unchanged"
      ),
      Status = factor(Status, levels = c("Significantly Increased", "Unchanged", "Significantly Decreased"))
    )
}

it_h_core <- left_join(it_h, kegg_long, by = "KO") %>%
  filter(PathwayType == "Core Metabolic") %>%
  classify_functional_status()

st_h_core <- left_join(st_h, kegg_long, by = "KO") %>%
  filter(PathwayType == "Core Metabolic") %>%
  classify_functional_status()

# 4. 改进的统计分析与可视化 ----------------------------
# 4.1 计算各通路状态比例
# 修正后的函数
calculate_pathway_stats <- function(df, group_name) {
  df %>%
    filter(!is.na(Level2), !is.na(Status)) %>%  # 修改这里
    group_by(Level2, Status) %>%
    summarise(KO_num = n(), .groups = "drop") %>%
    complete(Level2, Status, fill = list(KO_num = 0)) %>%
    group_by(Level2) %>%
    mutate(
      total = sum(KO_num),
      Percent = round(100 * KO_num / total, 1),
      Group = group_name
    ) %>%
    ungroup()
}
stat_it <- calculate_pathway_stats(it_h_core, "IT-H")
stat_st <- calculate_pathway_stats(st_h_core, "ST-H")
stat_all <- bind_rows(stat_it, stat_st) %>%
  mutate(Level2 = factor(Level2, levels = unique(Level2)))

# 4.2 使用混合效应模型进行更稳健的统计检验
# 准备建模数据
# 修正后的代码
model_data <- stat_all %>%
  dplyr::select(Level2, Group, Status, KO_num) %>%
  tidyr::pivot_wider(
    names_from = Status, 
    values_from = KO_num, 
    values_fill = 0
  ) %>%
  dplyr::mutate(
    total = `Significantly Increased` + Unchanged + `Significantly Decreased`,
        retention_rate = 100 * (Unchanged ) / total,
    Group = factor(Group, levels = c("IT-H", "ST-H"))
  )
model_data <- stat_all %>%
  dplyr::select(Level2, Group, Status, KO_num) %>%
  tidyr::pivot_wider(
    names_from = Status, 
    values_from = KO_num, 
    values_fill = 0
  ) %>%
  dplyr::mutate(
    total = `Significantly Increased` + Unchanged + `Significantly Decreased`,
    retention_rate = 100 * (Unchanged ) / total,
    Group = factor(Group, levels = c("IT-H", "ST-H"))
  )
# 拟合混合效应模型
retention_model <- lmer(
  retention_rate ~ Group + (1|Level2), 
  data = model_data
)

# 多重比较校正
retention_contrasts <- glht(
  retention_model, 
  linfct = mcp(Group = "Tukey")
)
retention_summary <- summary(retention_contrasts, test = adjusted("fdr"))

head(retention_summary)
# 提取显著性结果
sig_results <- broom::tidy(retention_summary) %>%
  mutate(
    sig_label = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

head(sig_results) 

# 5. 改进的可视化 ----------------------------
library(ggplot2)
library(ggrepel)

# 5.1 保留率对比图 (主图)
retention_plot_data <- model_data %>%
  group_by(Level2) %>%
  mutate(
    max_rate = max(retention_rate),
    y_position = max_rate + 5
  ) %>%
  # 修正匹配逻辑 - 直接添加显著性标记
  mutate(
    sig_label = ifelse(Group == "IT-H", 
                       unique(sig_results$sig_label), 
                       "")
  )


ggplot(retention_plot_data, aes(x = Level2, y = retention_rate, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
  geom_text(
    aes(y = y_position, label = sig_label),
    position = position_dodge(width = 0.8),
    vjust = 0, size = 6, color = "black"
  ) +
  scale_fill_manual(
    values = c("IT-H" = "#E64B35", "ST-H" = "#4DBBD5"),
    name = "Treatment"
  ) +
  labs(
    x = "Core Metabolic Pathway",
    y = "Functional Retention Rate (%)",
    title = "Core Metabolic Pathway Retention Under Heat Stress",
    subtitle = "Mixed-effects model with FDR correction: *p<0.05, **p<0.01, ***p<0.001"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.1)))

# 5.2 功能状态组成图 (补充图)
library(ggplot2)
library(ggridges)  # 用于河流图
library(patchwork) # 用于图形组合
library(viridis)   # 用于科学配色

## 数据准备 --------------------------------
# 确保stat_all包含正确的百分比数据
plot_data <- stat_all %>%
  mutate(
    Level2 = factor(Level2, levels = strict_core_pathways),
    Group = factor(Group, levels = c("IT-H", "ST-H")),
    Status = factor(Status, levels = c("Significantly Increased", "Unchanged", "Significantly Decreased"))
  )

## ISME期刊配色方案 ------------------------
isme_colors <- c(
  "Significantly Increased" = "#E64B35",  # 海绿色
  "Unchanged" = "#F5F5F5",               # 极浅灰
  "Significantly Decreased" = "#4DBBD5"   # 深红色
)

## 1. 百分比堆叠柱状图 (精细优化版) --------
stacked_bar <- ggplot(plot_data, aes(x = Group, y = Percent, fill = Status)) +
  geom_col(position = "fill", width = 0.7, color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(Percent > 5, paste0(round(Percent), "%"), "")),
    position = position_fill(vjust = 0.5),
    size = 3, color = "black"
  ) +
  facet_wrap(~ Level2, nrow = 1) +
  scale_fill_manual(values = isme_colors) +
  scale_y_continuous(
    labels = scales::percent,
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Percentage of Genes",
    title = "Functional State Composition by Metabolic Pathway",
    fill = "Expression Change"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(stacked_bar)

# 6. 结果导出 ----------------------------
# 保存统计结果
write_csv(model_data, "core_pathway_retention_stats.csv")
write_csv(sig_results, "statistical_test_results.csv")

# 保存模型摘要
sink("model_summary.txt")
print(summary(retention_model))
print(retention_summary)
sink()






















###################################################
#5.0
#【1. 数据准备与Level2功能归类】
library(tidyverse)
library(readr)
dir()
# 读取数据
kegg_paths <- read_csv("KEGG_paths.csv")
it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C-3.csv")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C-3.csv")

kegg_paths$KO <- str_extract(kegg_paths$KO, "K\\d+")
it_h$KO <- str_extract(it_h$KO, "K\\d+")
st_h$KO <- str_extract(st_h$KO, "K\\d+")

# 拆分Level2
kegg_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(
    Level2 = str_split(paths, "; ", simplify = TRUE)[,2]
  ) %>%
  filter(!is.na(Level2))

# 推荐Level2主功能类
key_level2 <- c(
  "Carbohydrate metabolism", "Energy metabolism", "Amino acid metabolism",
  "Lipid metabolism", "Nucleotide metabolism", "Environmental adaptation",
  "Signal transduction", "Membrane transport", "Replication and repair",
  "Xenobiotics biodegradation and metabolism"
)

core_ko_level2 <- kegg_long %>%
  filter(Level2 %in% key_level2) %>%
  select(KO, CorePathway = Level2) %>%
  distinct()

#【2. 保留/丢失判定与统计】
# 阈值判定
# Retained: 其它；Lost: PValue < 0.05 且 abs(logFC) > 1 且 logFC < 0
it_h_core <- left_join(it_h, core_ko_level2, by = "KO") %>%
  mutate(Status = case_when(
    is.na(CorePathway) ~ NA_character_,
    PValue < 0.05 & abs(logFC) > 1 & logFC < 0 ~ "Lost",
    TRUE ~ "Retained"
  ))

st_h_core <- left_join(st_h, core_ko_level2, by = "KO") %>%
  mutate(Status = case_when(
    is.na(CorePathway) ~ NA_character_,
    PValue < 0.05 & abs(logFC) > 1 & logFC < 0 ~ "Lost",
    TRUE ~ "Retained"
  ))

# 汇总
stat_it <- it_h_core %>%
  filter(!is.na(CorePathway), !is.na(Status)) %>%
  group_by(CorePathway, Status) %>%
  summarise(KO_num = n(), .groups = "drop") %>%
  mutate(Group = "Intertidal (IT-H)")

stat_st <- st_h_core %>%
  filter(!is.na(CorePathway), !is.na(Status)) %>%
  group_by(CorePathway, Status) %>%
  summarise(KO_num = n(), .groups = "drop") %>%
  mutate(Group = "Subtidal (ST-H)")

stat_all <- bind_rows(stat_it, stat_st)

#【3. 计算保留率、显著性检验（Fisher's Exact）和可视化表】
# 统计表（宽表形式）
stat_wide <- stat_all %>%
  mutate(Group = recode(Group, 
                        "Intertidal (IT-H)" = "IT_H", 
                        "Subtidal (ST-H)" = "ST_H")) %>%
  pivot_wider(names_from = Status, values_from = KO_num, values_fill = 0) %>%
  mutate(
    total = Retained + Lost,
    retention_rate = round(100 * Retained / total, 1),
    lost_rate = round(100 * Lost / total, 1)
  )

# 差异检验前，再次pivot_wider
stat_compare <- stat_wide %>%
  select(CorePathway, Group, Retained, Lost, total, retention_rate) %>%
  pivot_wider(names_from = Group, values_from = c(Retained, Lost, total, retention_rate), values_fill = 0) %>%
  rowwise() %>%
  mutate(
    p_value = {
      mat <- matrix(c(Retained_IT_H, Lost_IT_H, Retained_ST_H, Lost_ST_H), nrow = 2, byrow = TRUE)
      if (all(!is.na(mat)) && sum(mat) > 0 && all(mat >= 0)) {
        tryCatch(fisher.test(mat)$p.value, error = function(e) NA_real_)
      } else NA_real_
    },
    sig_label = case_when(
      is.na(p_value) ~ "",
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  ungroup()

# 长表供绘图
plot_data <- stat_all %>%
  left_join(stat_all %>%
              group_by(CorePathway, Group) %>%
              summarise(total = sum(KO_num), .groups = "drop"),
            by = c("CorePathway", "Group")) %>%
  mutate(
    Percent = round(100 * KO_num / total, 1),
    CorePathway = factor(CorePathway, levels = key_level2),
    Status = factor(Status, levels = c("Lost", "Retained"))
  )


#【4. 顶刊风格堆叠柱状图 + 显著性注释】
library(ggplot2)
library(ggrepel)

# 为可视化做整理：每个level2一行，IT-H和ST-H分列
# 我们只画“Retained”保留率，并在两组之间连线/星号标注
retention_plot <- stat_compare %>%
  select(CorePathway, retention_rate_IT_H = retention_rate_IT_H, retention_rate_ST_H = retention_rate_ST_H, sig_label) %>%
  pivot_longer(
    cols = starts_with("retention_rate_"),
    names_to = "Group",
    names_prefix = "retention_rate_",
    values_to = "Retention"
  ) %>%
  mutate(
    Group = recode(Group, "IT_H" = "IT-H", "ST_H" = "ST-H"),
    CorePathway = factor(CorePathway, levels = key_level2),
    Group = factor(Group, levels = c("IT-H", "ST-H"))
  )

# 计算星号标注位置（画在两柱顶上方的中点）
stat_compare$star_y <- apply(stat_compare[,c("retention_rate_IT_H","retention_rate_ST_H")], 1, function(x) max(x, na.rm = TRUE) + 5)

# 合并显著性注释（x为CorePathway, y为中点）
star_anno <- stat_compare %>%
  mutate(
    CorePathway = factor(CorePathway, levels = key_level2),
    x = CorePathway,
    y = star_y,
    label = sig_label
  ) %>%
  filter(label != "")

# 并列柱状图+星号标注
ggplot(retention_plot, aes(x = CorePathway, y = Retention, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6, color = "black") +
  # 星号标注在两柱之间
  geom_text(
    data = star_anno,
    aes(x = CorePathway, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 0, size = 7, color = "#E41A1C", fontface = "bold"
  ) +
  scale_fill_manual(values = c("IT-H" = "#E41A1C", "ST-H" = "#377EB8")) +
  labs(
    x = "KEGG Level2 pathway",
    y = "KO retention rate (%)",
    title = "Differential retention of core metabolic functions",
    subtitle = "*p<0.05, **p<0.01, ***p<0.001"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 14, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(limits = c(0, 120), expand = expansion(mult = c(0, 0.1)))

##########################################################################



#【5. 输出统计与补充表（可选）
write_csv(stat_compare, "KEGG_Level2_Retention_Fisher.csv")
write_csv(plot_data, "KEGG_Level2_PlotData.csv")























library(tidyverse)
# 假设你的KEGG_paths.csv格式为 KO, paths

kegg_paths <- read_csv("KEGG_paths.csv")
kegg_paths$KO <- str_extract(kegg_paths$KO, "K\\d+")

# 拆分所有Level3通路
kegg_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(
    Level1 = str_split(paths, "; ", simplify = TRUE)[,1],
    Level2 = str_split(paths, "; ", simplify = TRUE)[,2],
    Level3 = str_split(paths, "; ", simplify = TRUE)[,3]
  )

# 查看都有哪些Level2/Level3（辅助你挑重要通路）
table(kegg_long$Level2)
table(kegg_long$Level3)















library(tidyverse)
library(readr)

# ===== 1. 数据读取 =====
core_ko <- read_csv("Core_KO.csv") # KO, CorePathway
kegg_paths <- read_csv("KEGG_paths.csv") # KO, paths

it_h <- read_csv("KEGG_diff_results_IT-H_vs_IT-C-1.csv")
st_h <- read_csv("KEGG_diff_results_ST-H_vs_ST-C-1.csv")

# ===== 2. KO-通路归类（多对多合并） =====
# 确保KO格式统一
core_ko$KO <- str_extract(core_ko$KO, "K\\d+")
it_h$KO <- str_extract(it_h$KO, "K\\d+")
st_h$KO <- str_extract(st_h$KO, "K\\d+")

# ===== 3. 合并通路分类到差异表 =====
it_h_core <- left_join(it_h, core_ko, by = "KO")
st_h_core <- left_join(st_h, core_ko, by = "KO")

# ===== 4. 标记“丢失”/“保留”KO（以Significance为准）=====
# 这里假定Significance包含“IT_C enriched”或“ST_C enriched”为丢失，其余为保留
it_h_core <- it_h_core %>%
  mutate(Status = if_else(str_detect(Significance, "IT_C enriched"), "Lost", "Retained"))
st_h_core <- st_h_core %>%
  mutate(Status = if_else(str_detect(Significance, "ST_C enriched"), "Lost", "Retained"))

# ===== 5. 汇总每个核心通路的保留/丢失KO数 =====
stat_it <- it_h_core %>%
  filter(!is.na(CorePathway)) %>%
  group_by(CorePathway, Status) %>%
  summarise(KO_num = n(), .groups = "drop") %>%
  mutate(Group = "Intertidal (IT-H)")

stat_st <- st_h_core %>%
  filter(!is.na(CorePathway)) %>%
  group_by(CorePathway, Status) %>%
  summarise(KO_num = n(), .groups = "drop") %>%
  mutate(Group = "Subtidal (ST-H)")

# 合并
stat_all <- bind_rows(stat_it, stat_st)

# ===== 6. 计算每个通路的总KO数、保留率、丢失率 =====
stat_wide <- stat_all %>%
  pivot_wider(names_from = Status, values_from = KO_num, values_fill = 0) %>%
  mutate(
    total = Retained + Lost,
    retention_rate = round(100 * Retained / total, 1),
    lost_rate = round(100 * Lost / total, 1)
  )

# ===== 7. 可视化数据准备（堆叠格式） =====
plot_data <- stat_all %>%
  left_join(stat_all %>%
              group_by(CorePathway, Group) %>%
              summarise(total = sum(KO_num), .groups = "drop"),
            by = c("CorePathway", "Group")) %>%
  mutate(Percent = round(100 * KO_num / total, 1))

# 只展示主分析重点通路，可修改下面的通路列表
core_pathways_show <- c("TCA cycle", "Glycolysis", "Oxidative phosphorylation", 
                        "Fatty acid metabolism", "Amino acid metabolism", "Pyruvate metabolism", 
                        "HSP response", "ROS response")

plot_data <- plot_data %>%
  filter(CorePathway %in% core_pathways_show) %>%
  mutate(
    CorePathway = factor(CorePathway, levels = core_pathways_show),
    Status = factor(Status, levels = c("Lost", "Retained"))
  )

# ===== 8. 顶刊风格堆叠柱状图 =====
library(ggplot2)

ggplot(plot_data, aes(x = CorePathway, y = Percent, fill = Status)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +
  facet_wrap(~Group, nrow = 1) +
  scale_fill_manual(values = c("Lost" = "grey80", "Retained" = "#E41A1C")) +
  labs(
    x = "Core metabolic pathway",
    y = "KO retention rate (%)",
    title = "Retention of core metabolic pathway KOs under heat stress",
    subtitle = "Red = Retained; Grey = Lost; Pathway assignment based on KEGG/curated"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 15, face = "bold"),
    strip.text = element_text(size = 15, face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "top"
  ) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0.05)))

head(plot_data)
write.csv(plot_data, "plot_data.csv", row.names = FALSE)
write.csv(stat_it, "stat_it.csv", row.names = FALSE)
write.csv(stat_st, "stat_st.csv", row.names = FALSE)


# 输出归类比例
cat("归类KO比例 - IT-H:", mean(!is.na(it_h_core$CorePathway)), "\n")
cat("归类KO比例 - ST-H:", mean(!is.na(st_h_core$CorePathway)), "\n")
# 查看丢失/保留的总数
stat_all %>% group_by(Group, Status) %>% summarise(n = sum(KO_num))







###############################################################################
library(tidyverse)

# 1. 读取通路注释
kegg_paths <- read.csv("KEGG_paths.csv", stringsAsFactors = FALSE)
kegg_paths$KO <- str_extract(kegg_paths$KO, "K\\d+")

# 提取Level3路径
kegg_paths_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(Level3 = str_split(paths, "; ", simplify = TRUE)[, 3]) %>%
  filter(!is.na(Level3))

ROS_KOs <- c(
  "K04564", # SOD, Fe-Mn型
  "K04565", # SOD, Fe-Mn型
  "K00086", # SOD, Cu/Zn型
  "K03781", # Catalase
  "K00432", # Glutathione peroxidase
  "K00383", # Glutathione reductase
  "K00799", # Glutathione S-transferase
  "K03386", # Peroxiredoxin
  "K03671", # Thioredoxin peroxidase
  "K08348", # NADPH oxidase
  "K03385", # Alkyl hydroperoxide reductase
  "K05532", # Glutaredoxin
  "K01920", # Glutathione synthase
  "K01919"  # Glutamate-cysteine ligase
)

HSP_KOs <- c(
  "K04043", # DnaK (HSP70)
  "K04468", # GrpE (HSP70辅助因子)
  "K03686", # DnaJ (HSP40)
  "K04077", # GroEL (HSP60)
  "K04078", # GroES (HSP60辅助因子)
  "K04079", # HtpG (HSP90)
  "K04080", # HtpG (HSP90)
  "K13993", # sHSP/IbpA/IbpB (小分子热休克蛋白)
  "K03695", # ClpB (HSP100)
  "K03692", # ClpC/X/P (HSP100家族, 多个蛋白)
  "K06761", # HslO (HSP33)
"K13993"
  
  )



# 2. 定义手动关键词映射（生物学高标准）
assign_core_pathway <- function(level3, ko) {
  case_when(
    ko %in% HSP_KOs ~ "HSP response",
    ko %in% ROS_KOs ~ "ROS response",
    str_detect(level3, regex("Citrate cycle \\(TCA cycle\\)", ignore_case = TRUE)) ~ "TCA cycle",
    str_detect(level3, regex("Pyruvate metabolism", ignore_case = TRUE)) ~ "Pyruvate metabolism",
    str_detect(level3, regex("Glycolysis", ignore_case = TRUE)) ~ "Glycolysis",
    str_detect(level3, regex("Gluconeogenesis", ignore_case = TRUE)) ~ "Glycolysis",
    str_detect(level3, regex("Oxidative phosphorylation", ignore_case = TRUE)) ~ "Oxidative phosphorylation",
    str_detect(level3, regex("Electron transport chain", ignore_case = TRUE)) ~ "Oxidative phosphorylation",
    str_detect(level3, regex("Fatty acid", ignore_case = TRUE)) ~ "Fatty acid metabolism",
    str_detect(level3, regex("Beta-oxidation", ignore_case = TRUE)) ~ "Fatty acid metabolism",
    str_detect(level3, regex("Amino acid metabolism", ignore_case = TRUE)) ~ "Amino acid metabolism",
    str_detect(level3, regex("Valine, leucine and isoleucine degradation", ignore_case = TRUE)) ~ "Amino acid metabolism",
    str_detect(level3, regex("Glutamine and glutamate metabolism", ignore_case = TRUE)) ~ "Amino acid metabolism",
    TRUE ~ NA_character_
  )
}

kegg_paths_long <- kegg_paths_long %>%
  mutate(CorePathway = assign_core_pathway(Level3, KO))


# 3. 提取“核心通路KO”清单（只保留有明确分类的KO）
core_ko_table <- kegg_paths_long %>%
  filter(!is.na(CorePathway)) %>%
  select(KO, CorePathway) %>%
  distinct()

# 你可以用如下方式，直接获得TCA等功能类别及其KO列表：
# eg:
TCA_KOs <- core_ko_table %>% filter(CorePathway == "TCA cycle") %>% pull(KO)
Glycolysis_KOs <- core_ko_table %>% filter(CorePathway == "Glycolysis") %>% pull(KO)
OxPhos_KOs <- core_ko_table %>% filter(CorePathway == "Oxidative phosphorylation") %>% pull(KO)

# 如果需要全表导出，方便人工核查
write.csv(core_ko_table, "Core-KO.csv", row.names = FALSE)




########################################################
########################################################
library(tidyverse)
library(clusterProfiler)

# ==== 1. 读取和预处理数据 ====
it_h_data <- read.csv("KEGG_diff_results_IT-H_vs_IT-C.csv", stringsAsFactors = FALSE)
st_h_data <- read.csv("KEGG_diff_results_ST-H_vs_ST-C.csv", stringsAsFactors = FALSE)
kegg_paths <- read.csv("KEGG_paths.csv", stringsAsFactors = FALSE)

# 统一KO格式
it_h_data$KO <- str_extract(it_h_data$KO, "K\\d+")
st_h_data$KO <- str_extract(st_h_data$KO, "K\\d+")
kegg_paths$KO <- str_extract(kegg_paths$KO, "K\\d+")

# ==== 2. 筛选显著下调基因 ====
it_h_lost <- it_h_data %>% filter(grepl("IT_C enriched", Significance)) %>% select(KO) %>% distinct()
st_h_lost <- st_h_data %>% filter(grepl("ST_C enriched", Significance)) %>% select(KO) %>% distinct()

# ==== 3. KEGG二级通路分解 ====
kegg_paths_long <- kegg_paths %>%
  separate_rows(paths, sep = " \\| ") %>%
  mutate(level2 = str_extract(paths, "(?<=; )[^;]+(?=; )")) %>%
  select(KO, level2) %>%
  filter(!is.na(level2))

# ==== 4. 统计保留率 ====
# 分别处理IT-H和ST-H
get_retention <- function(all_data, lost_kos, group_name, kegg_paths_long) {
  # 合并注释
  data_annot <- all_data %>%
    select(KO) %>%
    left_join(kegg_paths_long, by = "KO") %>%
    filter(!is.na(level2))
  
  # 统计通路内KO总数
  total_tab <- data_annot %>% group_by(level2) %>% summarise(total_genes = n_distinct(KO))
  
  # 统计通路内丢失KO数量
  lost_in_pathway <- data_annot %>% filter(KO %in% lost_kos$KO) %>% group_by(level2) %>% summarise(lost_genes = n_distinct(KO))
  
  # 合并并计算保留率
  retention_tab <- total_tab %>%
    left_join(lost_in_pathway, by = "level2") %>%
    mutate(lost_genes = replace_na(lost_genes, 0),
           retained_genes = total_genes - lost_genes,
           retention_rate = retained_genes / total_genes * 100,
           group = group_name)
  
  return(retention_tab)
}

it_h_retention <- get_retention(it_h_data, it_h_lost, "IT-H", kegg_paths_long)
st_h_retention <- get_retention(st_h_data, st_h_lost, "ST-H", kegg_paths_long)

# ==== 5. Fisher显著性检验（和结果合并） ====
# 合并两组
combined <- full_join(it_h_retention, st_h_retention, by = "level2", suffix = c("_it", "_st"))

# 对每个通路做Fisher检验
library(purrr)

get_fisher_p <- function(retained_genes_it, lost_genes_it, retained_genes_st, lost_genes_st, ...) {
  p_value <- NA_real_
  if(!any(is.na(c(retained_genes_it, lost_genes_it, retained_genes_st, lost_genes_st)))) {
    matrix_test <- matrix(
      c(retained_genes_it, lost_genes_it, retained_genes_st, lost_genes_st),
      nrow = 2, byrow = TRUE
    )
    if(all(matrix_test >= 0) & sum(matrix_test) > 0) {
      ft <- fisher.test(matrix_test)
      p_value <- ft$p.value
    }
  }
  return(p_value)
}

combined$p_value <- pmap_dbl(
  combined[, c("retained_genes_it", "lost_genes_it", "retained_genes_st", "lost_genes_st")], 
  get_fisher_p
)

# ==== 6. 筛选并整理最终结果 ====
final_tab <- combined %>%
  filter(!is.na(p_value), (total_genes_it + total_genes_st) > 5) %>%
  filter(p_value < 0.05) %>%
  arrange(desc(retention_rate_it))

# ==== 7. 可视化 ====
library(ggplot2)
plot_data <- final_tab %>%
  select(level2, retention_rate_it, retention_rate_st, p_value) %>%
  pivot_longer(cols = starts_with("retention_rate"), 
               names_to = "group", 
               values_to = "retention_rate") %>%
  mutate(group = recode(group, retention_rate_it = "IT-H", retention_rate_st = "ST-H"))

ggplot(plot_data, aes(x = reorder(level2, -retention_rate), y = retention_rate, fill = group)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(
    data = plot_data %>% filter(group == "IT-H"),  # 只在一组上标注即可
    aes(label = ifelse(final_tab$p_value < 0.05, ifelse(final_tab$p_value < 0.01, "**", "*"), "")),
    position = position_dodge(width = 0.7), 
    vjust = -0.5, size = 5
  ) +
  scale_fill_manual(values = c("IT-H" = "#E41A1C", "ST-H" = "#377EB8")) +
  labs(
    x = "KEGG二级通路",
    y = "通路保留率 (%)",
    title = "热胁迫下核心代谢通路的保留率差异",
    subtitle = "显著性标注：*p<0.05, **p<0.01"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top"
  ) +
  scale_y_continuous(limits = c(0, 100))

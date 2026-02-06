
#####################################################################
#####################################################################
# ────────────────── 0. 载入包 ────────────────────────────────
library(tidyverse)   # 数据整理 + ggplot2
library(vegan)       # vegdist
library(FD)          # FRic / FDiv
library(ggpubr)      # stat_compare_means
library(ggsci)       # 色板
dir()

# ────────────────── 1. 读入数据 ──────────────────────────────
# a) KO × Sample (CSS 标准化) —— 用于 FNB & FDiv
abun <- read.csv("CSS_relative_KEGG_bases_insitu.csv",
                 header = TRUE, row.names = 1, check.names = FALSE)
meta <- read.csv("metadata.csv", stringsAsFactors = FALSE)

# b) KEGG -> 功能性状 —— 仅 FDiv 需要
paths <- read.csv("KEGG_paths.csv", stringsAsFactors = FALSE)

# ————— 1a. 提取组别样本向量
it_samples <- meta$Sample[meta$Group == "IT"]
st_samples <- meta$Sample[meta$Group == "ST"]

# ────────────────── 2. 计算 KO 级 FNB (标准化) ────────────────
## 拆分丰度矩阵
it_kegg <- abun[, it_samples]
st_kegg <- abun[, st_samples]

## FNB 函数
calc_fnb <- function(mat){
  apply(mat, 1, function(v){
    if(sum(v) == 0) return(NA_real_)   # 全零 KO
    p <- v / sum(v)
    1 / sum(p^2)
  })
}

it_fnb <- calc_fnb(it_kegg) / length(it_samples)   # 标准化
st_fnb <- calc_fnb(st_kegg) / length(st_samples)

fnb_long <- tibble(
  KO     = rownames(abun),
  IT     = it_fnb,
  ST     = st_fnb
) %>%
  pivot_longer(cols = c(IT, ST), names_to = "Group", values_to = "Value") %>%
  drop_na() %>%                        # 去掉全零 KO
  mutate(Index = "FNB")                # 标记指标

# ────────────────── 3. 计算样本级 FDiv ───────────────────────
## — 3a. 构建性状矩阵 (Heat_Score / Basic_Score 作为示例) —
heat_patterns  <- c("Oxidative phosphorylation","Glutathione metabolism",
                    "DMSP","Antioxidant","Heat shock")
basic_patterns <- c("Glycolysis / Gluconeogenesis","TCA",
                    "Lipid metabolism","Amino acid metabolism",
                    "Nucleotide metabolism")

trait <- paths %>%
  mutate(
    Heat_Score  = if_else(str_detect(paths, str_c(heat_patterns ,collapse="|")), 2, 0),
    Basic_Score = if_else(str_detect(paths, str_c(basic_patterns,collapse="|")), 1, 0)
  ) %>%
  filter(Heat_Score > 0 | Basic_Score > 0) %>%
  select(KO, Heat_Score, Basic_Score) %>%
  column_to_rownames("KO") %>% 
  as.matrix()

## — 3b. 留交集 KO 并计算 FDiv —
common_KO  <- intersect(rownames(trait), rownames(abun))
trait_sub  <- trait[common_KO, ]
abun_sub   <- abun[common_KO, ]

fd <- dbFD(
  x        = trait_sub,
  a        = t(abun_sub),
  stand.x  = TRUE,
  corr     = "cailliez",
  calc.FDiv = TRUE,
  calc.FRic = TRUE,
  messages  = FALSE
)

fdiv_long <- tibble(
  Sample = names(fd$FDiv),
  Value  = fd$FDiv
) %>%
  left_join(meta, by = "Sample") %>%
  select(Group, Value) %>%
  mutate(Index = "FDiv")
write_csv(fdiv_long, "fdiv.csv")
# ────────────────── 4. 合并两指标 & 统计检验 ────────────────
plot_df <- bind_rows(
  fnb_long %>% select(Group, Value, Index),
  fdiv_long
)

# 统计表
stat_table <- plot_df %>%
  group_by(Index) %>%
  summarise(
    P_value = wilcox.test(Value ~ Group)$p.value,
    .groups = "drop"
  )
print(stat_table)

# ────────────────── 5. 可视化 ───────────────────────────────
ggplot(plot_df, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.85, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
  facet_wrap(~ Index, scales = "free_y") +
  stat_compare_means(
    method = "wilcox.test",
    label  = "p.signif",
    label.y.npc = "top"
  ) +
  scale_fill_manual(values = c(IT = "#E64B35", ST = "#4DBBD5")) +scale_color_manual(values = c(IT = "#E64B35", ST = "#4DBBD5")) +
  labs(
    x = "",
    y = "Index value",
    title = "IT vs ST : Functional Niche Breadth (FNB) & Functional Divergence (FDiv)"
  ) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "bold")
  )
















# 加载库
library(tidyverse)
library(vegan)

# 1. 读取数据
kegg_data <- read.csv("CSS_relative_KEGG_bases_lab.csv", row.names = "KO", check.names = FALSE)
metadata <- read.csv("metadata.csv")

# 2. 按组拆分数据
it_samples <- metadata$Sample[metadata$Group == "IT"]
st_samples <- metadata$Sample[metadata$Group == "ST"]

it_kegg <- kegg_data[, it_samples]
st_kegg <- kegg_data[, st_samples]


sum(rowSums(it_kegg) == 0)
sum(rowSums(st_kegg) == 0)

# 3. 定义FNB计算函数（基于Levins' B）
calculate_fnb <- function(data_matrix) {
  apply(data_matrix, 1, function(x) {
    p <- x / sum(x)  # 包含所有样本（含零值）
    if (sum(x) == 0) return(NA)  # 仅排除全零模块
    1 / sum(p^2)
  })
}

# 4. 计算各组FNB并标准化
it_fnb <- calculate_fnb(it_kegg) / length(it_samples)  # 标准化除以样本数
st_fnb <- calculate_fnb(st_kegg) / length(st_samples)

# 合并结果
fnb_df <- data.frame(
  KO = rownames(kegg_data),
  IT_FNB = it_fnb,
  ST_FNB = st_fnb
) %>% drop_na()  # 移除全零模块

# 5. 组间比较（Wilcoxon检验）
wilcox_test <- wilcox.test(fnb_df$IT_FNB, fnb_df$ST_FNB, paired = FALSE)
print(paste("Wilcoxon p-value:", wilcox_test$p.value))

# 6. 可视化（箱线图）
fnb_long <- fnb_df %>% 
  pivot_longer(cols = -KO, names_to = "Group", values_to = "FNB")

ggplot(fnb_long, aes(x = Group, y = FNB, fill = Group)) +
  geom_boxplot() +
  labs(title = "Standardized Functional Niche Breadth (FNB)",
       y = "Standardized FNB (FNB/n)",
       x = "") +
  theme_minimal()

# 1. 合并所有样本
all_kegg <- kegg_data[, metadata$Sample]
all_groups <- metadata$Group

# 2. 定义置换检验函数
permute_fnb_diff <- function(data, groups, n_perm = 1000) {
  real_diff <- mean(calculate_fnb(data[, groups == "IT"]), na.rm = TRUE) - 
    mean(calculate_fnb(data[, groups == "ST"]), na.rm = TRUE)
  
  perm_diffs <- replicate(n_perm, {
    shuffled_groups <- sample(groups)
    mean(calculate_fnb(data[, shuffled_groups == "IT"]), na.rm = TRUE) - 
      mean(calculate_fnb(data[, shuffled_groups == "ST"]), na.rm = TRUE)
  })
  
  p_value <- sum(abs(perm_diffs) >= abs(real_diff)) / n_perm
  list(real_diff = real_diff, p_value = p_value, perm_diffs = perm_diffs)
}

# 3. 执行置换检验
set.seed(123)
perm_result <- permute_fnb_diff(all_kegg, all_groups)

# 输出结果
print(paste("Observed FNB difference (IT-ST):", perm_result$real_diff))
print(paste("Permutation p-value:", perm_result$p_value))

# 可视化置换分布
hist(perm_result$perm_diffs, main = "Permutation Distribution of FNB Differences",
     xlab = "IT-ST FNB Difference", breaks = 30)
abline(v = perm_result$real_diff, col = "red", lwd = 2)

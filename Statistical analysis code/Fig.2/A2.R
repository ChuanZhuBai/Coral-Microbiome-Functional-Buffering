cat("=== STEP 1: 初始化环境 ===\n")
rm(list = ls())
library(ggplot2)
library(dplyr)
library(purrr)
library(ggrepel)
dir()
# 加载数据
# 数据过滤（更严格）
abundance_base <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv") %>% 
  filter(Sample %in% colnames(abundance_base))
# 过滤低丰度菌（0.1%）
abundance_df <- abundance_base[rowSums(abundance_base) > 0.001, ]
prev <- rowMeans(abundance_base > 0)
#abundance_df <- abundance_base[prev >= 0.2, ]

# 数据验证
cat("\n数据维度验证:\n")
cat("物种表 - 样本数:", ncol(abundance_df), "物种数:", nrow(abundance_df), "\n")
cat("元数据 - 样本数:", nrow(metadata), "\n")

# 样本匹配
metadata <- metadata %>% filter(Sample %in% colnames(abundance_df))
cat("匹配后的有效样本数:", nrow(metadata), "\n")

# ----------------------
# 3. 计算标准化B值
# 3. 计算标准化B值（优化版）
# ----------------------
cat("\n=== STEP 3: 计算标准化B值（优化版） ===\n")

calculate_standardized_B <- function(abundance, group_samples, group_name) {
  B_values <- apply(abundance[, group_samples, drop = FALSE], 1, function(x) {
    if (sum(x) <= 0) return(NA_real_)   # << 新增
    # 关键改进1：包含零丰度样本（反映真实分布）
    p <- x / sum(x)  # 不过滤x > 0，保留零值
    
    # 计算原始B值
    B_raw <- 1 / sum(p^2)
    
    # 关键改进2：Hurlbert标准化（基于生境数n）
    n <- length(group_samples)
    B_standardized <- (B_raw - 1) / (n - 1)  # 标准化到[0,1]
    
    return(B_standardized)
  })
  
  result <- data.frame(
    species = names(B_values),
    B_standardized = unname(B_values),
    Group = group_name,
    stringsAsFactors = FALSE
  ) %>% 
    filter(!is.na(B_standardized)) %>% 
    arrange(B_standardized)
  
  cat(group_name, "组计算完成，有效物种数:", nrow(result), "\n")
  return(result)
}

# 分组样本（保持不变）
IT_samples <- metadata %>% filter(Group == "IT") %>% pull(Sample)
ST_samples <- metadata %>% filter(Group == "ST") %>% pull(Sample)

# 计算B值
B_IT <- calculate_standardized_B(abundance_df, IT_samples, "IT")
B_ST <- calculate_standardized_B(abundance_df, ST_samples, "ST")

# 合并数据
B_data <- bind_rows(B_IT, B_ST) %>% 
  mutate(Group = factor(Group, levels = c("IT", "ST")))

cat("\nB值统计摘要:\n")
print(B_data %>% group_by(Group) %>% 
        summarise(Median = median(B_standardized),
                  Mean = mean(B_standardized),
                  SD = sd(B_standardized)))


library(tidyverse)
library(irr)        # kappa2
library(broom)      # 整理统计检验结果

# ---- 1.1 阈值方案定义 ----
threshold_schemes <- list(
  q25_75 = list(lower = 0.25, upper = 0.75, type = "quantile"),
  q20_80 = list(lower = 0.20, upper = 0.80, type = "quantile"),
  q10_90 = list(lower = 0.10, upper = 0.90, type = "quantile"),
  tukey  = list(type = "tukey")  # Q1-1.5*IQR, Q3+1.5*IQR
)

# ---- 1.2 辅助函数：按方案打标签 ----
label_by_scheme <- function(B_df, scheme){
  B_all <- B_df$B_standardized
  if (scheme$type == "quantile") {
    lo <- quantile(B_all, probs = scheme$lower, na.rm = TRUE)
    hi <- quantile(B_all, probs = scheme$upper, na.rm = TRUE)
  } else if (scheme$type == "tukey") {
    Q1  <- quantile(B_all, 0.25, na.rm = TRUE)
    Q3  <- quantile(B_all, 0.75, na.rm = TRUE)
    IQRv <- Q3 - Q1
    lo <- max(0, Q1 - 1.5 * IQRv)
    hi <- min(1, Q3 + 1.5 * IQRv)
    # 若 IQR≈0 导致 lo≈hi，标注说明（后续会在 analyze 里处理跳过 Fisher）
    if (!is.finite(IQRv) || IQRv == 0) {
      message("[label_by_scheme] Tukey: IQR=0; thresholds collapse (lo==hi).")
    }
  } else stop("Unknown scheme type")
  
  lab <- dplyr::case_when(
    B_df$B_standardized <= lo ~ "Specialist",
    B_df$B_standardized >= hi ~ "Generalist",
    TRUE ~ "Intermediate"
  )
  
  B_df %>%
    dplyr::mutate(Label = factor(lab, levels = c("Specialist","Intermediate","Generalist")),
                  lo_thr = lo, hi_thr = hi)
}


# ---- 1.3 统计对比：IT vs ST 的 generalist 占比 & B 分布差异 ----
analyze_scheme <- function(B_df_labeled, scheme_name){
  
  # 1) Generalist 占比
  tab <- B_df_labeled %>%
    mutate(G_is1 = as.integer(Label == "Generalist")) %>%
    group_by(Group) %>%
    summarise(Generalist_n = sum(G_is1, na.rm = TRUE),
              Total = n(), .groups = "drop") %>%
    mutate(Generalist_prop = Generalist_n / Total)
  
  # 2) 构造 2×2 列联表（确保两行两列都存在，即使为 0）
  ct_tbl <- B_df_labeled %>%
    mutate(isG = factor(Label == "Generalist", levels = c(FALSE, TRUE))) %>%
    count(Group, isG) %>%
    tidyr::complete(Group = levels(B_df_labeled$Group),
                    isG   = factor(c(FALSE, TRUE), levels = c(FALSE, TRUE)),
                    fill  = list(n = 0)) %>%
    arrange(factor(Group, levels = levels(B_df_labeled$Group)))
  
  ct <- ct_tbl %>%
    tidyr::pivot_wider(names_from = isG, values_from = n) %>%
    select(`FALSE`, `TRUE`) %>% as.matrix()
  
  # 检查是否至少有一个 Generalist & 非 Generalist
  has_two_cols <- (ncol(ct) == 2L)
  has_var      <- all(colSums(ct) > 0)  # 两列均有计数才适合做 Fisher
  
  if (has_two_cols && has_var) {
    fisher_res <- fisher.test(ct)
    fisher_df <- tibble(
      scheme    = scheme_name,
      test      = "Fisher_generalist_prop",
      p_value   = fisher_res$p.value,
      odds_ratio= unname(fisher_res$estimate)
    )
  } else {
    fisher_df <- tibble(
      scheme    = scheme_name,
      test      = "Fisher_generalist_prop",
      p_value   = NA_real_,
      odds_ratio= NA_real_
    )
    message(sprintf("[analyze_scheme] '%s': no both categories for Fisher; returning NA.", scheme_name))
  }
  
  # 3) B 分布差异：K–S + Wilcoxon（ties 多时补充）
  B_IT <- B_df_labeled$B_standardized[B_df_labeled$Group == "IT"]
  B_ST <- B_df_labeled$B_standardized[B_df_labeled$Group == "ST"]
  
  ks_ok <- (length(B_IT) > 0 && length(B_ST) > 0)
  if (ks_ok) {
    ks_res <- suppressWarnings(ks.test(B_IT, B_ST))
    ks_df <- tibble(scheme = scheme_name, test = "KS_B_distribution",
                    p_value = ks_res$p.value, statistic = as.numeric(ks_res$statistic))
    wil_res <- suppressWarnings(wilcox.test(B_IT, B_ST, exact = FALSE))
    wil_df <- tibble(scheme = scheme_name, test = "Wilcoxon_B_distribution",
                     p_value = wil_res$p.value, statistic = as.numeric(wil_res$statistic))
  } else {
    ks_df  <- tibble(scheme = scheme_name, test = "KS_B_distribution",
                     p_value = NA_real_, statistic = NA_real_)
    wil_df <- tibble(scheme = scheme_name, test = "Wilcoxon_B_distribution",
                     p_value = NA_real_, statistic = NA_real_)
  }
  
  list(summary_props = tab,
       fisher = fisher_df,
       ks = ks_df,
       wilcox = wil_df,
       thresholds = B_df_labeled %>% distinct(lo_thr, hi_thr) %>% mutate(scheme = scheme_name))
}

# ---- 1.4 运行所有方案 ----
sens_results <- purrr::map2(threshold_schemes, names(threshold_schemes), ~{
  B_labeled <- label_by_scheme(B_data, .x)
  analyze_scheme(B_labeled, .y)
})

sens_prop <- bind_rows(map(sens_results, "summary_props"), .id="scheme_id")
sens_fish <- bind_rows(map(sens_results, "fisher"))
sens_ks   <- bind_rows(map(sens_results, "ks"))
sens_wil  <- bind_rows(map(sens_results, "wilcox"))
sens_thr  <- bind_rows(map(sens_results, "thresholds"))

print(sens_prop)
print(sens_thr)
print(sens_fish %>% arrange(p_value))
print(sens_ks   %>% arrange(p_value))
print(sens_wil  %>% arrange(p_value))





################################################
# ---- 2.1 计算占据度与 Shannon 广度（按组）----
safe_shannon <- function(x){
  p <- x / sum(x)
  p <- p[p > 0]
  if (length(p) == 0) return(NA_real_)
  -sum(p * log(p))
}

calc_occ_H <- function(mat, samples, group_name){
  sub <- mat[, samples, drop=FALSE]
  occ <- rowMeans(sub > 0)  # 占据度（0~1）
  H   <- apply(sub, 1, safe_shannon)  # 未标准化的 Shannon
  Hstd <- H / log(ncol(sub))          # 标准化到[0,1]
  tibble(species = rownames(sub),
         Group = group_name,
         occ = occ,
         H = H,
         Hstd = Hstd)
}

occH_IT <- calc_occ_H(abundance_df, IT_samples, "IT")
occH_ST <- calc_occ_H(abundance_df, ST_samples, "ST")
occH_all <- bind_rows(occH_IT, occH_ST)

# 与 B_data 合并
B_occH <- B_data %>%
  left_join(occH_all, by = c("species","Group"))

# ---- 2.2 基于占据度/标准化Shannon 打“广/专适”标签（25/75）----
label_by_quantile_vec <- function(v, lower=0.25, upper=0.75){
  lo <- quantile(v, lower, na.rm=TRUE)
  hi <- quantile(v, upper, na.rm=TRUE)
  lab <- case_when(
    v <= lo ~ "Specialist",
    v >= hi ~ "Generalist",
    TRUE ~ "Intermediate"
  )
  list(label = factor(lab, levels=c("Specialist","Intermediate","Generalist")),
       lo=lo, hi=hi)
}

add_alt_labels <- function(df){
  # 分组（IT/ST）分别求阈值，避免样本量差异偏倚
  df %>% group_by(Group) %>% group_modify(~{
    occ_lab <- label_by_quantile_vec(.x$occ)
    H_lab   <- label_by_quantile_vec(.x$Hstd)
    .x %>% mutate(
      Label_B = NA_character_,  # 占位，稍后回填
      Label_occ = occ_lab$label,
      occ_lo = occ_lab$lo, occ_hi = occ_lab$hi,
      Label_H   = H_lab$label,
      H_lo = H_lab$lo, H_hi = H_lab$hi
    )
  }) %>% ungroup()
}

# 1) 不要在 add_alt_labels() 里创建 Label_B 占位列（如果创建了，删掉）
# 这里直接把已生成的 B_occH_labeled 里的 Label_B（若存在）先移除
B_occH_labeled <- B_occH_labeled %>%
  dplyr::select(-dplyr::any_of("Label_B"))

# 2) 从 25/75 方案生成基于 B 的标签，并命名为 Label_B
B_25_75 <- label_by_scheme(B_data, threshold_schemes$q25_75) %>%
  dplyr::select(species, Group, Label) %>%
  dplyr::rename(Label_B = Label)

# 3) 左连接，把 Label_B 干净并唯一地合并到 B_occH_labeled
B_occH_labeled <- B_occH_labeled %>%
  dplyr::left_join(B_25_75, by = c("species","Group"))

# 快速自检
dplyr::glimpse(B_occH_labeled)
# 确认现在有：species, Group, B_standardized, occ, Hstd, Label_occ, Label_H, Label_B 等

# ---- 2.4 一致性（Cohen’s kappa）与相关性（Spearman）----
compute_kappa_cor <- function(df_grp){
  # 仅保留三类标签一致的水平
  k_occ <- kappa2(df_grp %>% select(Label_B, Label_occ))
  k_H   <- kappa2(df_grp %>% select(Label_B, Label_H))
  cor_occ <- cor(df_grp$B_standardized, df_grp$occ, method="spearman", use="complete.obs")
  cor_H   <- cor(df_grp$B_standardized, df_grp$Hstd, method="spearman", use="complete.obs")
  tibble(
    Group = unique(df_grp$Group),
    kappa_B_vs_occ = unname(k_occ$value),
    kappa_B_vs_H   = unname(k_H$value),
    spearman_B_occ = cor_occ,
    spearman_B_H   = cor_H
  )
}

kappa_summary <- B_occH_labeled %>% group_by(Group) %>% group_modify(~compute_kappa_cor(.x)) %>% ungroup()
print(kappa_summary)

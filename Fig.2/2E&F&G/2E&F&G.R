# =============================================================================
# Functional Redundancy Index (FRI) and its drivers, IT vs ST

# =============================================================================
library(tidyverse); library(ggpubr); library(patchwork); library(rstatix)

hab_cols <- c(Intertidal = "#E64B35FF", Subtidal = "#4DBBD5FF")
to_hab   <- function(g) factor(c(IT="Intertidal", ST="Subtidal")[g], levels = c("Intertidal","Subtidal"))

read_auto <- function(file) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}

stress_paths <- c(
  "Chaperones and folding catalysts", 
  "Base excision repair", "Mismatch repair", "Homologous recombination",
  "Oxidative phosphorylation", "Pentose phosphate pathway",
  "Glutathione metabolism", "Biosynthesis of unsaturated fatty acids",
   "Carotenoid biosynthesis")

# ---- load -------------------------------------------------------------------
genus_abun <- read_auto("genus.csv"); rownames(genus_abun) <- genus_abun[[1]]; genus_abun[[1]] <- NULL
genus_abun <- as.matrix(genus_abun); storage.mode(genus_abun) <- "numeric"
taxon_kegg <- read_auto("extracted_taxon_kegg_bacteria.csv")     # bacteria/archaea only: genus, KO
kegg_paths <- read_auto("KEGG_paths.csv")
metadata   <- read_auto("metadata.csv")
samp_col   <- intersect(c("Sample","SampleID"), names(metadata))[1]
metadata$Sample <- as.character(metadata[[samp_col]])
niche <- read_auto("Figure_2D_niche_breadth_values.csv")        # genus, B, Group, Class
meta_sg <- metadata[, c("Sample","Group")]

genus_abun <- genus_abun[, colnames(genus_abun) %in% metadata$Sample, drop = FALSE]

# ---- KO -> KEGG Level-3 pathway map -----------------------------------------
pcol <- setdiff(names(kegg_paths), "KO")[1]
kegg_long <- kegg_paths %>%
  transmute(KO, paths = gsub("\t", " ", .data[[pcol]])) %>%
  separate_rows(paths, sep = "\\|") %>% mutate(paths = str_trim(paths)) %>%
  filter(paths != "") %>%
  mutate(level3 = vapply(str_split(paths, ";"), function(p) {
    p <- str_trim(p[p != ""]); if (length(p) >= 3) p[3] else p[length(p)] }, character(1))) %>%
  filter(level3 != "") %>% distinct(KO, path = level3)
all_paths <- unique(kegg_long$path)

miss <- setdiff(stress_paths, all_paths)
if (length(miss)) message("WARNING: stress categories not in mapping: ", paste(miss, collapse="; "))

# ---- detected genera (relative abundance) per sample ------------------------
present_long <- as.data.frame(sweep(genus_abun, 2, colSums(genus_abun), "/")) %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "ab") %>%
  filter(ab > 0)

# =============================================================================
# KO-level redundancy — BOTH metrics computed ONCE, then pathway-averaged
#   KO_count = # distinct genera encoding the KO
#   KO_eff   = Hill q=1 effective # genera (abundance-weighted)
# =============================================================================
ko_red <- present_long %>%
  inner_join(taxon_kegg, by = "genus", relationship = "many-to-many") %>%
  group_by(Sample, KO) %>%
  summarise(KO_count = n_distinct(genus),
            KO_eff   = { p <- ab/sum(ab); exp(-sum(p*log(p))) }, .groups = "drop")

# pathway-level FRI for both metrics (long over Metric)
fri_path <- ko_red %>%
  inner_join(kegg_long, by = "KO", relationship = "many-to-many") %>%
  group_by(Sample, path) %>%
  summarise(FRI_count = mean(KO_count), FRI_eff = mean(KO_eff), .groups = "drop")

# whole-community FRI per sample (zero-filled over ALL pathways), both metrics
global_fri <- expand_grid(Sample = unique(metadata$Sample), path = all_paths) %>%
  left_join(fri_path, by = c("Sample","path")) %>%
  mutate(across(c(FRI_count, FRI_eff), ~replace_na(., 0))) %>%
  group_by(Sample) %>%
  summarise(Global_FRI = mean(FRI_count), Global_FRI_eff = mean(FRI_eff), .groups = "drop") %>%
  left_join(meta_sg, by = "Sample") %>% mutate(Habitat = to_hab(Group))

# =============================================================================
# 2E — whole-community FRI (main metric), IT vs ST
# =============================================================================
message(sprintf("2H FRI_count: IT median=%.2f ST median=%.2f | Wilcoxon p=%.3g",
                median(global_fri$Global_FRI[global_fri$Group=="IT"]),
                median(global_fri$Global_FRI[global_fri$Group=="ST"]),
                wilcox.test(Global_FRI ~ Group, data = global_fri)$p.value))
message(sprintf("[Check] FRI_eff (Hill q=1): IT median=%.3f ST median=%.3f | Wilcoxon p=%.3g",
                median(global_fri$Global_FRI_eff[global_fri$Group=="IT"]),
                median(global_fri$Global_FRI_eff[global_fri$Group=="ST"]),
                wilcox.test(Global_FRI_eff ~ Group, data = global_fri)$p.value))

# anti-cherry-picking: of ALL pathways, how many show IT>ST (main metric)?
path_dir <- fri_path %>% left_join(meta_sg, by = "Sample") %>% group_by(path) %>%
  summarise(medIT = median(FRI_count[Group=="IT"]), medST = median(FRI_count[Group=="ST"]),
            p = tryCatch(wilcox.test(FRI_count ~ Group, exact = FALSE)$p.value,
                         error = function(e) NA_real_), .groups = "drop")
message(sprintf("Across ALL %d pathways: %d IT>ST median (%d of them p<0.05)",
                nrow(path_dir), sum(path_dir$medIT > path_dir$medST, na.rm=TRUE),
                sum(path_dir$medIT > path_dir$medST & path_dir$p < 0.05, na.rm=TRUE)))

p_2h <- ggplot(global_fri, aes(Habitat, Global_FRI, fill = Habitat)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.5) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("Intertidal","Subtidal")),
                     label = "p.signif", size = 5) +
  scale_fill_manual(values = hab_cols) +
  labs(x = NULL, y = "Mean Functional Redundancy Index") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text = element_text(color = "black"))

# =============================================================================
# 2F — 9 a priori stress categories (main metric), Wilcoxon + BH
# =============================================================================
fri_stress <- expand_grid(Sample = unique(metadata$Sample), path = stress_paths) %>%
  left_join(fri_path, by = c("Sample","path")) %>%
  mutate(across(c(FRI_count, FRI_eff), ~replace_na(., 0))) %>%
  left_join(meta_sg, by = "Sample")

stat_2i <- fri_stress %>% group_by(path) %>%
  summarise(p = if (sd(FRI_count) == 0) NA_real_
            else tryCatch(wilcox.test(FRI_count ~ Group, exact = FALSE)$p.value,
                          error = function(e) NA_real_), .groups = "drop") %>%
  mutate(p.adj = p.adjust(p, "BH"),
         p.adj.signif = case_when(is.na(p.adj)~"n.a.", p.adj<1e-4~"****", p.adj<1e-3~"***",
                                  p.adj<1e-2~"**", p.adj<0.05~"*", TRUE~"ns"))

summ_2i <- fri_stress %>% group_by(path, Group) %>%
  summarise(mean = mean(FRI_count), se = sd(FRI_count)/sqrt(n()), .groups = "drop")
ord <- summ_2i %>% group_by(path) %>% summarise(m = mean(mean)) %>% arrange(m) %>% pull(path)
summ_2i$path    <- factor(summ_2i$path, levels = ord)
summ_2i$Habitat <- factor(to_hab(summ_2i$Group), levels = c("Subtidal", "Intertidal")) 

label_pos <- summ_2i %>% group_by(path) %>%
  summarise(label_y = max(mean + se) + max(mean + se) * 0.04, .groups = "drop")
stat_2i_plot <- stat_2i %>% left_join(label_pos, by = "path") %>%
  mutate(path = factor(path, levels = ord))
xlim_max <- max(label_pos$label_y) * 1.10

p_2i <- ggplot(summ_2i, aes(path, mean, fill = Habitat)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se, group = Habitat),
                position = position_dodge(width = 0.75), width = 0.25) +
  geom_text(data = stat_2i_plot, aes(x = path, y = label_y, label = p.adj.signif),
            inherit.aes = FALSE, size = 3.8, hjust = 0) +
  scale_fill_manual(values = hab_cols, breaks = c("Intertidal","Subtidal"),
                    limits = c("Intertidal","Subtidal")) +
  scale_y_continuous(limits = c(0, xlim_max), expand = expansion(mult = c(0, 0.02))) +
  coord_flip() +
  labs(x = NULL, y = "Mean Functional Redundancy Index") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", legend.title = element_blank(),
        axis.text = element_text(color = "black"))

# =============================================================================
# 2G — pooled regression on global generalist richness (+ within-habitat check)
# =============================================================================

if (!"Class" %in% names(niche))
  stop("niche file lacks 'Class'; regenerate Figure_2D_niche_breadth_values.csv")
gen_global  <- niche %>% filter(Class == "Generalist") %>% pull(genus) %>% unique()
B_threshold <- as.numeric(quantile(niche$B, 0.75, na.rm = TRUE))   # reported only
message(sprintf("Generalist genera from Fig 2D classification: %d (75th-pct pooled B = %.3f)",
                length(gen_global), B_threshold))

gen_rich <- present_long %>% filter(genus %in% gen_global) %>%
  group_by(Sample) %>% summarise(Generalist_Richness = n_distinct(genus), .groups = "drop")
reg_df <- global_fri %>% left_join(gen_rich, by = "Sample") %>%
  mutate(Generalist_Richness = replace_na(Generalist_Richness, 0))

fit_all <- summary(lm(Global_FRI ~ Generalist_Richness, data = reg_df))
message(sprintf("2J pooled regression: R2=%.3f, p=%.3g", fit_all$r.squared, fit_all$coefficients[2,4]))
for (g in c("IT","ST")) {
  d <- reg_df[reg_df$Group == g, ]
  ct <- cor.test(d$Generalist_Richness, d$Global_FRI)
  message(sprintf("   within %s: r=%.3f, p=%.3g (n=%d)", g, ct$estimate, ct$p.value, nrow(d)))
}

p_2j <- ggplot(reg_df, aes(Generalist_Richness, Global_FRI)) +
  geom_smooth(method = "lm", color = "black", fill = "grey85", alpha = 0.6, linewidth = 1) +
  geom_point(aes(fill = Habitat), shape = 21, size = 4, stroke = 0.4, alpha = 0.85) +
  annotate("text", x = min(reg_df$Generalist_Richness), y = max(reg_df$Global_FRI),
           label = sprintf("R² = %.2f, p = %.3g", fit_all$r.squared, fit_all$coefficients[2,4]),
           hjust = 0, vjust = 1, size = 4, fontface = "bold") +
  scale_fill_manual(values = hab_cols) +
  labs(x = "Richness of generalist genera (per sample)", y = "Mean Functional Redundancy Index") +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"), legend.position = "top", legend.title = element_blank())

# =============================================================================
# Table S7 — 9 stress categories x {FRI_count, FRI_eff}: Mean+-SE + Wilcoxon
# =============================================================================
tabS5 <- fri_stress %>%
  pivot_longer(c(FRI_count, FRI_eff), names_to = "Metric", values_to = "value") %>%
  group_by(path, Metric) %>%
  group_modify(~{
    d <- .x
    ms <- d %>% group_by(Group) %>% summarise(m = mean(value), se = sd(value)/sqrt(n()), .groups="drop")
    itv <- ms[ms$Group=="IT",]; stv <- ms[ms$Group=="ST",]
    if (sd(d$value) == 0) {
      tibble(IT_mean_SE = sprintf("%.3f +/- %.3f", itv$m, itv$se),
             ST_mean_SE = sprintf("%.3f +/- %.3f", stv$m, stv$se),
             statistic = NA, p = NA, effsize = NA)
    } else {
      wt <- wilcox.test(value ~ Group, data = d, exact = FALSE)
      es <- d %>% wilcox_effsize(value ~ Group)
      tibble(IT_mean_SE = sprintf("%.3f +/- %.3f", itv$m, itv$se),
             ST_mean_SE = sprintf("%.3f +/- %.3f", stv$m, stv$se),
             statistic = as.numeric(wt$statistic), p = wt$p.value, effsize = es$effsize)
    }
  }) %>% ungroup() %>%
  mutate(p.adj = p.adjust(p, "BH"),
         p.adj.signif = case_when(is.na(p.adj)~"n.a.", p.adj<1e-4~"****", p.adj<1e-3~"***",
                                  p.adj<1e-2~"**", p.adj<0.05~"*", TRUE~"ns"),
         magnitude = case_when(is.na(effsize)~NA_character_, effsize<0.1~"negligible",
                               effsize<0.3~"small", effsize<0.5~"moderate", TRUE~"large")) %>%
  arrange(Metric, p.adj)

# =============================================================================
# SAVE — figures + essential source data only
# =============================================================================
ggsave("Figure_2H_global_FRI.pdf", p_2h, width = 4, height = 6)
ggsave("Figure_2I_stress_FRI.pdf",  p_2i, width = 10, height = 6)
ggsave("Figure_2J_regression.pdf",  p_2j, width = 8,   height = 7)

write.csv(reg_df %>% select(Sample, Group, Habitat, Generalist_Richness,
                            Global_FRI, Global_FRI_eff),
          "SourceData_2HJ_per_sample.csv", row.names = FALSE)
write.csv(tabS5, "TableS5_FRI_stress_stats.csv", row.names = FALSE)

print(p_2h); print(p_2i); print(p_2j)
print(tabS5 %>% filter(Metric == "FRI_count"))


# =============================================================================
# =============================================================================
path_dir_full <- fri_path %>%
  left_join(meta_sg, by = "Sample") %>%
  group_by(path) %>%
  summarise(
    n_IT = sum(Group == "IT"), n_ST = sum(Group == "ST"),   
    
    median_IT = median(FRI_count[Group == "IT"]),
    median_ST = median(FRI_count[Group == "ST"]),
    p = tryCatch(wilcox.test(FRI_count ~ Group, exact = FALSE)$p.value,
                 error = function(e) NA_real_),
    
    median_IT_eff = median(FRI_eff[Group == "IT"]),
    median_ST_eff = median(FRI_eff[Group == "ST"]),
    p_eff = tryCatch(wilcox.test(FRI_eff ~ Group, exact = FALSE)$p.value,
                     error = function(e) NA_real_),
    .groups = "drop") %>%
  mutate(
    diff_median = median_IT - median_ST,
    direction   = case_when(median_IT > median_ST ~ "IT>ST",
                            median_IT < median_ST ~ "ST>IT", TRUE ~ "tie"),
    p.adj       = p.adjust(p, "BH"),
    p_eff.adj   = p.adjust(p_eff, "BH"),
    signif_BH   = case_when(is.na(p.adj) ~ "n.a.", p.adj < 1e-4 ~ "****",
                            p.adj < 1e-3 ~ "***", p.adj < 1e-2 ~ "**",
                            p.adj < 0.05 ~ "*", TRUE ~ "ns")) %>%
  arrange(p.adj, desc(diff_median))

write.csv(path_dir_full, "TableSx_all_pathways_IT_vs_ST.csv", row.names = FALSE)


message(sprintf("Exported %d pathways | IT>ST: %d | raw p<0.05: %d | BH q<0.05: %d",
                nrow(path_dir_full),
                sum(path_dir_full$direction == "IT>ST"),
                sum(path_dir_full$direction == "IT>ST" & path_dir_full$p < 0.05, na.rm = TRUE),
                sum(path_dir_full$direction == "IT>ST" & path_dir_full$p.adj < 0.05, na.rm = TRUE)))
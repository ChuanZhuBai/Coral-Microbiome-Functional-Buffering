# =============================================================================
# functional retention rate (FRR) under heat, intertidal vs subtidal
# =============================================================================
library(tidyverse); library(rstatix)

LOG2FC_RETAIN <- 1
MIN_KO_PATH   <- 5

EXCLUDE_L1 <- c("Human Diseases", "Organismal Systems")
set.seed(1)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.table(file, header = TRUE, sep = sep, row.names = row.names,
             check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
}
canonical_id <- function(x) {
  x %>%
    gsub("^(IT|ST)[-]?C([0-9]+)$", "\\1-C\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)C$",     "\\1-C\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)A$",     "\\1-C\\2", .) %>%
    gsub("^(IT|ST)[-]?H([0-9]+)$", "\\1-H\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)H$",     "\\1-H\\2", .)
}

read_paths_L2 <- function(file) {
  raw <- read_auto(file); pcol <- setdiff(names(raw), "KO")[1]
  raw %>%
    transmute(KO = str_extract(KO, "K\\d+"), paths = .data[[pcol]]) %>%
    separate_rows(paths, sep = "\\s*\\|\\s*") %>%
    mutate(parts = str_split(paths, ";\\s*"),
           L1 = map_chr(parts, ~ if (length(.x) >= 1) str_trim(.x[1]) else NA_character_),
           L2 = map_chr(parts, ~ if (length(.x) >= 2) str_trim(.x[2]) else NA_character_)) %>%
    filter(!is.na(L2), L2 != "") %>%
    filter(!str_detect(L1, regex("Brite|Not Included", ignore_case = TRUE))) %>%
    filter(!L1 %in% EXCLUDE_L1) %>%               
    distinct(KO, L1, L2)
}


counts <- as.matrix(read_auto("ko_abund.csv", row.names = 1)); storage.mode(counts) <- "numeric"
colnames(counts) <- canonical_id(colnames(counts)); counts <- round(counts)
meta <- read_auto("metadata.csv"); colnames(meta)[1:2] <- c("Sample","Group")
meta$Sample <- canonical_id(meta$Sample)
meta <- meta[meta$Group %in% c("IT-C","IT-H","ST-C","ST-H"), ]

exp_samples <- intersect(meta$Sample, colnames(counts))
message(sprintf("Experimental samples used: %d / %d", length(exp_samples), nrow(meta)))
counts <- counts[, exp_samples, drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
rel <- sweep(counts, 2, colSums(counts), "/")

kegg_long <- read_paths_L2("KEGG_paths.csv")
message(sprintf("KEGG Level-2 categories retained (bacteria-relevant): %d",
                n_distinct(kegg_long$L2)))

# ---- pair Heat with Control -------------------------------------------------
heated <- meta$Sample[grepl("-H[0-9]+$", meta$Sample)]
pairs <- tibble(heat = heated, ctrl = sub("-H", "-C", heated)) %>%
  filter(heat %in% colnames(rel), ctrl %in% colnames(rel)) %>%
  mutate(Habitat = ifelse(grepl("^IT", heat), "Intertidal", "Subtidal"))
message(sprintf("Pairs: %d (IT=%d, ST=%d)", nrow(pairs),
                sum(pairs$Habitat=="Intertidal"), sum(pairs$Habitat=="Subtidal")))

pc <- apply(rel, 1, function(x){ nz <- x[x>0]; if(length(nz)) min(nz)/2 else 1e-9 })

# per-KO retained flag for each pair
retain_by_pair <- pmap(pairs, function(heat, ctrl, Habitat){
  l2fc <- log2((rel[, heat] + pc) / (rel[, ctrl] + pc))
  tibble(KO = rownames(rel), retained = abs(l2fc) <= LOG2FC_RETAIN, Sample = heat, Habitat = Habitat)
})

# ============================================================================
# (1) PRIMARY: OVERALL retention per pair (all bacteria-relevant KOs), IT vs ST
# ============================================================================
bact_kos <- unique(kegg_long$KO)
overall <- map_dfr(retain_by_pair, function(d){
  dd <- d %>% filter(KO %in% bact_kos)
  tibble(Sample = dd$Sample[1], Habitat = dd$Habitat[1],
         Overall_FRR = mean(dd$retained) * 100)
})
overall$Habitat <- factor(overall$Habitat, levels = c("Intertidal","Subtidal"))
wt_overall <- overall %>% wilcox_test(Overall_FRR ~ Habitat)
eff_overall <- overall %>% wilcox_effsize(Overall_FRR ~ Habitat)
message(sprintf("[OVERALL] IT mean=%.1f%% ST mean=%.1f%% | Wilcoxon p=%.3g | r=%.3f (%s)",
                mean(overall$Overall_FRR[overall$Habitat=="Intertidal"]),
                mean(overall$Overall_FRR[overall$Habitat=="Subtidal"]),
                wt_overall$p, eff_overall$effsize, eff_overall$magnitude))
write.csv(overall, "SourceData_3E_overall_FRR.csv", row.names = FALSE)

# ============================================================================
# (2) per Level-2 pathway FRR + direction-consistency sign test
# ============================================================================
frr_long <- map_dfr(retain_by_pair, function(d){
  d %>% inner_join(kegg_long, by = "KO") %>%
    group_by(L1, L2) %>%
    summarise(n_ko = n_distinct(KO), FRR = mean(retained) * 100, .groups = "drop") %>%
    filter(n_ko >= MIN_KO_PATH) %>%
    mutate(Sample = d$Sample[1], Habitat = d$Habitat[1])
})
stat_path <- frr_long %>%
  group_by(L1, L2) %>% filter(n_distinct(Habitat) == 2) %>%
  summarise(IT_mean = mean(FRR[Habitat=="Intertidal"]),
            ST_mean = mean(FRR[Habitat=="Subtidal"]),
            IT_se = sd(FRR[Habitat=="Intertidal"])/sqrt(sum(Habitat=="Intertidal")),
            ST_se = sd(FRR[Habitat=="Subtidal"])/sqrt(sum(Habitat=="Subtidal")),
            p = tryCatch(wilcox.test(FRR ~ Habitat, exact=FALSE)$p.value, error=function(e) NA_real_),
            .groups="drop") %>%
  mutate(p.adj = p.adjust(p, "BH"),
         signif = case_when(is.na(p.adj)~"n.a.", p.adj<0.001~"***", p.adj<0.01~"**",
                            p.adj<0.05~"*", TRUE~"ns"),
         direction = ifelse(IT_mean > ST_mean, "IT>ST", "ST>=IT"))

n_tot <- nrow(stat_path); n_it <- sum(stat_path$IT_mean > stat_path$ST_mean)
sign_p <- binom.test(n_it, n_tot, p = 0.5, alternative = "two.sided")$p.value
message(sprintf("[SIGN TEST] %d of %d Level-2 pathways IT>ST | binomial p=%.3g",
                n_it, n_tot, sign_p))

message(sprintf("Distinct L2 names: %d ; rows in stat_path: %d", n_distinct(stat_path$L2), nrow(stat_path)))
if (n_distinct(stat_path$L2) != nrow(stat_path))
  message("NOTE: some L2 names appear under >1 L1 super-class; each (L1,L2) pair is a separate row.")
print(stat_path %>% arrange(IT_mean - ST_mean) %>%
        select(L1, L2, IT_mean, ST_mean, p.adj, signif, direction), n = 100)
write.csv(stat_path, "Figure_3E_FRR_perpathway_stats.csv", row.names = FALSE)
write.csv(frr_long,  "SourceData_3E_FRR_per_pair.csv", row.names = FALSE)


library(ggpubr)
hab_cols <- c(Intertidal = "#E64B35FF", Subtidal = "#4DBBD5FF")
wt_pl <- overall %>% wilcox_test(Overall_FRR ~ Habitat) %>% add_xy_position()
p_overall <- ggplot(overall, aes(Habitat, Overall_FRR, fill = Habitat)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +
  stat_pvalue_manual(wt_pl, label = "p = {p}", tip.length = 0.01) +
  scale_fill_manual(values = hab_cols) +
  labs(x = NULL, y = "Overall functional retention rate (%)",
       subtitle = sprintf("Per-pair Wilcoxon p = %.3f, r = %.2f (%s); direction: %d/%d pathways IT>ST, sign test p = %.2g",
                          wt_overall$p, eff_overall$effsize, eff_overall$magnitude, n_it, n_tot, sign_p),
       title = "Overall functional retention under heat") +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face="bold", hjust=0.5),
        plot.subtitle = element_text(hjust=0.5, size=8), axis.text = element_text(color="black"))

# ---- DIRECTION-CONSISTENCY plot: dumbbell (ALL pathways, IT vs ST) -----------
dumbbell_df <- stat_path %>%
  mutate(lean = ifelse(IT_mean > ST_mean, "IT higher", "ST higher"),
         diff = IT_mean - ST_mean) %>%
  arrange(diff) %>%
  mutate(L2 = factor(L2, levels = L2))            
stopifnot(nrow(dumbbell_df) == n_tot)             

p_dumbbell <- ggplot(dumbbell_df, aes(y = L2)) +
  geom_segment(aes(x = ST_mean, xend = IT_mean, yend = L2, color = lean), linewidth = 1) +
  geom_point(aes(x = ST_mean, shape = "Subtidal"), color = hab_cols["Subtidal"], size = 2.8) +
  geom_point(aes(x = IT_mean, shape = "Intertidal"), color = hab_cols["Intertidal"], size = 2.8) +
  scale_color_manual(name = "Direction",
                     values = c("IT higher" = "#E64B3577", "ST higher" = "#4DBBD5")) +
  scale_shape_manual(name = "Habitat", values = c("Intertidal" = 16, "Subtidal" = 16)) +
  guides(shape = guide_legend(override.aes = list(
    color = c(hab_cols["Intertidal"], hab_cols["Subtidal"])))) +
  labs(x = "Functional retention rate (%)", y = NULL,
       title = "Direction consistency across functional categories",
       subtitle = sprintf("%d of %d bacteria-relevant Level-2 categories show IT > ST (sign test p = %.2g)",
                          n_it, n_tot, sign_p)) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face="bold", hjust=0.5, size=12),
        plot.subtitle = element_text(hjust=0.5, size=9),
        axis.text = element_text(color="black"),
        axis.text.y = element_text(size = 8),
        legend.position = "top")
ggsave("Figure_3E_direction_dumbbell.pdf", p_dumbbell, width = 18, height = 10)


plot_df <- frr_long %>%
  mutate(key = paste(L1, L2, sep = " | ")) %>%
  group_by(L1, L2, key, Habitat) %>%
  summarise(mean = mean(FRR), se = sd(FRR)/sqrt(n()), .groups = "drop") %>%
  mutate(Habitat = factor(Habitat, levels = c("Intertidal","Subtidal")))
ord <- plot_df %>% group_by(key) %>% summarise(m = mean(mean), .groups="drop") %>%
  arrange(m) %>% pull(key)
plot_df$key <- factor(plot_df$key, levels = ord)
lab_df <- stat_path %>% mutate(key = paste(L1, L2, sep = " | "),
                               key = factor(key, levels = ord),
                               y = pmax(IT_mean + IT_se, ST_mean + ST_se, na.rm = TRUE))
stopifnot(nlevels(plot_df$key) == n_tot)        # guard: all pathways plotted
p_path <- ggplot(plot_df, aes(key, mean, fill = Habitat)) +
  geom_col(position = position_dodge(width=0.75), width=0.65) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), position=position_dodge(width=0.75), width=0.25) +
  geom_text(data = lab_df, aes(x=key, y=y+2.5, label=signif), inherit.aes=FALSE, size=3, hjust=0) +
  scale_fill_manual(values = hab_cols) +
  scale_x_discrete(labels = function(k) sub("^.* \\| ", "", k)) +   # show only the L2 name
  facet_grid(L1 ~ ., scales="free_y", space="free_y") +
  coord_flip(ylim = c(0, 105)) +
  labs(x = NULL, y = "Functional retention rate (%)  (|log2FC| <= 1, Heat vs Control)",
       title = "Per-pathway functional retention (bacteria-relevant Level-2)") +
  theme_classic(base_size = 10) +
  theme(legend.position="top", legend.title=element_blank(), axis.text=element_text(color="black"),
        strip.text.y = element_text(angle=0, face="bold"),
        strip.background = element_rect(fill="#f0f0f0", color=NA),
        plot.title = element_text(face="bold", hjust=0.5, size=12))

ggsave("Figure_3E_overall_retention.pdf", p_overall, width = 4, height = 6)
ggsave("Figure_3E_perpathway_retention.pdf", p_path, width = 15, height = 8)
print(p_overall); print(p_dumbbell); print(p_path)

# =============================================================================
# THRESHOLD SENSITIVITY ANALYSIS
# The "retained" criterion (|log2FC| <= 1, i.e. within 2-fold) is a magnitude
# threshold. To show the conclusion (IT retains more function than ST) does not
# depend on this specific cutoff, the whole analysis is repeated at 1.5-fold
# (|log2FC| <= 0.585), 2-fold (<= 1, main), and 3-fold (<= 1.585). For each
# threshold we report: overall FRR per habitat, the per-pair Wilcoxon, the
# effect size, and the across-pathway direction-consistency sign test.
# =============================================================================
sensitivity_at <- function(cut) {
  # per-KO retained flag per pair at this threshold
  rbp <- pmap(pairs, function(heat, ctrl, Habitat){
    l2fc <- log2((rel[, heat] + pc) / (rel[, ctrl] + pc))
    tibble(KO = rownames(rel), retained = abs(l2fc) <= cut,
           Sample = heat, Habitat = Habitat)
  })
  # (a) overall FRR per pair (bacteria-relevant KOs), IT vs ST
  ov <- map_dfr(rbp, function(d){
    dd <- d %>% filter(KO %in% bact_kos)
    tibble(Habitat = dd$Habitat[1], Overall_FRR = mean(dd$retained) * 100)
  })
  ov$Habitat <- factor(ov$Habitat, levels = c("Intertidal","Subtidal"))
  wt  <- ov %>% wilcox_test(Overall_FRR ~ Habitat)
  eff <- ov %>% wilcox_effsize(Overall_FRR ~ Habitat)
  # (b) per-pathway FRR -> direction-consistency sign test
  fl <- map_dfr(rbp, function(d){
    d %>% inner_join(kegg_long, by = "KO") %>%
      group_by(L1, L2) %>%
      summarise(n_ko = n_distinct(KO), FRR = mean(retained) * 100, .groups = "drop") %>%
      filter(n_ko >= MIN_KO_PATH) %>% mutate(Habitat = d$Habitat[1])
  })
  sp <- fl %>% group_by(L1, L2) %>% filter(n_distinct(Habitat) == 2) %>%
    summarise(IT = mean(FRR[Habitat=="Intertidal"]),
              ST = mean(FRR[Habitat=="Subtidal"]), .groups = "drop")
  n_t <- nrow(sp); n_i <- sum(sp$IT > sp$ST)
  signp <- binom.test(n_i, n_t, 0.5)$p.value
  tibble(fold = sprintf("%.1f-fold", 2^cut),
         log2FC_cut = cut,
         IT_overall = mean(ov$Overall_FRR[ov$Habitat=="Intertidal"]),
         ST_overall = mean(ov$Overall_FRR[ov$Habitat=="Subtidal"]),
         overall_p = wt$p, overall_r = eff$effsize, magnitude = eff$magnitude,
         pathways_IT_gt_ST = sprintf("%d/%d", n_i, n_t), sign_test_p = signp)
}

sens <- bind_rows(lapply(c(0.585, 1, 1.585), sensitivity_at))
message("\n===== Threshold sensitivity analysis =====")
print(as.data.frame(sens), digits = 3)
write.csv(sens, "Figure_3E_threshold_sensitivity.csv", row.names = FALSE)












plot_df <- stat_path %>%
  filter(L2 != "Cellular community - eukaryotes") %>% 
  mutate(lean = ifelse(IT_mean > ST_mean, "IT higher", "ST higher"),
         diff = IT_mean - ST_mean,
         abs_diff = abs(diff)) %>%      
  arrange(diff) %>%                     
  mutate(L2 = factor(L2, levels = L2))            


n_total_filtered <- nrow(plot_df)
n_it_filtered <- sum(plot_df$diff > 0)


p_bar_sameside <- ggplot(plot_df, aes(x = abs_diff, y = L2, fill = lean)) +
  geom_col(width = 0.7, alpha = 0.9) +
  
  geom_text(aes(label = sprintf("%.1f%%", abs_diff)),
            hjust = -0.15, size = 3.5, color = "black", fontface = "bold") +
  
  scale_fill_manual(values = c("IT higher" = "#E64B35", "ST higher" = "#4DBBD5")) +
  
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(x = "Absolute Retention Rate Difference (|IT - ST|, %)", 
       y = NULL,
       title = "Relative Functional Retention Difference",
       subtitle = sprintf("%d of %d prokaryotic Level-2 categories show IT > ST", 
                          n_it_filtered, n_total_filtered)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"), 
        axis.text.x = element_text(color = "black"),
        axis.line.y = element_line(color = "black", linewidth = 0.8), 
        axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"))


ggsave("Figure_3E_SameSide_BarPlot.pdf", p_bar_sameside, width = 12, height = 6)
print(p_bar_sameside)

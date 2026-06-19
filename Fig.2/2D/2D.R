# =============================================================================
# Fig. 2D + 2E — Standardized Levins niche breadth (B) of bacterial genera, IT vs ST
# =============================================================================
library(spaa)     # niche.width
library(irr)      # kappa2 
library(dplyr)
library(ggplot2)
library(ggpubr)

col_it <- "#E64B35"; col_st <- "#4DBBD5"
habitat_colors <- c(IT = col_it, ST = col_st)
habitat_labs   <- c(IT = "Intertidal", ST = "Subtidal")
PREV <- 0.20      
set.seed(1)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}

# ---- load relative-abundance table + metadata -------------------------------
abund <- as.matrix(read_auto("genus.csv", row.names = 1))
storage.mode(abund) <- "numeric"
if (anyNA(abund) || any(abund < 0)) stop("Invalid values in genus.csv")

meta <- read_auto("metadata.csv")
meta <- meta[meta$Sample %in% colnames(abund), ]
abund <- abund[, meta$Sample, drop = FALSE]

# ---- 20% prevalence cutoff across ALL samples ---------------------
abund <- abund[rowMeans(abund > 0) >= PREV, , drop = FALSE]
message(sprintf("Genera retained at %.0f%% prevalence: %d", PREV * 100, nrow(abund)))

IT <- meta$Sample[meta$Group == "IT"]
ST <- meta$Sample[meta$Group == "ST"]


levins_std <- function(mat, samples, group) {
  sub <- mat[, samples, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]          
  bw  <- spaa::niche.width(t(sub), method = "levins")   # B0 = 1/sum(p_i^2)
  B0  <- setNames(as.numeric(bw[1, ]), colnames(bw))
  n   <- length(samples)
  data.frame(genus = names(B0), B = (B0 - 1) / (n - 1),
             Group = group, stringsAsFactors = FALSE)
}
B_data <- rbind(levins_std(abund, IT, "IT"), levins_std(abund, ST, "ST"))
B_data$Group <- factor(B_data$Group, levels = c("IT", "ST"))

med <- aggregate(B ~ Group, B_data, median)
message(sprintf("Median B: IT=%.3f ST=%.3f",
                med$B[med$Group == "IT"], med$B[med$Group == "ST"]))


bi <- B_data$B[B_data$Group == "IT"]; bs <- B_data$B[B_data$Group == "ST"]
ks  <- suppressWarnings(ks.test(bi, bs))               
wil <- suppressWarnings(wilcox.test(bi, bs))           
message(sprintf("KS D=%.3f p=%.3g | Wilcoxon p=%.3g", ks$statistic, ks$p.value, wil$p.value))

classify <- function(B, lo, hi)
  factor(ifelse(B <= lo, "Specialist", ifelse(B >= hi, "Generalist", "Intermediate")),
         levels = c("Specialist", "Intermediate", "Generalist"))

schemes <- list(q25_75 = c(.25, .75), q20_80 = c(.20, .80), q10_90 = c(.10, .90))
sens <- lapply(names(schemes), function(nm) {
  pr <- schemes[[nm]]
  lo_s <- quantile(B_data$B, pr[1]); hi_s <- quantile(B_data$B, pr[2])
  lab <- classify(B_data$B, lo_s, hi_s)
  tab <- table(B_data$Group, lab == "Generalist")
  ft  <- fisher.test(tab)
  prop <- tapply(lab == "Generalist", B_data$Group, mean)
  data.frame(scheme = nm, lo = lo_s, hi = hi_s,
             IT_generalist_pct = round(100 * prop["IT"], 1),
             ST_generalist_pct = round(100 * prop["ST"], 1),
             fisher_p = ft$p.value, odds_ratio = unname(ft$estimate))
})
sens <- do.call(rbind, sens)
print(sens)
write.csv(sens, "Figure_2D_threshold_sensitivity_TableS2.csv", row.names = FALSE)

# ---- MAIN definition: upper/lower quartile (25%/75%) — used for BOTH the
#      classification AND the figure (single source of truth) -----------------
lo25 <- as.numeric(quantile(B_data$B, .25))
hi75 <- as.numeric(quantile(B_data$B, .75))
B_data$Class <- classify(B_data$B, lo25, hi75)
write.csv(B_data, "Figure_2D_niche_breadth_values.csv", row.names = FALSE)
message(sprintf("Main thresholds (25%%/75%%): lo=%.4f hi=%.4f", lo25, hi75))

base_theme <- theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(linewidth = 0.4),
        legend.position = "none",
        legend.title = element_blank())

dens_max <- max(sapply(split(B_data$B, B_data$Group),
                       function(v) max(density(v, from = 0, to = 1)$y)))
ecdf_df <- B_data %>% group_by(Group) %>% arrange(B) %>%
  mutate(cum = (row_number() / n()) * dens_max) %>% ungroup()
ytop <- dens_max * 1.04
ks_lab <- sprintf("K-S test: D = %.2f\np %s", ks$statistic,
                  if (ks$p.value < 0.001) "< 0.001" else sprintf("= %.3f", ks$p.value))

p2d <- ggplot() +
  geom_density(data = B_data, aes(B, fill = Group), alpha = 0.35, color = NA) +
  geom_step(data = ecdf_df, aes(B, cum, color = Group), linewidth = 0.7) +
  geom_vline(xintercept = c(lo25, hi75), linetype = "dashed", color = "grey25", linewidth = 0.4) +
  geom_vline(data = med, aes(xintercept = B, color = Group),
             linetype = "dashed", linewidth = 0.7) +
  annotate("text", x = lo25, y = ytop, label = sprintf("Specialists\nB < %.2f", lo25),
           hjust = 1.05, vjust = 1, size = 3, color = "grey25") +
  annotate("text", x = hi75, y = ytop, label = sprintf("Generalists\nB > %.2f", hi75),
           hjust = -0.05, vjust = 1, size = 3, color = "grey25") +
  annotate("text", x = med$B[med$Group == "IT"], y = dens_max * 0.55,
           label = sprintf("Median = %.2f", med$B[med$Group == "IT"]),
           angle = 90, vjust = -0.4, size = 3, color = col_it) +
  annotate("text", x = med$B[med$Group == "ST"], y = dens_max * 0.55,
           label = sprintf("Median = %.2f", med$B[med$Group == "ST"]),
           angle = 90, vjust = 1.2, size = 3, color = col_st) +
  annotate("text", x = 0.98, y = dens_max * 0.30, label = ks_lab,
           hjust = 1, vjust = 1, size = 3.2, color = "black") +
  scale_fill_manual(values = habitat_colors, labels = habitat_labs) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  scale_x_continuous("Niche breadth (B)", limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous("Density", expand = expansion(mult = c(0, 0.10)),
                     sec.axis = sec_axis(~ . / dens_max, name = "Cumulative proportion",
                                         labels = scales::percent)) +
  base_theme +
  theme(legend.position = c(0.5, 0.88), legend.direction = "horizontal",
        legend.background = element_blank())

print(p2d)
ggsave("Figure_2D_niche_breadth.pdf", p2d, width = 7, height = 4)
ggsave("Figure_2D_niche_breadth.png", p2d, width = 7, height = 4, dpi = 600)


sub2e <- B_data[B_data$Class %in% c("Specialist", "Generalist"), ]
sub2e$Class <- factor(sub2e$Class, levels = c("Specialist", "Generalist"))
sub2e$Class <- factor(sub2e$Class,
                      labels = c("Specialists", "Generalists"))  # facet titles

# per-facet stats: Wilcoxon p, rank-biserial r, and n per habitat -------------
library(rstatix)
stat2e <- sub2e %>% group_by(Class) %>%
  wilcox_test(B ~ Group) %>% ungroup()
eff2e <- sub2e %>% group_by(Class) %>%
  wilcox_effsize(B ~ Group) %>% ungroup()
n2e <- sub2e %>% count(Class, Group)

# label data: place p/r near the top of each facet (uses free-y, so per panel)
lab2e <- sub2e %>% group_by(Class) %>%
  summarise(ymax = max(B), ymin = min(B), .groups = "drop") %>%
  left_join(stat2e %>% select(Class, p), by = "Class") %>%
  left_join(eff2e  %>% select(Class, effsize), by = "Class") %>%
  mutate(lab = sprintf("r = %.2f\np = %s", effsize,
                       ifelse(p < 0.001, formatC(p, format = "e", digits = 2),
                              sprintf("%.3f", p))),
         ylab = ymax + (ymax - ymin) * 0.10)
n_lab <- n2e %>% left_join(sub2e %>% group_by(Class) %>%
                             summarise(ymin = min(B), .groups="drop"), by = "Class") %>%
  mutate(x = ifelse(Group == "IT", 1, 2),
         ylab = ymin - 0.0)   # n labels near bottom

p2e <- ggplot(sub2e, aes(Group, B, fill = Group)) +
  geom_violin(aes(color = Group), width = 0.85, alpha = 0.45,
              linewidth = 0.3, trim = FALSE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.95,
               linewidth = 0.3, color = "grey20") +
  geom_jitter(aes(color = Group), width = 0.08, size = 0.6, alpha = 0.4) +
  geom_text(data = lab2e, aes(x = 1.5, y = ylab, label = lab),
            inherit.aes = FALSE, size = 3.4, vjust = 1) +
  geom_text(data = n_lab, aes(x = x, y = ylab, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, vjust = 1.8) +
  facet_wrap(~ Class, scales = "free_y") +                 # independent y per panel
  scale_fill_manual(values = habitat_colors, labels = habitat_labs) +
  scale_color_manual(values = habitat_colors, labels = habitat_labs) +
  scale_x_discrete(labels = habitat_labs) +
  scale_y_continuous("Standardized niche breadth (B)",
                     expand = expansion(mult = c(0.05, 0.15))) +
  labs(x = NULL) +
  base_theme +
  theme(legend.position = c(0.92, 0.92), legend.direction = "vertical",
        legend.background = element_blank(),
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text = element_text(face = "bold", size = 12),
        panel.spacing = unit(1, "lines"))

print(p2e)
ggsave("Figure_2E_specialist_generalist.pdf", p2e, width = 3, height = 6)
ggsave("Figure_2E_specialist_generalist.png", p2e, width = 3, height = 6, dpi = 600)

message(sprintf("Median B: IT=%.3f ST=%.3f | quartiles lo=%.3f hi=%.3f | KS p=%.3g",
                med$B[med$Group == "IT"], med$B[med$Group == "ST"], lo25, hi75, ks$p.value))

occ_shannon <- function(mat, samples, group) {
  sub <- mat[, samples, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  Hb  <- spaa::niche.width(t(sub), method = "shannon")        
  H   <- setNames(as.numeric(Hb[1, ]), colnames(Hb))
  data.frame(genus = rownames(sub), Group = group,
             occ = rowMeans(sub > 0),                          
             Hstd = H / log(length(samples)),                  
             stringsAsFactors = FALSE)
}
val <- rbind(occ_shannon(abund, IT, "IT"), occ_shannon(abund, ST, "ST"))
val <- merge(B_data[, c("genus", "Group", "B")], val, by = c("genus", "Group"))

validation <- val %>% group_by(Group) %>% summarise(
  rho_B_occ  = cor(B, occ,  method = "spearman"),
  rho_B_Hstd = cor(B, Hstd, method = "spearman"),
  kappa_B_occ  = kappa2(data.frame(classify(B, quantile(B,.25), quantile(B,.75)),
                                   classify(occ, quantile(occ,.25), quantile(occ,.75))))$value,
  kappa_B_Hstd = kappa2(data.frame(classify(B, quantile(B,.25), quantile(B,.75)),
                                   classify(Hstd, quantile(Hstd,.25), quantile(Hstd,.75))))$value,
  .groups = "drop")
print(validation)
write.csv(validation, "Figure_2D_validation_TableS3.csv", row.names = FALSE)
# =============================================================================
# Prevalence-threshold sensitivity for niche breadth 
# =============================================================================
library(spaa)
library(dplyr)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}

abund0 <- as.matrix(read_auto("genus.csv", row.names = 1))   # relative abundance
storage.mode(abund0) <- "numeric"
meta <- read_auto("metadata.csv")
meta <- meta[meta$Sample %in% colnames(abund0), ]
abund0 <- abund0[, meta$Sample, drop = FALSE]
IT <- meta$Sample[meta$Group == "IT"]; ST <- meta$Sample[meta$Group == "ST"]

levins_std <- function(mat, samples) {
  sub <- mat[, samples, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  bw  <- spaa::niche.width(t(sub), method = "levins")
  B0  <- as.numeric(bw[1, ])
  (B0 - 1) / (length(samples) - 1)
}

PREV_GRID <- c(0, 0.05, 0.10, 0.20, 0.30, 0.50)

run_prev <- function(prev) {
  a  <- abund0[rowMeans(abund0 > 0) >= prev, , drop = FALSE]
  bi <- levins_std(a, IT); bs <- levins_std(a, ST)
  B  <- c(bi, bs); grp <- rep(c("IT", "ST"), c(length(bi), length(bs)))
  lo <- quantile(B, .25); hi <- quantile(B, .75)
  isG <- B >= hi
  data.frame(
    prevalence  = prev,
    n_genera    = nrow(a),
    median_IT   = round(median(bi), 3),
    median_ST   = round(median(bs), 3),
    KS_p        = signif(suppressWarnings(ks.test(bi, bs)$p.value), 3),
    Wilcoxon_p  = signif(suppressWarnings(wilcox.test(bi, bs)$p.value), 3),
    IT_gen_pct  = round(100 * mean(isG[grp == "IT"]), 1),
    ST_gen_pct  = round(100 * mean(isG[grp == "ST"]), 1),
    fisher_p    = signif(fisher.test(table(grp, isG))$p.value, 3)
  )
}

prev_sens <- do.call(rbind, lapply(PREV_GRID, run_prev))
print(prev_sens)
write.csv(prev_sens, "Figure_2D_prevalence_sensitivity.csv", row.names = FALSE)


# =============================================================================
# Validation of the niche-breadth filter with MicroNiche (Finn et al. 2020,
# FEMS Microbiol Ecol). 
# =============================================================================
library(MicroNiche)
library(spaa)
library(dplyr)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}


counts <- as.matrix(read_auto("genus_abund.csv", row.names = 1))
storage.mode(counts) <- "numeric"
if (any(abs(counts - round(counts)) > 1e-6))
  stop("MicroNiche/LOQ needs integer counts; genus_abund.csv is not count data.")
counts <- round(counts)

meta <- read_auto("metadata.csv")
samp_col <- intersect(c("Sample", "SampleID"), names(meta))[1]
meta$Sample <- as.character(meta[[samp_col]])
missing <- setdiff(meta$Sample, colnames(counts))
if (length(missing)) stop("Missing field samples in genus_abund.csv: ",
                          paste(missing, collapse = ", "))
counts <- counts[, meta$Sample, drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]

IT <- meta$Sample[meta$Group == "IT"]; ST <- meta$Sample[meta$Group == "ST"]


df <- data.frame(Taxon = rownames(counts), counts, check.names = FALSE)
R  <- ncol(counts)
sampleInfo <- colnames(counts)                 
micro <- levins.Bn(df, R, sampleInfo)          


micro$Taxon <- rownames(micro)                      
below_loq <- micro$Below.LOQ == "Y"                 
keep_taxa <- micro$Taxon[!below_loq]                
message(sprintf("Taxa below LOQ (dropped): %d | retained: %d of %d",
                sum(below_loq), length(keep_taxa), nrow(micro)))


prev_all <- rowMeans(counts > 0)
message(sprintf("LOQ-retained taxa prevalence: min=%.1f%% median=%.1f%% | 20%% rule keeps %d taxa",
                100 * min(prev_all[rownames(counts) %in% keep_taxa]),
                100 * median(prev_all[rownames(counts) %in% keep_taxa]),
                sum(prev_all >= 0.20)))


rel <- sweep(counts, 2, colSums(counts), "/")
rel <- rel[rownames(rel) %in% keep_taxa, , drop = FALSE]

levins_std <- function(mat, samples, group) {
  sub <- mat[, samples, drop = FALSE]; sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  bw  <- spaa::niche.width(t(sub), method = "levins")
  data.frame(species = colnames(bw),
             B = (as.numeric(bw[1, ]) - 1) / (length(samples) - 1),
             Group = group, stringsAsFactors = FALSE)
}
B <- rbind(levins_std(rel, IT, "IT"), levins_std(rel, ST, "ST"))

bi <- B$B[B$Group == "IT"]; bs <- B$B[B$Group == "ST"]
lo <- quantile(B$B, .25); hi <- quantile(B$B, .75)
isG <- B$B >= hi
ks  <- suppressWarnings(ks.test(bi, bs))
wil <- suppressWarnings(wilcox.test(bi, bs))
ft  <- fisher.test(table(B$Group, isG))

cat("\n=== LOQ-filtered niche-breadth comparison (validation) ===\n")
cat(sprintf("Median B: IT=%.3f  ST=%.3f\n", median(bi), median(bs)))
cat(sprintf("KS D=%.3f p=%.3g | Wilcoxon p=%.3g\n", ks$statistic, ks$p.value, wil$p.value))
cat(sprintf("Generalist%% (25/75): IT=%.1f  ST=%.1f  | Fisher p=%.3g\n",
            100 * mean(isG[B$Group == "IT"]), 100 * mean(isG[B$Group == "ST"]), ft$p.value))


padj_col <- names(micro)[grepl("adj", names(micro), ignore.case = TRUE)][1]
if (!is.na(padj_col)) {
  sig <- micro$P.adj < 0.05
  micro$class <- ifelse(sig & micro$Bn >= median(micro$Bn), "Generalist",
                        ifelse(sig & micro$Bn <  median(micro$Bn), "Specialist", "Intermediate"))
  write.csv(micro, "Niche_MicroNiche_levinsBn_nullmodel.csv", row.names = TRUE)
  print(table(micro$class))
}
# =============================================================================
# Fig. 2F — KO-level functional diversity, IT vs ST
# =============================================================================
library(vegan)
library(metagenomeSeq)
library(spaa)
library(ggplot2)
library(ggpubr)
library(patchwork)

col_it <- "#E64B35"; col_st <- "#4DBBD5"
habitat_labs <- c(IT = "Intertidal", ST = "Subtidal")
PREV <- 0.20
N_RAREFY <- 100
set.seed(1)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}

counts <- as.matrix(read_auto("kegg_abund.csv", row.names = 1)); storage.mode(counts) <- "numeric"
if (any(abs(counts - round(counts)) > 1e-6)) stop("kegg_abund.csv must be integer counts.")
counts <- round(counts)
meta <- read_auto("metadata.csv")
samp_col <- intersect(c("Sample", "SampleID"), names(meta))[1]
meta$Sample <- as.character(meta[[samp_col]])
missing <- setdiff(meta$Sample, colnames(counts))
if (length(missing)) stop("Field samples missing from kegg_abund.csv: ", paste(missing, collapse = ", "))
counts <- counts[, meta$Sample, drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
message(sprintf("KO table: %d features x %d samples", nrow(counts), ncol(counts)))
IT <- meta$Sample[meta$Group == "IT"]; ST <- meta$Sample[meta$Group == "ST"]

# ============================================================================
# (1) KO Shannon — mean over N rarefactions to the minimum depth (no extrapolation)
# ============================================================================
depth <- colSums(counts); dmin <- min(depth)
message(sprintf("Rarefy depth = %d ; averaging Shannon over %d rarefactions", dmin, N_RAREFY))
ct <- t(counts)                                   # samples x KOs
sh_mat <- replicate(N_RAREFY, diversity(rrarefy(ct, dmin), index = "shannon"))
shannon <- rowMeans(sh_mat)
sh_df <- data.frame(Sample = names(shannon), Shannon = as.numeric(shannon))
sh_df <- merge(sh_df, meta[, c("Sample", "Group")], by = "Sample")
sh_df$Habitat <- factor(habitat_labs[sh_df$Group], levels = habitat_labs)

print(tapply(sh_df$Shannon, sh_df$Group, summary))
sh_wt <- wilcox.test(Shannon ~ Group, data = sh_df)
message(sprintf("KO Shannon (mean of %d rarefactions): IT med=%.3f ST med=%.3f | Wilcoxon p=%.3g",
                N_RAREFY, median(sh_df$Shannon[sh_df$Group == "IT"]),
                median(sh_df$Shannon[sh_df$Group == "ST"]), sh_wt$p.value))

# ============================================================================
# (2) KO standardized Levins' B — CSS-normalized (metagenomeSeq), UNCHANGED
# ============================================================================
mr  <- newMRexperiment(counts)
mr  <- cumNorm(mr, p = cumNormStatFast(mr))
css <- MRcounts(mr, norm = TRUE, log = FALSE)
css <- css[rowMeans(counts > 0) >= PREV, , drop = FALSE]

levins_std <- function(mat, samples, group) {
  sub <- mat[, samples, drop = FALSE]; sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  bw  <- spaa::niche.width(t(sub), method = "levins")
  data.frame(KO = colnames(bw),
             B = (as.numeric(bw[1, ]) - 1) / (length(samples) - 1),
             Group = group, stringsAsFactors = FALSE)
}
B_df <- rbind(levins_std(css, IT, "IT"), levins_std(css, ST, "ST"))
B_df$Habitat <- factor(habitat_labs[B_df$Group], levels = habitat_labs)
B_wt <- wilcox.test(B ~ Group, data = B_df)
B_ks <- suppressWarnings(ks.test(B_df$B[B_df$Group == "IT"], B_df$B[B_df$Group == "ST"]))
message(sprintf("KO Levins' B: IT med=%.3f ST med=%.3f | Wilcoxon p=%.3g KS p=%.3g",
                median(B_df$B[B_df$Group == "IT"]), median(B_df$B[B_df$Group == "ST"]),
                B_wt$p.value, B_ks$p.value))

write.csv(sh_df, "Figure_2F_KO_shannon.csv", row.names = FALSE)
write.csv(B_df,  "Figure_2F_KO_levinsB_css.csv", row.names = FALSE)

# ============================================================================
# (3) Plot — two stacked panels
# ============================================================================
base_theme <- theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black"), axis.line = element_line(linewidth = 0.4),
        legend.position = "none")
fill_sc <- scale_fill_manual(values = setNames(c(col_it, col_st), habitat_labs))
col_sc  <- scale_color_manual(values = setNames(c(col_it, col_st), habitat_labs))
cmp     <- list(c("Intertidal", "Subtidal"))

p_sh <- ggplot(sh_df, aes(Habitat, Shannon, fill = Habitat)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85, color = "grey20", linewidth = 0.3) +
  geom_jitter(aes(color = Habitat), width = 0.12, size = 1.4, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", comparisons = cmp, label = "p.signif", size = 5) +
  fill_sc + col_sc + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "KO-level Shannon diversity") + base_theme

p_B <- ggplot(B_df, aes(Habitat, B, fill = Habitat)) +
  geom_violin(aes(color = Habitat), width = 0.85, alpha = 0.45, linewidth = 0.3, trim = FALSE) +
  geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.95, color = "grey20", linewidth = 0.3) +
  stat_compare_means(method = "wilcox.test", comparisons = cmp, label = "p.signif", size = 5) +
  fill_sc + col_sc + scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "KO-level Levins' B") + base_theme

fig2f <- p_sh / p_B
ggsave("Figure_2F_KO_diversity.pdf", fig2f, width = 4, height = 6)
ggsave("Figure_2F_KO_diversity.png", fig2f, width = 4, height = 6, dpi = 600)
print(fig2f)

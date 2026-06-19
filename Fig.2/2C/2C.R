# =============================================================================
# Fig. 2C — Genus-level Shannon alpha diversity: Intertidal vs Subtidal
# =============================================================================
library(vegan)
library(ggplot2)
library(ggpubr)

habitat_colors <- c("Intertidal" = "#E64B35FF", "Subtidal" = "#4DBBD5FF")
SEED <- 1
set.seed(SEED)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names,
           check.names = FALSE)
}


counts <- as.matrix(read_auto("genus_abund.csv", row.names = 1))  # genus READ COUNTS
storage.mode(counts) <- "numeric"
if (anyNA(counts) || any(counts < 0)) stop("Invalid values in genus_abund.csv")
if (any(abs(counts - round(counts)) > 1e-6))
  stop("Rarefaction needs integer counts; genus_abund.csv is not count data.")
counts <- round(counts)

meta <- read_auto("metadata.csv")
stopifnot(all(c("Sample", "Group") %in% colnames(meta)))
meta$Sample  <- as.character(meta$Sample)
meta$Habitat <- ifelse(meta$Group == "IT", "Intertidal", "Subtidal")

missing <- setdiff(meta$Sample, colnames(counts))
if (length(missing))
  stop("Field samples missing from genus_abund.csv: ", paste(missing, collapse = ", "))
counts <- counts[, meta$Sample, drop = FALSE]          


depth <- colSums(counts)
depth_p <- wilcox.test(depth ~ meta$Group)$p.value
message(sprintf("Sequencing depth (reads): min=%d max=%d | Wilcoxon depth diff p=%.6f",
                min(depth), max(depth), depth_p))


counts_rare <- rrarefy(t(counts), sample = min(depth))  
shannon <- diversity(counts_rare, index = "shannon")

alpha_df <- data.frame(Sample = names(shannon), Shannon = as.numeric(shannon))
alpha_df <- merge(alpha_df, meta[, c("Sample", "Group", "Habitat")], by = "Sample")
alpha_df$Habitat <- factor(alpha_df$Habitat, levels = c("Intertidal", "Subtidal"))


wt <- wilcox.test(Shannon ~ Habitat, data = alpha_df)
med <- aggregate(Shannon ~ Habitat, alpha_df, median)
message(sprintf("Shannon medians: Intertidal=%.3f Subtidal=%.3f | Wilcoxon p=%.4g",
                med$Shannon[med$Habitat == "Intertidal"],
                med$Shannon[med$Habitat == "Subtidal"], wt$p.value))

write.csv(alpha_df, "Figure_2C_shannon_genus.csv", row.names = FALSE)

# ---- plot -------------------------------------------------------------------
p <- ggplot(alpha_df, aes(Habitat, Shannon, fill = Habitat)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.85, color = "black") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.6) +
  scale_fill_manual(values = habitat_colors) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("Intertidal", "Subtidal")),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols   = c("***", "**", "*", "ns")),
                     size = 6) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = NULL, y = "Shannon index (genus level)") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

print(p)
ggsave("Figure_2C_genus_Shannon.pdf", p, width = 4, height = 4.5)
ggsave("Figure_2C_genus_Shannon.png", p, width = 4, height = 4.5, dpi = 600)

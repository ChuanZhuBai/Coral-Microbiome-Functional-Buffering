# =============================================================================
# Fig. 2A — Genus-level PCoA (Bray-Curtis) of field communities: Intertidal vs Subtidal
# =============================================================================
library(vegan)
library(ggplot2)

habitat_colors <- c("Intertidal" = "#E64B35FF", "Subtidal" = "#4DBBD5FF")
habitat_shapes <- c("Intertidal" = 16, "Subtidal" = 17)


read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names,
           check.names = FALSE)
}


abund <- as.matrix(read_auto("genus.csv", row.names = 1))
storage.mode(abund) <- "numeric"
if (anyNA(abund)) stop("Missing/non-numeric values in genus.csv")
if (any(abund < 0)) stop("Negative values in genus.csv")


cs <- colSums(abund)
if (max(abs(cs - 1)) > 0.01) {
  message("Columns do not sum to 1; converting to relative abundance.")
  abund <- sweep(abund, 2, cs, "/")
}
abund <- abund[rowSums(abund) > 0, , drop = FALSE]   


meta <- read_auto("metadata.csv")
if (!setequal(colnames(abund), meta$SampleID))
  stop("Sample names in genus.csv and metadata.csv do not match.")
meta <- meta[match(colnames(abund), meta$SampleID), ]   
meta$Habitat <- ifelse(meta$Group == "IT", "Intertidal", "Subtidal")


bray <- vegdist(t(abund), method = "bray")

set.seed(123)
per <- adonis2(bray ~ Group, data = meta, permutations = 9999)  
R2  <- per$R2[1]
pval <- per$`Pr(>F)`[1]


bd <- betadisper(bray, meta$Group)
bd_p <- permutest(bd, permutations = 9999)$tab$`Pr(>F)`[1]
message(sprintf("PERMANOVA R2=%.3f p=%.4f | betadisper p=%.4f", R2, pval, bd_p))

pc <- cmdscale(bray, k = 2, eig = TRUE)
eig <- pc$eig
pc1 <- round(eig[1] / sum(eig) * 100, 1)
pc2 <- round(eig[2] / sum(eig) * 100, 1)

scores <- data.frame(SampleID = rownames(pc$points),
                     PCoA1 = pc$points[, 1], PCoA2 = pc$points[, 2])
scores <- merge(scores, meta[, c("SampleID", "Habitat")], by = "SampleID")

# ---- plot -------------------------------------------------------------------
anno <- sprintf("PERMANOVA\nR\u00B2 = %.2f\np %s", R2,
                if (pval < 0.001) "< 0.001" else sprintf("= %.3f", pval))

p <- ggplot(scores, aes(PCoA1, PCoA2, color = Habitat, fill = Habitat, shape = Habitat)) +
  stat_ellipse(geom = "polygon", alpha = 0.15, linetype = 1, show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.85) +
  scale_color_manual(name = NULL, values = habitat_colors) +
  scale_fill_manual(name = NULL, values = habitat_colors) +
  scale_shape_manual(name = NULL, values = habitat_shapes) +
  labs(x = sprintf("PCoA1 (%.1f%%)", pc1), y = sprintf("PCoA2 (%.1f%%)", pc2)) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3,
           label = anno, size = 4.5) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = "black"),
    legend.position = c(0.85, 0.9),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.border = element_rect(linewidth = 1.2, color = "black")
  )

print(p)
ggsave("Figure_2A_genus_PCoA.pdf", p, width = 6, height = 5.5)
ggsave("Figure_2A_genus_PCoA.png", p, width = 6, height = 5.5, dpi = 600)
print(tapply(bd$distances, meta$Group, mean))  

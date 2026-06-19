# =============================================================================
# habitat comparison of microbiome stability under heat stress
# =============================================================================
library(vegan); library(tidyverse); library(rstatix); library(ggpubr)
library(patchwork); library(metagenomeSeq)

hab_cols <- c(Intertidal = "#E64B35FF", Subtidal = "#4DBBD5FF")
set.seed(1)

read_auto <- function(file, row.names = 1) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}
canonical_id <- function(x) {
  x %>%
    gsub("^(IT|ST)[-]?C([0-9]+)$", "\\1-C\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)C$",     "\\1-C\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)A$",     "\\1-C\\2", .) %>%
    gsub("^(IT|ST)[-]?H([0-9]+)$", "\\1-H\\2", .) %>%
    gsub("^(IT|ST)([0-9]+)H$",     "\\1-H\\2", .)
}

# ---- metadata: Habitat + Treatment ------------------------------------------
metadata <- read_auto("metadata.csv", row.names = NULL)
colnames(metadata)[1:2] <- c("Sample","Group")
metadata$Sample <- canonical_id(as.character(metadata$Sample))
metadata <- metadata %>%
  filter(Group %in% c("IT-C","IT-H","ST-C","ST-H")) %>%
  mutate(Habitat   = factor(ifelse(grepl("^IT", Group), "Intertidal", "Subtidal"),
                            levels = c("Intertidal","Subtidal")),
         Treatment = ifelse(grepl("-C$", Group), "Control", "Heat"))

# ---- per-sample stability = mean BC distance to the opposite treatment -------
# (within the same habitat), computed on a properly normalised table.
per_sample_stability <- function(file, norm = c("relabund","css")) {
  norm <- match.arg(norm)
  mat <- as.matrix(read_auto(file)); storage.mode(mat) <- "numeric"
  colnames(mat) <- canonical_id(colnames(mat))
  common <- intersect(metadata$Sample, colnames(mat))
  miss <- setdiff(metadata$Sample, colnames(mat))
  if (length(miss)) message("WARNING [", file, "] not found: ", paste(miss, collapse=", "))
  message(sprintf("[%s] matched %d / %d samples", file, length(common), nrow(metadata)))
  mat <- mat[, common, drop = FALSE]
  
  if (norm == "relabund") {
    mat <- sweep(mat, 2, colSums(mat), "/")
  } else {
    obj <- newMRexperiment(mat); obj <- cumNorm(obj, p = cumNormStatFast(obj))
    mat <- MRcounts(obj, norm = TRUE, log = FALSE)
    mat <- sweep(mat, 2, colSums(mat), "/")
  }
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  
  D <- as.matrix(vegdist(t(mat), method = "bray"))
  md <- metadata[match(common, metadata$Sample), ]
  
  # for each sample: mean distance to opposite-treatment samples of same habitat
  res <- map_dfr(seq_along(common), function(i) {
    s <- common[i]; hab <- md$Habitat[i]; trt <- md$Treatment[i]
    opp <- md$Sample[md$Habitat == hab & md$Treatment != trt]
    tibble(Sample = s, Habitat = hab, Treatment = trt,
           Stability_distance = mean(D[s, opp]))
  })
  res
}

# ---- analysis + plot --------------------------------------------------------
analyse_panel <- function(file, norm, title, ylab) {
  d <- per_sample_stability(file, norm)
  d$Habitat <- factor(d$Habitat, levels = c("Intertidal","Subtidal"))
  
  wt  <- d %>% wilcox_test(Stability_distance ~ Habitat) %>%
    add_significance() %>% add_y_position()
  eff <- d %>% wilcox_effsize(Stability_distance ~ Habitat)
  message(sprintf("[%s] IT median=%.3f ST median=%.3f | Wilcoxon p=%.3g | r=%.3f (%s)",
                  title,
                  median(d$Stability_distance[d$Habitat=="Intertidal"]),
                  median(d$Stability_distance[d$Habitat=="Subtidal"]),
                  wt$p, eff$effsize, eff$magnitude))
  
  p <- ggplot(d, aes(Habitat, Stability_distance)) +
    geom_violin(aes(fill = Habitat), alpha = 0.2, color = NA) +
    geom_boxplot(aes(fill = Habitat), width = 0.2, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.1, alpha = 0.5, size = 1.6) +
    scale_fill_manual(values = hab_cols) +
    stat_pvalue_manual(wt, label = "p = {p}{p.signif}", tip.length = 0.01) +
    labs(title = title, subtitle = sprintf("Effect size r = %.3f (%s)", eff$effsize, eff$magnitude),
         x = NULL, y = ylab) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face="bold", hjust=0.5),
          plot.subtitle = element_text(hjust=0.5, size=10),
          legend.position = "none", axis.text = element_text(color="black"))
  list(plot = p, data = d, stat = wt, eff = eff)
}

res_C <- analyse_panel("genus_abund.csv", "relabund", "Taxonomic Stability",
                       "Mean Bray-Curtis distance\nto opposite treatment (Control vs Heat)")
res_D <- analyse_panel("ko_abund.csv",    "css",      "Functional Stability",
                       "Mean Bray-Curtis distance\nto opposite treatment (Control vs Heat)")

# ---- export -----------------------------------------------------------------
write.csv(res_C$data, "SourceData_3C_stability_genus.csv", row.names = FALSE)
write.csv(res_D$data, "SourceData_3D_stability_KO.csv",    row.names = FALSE)


stat_C <- res_C$stat %>%
  left_join(res_C$eff, by = c(".y.", "group1", "group2")) %>%
  mutate(Panel = "Fig.3C (Taxonomic)")


stat_D <- res_D$stat %>%
  left_join(res_D$eff, by = c(".y.", "group1", "group2")) %>%
  mutate(Panel = "Fig.3D (Functional)")


stat_out <- bind_rows(stat_C, stat_D) %>%
  select(Panel, everything()) 

stat_out_clean <- stat_out %>%
  mutate(across(where(is.list), ~ sapply(.x, paste, collapse = " vs ")))


write.csv(stat_out_clean, "Stats_3C_3D_stability.csv", row.names = FALSE)


fig <- res_C$plot + res_D$plot + plot_annotation(tag_levels = list(c("C", "D")))
ggsave("Figure_3C_3D_stability.pdf", fig, width = 8, height = 5)
print(fig)
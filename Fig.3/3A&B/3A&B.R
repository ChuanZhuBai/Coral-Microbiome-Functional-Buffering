# =============================================================================
# CPCoA (db-RDA) of the in-lab heat-stress experiment:
# =============================================================================
library(tidyverse); library(vegan); library(metagenomeSeq); library(patchwork)

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


metadata <- read_auto("metadata.csv", row.names = NULL)
samp_col <- intersect(c("Sample","SampleID"), names(metadata))[1]
metadata$Sample <- as.character(metadata[[samp_col]])
metadata$Group  <- factor(metadata$Group, levels = c("IT-C","IT-H","ST-C","ST-H"))
metadata <- metadata[!is.na(metadata$Group), ]
rownames(metadata) <- metadata$Sample


prep_table <- function(file) {
  mat <- as.matrix(read_auto(file)); storage.mode(mat) <- "numeric"
  colnames(mat) <- canonical_id(colnames(mat))
  if (any(duplicated(colnames(mat))))
    stop("Duplicated sample names after canonicalisation in ", file, ": ",
         paste(unique(colnames(mat)[duplicated(colnames(mat))]), collapse=", "))
  common <- intersect(metadata$Sample, colnames(mat))
  miss <- setdiff(metadata$Sample, colnames(mat))
  if (length(miss)) message("WARNING [", file, "] metadata samples not found: ",
                            paste(miss, collapse = ", "))
  message(sprintf("[%s] matched %d / %d experimental samples", file, length(common), nrow(metadata)))
  mat[, common, drop = FALSE]
}


run_permanova <- function(dist_mat, env) {
  set.seed(123)
  glob <- adonis2(dist_mat ~ Group, data = env, permutations = 999)
  contrasts <- list(c("IT-C","IT-H"), c("ST-C","ST-H"))
  labels    <- c("Intertidal (Control vs Heat)", "Subtidal (Control vs Heat)")
  pw <- map_dfr(seq_along(contrasts), function(i) {
    g <- contrasts[[i]]
    se <- env[env$Group %in% g, ]; se$Group <- droplevels(se$Group)
    if (nlevels(se$Group) < 2) return(tibble(Comparison = labels[i], R2 = NA, P_raw = NA))
    sd <- as.dist(as.matrix(dist_mat)[rownames(se), rownames(se)])
    set.seed(123 + i)
    r <- adonis2(sd ~ Group, data = se, permutations = 999)
    tibble(Comparison = labels[i], R2 = r$R2[1], P_raw = r$`Pr(>F)`[1])
  })
  pw$P_adj_BH <- p.adjust(pw$P_raw, method = "BH")
  list(global = glob, pairwise = pw)
}


col_pal <- c("IT-C"="#FC8D62","IT-H"="#E41A1C","ST-C"="#4DBBD5","ST-H"="#377EB8")
make_cpcoa_plot <- function(mod, env, title, pairwise) {
  sc <- as.data.frame(scores(mod, display = "sites", choices = 1:2))
  colnames(sc) <- c("CAP1","CAP2")
  sc$Group <- env$Group[match(rownames(sc), env$Sample)]
  
 
  cap_eig <- mod$CCA$eig
  cap1 <- round(cap_eig[1] / sum(cap_eig) * 100, 1)
  cap2 <- round(cap_eig[2] / sum(cap_eig) * 100, 1)
  
  tot  <- round(mod$CCA$tot.chi / mod$tot.chi * 100, 1)
  pval <- anova.cca(mod, permutations = 999)$`Pr(>F)`[1]
  
  
  pw_txt <- paste0("Pairwise PERMANOVA:\n",
                   paste(sprintf("%s: R\u00b2 = %.2f, p = %.3g",
                                 pairwise$Comparison, pairwise$R2, pairwise$P_raw),
                         collapse = "\n"))
  
  ggplot(sc, aes(CAP1, CAP2, color = Group, fill = Group)) +
    stat_ellipse(aes(group = Group), geom = "polygon", alpha = 0.1, linetype = 2, show.legend = FALSE) +
    geom_point(aes(shape = Group), size = 3.5, stroke = 0.6) +
    annotate("text", x = -Inf, y = -Inf, hjust = -0.05, vjust = -0.3,
             label = pw_txt, size = 2.9, lineheight = 0.95) +
    scale_color_manual(values = col_pal) + scale_fill_manual(values = col_pal) +
    scale_shape_manual(values = c("IT-C"=16,"IT-H"=16,"ST-C"=17,"ST-H"=17)) +
    labs(title = title,
         subtitle = sprintf("Constrained variance: %.1f%%, p = %s",
                            tot, if (pval < 0.001) "< 0.001" else sprintf("%.3f", pval)),
         x = sprintf("CPCoA 1 (%.1f%%)", cap1), y = sprintf("CPCoA 2 (%.1f%%)", cap2)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), plot.title = element_text(face="bold", hjust=0.5),
          plot.subtitle = element_text(hjust=0.5, size=9), axis.text = element_text(color="black"),
          legend.title = element_blank())
}


genus <- prep_table("genus_abund.csv")
env_g <- metadata[colnames(genus), ]
genus_rel <- sweep(genus, 2, colSums(genus), "/")
genus_rel <- genus_rel[rowSums(genus_rel) > 0, , drop = FALSE]
dist_g <- vegdist(t(genus_rel), method = "bray")
stat_g <- run_permanova(dist_g, env_g)
cpcoa_g <- capscale(dist_g ~ Group, data = env_g)
p_3A <- make_cpcoa_plot(cpcoa_g, env_g, "Community Structure (Genus)", stat_g$pairwise)
message(sprintf("genus: global R2=%.3f, p=%.3g", stat_g$global$R2[1], stat_g$global$`Pr(>F)`[1]))
print(stat_g$pairwise)


ko <- prep_table("ko_abund.csv")
ko <- ko[!tolower(rownames(ko)) %in% c("unclassified","unassigned"), , drop = FALSE]
env_k <- metadata[colnames(ko), ]
obj <- newMRexperiment(ko); obj <- cumNorm(obj, p = cumNormStatFast(obj))
ko_css <- MRcounts(obj, norm = TRUE, log = FALSE)
ko_css <- ko_css[rowSums(ko_css) > 0, , drop = FALSE]
dist_k <- vegdist(t(ko_css), method = "bray")
stat_k <- run_permanova(dist_k, env_k)
cpcoa_k <- capscale(dist_k ~ Group, data = env_k)
p_3B <- make_cpcoa_plot(cpcoa_k, env_k, "Functional Structure (KOs)", stat_k$pairwise)
message(sprintf("KO: global R2=%.3f, p=%.3g", stat_k$global$R2[1], stat_k$global$`Pr(>F)`[1]))
print(stat_k$pairwise)


write.csv(stat_g$pairwise, "TableS_pairwise_PERMANOVA_genus.csv", row.names = FALSE)
write.csv(stat_k$pairwise, "TableS_pairwise_PERMANOVA_KO.csv", row.names = FALSE)


sc_g <- as.data.frame(scores(cpcoa_g, display="sites", choices=1:2)) %>%
  rownames_to_column("Sample") %>% left_join(env_g[,c("Sample","Group")], by="Sample") %>%
  mutate(Level = "Genus")
sc_k <- as.data.frame(scores(cpcoa_k, display="sites", choices=1:2)) %>%
  rownames_to_column("Sample") %>% left_join(env_k[,c("Sample","Group")], by="Sample") %>%
  mutate(Level = "KO")
write.csv(bind_rows(sc_g, sc_k), "SourceData_3AB_CPCoA_scores.csv", row.names = FALSE)


summary_tab <- tibble(
  Panel = c("3A Genus","3B KO"),
  Constrained_variance_pct = c(round(cpcoa_g$CCA$tot.chi/cpcoa_g$tot.chi*100,1),
                               round(cpcoa_k$CCA$tot.chi/cpcoa_k$tot.chi*100,1)),
  Global_R2 = c(stat_g$global$R2[1], stat_k$global$R2[1]),
  Global_p  = c(stat_g$global$`Pr(>F)`[1], stat_k$global$`Pr(>F)`[1]))
write.csv(summary_tab, "SourceData_3AB_global_summary.csv", row.names = FALSE)
print(summary_tab)

fig3AB <- p_3A + p_3B + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave("Figure_3A_3B_CPCoA.pdf", fig3AB, width = 10, height = 5)
print(p_3A); print(p_3B)
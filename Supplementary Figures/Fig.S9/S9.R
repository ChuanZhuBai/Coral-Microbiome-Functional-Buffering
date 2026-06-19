# =============================================================================
# Experimental alpha diversity (Shannon) — genus and KO level
# =============================================================================
library(tidyverse); library(vegan); library(rstatix); library(patchwork)

set.seed(123)
N_RAREFY <- 100                                  

col_pal <- c("IT-C"="#FC8D62","IT-H"="#E41A1C","ST-C"="#4DBBD5","ST-H"="#377EB8")

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


shannon_rarefied <- function(count_file) {
  mat <- as.matrix(read_auto(count_file)); storage.mode(mat) <- "numeric"
  colnames(mat) <- canonical_id(colnames(mat))
  common <- intersect(metadata$Sample, colnames(mat))
  miss <- setdiff(metadata$Sample, colnames(mat))
  if (length(miss)) message("WARNING [", count_file, "] not found: ", paste(miss, collapse=", "))
  message(sprintf("[%s] matched %d / %d samples", count_file, length(common), nrow(metadata)))
  mat <- round(mat[, common, drop = FALSE])
  depth <- min(colSums(mat))
  message(sprintf("   rarefy depth = %d; averaging %d rarefactions", depth, N_RAREFY))
  
  sh <- replicate(N_RAREFY, {
    r <- rrarefy(t(mat), sample = depth)          
    diversity(r, index = "shannon")
  })
  tibble(Sample = common, Shannon = rowMeans(sh))
}


analyse_panel <- function(count_file, ylab, title) {
  d <- shannon_rarefied(count_file) %>%
    left_join(metadata[, c("Sample","Group")], by = "Sample")
  kw <- d %>% kruskal_test(Shannon ~ Group)
  
  p <- ggplot(d, aes(Group, Shannon, fill = Group)) +
    geom_violin(alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.85, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.45, size = 1.5) +
    scale_fill_manual(values = col_pal) +
    labs(title = title, x = NULL, y = ylab,
         subtitle = sprintf("Kruskal-Wallis p = %.3g", kw$p)) +
    theme_bw(base_size = 12) +
    theme(panel.grid = element_blank(), legend.position = "none",
          plot.title = element_text(face="bold", hjust=0.5),
          axis.text.x = element_text(face="bold"), axis.text = element_text(color="black"))
  
  list(plot = p, data = d, kw = kw)
}


res_gen <- analyse_panel("genus_abund.csv", "Shannon index (genus)", "Taxonomic alpha diversity")
res_ko  <- analyse_panel("ko_abund.csv",    "Shannon index (KO)",    "Functional alpha diversity")
message(sprintf("Genus KW p=%.3g | KO KW p=%.3g", res_gen$kw$p, res_ko$kw$p))


write.csv(bind_rows(res_gen$data %>% mutate(Level="Genus"),
                    res_ko$data  %>% mutate(Level="KO")),
          "SourceData_experimental_Shannon.csv", row.names = FALSE)
kw_tab <- bind_rows(res_gen$kw %>% mutate(Level="Genus"),
                    res_ko$kw  %>% mutate(Level="KO"))
write.csv(kw_tab, "Stats_KW_experimental_Shannon.csv", row.names = FALSE)
print(kw_tab)

fig <- res_gen$plot + res_ko$plot + plot_annotation(tag_levels = "A")
ggsave("Figure_S_experimental_Shannon.pdf", fig, width = 10, height = 5.5)
print(fig)
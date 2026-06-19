# =============================================================================
# KEGG pathway over-representation (enrichment) bubble plot, IT vs ST
# =============================================================================
library(edgeR)
library(clusterProfiler)
library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(forcats)
library(patchwork)
library(tidytext)
library(ggrepel) 
library(ggsci)
library(ggforce)

FDR_KO  <- 0.05    # significance for differential KOs
LFC_KO  <- 1
MIN_PREV <- 0.1   # prevalence filter for edgeR 
FDR_PATH <- 0.05   # significance for enriched pathways
TOPN <- 10         # pathways shown per direction
set.seed(1)

read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, row.names = row.names, check.names = FALSE)
}

# ---- (1) differential KOs via edgeR  -------
counts <- as.matrix(read_auto("kegg_abund.csv", row.names = 1)); storage.mode(counts) <- "numeric"
counts <- round(counts)
meta <- read_auto("metadata.csv")
samp_col <- intersect(c("Sample", "SampleID"), names(meta))[1]
meta$Sample <- as.character(meta[[samp_col]])
missing <- setdiff(meta$Sample, colnames(counts))
if (length(missing)) stop("Field samples missing from kegg_abund.csv: ", paste(missing, collapse = ", "))
counts <- counts[, meta$Sample, drop = FALSE]

group <- factor(meta$Group, levels = c("ST", "IT"))     
keep  <- rowMeans(counts > 0) >= MIN_PREV               
y <- DGEList(counts = counts[keep, , drop = FALSE], group = group)
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = 2)                        
res <- topTags(qlf, n = Inf)$table
res$KO <- rownames(res)

# ---- (2) KO -> KEGG Level-3 pathway map (TERM2GENE) --------------------------
kegg_paths <- read_auto("KEGG_paths.csv")
pcol <- setdiff(names(kegg_paths), "KO")[1]              
t2g <- kegg_paths %>%
  transmute(KO, paths = gsub("\t", " ", .data[[pcol]])) %>%
  separate_rows(paths, sep = "\\|") %>%                  
  mutate(paths = str_trim(paths)) %>%
  filter(paths != "") %>%
  mutate(Level3 = vapply(str_split(paths, ";"), function(p) {
    p <- str_trim(p[p != ""]); if (length(p) >= 3) p[3] else p[length(p)]
  }, character(1))) %>%
  filter(!is.na(Level3), Level3 != "") %>%
  distinct(Level3, KO) %>%
  select(term = Level3, gene = KO)                       # TERM2GENE: term, gene

# ---- (3) define foreground KO sets + annotated background --------------------
universe <- intersect(res$KO, unique(t2g$gene))          
sig <- res[res$FDR < FDR_KO & abs(res$logFC) > LFC_KO, ]
IT_KO <- intersect(sig$KO[sig$logFC > 0], universe)      # enriched in intertidal
ST_KO <- intersect(sig$KO[sig$logFC < 0], universe)      # enriched in subtidal
message(sprintf("Differential KOs (annotated): %d IT-enriched, %d ST-enriched (background %d)",
                length(IT_KO), length(ST_KO), length(universe)))
run_ora <- function(genes, label) {
  e <- enricher(gene = genes, universe = universe, TERM2GENE = t2g,
                pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                minGSSize = 3, maxGSSize = 500)
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  df <- as.data.frame(e)
  df$RichFactor <- df$Count / as.numeric(sub("/\\d+", "", df$BgRatio))  # k / M
  df$Direction <- label
  df
}
ora <- bind_rows(run_ora(IT_KO, "Intertidal enriched"),
                 run_ora(ST_KO, "Subtidal enriched"))
write.csv(ora, "Figure_2G_KEGG_enrichment_full.csv", row.names = FALSE)


# ---- (4) bubble plot — two independent panels (Top-Bottom Facets) ---------
if (is.null(ora) || nrow(ora) == 0) {
  stop("No ORA result returned. Check KO thresholds, TERM2GENE mapping and universe.")
}

plot_df <- ora %>%
  filter(!is.na(p.adjust), p.adjust < FDR_PATH) %>%
  mutate(
    Direction = factor(Direction, 
                       levels = c("Intertidal enriched", "Subtidal enriched"),
                       labels = c("Intertidal Enriched", "Subtidal Enriched")),
    neglogFDR = -log10(pmax(p.adjust, .Machine$double.xmin)),
    Description = str_squish(Description)
  ) %>%
  group_by(Direction) %>%
  arrange(p.adjust, desc(Count), .by_group = TRUE) %>%
  slice_head(n = TOPN) %>%
  ungroup()

if (nrow(plot_df) == 0) {
  stop(sprintf("No enriched pathways passed FDR < %.2f.", FDR_PATH))
}


plot_df <- plot_df %>%
  mutate(Description = reorder_within(Description, RichFactor, Direction))

xmax <- max(plot_df$RichFactor, na.rm = TRUE) * 1.15


fig2g <- ggplot(plot_df, aes(x = RichFactor, y = Description)) +
  
  geom_point(aes(size = Count, color = neglogFDR), alpha = 0.9) +
  
  geom_text_repel(aes(label = sprintf("%.2f", RichFactor)),
                  color = "grey20", size = 3.2, direction = "x",
                  nudge_x = 0.03, segment.color = NA, show.legend = FALSE, seed = 42) +         
  
  
  facet_col(vars(Direction), scales = "free_y", space = "free") +
  
  
  scale_y_reordered() +
  
  scale_color_gradient(low = "#4DBBD5", high = "#E64B35",
                       name = expression(-log[10]~"(FDR)")) +
  
  scale_size_continuous(range = c(2.5, 7.5), name = "Number",
                        breaks = cnt_breaks, limits = range(plot_df$Count)) +
  
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "Rich factor", y = NULL) +
  
  theme_bw(base_size = 11, base_family = "sans") +
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "black", linewidth = 0.8), 
    
   
    strip.background = element_rect(fill = "grey85", color = "black", linewidth = 0.8),
    strip.text = element_text(color = "black", size = 11, face = "bold", margin = margin(t=4, b=4)),
    
    axis.text = element_text(color = "black"),
    axis.text.y = element_text(size = 9.5),
    axis.title.x = element_text(margin = margin(t = 8)),
    
    legend.position = "right",
    legend.key = element_blank()
  )


fig_h <- max(5.0, 2.0 + 0.28 * nrow(plot_df))


ggsave("Figure_2G_KEGG_enrichment_FacetStyle.pdf", fig2g,
       width = 6.8, height = fig_h, device = grDevices::cairo_pdf)

ggsave("Figure_2G_KEGG_enrichment_FacetStyle.png", fig2g,
       width = 6.8, height = fig_h, dpi = 600)

print(fig2g)



# =============================================================================
# Fig. S4 — Volcano plot of differentially abundant KOs (field: Intertidal vs Subtidal)
# =============================================================================



# ---- packages ----
library(edgeR)      # DGEList, calcNormFactors(TMM), glmQLFit/Test
library(ggplot2)
library(scales)


LFC_CUT       <- 1          # |log2FC| threshold (Methods: > 1)
FDR_CUT       <- 0.05       # FDR threshold (Methods: < 0.05)
FILTER_METHOD <- "prevalence"  # "prevalence": present in >= MIN_PREV_SAMP of samples

MIN_PREV_SAMP <- 0.10       # 10% prevalence across ALL samples 

TEST          <- "QL"       # "QL" (glmQLFTest, recommended) or "exact" (exactTest)
SEED          <- 1          # edgeR QL is deterministic; kept only for full reproducibility
set.seed(SEED)

habitat_colors <- c(
  "Intertidal-enriched" = "#E64B35FF",   # up in IT
  "Subtidal-enriched"   = "#4DBBD5FF",   # up in ST
  "Not significant"     = "grey80"
)


script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
if (length(script_file) == 1) {
  setwd(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
}


read_auto <- function(file, row.names = NULL) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.table(file, header = TRUE, sep = sep, row.names = row.names,
             check.names = FALSE, quote = "", comment.char = "",
             stringsAsFactors = FALSE)
}


read_paths <- function(file) {
  lines <- readLines(file, warn = FALSE)
  lines <- lines[nzchar(lines)]
  sep   <- if (grepl("\t", lines[1])) "\t" else ","
  body  <- lines[-1]                                   
  pos   <- regexpr(sep, body, fixed = TRUE)            
  ko    <- ifelse(pos > 0, substr(body, 1, pos - 1), body)
  pth   <- ifelse(pos > 0, substr(body, pos + 1, nchar(body)), NA_character_)
  data.frame(KO = trimws(ko), paths = trimws(pth),
             check.names = FALSE, stringsAsFactors = FALSE)
}


counts <- as.matrix(read_auto("kegg_abund.csv", row.names = 1))
storage.mode(counts) <- "numeric"
if (anyNA(counts))                stop("Missing/non-numeric values in kegg_abund.csv")
if (any(counts < 0))              stop("Negative counts in kegg_abund.csv")
if (any(abs(counts - round(counts)) > 1e-6))
  stop("edgeR expects raw integer counts, not normalized abundances.")
counts <- round(counts)

meta  <- read_auto("metadata.csv")
paths <- read_paths("KEGG_paths.csv")
colnames(meta)[1:2] <- c("Sample", "Group")


if (!setequal(colnames(counts), meta$Sample))
  stop("Sample names in kegg_abund.csv and metadata.csv do not match.")
counts <- counts[, meta$Sample, drop = FALSE]


group <- factor(meta$Group, levels = c("ST", "IT"))
message(sprintf("Samples: %d IT, %d ST | KOs (raw): %d",
                sum(group == "IT"), sum(group == "ST"), nrow(counts)))


y <- DGEList(counts = counts, group = group)

if (FILTER_METHOD == "prevalence") {
  keep <- rowMeans(counts > 0) >= MIN_PREV_SAMP    
} else {
  keep <- filterByExpr(y, group = group)           
}
y <- y[keep, , keep.lib.sizes = FALSE]
message(sprintf("KOs retained: %d of %d  (%s filter)",
                sum(keep), length(keep), FILTER_METHOD))

y <- calcNormFactors(y, method = "TMM")          
design <- model.matrix(~ group)                 

if (TEST == "QL") {
  y   <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  res <- topTags(glmQLFTest(fit, coef = 2), n = Inf, sort.by = "none")$table
} else {                                          
  y   <- estimateDisp(y)
  res <- topTags(exactTest(y, pair = c("ST", "IT")), n = Inf, sort.by = "none")$table
}

res$KO <- rownames(res)

# =============================================================================
# Classify + annotate
# =============================================================================
res$direction <- "Not significant"
res$direction[res$FDR < FDR_CUT & res$logFC >  LFC_CUT] <- "Intertidal-enriched"
res$direction[res$FDR < FDR_CUT & res$logFC < -LFC_CUT] <- "Subtidal-enriched"
res$direction <- factor(res$direction, levels = names(habitat_colors))

res <- merge(res, paths[, c("KO", "paths")], by = "KO", all.x = TRUE)

summary_tbl <- data.frame(
  total_tested        = nrow(res),
  intertidal_enriched = sum(res$direction == "Intertidal-enriched"),
  subtidal_enriched   = sum(res$direction == "Subtidal-enriched"),
  significant_total   = sum(res$direction != "Not significant"),
  not_significant     = sum(res$direction == "Not significant")
)
print(summary_tbl)

write.csv(res[order(res$FDR), ], "Figure_S4_DA_results.csv", row.names = FALSE)
write.csv(summary_tbl, "Figure_S4_summary.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Figure_S4_sessionInfo.txt")

# =============================================================================

res$negLog10FDR <- -log10(res$FDR)
fin <- is.finite(res$negLog10FDR)
if (any(!fin)) res$negLog10FDR[!fin] <- max(res$negLog10FDR[fin]) * 1.05


res <- res[order(res$direction != "Not significant"), ]


n_by <- table(res$direction)
legend_labels <- setNames(
  sprintf("%s (%d)", names(habitat_colors), n_by[names(habitat_colors)]),
  names(habitat_colors)
)

p <- ggplot(res, aes(logFC, negLog10FDR, color = direction)) +
  geom_point(size = 0.9, alpha = 0.7) +
  scale_color_manual(values = habitat_colors, labels = legend_labels,
                     drop = FALSE, name = NULL) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = -log10(FDR_CUT), linetype = "dashed",
             color = "grey40", linewidth = 0.3) +
  labs(
    x = expression(log[2]~"fold change (Intertidal / Subtidal)"),
    y = expression(-log[10]~"(FDR)")
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

print(p)
ggsave("Figure_S4_Volcano.png", p, width = 7, height = 6, dpi = 300)
ggsave("Figure_S4_Volcano.pdf", p, width = 7, height = 6)


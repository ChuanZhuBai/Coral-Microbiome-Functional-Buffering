# =============================================================================
# exploratory co-occurrence networks of KO functions  (IT-H vs ST-H, n=6)
# =============================================================================
suppressMessages({ library(tidyverse); library(igraph) })

PREV_IN_GROUP <- 0.5
RHO_CUT       <- 0.60
P_CUT         <- 0.05
USE_FDR       <- FALSE
TOP_N_NODES   <- 50
TOP_N         <- 20
VERIFY_TOP    <- 150     
set.seed(1)


PANEL <- c("eco","bsu","sau","pae","mtu","syn","bth","hpy","vch","cje","son",
           "gsu","dra","spn","ype","nme","atu","sco","cgl","ctr","tma","bbu",
           "rsp","lpl"                                  
           )    

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

read_paths_annot <- function(file) {
  raw <- read_auto(file); pcol <- setdiff(names(raw), "KO")[1]
  raw %>%
    transmute(KO = str_extract(KO, "K\\d+"), paths = .data[[pcol]]) %>%
    separate_rows(paths, sep = "\\s*\\|\\s*") %>%
    mutate(parts = str_split(paths, ";\\s*"),
           L2 = map_chr(parts, ~ if (length(.x) >= 2) str_trim(.x[2]) else NA_character_),
           L3 = map_chr(parts, ~ if (length(.x) >= 3) str_trim(.x[3]) else NA_character_)) %>%
    filter(!is.na(L2)) %>% distinct(KO, L2, L3) %>%
    group_by(KO) %>% summarise(L2 = paste(unique(L2), collapse="; "),
                               L3 = paste(unique(L3[!is.na(L3)]), collapse="; "), .groups="drop")
}

is_unannotated <- function(x) grepl("uncl|unass|unknown|^k$|^-$|^$", x, ignore.case = TRUE)


counts <- as.matrix(read_auto("ko_abund.csv", row.names = 1)); storage.mode(counts) <- "numeric"
colnames(counts) <- canonical_id(colnames(counts)); counts <- round(counts)
meta <- read_auto("metadata.csv"); colnames(meta)[1:2] <- c("Sample","Group")
meta$Sample <- canonical_id(meta$Sample)
annot <- read_paths_annot("KEGG_paths.csv")

if (!requireNamespace("vegan", quietly = TRUE)) stop("install.packages('vegan') needed.")


exp_all <- intersect(meta$Sample[meta$Group %in% c("IT-H","ST-H")], colnames(counts))
depth_min <- min(colSums(counts[, exp_all, drop = FALSE]))
message(sprintf("Rarefying experimental samples to common depth = %d (full table, incl. Unclassified)", depth_min))
set.seed(1)
counts[, exp_all] <- t(vegan::rrarefy(t(counts[, exp_all, drop = FALSE]), sample = depth_min))


counts <- counts[!is_unannotated(rownames(counts)), , drop = FALSE]

# =============================================================================
# EUKARYOTE REMOVAL via a KEGG prokaryote KO reference set (keggLink panel)
# =============================================================================
if (!requireNamespace("KEGGREST", quietly = TRUE))
  stop("KEGGREST needed: install.packages('BiocManager'); BiocManager::install('KEGGREST')")

PROK_SET_TXT   <- "KEGG_prok_KO_set.txt"
PROK_CHECK_CSV <- "Figure_3F_KO_prokaryote_check.csv"


build_prok_set <- function() {
  if (file.exists(PROK_SET_TXT)) {
    P <- readLines(PROK_SET_TXT)
    message(sprintf("Loaded cached prokaryote KO set: %d KOs (delete %s to rebuild)",
                    length(P), PROK_SET_TXT))
    return(P)
  }
  got <- character(0); ok <- 0L
  for (o in PANEL) {
    v <- tryCatch(unique(sub("ko:", "", KEGGREST::keggLink("ko", o))),
                  error = function(e) character(0))
    if (length(v)) { got <- union(got, v); ok <- ok + 1L
    message(sprintf("   %s: %d KOs (cumulative %d)", o, length(v), length(got))) }
    Sys.sleep(0.1)
  }
  message(sprintf("Prokaryote panel: %d/%d organisms OK -> %d KOs in reference set",
                  ok, length(PANEL), length(got)))
  if (length(got) < 3000)
    warning("Prokaryote reference set looks small (<3000 KOs); consider expanding PANEL.")
  writeLines(got, PROK_SET_TXT)
  got
}
P <- build_prok_set()
stopifnot(length(P) > 1000)


top_abund_kos <- function(grp, K) {
  s   <- intersect(meta$Sample[meta$Group == grp], colnames(counts))
  cnt <- counts[, s, drop = FALSE]
  cnt <- cnt[rowSums(cnt > 0) >= ceiling(PREV_IN_GROUP * length(s)), , drop = FALSE]
  rel <- sweep(cnt, 2, colSums(cnt), "/")
  names(sort(rowMeans(rel), decreasing = TRUE))[seq_len(min(K, nrow(cnt)))]
}
cand <- union(top_abund_kos("IT-H", VERIFY_TOP), top_abund_kos("ST-H", VERIFY_TOP))
in_P <- cand %in% P

KEEP_KOS   <- c("K11180")               
flagged    <- cand[!in_P]
EUK_REMOVE <- setdiff(flagged, KEEP_KOS)
message(sprintf("Candidates: %d | flagged not-in-panel: %d | kept by override: %d | removed: %d",
                length(cand), length(flagged), length(intersect(flagged, KEEP_KOS)), length(EUK_REMOVE)))


enrich <- function(kos) {
  empty <- tibble(KO = character(), symbol = character(), name = character())
  if (!length(kos)) return(empty)
  fetch <- function(ids) tryCatch(KEGGREST::keggGet(ids), error = function(e) NULL)
  recs <- list()
  for (g in split(kos, ceiling(seq_along(kos) / 10))) {
    rr <- fetch(g)
    if (is.null(rr)) rr <- unlist(lapply(g, fetch), recursive = FALSE)  # 1-by-1 on chunk failure
    if (!is.null(rr)) for (r in rr) {
      if (is.null(r$ENTRY)) next
      ko <- sub("\\s.*$", "", unname(r$ENTRY))
      recs[[ko]] <- tibble(
        KO     = ko,
        symbol = if (!is.null(r$SYMBOL)) paste(r$SYMBOL, collapse = ", ") else NA_character_,
        name   = if (!is.null(r$NAME))   paste(r$NAME,   collapse = "; ")  else NA_character_)
    }
    Sys.sleep(0.1)
  }
  if (!length(recs)) empty else bind_rows(recs)   # always 3 columns, never 0-col
}
names_tab <- tryCatch(enrich(flagged),
                      error = function(e) tibble(KO = character(), symbol = character(), name = character()))
audit <- tibble(KO = cand, in_prokaryote = in_P, removed = cand %in% EUK_REMOVE) %>%
  left_join(names_tab, by = "KO") %>% arrange(desc(removed), in_prokaryote, KO)
write.csv(audit, PROK_CHECK_CSV, row.names = FALSE)
if (length(EUK_REMOVE)) message("  removed: ", paste(EUK_REMOVE, collapse = ", "))

counts <- counts[!rownames(counts) %in% EUK_REMOVE, , drop = FALSE]
stopifnot(!"K16196" %in% rownames(counts))          # canary: EIF2AK4/GCN2 must be gone
message(sprintf("KOs remaining after eukaryote removal: %d", nrow(counts)))

if (!requireNamespace("psych", quietly = TRUE)) stop("install.packages('psych') required.")


analyse_group <- function(grp) {
  samples <- intersect(meta$Sample[meta$Group == grp], colnames(counts))
  message(sprintf("[%s] %d samples", grp, length(samples)))
  cnt <- counts[, samples, drop = FALSE]
  cnt <- cnt[!is_unannotated(rownames(cnt)), , drop = FALSE]
  keep <- rowSums(cnt > 0) >= ceiling(PREV_IN_GROUP * length(samples))
  cnt <- cnt[keep, , drop = FALSE]
  
  rel_full <- sweep(cnt, 2, colSums(cnt), "/")
  ord <- order(rowMeans(rel_full), decreasing = TRUE)
  cnt <- cnt[head(ord, TOP_N_NODES), , drop = FALSE]
  message(sprintf("   prevalent KOs: %d -> top %d by abundance used for network",
                  length(ord), nrow(cnt)))
  rel <- sweep(cnt, 2, colSums(cnt), "/")
  mean_rel <- rowMeans(rel) * 100
  
  ct <- psych::corr.test(t(rel), method = "spearman", adjust = "none", ci = FALSE)
  R <- ct$r; P <- ct$p
  if (USE_FDR) {
    P[upper.tri(P)] <- p.adjust(P[upper.tri(P)], "BH")
    P[lower.tri(P)] <- t(P)[lower.tri(P)]
  }
  adj <- (abs(R) > RHO_CUT) & (P < P_CUT); diag(adj) <- FALSE
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
  if (ecount(g) > 0) {
    el <- as_edgelist(g, names = FALSE)
    E(g)$rho  <- R[el]
    E(g)$sign <- ifelse(E(g)$rho > 0, "positive", "negative")
    E(g)$correlation_type     <- E(g)$sign
    E(g)$absolute_correlation <- abs(E(g)$rho)
  }
  g <- delete_vertices(g, degree(g) == 0)
  message(sprintf("   network: %d nodes, %d edges (rho>%.2f, %s p<%.2f)",
                  vcount(g), ecount(g), RHO_CUT, if (USE_FDR) "BH-FDR" else "raw", P_CUT))
  
  hubs <- NULL
  if (vcount(g) > 0) {
    V(g)$degree   <- degree(g)
    V(g)$mean_rel <- mean_rel[V(g)$name]
    V(g)$L2 <- annot$L2[match(V(g)$name, annot$KO)]; V(g)$L2[is.na(V(g)$L2)] <- "Unannotated"
    V(g)$L3 <- annot$L3[match(V(g)$name, annot$KO)]; V(g)$L3[is.na(V(g)$L3)] <- ""
    hubs <- tibble(KO=V(g)$name, degree=V(g)$degree, mean_rel=V(g)$mean_rel,
                   L2=V(g)$L2, L3=V(g)$L3) %>% arrange(desc(degree)) %>% slice_head(n=TOP_N)
    write.csv(hubs, sprintf("Figure_3F_hubs_%s.csv", grp), row.names = FALSE)
    V(g)$id <- seq_len(vcount(g)); V(g)$label <- V(g)$name
    write_graph(g, sprintf("Figure_3F_network_%s.graphml", grp), format = "graphml")
    tryCatch(write_graph(g, sprintf("Figure_3F_network_%s.gml", grp), format = "gml"),
             error = function(e) message("   (GML skipped: ", conditionMessage(e),
                                         "); use .graphml"))
  }
  list(grp=grp, graph=g, hubs=hubs)
}

res_ITH <- analyse_group("IT-H")
res_STH <- analyse_group("ST-H")

# guard: every network node was verified prokaryote (else raise VERIFY_TOP)
allnodes <- unique(c(if(!is.null(res_ITH$graph)) V(res_ITH$graph)$name,
                     if(!is.null(res_STH$graph)) V(res_STH$graph)$name))
notP <- setdiff(allnodes, P)
if (length(notP))
  warning("These network nodes are not in the prokaryote set -- increase VERIFY_TOP: ",
          paste(notP, collapse = ", "))

if (!is.null(res_ITH$hubs)) { cat("\n== IT-H hubs ==\n"); print(res_ITH$hubs, n=TOP_N) }
if (!is.null(res_STH$hubs)) { cat("\n== ST-H hubs ==\n"); print(res_STH$hubs, n=TOP_N) }


plot_net <- function(res, base_col) {
  g <- res$graph
  if (is.null(g) || vcount(g)==0) { message("No network for ", res$grp); return(invisible()) }
  ramp <- colorRampPalette(c("grey85", base_col))(100)
  v <- V(g)$mean_rel; v[is.na(v)] <- 0
  V(g)$color <- ramp[as.integer(scales::rescale(rank(v), to=c(1,100)))]
  V(g)$size  <- scales::rescale(V(g)$degree, to=c(3,12))
  hub_kos <- res$hubs$KO[1:min(8, nrow(res$hubs))]
  V(g)$label <- ifelse(V(g)$name %in% hub_kos, V(g)$name, NA)
  E(g)$color <- ifelse(E(g)$sign=="positive", "#D6604D66", "#4393C366")
  pdf(sprintf("Figure_3F_network_%s_preview.pdf", res$grp), width=7, height=7)
  set.seed(1)
  plot(g, vertex.frame.color=NA, vertex.label.cex=0.6, vertex.label.color="black",
       layout=layout_with_fr,
       main=sprintf("%s (n=6, EXPLORATORY): %d nodes, %d edges | |rho|>%.2f, raw p<%.2f",
                    res$grp, vcount(g), ecount(g), RHO_CUT, P_CUT))
  dev.off()
}
plot_net(res_ITH, "#E41A1C"); plot_net(res_STH, "#377EB8")


topo <- function(res) {
  g <- res$graph
  if (is.null(g) || vcount(g)==0) return(tibble(Group=res$grp, nodes=0, edges=0,
                                                avg_degree=NA, modularity=NA, density=NA, pos_edge_frac=NA))
  tibble(Group=res$grp, nodes=vcount(g), edges=ecount(g), avg_degree=mean(degree(g)),
         modularity=modularity(cluster_fast_greedy(as.undirected(g))),
         density=edge_density(g), pos_edge_frac=mean(E(g)$sign=="positive"))
}
topo_tab <- bind_rows(topo(res_ITH), topo(res_STH))
write.csv(topo_tab, "Figure_3F_network_topology.csv", row.names = FALSE)
print(topo_tab)
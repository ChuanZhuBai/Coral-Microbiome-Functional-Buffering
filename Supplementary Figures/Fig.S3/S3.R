# =============================================================================
# Fig. S3 / Table S4 — Genus-level co-occurrence networks (field: IT vs ST)
# =============================================================================


library(igraph)  
library(Hmisc)    


MIN_ABUND        <- 1e-4      
                              
MIN_PREV         <- 0.20      
R_CUT            <- 0.6       
P_CUT            <- 0.01      
KEEP_SIGN        <- "both"    
WEIGHT_MODE      <- "abs"     
                              
COMMUNITY_METHOD <- "louvain" 
SEED             <- 1

set.seed(SEED)


read_abund <- function(file) {
  m <- read.csv(file, header = TRUE, row.names = 1, check.names = FALSE)
  m <- as.matrix(m)
  storage.mode(m) <- "numeric"
  if (anyNA(m))      stop("Missing/non-numeric values in ", file)
  if (any(m < 0))    stop("Negative abundances in ", file)
  
  cs <- colSums(m)
  message(sprintf("%s: %d taxa x %d samples | colSums range [%.3g, %.3g]",
                  file, nrow(m), ncol(m), min(cs), max(cs)))
  m
}

filter_abund_prev <- function(m, min_abund, min_prev) {
  prev <- rowMeans(m > 0)
  keep <- rowMeans(m) > min_abund & prev > min_prev   
  m[keep, , drop = FALSE]
}


build_network <- function(otu,
                          r_cut = R_CUT, p_cut = P_CUT,
                          keep_sign = KEEP_SIGN, weight_mode = WEIGHT_MODE,
                          community = COMMUNITY_METHOD) {
  taxa <- rownames(otu)
  n    <- length(taxa)

  sp <- rcorr(t(otu), type = "spearman")
  R  <- sp$r
  P  <- sp$P
  diag(R) <- 0

  
  padj <- matrix(NA_real_, n, n, dimnames = list(taxa, taxa))
  ut <- upper.tri(P)
  padj[ut] <- p.adjust(P[ut], method = "BH")       
  padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]

  
  rmag <- if (keep_sign == "positive") R else abs(R)
  sel  <- (rmag > r_cut) & (padj < p_cut)
  sel[is.na(sel)] <- FALSE          
  diag(sel) <- FALSE

  
  wmat <- if (weight_mode == "abs") abs(R) else R
  A <- matrix(0, n, n, dimnames = list(taxa, taxa))
  A[sel] <- wmat[sel]

  
  keep <- rowSums(sel) > 0
  A <- A[keep, keep, drop = FALSE]

  g <- graph_from_adjacency_matrix(A, mode = "undirected",
                                   weighted = TRUE, diag = FALSE)
  g <- simplify(g, edge.attr.comb = "first") 

  
  if (ecount(g) > 0) {
    ends_m <- ends(g, E(g))
    E(g)$rho <- mapply(function(a, b) R[a, b], ends_m[, 1], ends_m[, 2])
  }

  
  comm <- if (community == "louvain")
    cluster_louvain(g, weights = E(g)$weight)
  else
    cluster_walktrap(g, weights = E(g)$weight)

  V(g)$community <- membership(comm)
  V(g)$degree    <- degree(g)
  V(g)$abundance <- rowMeans(otu[V(g)$name, , drop = FALSE])
  attr(g, "comm") <- comm
  g
}


compute_topo <- function(g) {
  comm  <- attr(g, "comm")
  comps <- components(g)$csize
  data.frame(
    nodes         = vcount(g),
    edges         = ecount(g),
    avg_degree    = mean(degree(g)),
    density       = edge_density(g),
    modularity    = modularity(comm),
    clustering    = transitivity(g, type = "average"),
    avg_path_len  = mean_distance(g, directed = FALSE, weights = NA),  
    lcc_rel_size  = max(comps) / vcount(g),
    n_modules     = length(unique(membership(comm))),
    row.names     = NULL
  )
}


otu_it <- filter_abund_prev(read_abund("IT.csv"), MIN_ABUND, MIN_PREV)
otu_st <- filter_abund_prev(read_abund("ST.csv"), MIN_ABUND, MIN_PREV)

net_it <- build_network(otu_it)
net_st <- build_network(otu_st)

write_graph(net_it, "IT_network.gml", format = "gml")
write_graph(net_st, "ST_network.gml", format = "gml")

topo_df <- rbind(IT = compute_topo(net_it),
                 ST = compute_topo(net_st))
write.csv(topo_df, "network_topology_metrics.csv")
print(topo_df)

writeLines(capture.output(sessionInfo()), "Figure_S3_sessionInfo.txt")


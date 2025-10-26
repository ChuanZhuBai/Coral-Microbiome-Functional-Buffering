# === 0. 载入包 ===
library(igraph)
library(Hmisc)
library(ggplot2)
dir()
# === 1. 数据导入与过滤 ===
# 读入 IT 组或 ST 组相对丰度表（行：属，列：样本）
otu = read.csv("IT.csv", header = TRUE, row.names = 1, check.names = FALSE)
otu = as.matrix(otu)

# 过滤：平均相对丰度 > 0.01%
min_abund = 1e-4
otu = otu[rowMeans(otu) > min_abund, ]

# === 2. 相关性与 FDR 校正 ===
sp.cor = rcorr(t(otu), type = "spearman")
r.cor  = sp.cor$r
p.cor  = sp.cor$P

get_fdr_mat <- function(pval_mat) {
  idx <- upper.tri(pval_mat)
  pvals <- pval_mat[idx]  # 转换为向量
  padj_vals <- p.adjust(pvals, method = "fdr")
  padj <- matrix(0, nrow(pval_mat), ncol(pval_mat))
  padj[idx] <- padj_vals
  padj[lower.tri(padj)] <- t(padj)[lower.tri(padj)]
  diag(padj) <- 0
  padj
}
p.adj = get_fdr_mat(p.cor)

# === 3. 网络构建（只保留正相关） ===
r.cut = 0.6; p.cut = 0.05
r.mat = r.cor
r.mat[r.cor <= r.cut] = 0
r.mat[p.adj > p.cut]  = 0

# 去除孤立节点
keep = which(rowSums(r.mat != 0) > 0)
r.mat = r.mat[keep, keep]

# igraph 对象
net_pos = graph_from_adjacency_matrix(r.mat, mode = "undirected", weighted = TRUE, diag = FALSE)
net_pos = simplify(net_pos)

# 输出 gml
write_graph(net_pos, "IT_network_pos.gml", format = "gml")
write_graph(net_pos, "ST_network_pos.gml", format = "gml")
# 社区与节点属性
comm = cluster_walktrap(net_pos)
V(net_pos)$community = comm$membership
V(net_pos)$degree    = degree(net_pos)

# === 4. 随机攻击稳定性模拟函数 ===
simulate_attack <- function(graph, prop_seq = seq(0, 1, by = 0.05), nrep = 50) {
  N = vcount(graph)
  res = data.frame()
  for (p in prop_seq) {
    # 多次随机删除
    comps = replicate(nrep, {
      g2 = delete_vertices(graph, sample(V(graph), size = floor(p * N)))
      max(components(g2)$csize) / N
    })
    res = rbind(res,
                data.frame(
                  prop_removed  = p,
                  mean_lcc_size = mean(comps),
                  sd_lcc_size   = sd(comps)
                )
    )
  }
  res
}

# 分别对 IT 和 ST 网络模拟
it_attack = simulate_attack(net_pos) 
# 读 ST 数据并重复步骤 1–3，构建 net_pos_st，然后：
# otu_st = read.csv("ST.csv", ...); ...; net_pos_st = ...
st_attack = simulate_attack(net_pos)

it_attack$Group = "IT"
st_attack$Group = "ST"
stab_df = rbind(it_attack, st_attack)

# === 5. 绘制稳定性曲线图 ===
p <- ggplot(stab_df, aes(x = prop_removed, y = mean_lcc_size, color = Group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  geom_ribbon(aes(ymin = mean_lcc_size - sd_lcc_size,
                  ymax = mean_lcc_size + sd_lcc_size,
                  fill = Group),
              alpha = 0.2, color = NA) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "移除节点比例", y = "最大连通分量相对大小",
       title = "Panel C：微生物网络随机攻击稳定性",
       color = NULL, fill = NULL) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title  = element_text(size = 12),
    legend.position = c(0.8, 0.8)
  )
print(p)

# === 6. Gephi 可视化提示 ===
# - 在 Gephi 中导入 IT_network_pos.gml / ST_network_pos.gml
# - 节点颜色映射到 community，节点大小映射到 degree
# - 布局可选 ForceAtlas2，标签隐藏或只显示高 degree 节点

# End of Panel C























# =============================================================================
# Panel C: 微生物共现网络构建与随机/定向攻击模拟 — 全面优化脚本
# High‐Level SCI Standard Co‐occurrence Network & Robustness Analysis in R
# =============================================================================

# === 0. 载入必要 R 包 ===
library(igraph)       # network analysis
library(Hmisc)        # rcorr()
library(ggplot2)      # plotting
library(RColorBrewer) # color palettes
library(pracma)       # trapz() for AUC

# =============================================================================
# 1. 数据导入与预处理
# =============================================================================

# 1.1 读入 IT & ST 相对丰度表（行：属/OTU，列：样本）
otu_it <- read.csv("IT.csv", header = TRUE, row.names = 1, check.names = FALSE)
otu_st <- read.csv("ST.csv", header = TRUE, row.names = 1, check.names = FALSE)
otu_it <- as.matrix(otu_it)
otu_st <- as.matrix(otu_st)

# 1.2 联合过滤：平均丰度 & 检出率
filter_abund_prev <- function(otu_mat, min_abund = 1e-4, min_prev = 0.2) {
  prev <- rowMeans(otu_mat > 0)
  keep <- which(rowMeans(otu_mat) > min_abund & prev >= min_prev)
  otu_mat[keep, , drop = FALSE]
}

otu_it <- filter_abund_prev(otu_it, min_abund = 1e-4, min_prev = 0.2)
otu_st <- filter_abund_prev(otu_st, min_abund = 1e-4, min_prev = 0.2)

# =============================================================================
# 2. 共现网络构建函数（Spearman + FDR + 阈值 + Consensus Bootstrap）
# =============================================================================

build_consensus_network <- function(otu_mat,
                                    r.cut = 0.6, p.cut = 0.05,
                                    nboot = 500, cons_thresh = 0.8) {
  taxa <- rownames(otu_mat)
  ntax <- length(taxa)
  
  # 2.1 bootstrap 生成共识边计数矩阵
  edge_counts <- matrix(0, ntax, ntax, dimnames = list(taxa, taxa))
  
  for (b in 1:nboot) {
    samp_cols <- sample(ncol(otu_mat), replace = TRUE)
    mat_boot   <- otu_mat[, samp_cols]
    sp         <- rcorr(t(mat_boot), type = "spearman")
    Rb         <- sp$r
    Pb         <- sp$P
    
    # FDR 校正
    padjb <- matrix(0, ntax, ntax, dimnames = list(taxa, taxa))
    idx   <- upper.tri(Pb)
    padjb[idx] <- p.adjust(Pb[idx], method = "fdr")
    padjb[lower.tri(padjb)] <- t(padjb)[lower.tri(padjb)]
    
    # 应用阈值，提取强正相关边
    sel <- (Rb > r.cut) & (padjb < p.cut)
    edge_counts[sel] <- edge_counts[sel] + 1
  }
  
  # 2.2 构建共识矩阵（只保留 ≥ cons_thresh 的边）
  consensus <- (edge_counts / nboot) >= cons_thresh
  R0 <- matrix(0, ntax, ntax, dimnames = list(taxa, taxa))
  
  # 再计算全数据相关度，用以赋予权重
  sp_full <- rcorr(t(otu_mat), type = "spearman")
  Rfull   <- sp_full$r
  
  R0[consensus] <- Rfull[consensus]
  diag(R0) <- 0
  
  # 2.3 删除孤立节点
  keep <- which(rowSums(R0 != 0) > 0)
  R0 <- R0[keep, keep]
  
  # 2.4 转换为 igraph 对象
  g <- graph_from_adjacency_matrix(R0, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- simplify(g)
  
  # 2.5 社区检测与节点属性
  comm <- cluster_walktrap(g)
  V(g)$community <- comm$membership
  V(g)$degree    <- degree(g)
  V(g)$abundance <- rowMeans(otu_mat[V(g)$name, ])
  
  return(g)
}

# 构建 IT & ST 网络
net_it <- build_consensus_network(otu_it,
                                  r.cut = 0.6, p.cut = 0.05,
                                  nboot = 500, cons_thresh = 0.8)
net_st <- build_consensus_network(otu_st,
                                  r.cut = 0.6, p.cut = 0.05,
                                  nboot = 500, cons_thresh = 0.8)

# 导出 GML 供 Gephi 进一步美化
write_graph(net_it, "IT_consensus_network.gml", format = "gml")
write_graph(net_st, "ST_consensus_network.gml", format = "gml")

# =============================================================================
# 3. 网络拓扑指标计算（用于论文表格/结果）
# =============================================================================

compute_topo <- function(g) {
  comps <- components(g)$csize
  data.frame(
    nodes        = vcount(g),
    edges        = ecount(g),
    density      = edge_density(g),
    avg_degree   = mean(degree(g)),
    clustering   = transitivity(g, type = "average"),
    modularity   = modularity(cluster_walktrap(g)),
    lcc_rel_size = max(comps) / vcount(g)
  )
}

topo_it <- compute_topo(net_it)
topo_st <- compute_topo(net_st)
topo_df <- rbind(IT = topo_it, ST = topo_st)
write.csv(topo_df, "network_topology_metrics.csv")

# =============================================================================
# 4. 随机 vs 定向攻击稳定性模拟
# =============================================================================

simulate_attack <- function(graph, prop_seq = seq(0,1,by=0.05),
                            nrep = 100, targeted = FALSE) {
  N   <- vcount(graph)
  deg <- degree(graph)
  out <- data.frame()
  
  for (p in prop_seq) {
    lccs <- replicate(nrep, {
      if (targeted) {
        # 节点度越高 => 优先被删除（权重采样）
        rmv <- sample(V(graph), size = floor(p * N), prob = deg^2)
      } else {
        rmv <- sample(V(graph), size = floor(p * N))
      }
      g2 <- delete_vertices(graph, rmv)
      if (vcount(g2) == 0) return(0)
      max(components(g2)$csize) / N
    })
    out <- rbind(out, data.frame(
      prop_removed  = p,
      mean_lcc      = mean(lccs),
      sd_lcc        = sd(lccs),
      type          = ifelse(targeted, "Targeted", "Random")
    ))
  }
  out
}

# IT 网络模拟
it_rand   <- simulate_attack(net_it, targeted = FALSE)
it_target <- simulate_attack(net_it, targeted = TRUE)
it_all    <- rbind(it_rand, it_target); it_all$Group <- "IT"

# ST 网络模拟
st_rand   <- simulate_attack(net_st, targeted = FALSE)
st_target <- simulate_attack(net_st, targeted = TRUE)
st_all    <- rbind(st_rand, st_target); st_all$Group <- "ST"

stab_df <- rbind(it_all, st_all)
write.csv(stab_df, "network_stability_simulation.csv", row.names = FALSE)

# =============================================================================
# 5. 鲁棒性 AUC 计算与置换检验
# =============================================================================

calculate_auc <- function(df) {
  df_split <- split(df, df$type)
  sapply(df_split, function(d) trapz(d$prop_removed, d$mean_lcc))
}

auc_it    <- calculate_auc(subset(stab_df, Group=="IT"))
auc_st    <- calculate_auc(subset(stab_df, Group=="ST"))
# 可按 type 分别输出 IT_random, IT_targeted, ST_random, ST_targeted

# 置换检验示例：比较 IT vs ST 随机攻击 AUC 差异
it_vals <- subset(stab_df, Group=="IT" & type=="Random")$mean_lcc
st_vals <- subset(stab_df, Group=="ST" & type=="Random")$mean_lcc
diff0   <- trapz(unique(stab_df$prop_removed), it_vals) -
  trapz(unique(stab_df$prop_removed), st_vals)

combined_vals <- c(it_vals, st_vals)
nperm <- 999
perm_diffs <- replicate(nperm, {
  samp <- sample(combined_vals)
  auc1 <- trapz(unique(stab_df$prop_removed), samp[1:length(it_vals)])
  auc2 <- trapz(unique(stab_df$prop_removed), samp[(length(it_vals)+1):length(samp)])
  auc1 - auc2
})
pval <- (sum(abs(perm_diffs) >= abs(diff0)) + 1) / (nperm + 1)

# =============================================================================
# 6. 绘图：网络稳定性（随机 vs 定向，IT vs ST）
# =============================================================================

p <- ggplot(stab_df, aes(x = prop_removed, y = mean_lcc, 
                         color = Group, linetype = type, fill = Group)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = mean_lcc - sd_lcc, ymax = mean_lcc + sd_lcc),
              alpha = 0.2, colour = NA) +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  scale_x_continuous(labels = scales::percent, limits = c(0,1)) +
  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
  labs(x = "Node removal ratio",
       y = "Relative size of the largest connected component",
       title = "Stability of random attacks on microbial networks",
       subtitle = sprintf("IT_random AUC=%.3f, ST_random AUC=%.3f, P=%.3f",
                          auc_it["Random"], auc_st["Random"], pval)) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        plot.title    = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5))
print(p)

# 保存图形
ggsave("PanelC_network_stability.pdf", plot = p, width = 8, height = 6, dpi = 300)

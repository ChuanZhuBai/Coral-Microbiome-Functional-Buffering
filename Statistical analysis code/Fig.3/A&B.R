# 加载必要包
library(tidyverse)
library(vegan)
library(metagenomeSeq)
library(ggplot2)
library(patchwork)

# ==========================================
# 1. 辅助函数定义 (Helper Functions)
# ==========================================

# A. 自动 ID 转换函数
rename_samples <- function(names_vec) {
  names_vec %>%
    gsub("restI([0-9]+)A", "IT-C\\1", .) %>%
    gsub("restI([0-9]+)H", "IT-H\\1", .) %>%
    gsub("restS([0-9]+)A", "ST-C\\1", .) %>%
    gsub("restS([0-9]+)H", "ST-H\\1", .)
}

# B. 统计全套分析函数 (adonis2 + Pairwise)
run_stats_logic <- function(dist_mat, env) {
  # 1. Global PERMANOVA
  set.seed(123)
  adonis_glob <- adonis2(dist_mat ~ Group, data = env, permutations = 999)
  
  # 2. Pairwise PERMANOVA (手动对比，保证最严谨)
  groups <- levels(factor(env$Group))
  pairwise_res <- list()
  combn_matrix <- combn(groups, 2)
  
  for(i in 1:ncol(combn_matrix)) {
    g1 <- combn_matrix[1, i]
    g2 <- combn_matrix[2, i]
    sub_env <- env[env$Group %in% c(g1, g2), ]
    # 提取子距离矩阵
    sub_dist <- as.dist(as.matrix(dist_mat)[rownames(sub_env), rownames(sub_env)])
    res <- adonis2(sub_dist ~ Group, data = sub_env, permutations = 999)
    pairwise_res[[i]] <- data.frame(Comparison = paste(g1, "vs", g2), 
                                    R2 = res$R2[1], 
                                    P = res$`Pr(>F)`[1])
  }
  return(list(Global = adonis_glob, Pairwise = do.call(rbind, pairwise_res)))
}

# ==========================================
# 2. 数据载入与元数据准备
# ==========================================
metadata <- read.csv("metadata.csv")
metadata$Group <- factor(metadata$Group, levels = c("IT-C", "IT-H", "ST-C", "ST-H"))
rownames(metadata) <- metadata$Sample

# ==========================================
# 3. 属水平分析 (Relative Abundance 方案)
# ==========================================
genus_raw <- read.csv("genus.csv", row.names = 1, check.names = FALSE)
colnames(genus_raw) <- rename_samples(colnames(genus_raw))

# 筛选实验组并过滤
exp_genus <- genus_raw[, metadata$Sample]
genus_rel <- sweep(exp_genus, 2, colSums(exp_genus), "/")
genus_rel <- genus_rel[rowSums(genus_rel) > 0, ]
write.csv(genus_rel, "Experimental_Genus_RelAbund_Final.csv")

# 统计
dist_gen <- vegdist(t(genus_rel), method = "bray")
stats_gen <- run_stats_logic(dist_gen, metadata)
cpcoa_gen <- capscale(dist_gen ~ Group, data = metadata)

# ==========================================
# 4. KO 水平分析 (CSS Normalization 方案)
# ==========================================
ko_raw <- read.csv("ko_bases.csv", row.names = 1, check.names = FALSE)
# 预清洗
ko_cleaned <- ko_raw[!tolower(rownames(ko_raw)) %in% c("unclassified", "unassigned"), ]
colnames(ko_cleaned) <- rename_samples(colnames(ko_cleaned))
exp_ko <- ko_cleaned[, metadata$Sample]
exp_ko <- exp_ko[rowSums(exp_ko) > 0, ]

# CSS 归一化 (不进行相对丰度二次转化，保持 CSS 原始逻辑)
cat("正在执行 KO 层的 CSS 归一化...\n")
obj <- newMRexperiment(exp_ko)
p <- cumNormStatFast(obj) # 自动计算分位数
obj <- cumNorm(obj, p = p)
ko_css <- as.data.frame(MRcounts(obj, norm = TRUE, log = FALSE))
write.csv(ko_css, "Experimental_KO_CSS_Normalized_Final.csv")

# 统计
dist_ko <- vegdist(t(ko_css), method = "bray")
stats_ko <- run_stats_logic(dist_ko, metadata)
cpcoa_ko <- capscale(dist_ko ~ Group, data = metadata)

# ==========================================
# 5. 绘图与结果导出
# ==========================================

# 统一配色：IT红色系，ST蓝色系
color_pal <- c("IT-C" = "#FC8D62", "IT-H" = "#E41A1C", "ST-C" = "#80B1D3", "ST-H" = "#377EB8")

make_publication_plot <- function(mod, title, palette) {
  # 提取坐标
  plot_data <- as.data.frame(scores(mod, display = "sites"))
  plot_data$Group <- metadata$Group
  
  # 提取方差解释度
  eig <- summary(mod)$concont$importance
  cap1 <- round(eig[2,1] * 100, 1)
  cap2 <- round(eig[2,2] * 100, 1)
  total_expl <- round((eig[2,1] + eig[2,2]) * 100, 1)
  
  # 提取 anova.cca 的 P 值
  p_val <- anova.cca(mod)$`Pr(>F)`[1]
  
  ggplot(plot_data, aes(x = CAP1, y = CAP2, color = Group, fill = Group)) +
    stat_ellipse(geom = "polygon", alpha = 0.1, linetype = 2) +
    geom_point(aes(shape = Group), size = 4, stroke = 0.8) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    scale_shape_manual(values = c(16, 16, 17, 17)) +
    labs(title = title, 
         subtitle = paste0("Constrained variance: ", total_expl, "% (p = ", p_val, ")"),
         x = paste0("CPCoA 1 (", cap1, "%)"),
         y = paste0("CPCoA 2 (", cap2, "%)")) +
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          plot.title = element_text(face="bold", hjust=0.5, size=14),
          legend.title = element_text(face="bold"))
}

p_gen <- make_publication_plot(cpcoa_gen, "Bacterial Community Structure (Genus)", color_pal)
p_ko  <- make_publication_plot(cpcoa_ko, "Functional Structure (KO)", color_pal)

# 打印查看
print(p_gen)
print(p_ko)

# 导出所有报表
write.csv(stats_gen$Pairwise, "Table_S_Pairwise_PERMANOVA_Genus.csv", row.names = FALSE)
write.csv(stats_ko$Pairwise, "Table_S_Pairwise_PERMANOVA_KO.csv", row.names = FALSE)

# 导出最终汇总日志
sink("Full_Statistical_Summary_Report.txt")
cat("=== GENUS LEVEL GLOBAL PERMANOVA (adonis2) ===\n")
print(stats_gen$Global)
cat("\n\n=== KO LEVEL GLOBAL PERMANOVA (adonis2) ===\n")
print(stats_ko$Global)
sink()

# 拼图并保存
final_fig <- p_gen + p_ko + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
ggsave("Figure_3_Combined_CPCoA_Final.pdf", final_fig, width = 10, height = 5)

cat("所有分析圆满完成，结果已存至当前目录。\n")

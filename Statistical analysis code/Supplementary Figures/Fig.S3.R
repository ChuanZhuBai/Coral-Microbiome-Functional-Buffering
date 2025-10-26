keep <- filterByExpr(
  dge,
  group = dge$samples$group,
  min.count = 5,
  min.total.count = 20,
  large.n = 10,
  min.prop = 0.2
)
dge <- dge[keep, , keep.lib.sizes = FALSE]


##############################################
############################################################################
# ==================== 1. 加载所有必要包 ====================
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)  # 用于数据整理

# ==================== 2. 数据加载与预处理 ====================
# 读取原始counts数据
counts_data <- read.csv("KEGG_bases_insitu.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv", header = TRUE)

# 读取通路数据（确保正确分隔）
kegg_paths <- read.csv("KEGG_paths.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# 检查列名并处理（如果需要）
if("KO.paths" %in% colnames(kegg_paths)) {
  kegg_paths <- separate(kegg_paths, col = "KO.paths", into = c("KO", "paths"), sep = ",", extra = "merge")
}

# ==================== 3. 创建DGEList对象 ====================
dge <- DGEList(
  counts = counts_data,
  group = factor(metadata$Group, levels = c("ST", "IT")),  # ST为参照
  genes = rownames(counts_data)
)

# ==================== 4. 过滤低表达基因 ====================
keep <- rowSums(counts_data > 0) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat(sprintf("过滤后保留KO数量: %d (原数据: %d)\n", 
            nrow(dge), nrow(counts_data)))
# ==================== 5. 标准化与差异分析 ====================
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)
dge <- estimateDisp(dge, design, robust = TRUE, prior.df = 1)
fit <- glmQLFit(dge, design, robust = TRUE)
contrast <- makeContrasts(IT - ST, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# ==================== 6. 提取结果 ====================
results <- topTags(qlf, n = Inf, adjust.method = "BH")$table
results$KO <- rownames(results)

# ==================== 7. 提取通路第三级名称 ====================
extract_level3 <- function(path_string) {
  paths <- unlist(strsplit(path_string, " \\| "))
  level3 <- sapply(strsplit(paths, "; "), function(x) {
    if(length(x) >= 3) x[3] else NA
  })
  paste(unique(na.omit(level3)), collapse = " | ")
}

kegg_paths$Level3 <- sapply(kegg_paths$paths, extract_level3)

# 合并通路信息
results <- merge(results, kegg_paths[, c("KO", "Level3")], by = "KO", all.x = TRUE)

# ==================== 8. 准备标注数据 ====================
# 定义显著性分组
results <- results %>%
  mutate(
    Significance = case_when(
      FDR < 0.05 & logFC > 1 ~ "IT UP",
      FDR < 0.05 & logFC < -1 ~ "ST UP",
      TRUE ~ "Not significant"
    )
  )

# 获取IT组和ST组各自前5的KO
it_top5 <- results %>%
  filter(logFC > 0) %>%
  arrange(desc(abs(logFC))) %>%
  head(0)

st_top5 <- results %>%
  filter(logFC < 0) %>%
  arrange(desc(abs(logFC))) %>%
  head(0)

label_data <- rbind(it_top5, st_top5) %>%
  mutate(
    Label = paste0(KO, ": ", ifelse(is.na(Level3), "Unknown", Level3)),
    Direction = ifelse(logFC > 0, "IT", "ST")
  )

# ==================== 9. 设置颜色方案 ====================
colors <- c(
  "IT UP" = "#E64B35",  # 深红色"IT-H" = "#E64B35", "ST-H" = "#4DBBD5"
  "ST UP" = "#4DBBD5",  # 深蓝色
  "Not significant" = "grey80"                  # 灰色
)

# ==================== 10. 绘制顶刊级火山图 ====================
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(FDR))) +
  # 散点图层（按显著性分组）
  geom_point(
    aes(fill = Significance, size = Significance, alpha = Significance),
    shape = 21, color = "white", stroke = 0.3
  ) +
  scale_fill_manual(values = colors) +
  scale_size_manual(values = c(
    "IT UP" = 3.5,
    "ST UP" = 3.5,
        "Not significant" = 3.0
  )) +
  scale_alpha_manual(values = c(
    "IT UP" = 0.8,
    "ST UP" = 0.8,
    
    "Not significant" = 0.7
  )) +
  
  # 参考线
  geom_vline(xintercept = c( -1, 0, 1), 
             linetype = c(2, 1, 2), color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, color = "grey40", linewidth = 0.5) +
  
  # 智能标注
  geom_text_repel(
    data = label_data,
    aes(label = Label, color = Direction),
    size = 3,
    fontface = "italic",
    box.padding = 0.6,
    point.padding = 0.3,
    segment.color = "grey40",
    segment.size = 0.3,
    min.segment.length = 0.2,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("IT" = "#D53E4F", "ST" = "#3288BD")) +
  
  # 坐标轴和标题
  labs(
    x = expression(log[2]("Fold Change (IT/ST)")),
    y = expression(-log[10]("FDR")),
    title = "Differential KEGG Pathway Abundance",
    subtitle = sprintf("IT (n=%d) vs ST (n=%d) | Top 5 most significant KOs labeled", 
                       sum(metadata$Group == "IT"), 
                       sum(metadata$Group == "ST")),
    fill = "Significance"
  ) +
  
  # 顶刊主题设置
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey40", margin = margin(b = 20)),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(color = "black", size = 10),
    panel.grid.major = element_line(linewidth = 0.2, color = "grey90"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 9),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    aspect.ratio = 0.8
  )



# ==================== 11. 保存结果 ====================
# 保存高清PDF（适合投稿）
ggsave(
  "Volcano_plot_IT_vs_ST.pdf",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 600,
  device = cairo_pdf
)

# 保存分析结果
write.csv(
  results[order(results$FDR), ],
  "KEGG_diff_results_IT_vs_ST.csv",
  row.names = FALSE
)

# 打印统计摘要
cat("=== Analysis Summary ===\n")
cat("Total KOs analyzed:", nrow(results), "\n")
cat("Significant KOs (FDR<0.05 & |logFC|>1):", sum(results$FDR < 0.05 & abs(results$logFC) > 1), "\n")
cat("IT-enriched KOs (FDR<0.05 & LFC>1.5):", sum(results$FDR < 0.05 & results$logFC > 1), "\n")
cat("ST-enriched KOs (FDR<0.05 & LFC<-1.5):", sum(results$FDR < 0.05 & results$logFC < -1), "\n")

# 显示图形
print(volcano_plot)


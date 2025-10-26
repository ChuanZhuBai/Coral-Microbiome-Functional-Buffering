keep <- filterByExpr(
  dge,
  group = dge$samples$group,
  min.count = 5,
  min.total.count = 10,  # 因样本量较少，放宽过滤阈值
  large.n = 2,
  min.prop = 0.2
)
dge <- dge[keep, , keep.lib.sizes = FALSE]


##############################################
# 潮间带vs潮下带实验室对照组(IT-C vs ST-C) KEGG通路差异火山图
# 最终优化版本 - 仅保留基础差异分析
############################################################################

# ==================== 1. 加载所有必要包 ====================
library(edgeR)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
dir()
# ==================== 2. 数据加载与预处理 ====================
# 读取原始counts数据（确保包含IT-C和ST-C样本）
counts_data <- read.csv("KEGG_bases_lab.csv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv", header = TRUE) %>%
  filter(Group %in% c("IT_H", "ST_H"))  # 仅保留实验室对照组

# 确保counts数据与metadata样本顺序一致
counts_data <- counts_data[, metadata$Sample]

# 读取通路数据
kegg_paths <- read.csv("KEGG_paths.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

# 提取KEGG Level3通路名称
extract_level3 <- function(path_string) {
  paths <- unlist(strsplit(path_string, " \\| "))
  level3 <- sapply(strsplit(paths, "; "), function(x) {
    if(length(x) >= 3) x[3] else NA
  })
  paste(unique(na.omit(level3)), collapse = " | ")
}
kegg_paths$Level3 <- sapply(kegg_paths$paths, extract_level3)

# ==================== 3. 创建DGEList对象 ====================
dge <- DGEList(
  counts = counts_data,
  group = factor(metadata$Group, levels = c("ST_H", "IT_H")),  # ST-C为参照
  genes = rownames(counts_data)
)

# ==================== 4. 过滤低表达基因 ====================
keep <- rowSums(counts_data > 0) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]

cat(sprintf("过滤后保留KO数量: %d (原数据: %d)\n", nrow(dge), nrow(counts_data)))

# ==================== 5. 标准化与差异分析 ====================
dge <- calcNormFactors(dge, method = "TMM")
design <- model.matrix(~0 + group, data = dge$samples)
colnames(design) <- levels(dge$samples$group)
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
contrast <- makeContrasts(`IT_H` - `ST_H`, levels = design)
qlf <- glmQLFTest(fit, contrast = contrast)

# ==================== 6. 提取结果 ====================
results <- topTags(qlf, n = Inf, adjust.method = "BH")$table
results$KO <- rownames(results)
results <- merge(results, kegg_paths[, c("KO", "Level3")], by = "KO", all.x = TRUE)

# ==================== 7. 准备标注数据 ====================
# 定义显著性分组
results <- results %>%
  mutate(
    Significance = case_when(
      PValue < 0.05 & logFC > 1 ~ "IT_H UP (PValue<5%, LFC>1)",
      PValue < 0.05 & logFC < -1 ~ "ST_H UP (PValue<5%, LFC>1)",
      TRUE ~ "Not significant"
    )
  )

# 选择标注点：按显著性排序
label_data <- results %>%
  filter(PValue < 0.05 & abs(logFC) > 1) %>%
  arrange(PValue) %>%
  head(0) %>%  # 取最显著的20个点标注
  mutate(
    Label = paste0(KO, "\n", strtrim(Level3, 40)),  # 限制标签长度
    Direction = ifelse(logFC > 0, "IT_H", "ST_H")
  )

# ==================== 8. 设置颜色方案 ====================
colors <- c(
  "IT_H UP (PValue<5%, LFC>1)" = "#E64B35",  # 深红色"#E64B35", "ST-H" = "#4DBBD5"
  "ST_H UP (PValue<5%, LFC>1)" = "#4DBBD5",  # 深蓝色
  "Not significant" = "grey80"
)

# ==================== 9. 绘制实验室对照组火山图 ====================
volcano_plot <- ggplot(results, aes(x = logFC, y = -log10(PValue))) +
  # 散点图层
  geom_point(
    aes(fill = Significance, size = Significance, alpha = Significance),
    shape = 21, color = "white", stroke = 0.3
  ) +
  scale_fill_manual(values = colors) +
  scale_size_manual(values = c(
    "IT_H UP (PValue<5%, LFC>1)" = 4,
    "ST_H UP (PValue<5%, LFC>1)" = 4,
    "Not significant" = 3
  )) +
  scale_alpha_manual(values = c(
    "IT_H UP (PValue<5%, LFC>1)" = 0.9,
    "ST_H UP (PValue<5%, LFC>1)" = 0.9,
    "Not significant" = 0.7
  )) +
  
  # 智能标注
  geom_text_repel(
    data = label_data,
    aes(label = Label, color = Direction),
    size = 3.2,
    fontface = "italic",
    box.padding = 0.7,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.3,
    min.segment.length = 0.1,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("IT_H" = "#D73027", "ST_H" = "#4575B4")) +
  
  # 参考线
  geom_vline(xintercept = c(-1, 0, 1), 
             linetype = c(2, 1, 2), color = "grey40", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = 3, color = "grey40", linewidth = 0.5) +
  
  # 坐标轴和标题
  labs(
    x = expression(log[2]("Fold Change (IT_H/ST_H)")),
    y = expression(-log[10]("PValue")),
    title = "Baseline Functional Differences in Coral Microbiome",
    subtitle = "KEGG pathway abundance: IT_H (lab-control) vs ST_H (lab-control)",
    fill = "Significance"
  ) +
  
  # 顶刊主题优化
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 18, margin = margin(b = 8)),
    plot.subtitle = element_text(hjust = 0.5, size = 13, color = "grey30", margin = margin(b = 15)),
    axis.title = element_text(face = "bold", size = 13),
    axis.text = element_text(color = "black", size = 11),
    panel.grid.major = element_line(linewidth = 0.25, color = "grey85"),
    panel.grid.minor = element_blank(),
    legend.position = "right",
    legend.box.spacing = unit(0.5, "cm"),
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.margin = margin(1.2, 1.2, 1.2, 1.2, "cm"),
    aspect.ratio = 0.75
  )

# ==================== 10. 保存结果 ====================
ggsave(
  "Volcano_plot_IT-H_vs_ST-H.pdf",
  plot = volcano_plot,
  width = 12,
  height = 8,
  dpi = 600,
  device = cairo_pdf
)



write.csv(
  results[order(results$FDR), ],
  "KEGG_diff_results_IT-H_vs_ST-H_1.csv",
  row.names = FALSE
)

# 打印统计摘要
cat("=== IT-H vs ST-H Analysis Summary ===\n")
cat("Total KOs analyzed:", nrow(results), "\n")
cat("Significant KOs (PValue < 0.05 & |logFC|>1):", sum(results$PValue < 0.05 & abs(results$logFC) > 1), "\n")
cat("IT-H enriched KOs:", sum(results$PValue < 0.05 & results$logFC > 1), "\n")
cat("ST-H enriched KOs:", sum(results$PValue < 0.05 & results$logFC < -1), "\n")

# 显示图形
print(volcano_plot)

































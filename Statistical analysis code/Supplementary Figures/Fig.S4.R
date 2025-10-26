#############################################################
##############################################################
#############################################################
# 1. 加载包 ------------------------------------------------------------------
library(vegan)
library(ggplot2)
library(dplyr)
library(ggrepel)
cat("===== 所有包加载完成 =====\n")

# 2. 数据准备 ----------------------------------------------------------------
# 读取数据
otu_table <- read.csv("relative-species.csv", row.names = 1, check.names = FALSE)
meta_data <- read.csv("metadata.csv")

# 检查样本顺序是否一致
if(!all(colnames(otu_table) == meta_data$Sample)) {
  warning("样本顺序不一致，正在调整...")
  otu_table <- otu_table[, meta_data$Sample]
}

cat("===== 数据预览 =====\n")
cat("物种表维度:", dim(otu_table), "\n")
cat("元数据分组:\n")
print(table(meta_data$Group))

# 3. PCoA分析 ---------------------------------------------------------------
# 转置并计算Bray-Curtis距离
dist_matrix <- vegdist(t(otu_table), method = "bray")

# 执行PCoA
pcoa <- cmdscale(dist_matrix, k = 3, eig = TRUE)
eig_percent <- round(pcoa$eig / sum(pcoa$eig[pcoa$eig > 0]) * 100, 1)

# 合并坐标与元数据
pcoa_points <- as.data.frame(pcoa$points[, 1:2]) %>%
  setNames(c("PC1", "PC2")) %>%
  tibble::rownames_to_column("Sample")
plot_data <- merge(pcoa_points, meta_data, by = "Sample")

cat("\n===== PCoA结果 =====\n")
cat("主坐标解释度:", eig_percent[1:2], "%\n")

# 4. 统计检验 ---------------------------------------------------------------
# PERMANOVA分析分组差异
adonis_result <- adonis2(dist_matrix ~ Group, data = meta_data, permutations = 9999)
cat("\n===== PERMANOVA结果 =====\n")
print(adonis_result)

# 5. 可视化 -----------------------------------------------------------------
# 科学配色方案
group_colors <- c("IT-C" = "#E64B35", "IT-H" = "#E41A1C",
                  "ST-C" = "#4DBBD5", "ST-H" = "#377EB8")

# 按您要求的形状映射
group_shapes <- c("IT-C" = 16,  # 圆形
                  "IT-H" = 16,  # 三角形
                  "ST-C" = 17,  # 方形
                  "ST-H" = 17)  # 菱形

# 基础PCoA图
p <- ggplot(plot_data, aes(x = PC1, y = PC2)) +
  # 修改点图层：分开颜色和形状映射
  geom_point(aes(color = Group, shape = Group), size = 4, alpha = 0.8) +
  # 添加椭圆（按颜色分组）
  stat_ellipse(aes(color = Group), level = 0.8, linetype = 2, linewidth = 0.5) +
  # 手动设置颜色和形状
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = group_shapes) +
  labs(
    x = paste0("PCo1 (", eig_percent[1], "%)"),
    y = paste0("PCo2 (", eig_percent[2], "%)"),
    title = "PCoA of Microbial Community Structure",
    subtitle = paste("PERMANOVA R2 =", round(adonis_result$R2[1], 3), 
                     "p =", adonis_result$`Pr(>F)`[1]),
    caption = "Shapes: IT=circle/triangle, ST=square/diamond"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.key = element_blank()
  ) +
  guides(
    color = guide_legend(override.aes = list(shape = c(16, 17, 15, 18))),
    shape = "none"  # 隐藏形状图例避免重复
  )

# 添加质心标签（带形状）
centroids <- plot_data %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2)) %>%
  mutate(Shape = group_shapes[Group])

p <- p + 
  geom_label_repel(
    data = centroids,
    aes(x = PC1, y = PC2, label = Group, color = Group),
    size = 4, show.legend = FALSE
  )

print(p)

# 保存图形
ggsave("PCoA_plot.png", p, width = 8, height = 6, dpi = 300)
cat("\n===== 分析完成，图形已保存 =====\n")












# 6. 绘图设置 ----------------------------------------------------------------
# 科学配色方案（ColorBrewer Set1修改版）
group_colors <- c("IT-C" = "#E64B35",  # 红色
                  "IT-H" = "#E41A1C", # 橙色
                  "ST-C" = "#4DBBD5", # 蓝色
                  "ST-H" = "#377EB8")   # 绿色

# 形状方案
group_shapes <- c("IT-C" = 16,  # 圆形
                  "IT-H" = 16,   # 三角形
                  "ST-C" = 17,   # 方形
                  "ST-H" = 17)   # 菱形


library(tidyverse)
library(ggsci)
library(scales)
library(patchwork)
dir()
# ==========================================
# 1. 载入数据 (实验组数据)
# ==========================================
# 使用之前生成的实验组专用数据表
data_rel_lab <- read.csv("Phylum_Relative_Abundance_Lab.csv", row.names = 1, check.names = FALSE)
metadata_lab <- read.csv("Metadata_Lab_Final.csv")

# ==========================================
# 2. 数据清洗与合并
# ==========================================
top_n <- 10

# 计算平均丰度并排除 "Unclassified Bacteria"
avg_abundance <- rowMeans(data_rel_lab)
top_phyla <- names(sort(avg_abundance, decreasing = TRUE))
top_phyla <- top_phyla[top_phyla != "Unclassified Bacteria"][1:top_n]

# 转换长表
data_long <- data_rel_lab %>%
  rownames_to_column("Phylum") %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Phylum_Display = ifelse(Phylum %in% top_phyla, Phylum, "Others"))

# 排序逻辑：Others 在最上方，Top N 按丰度从下往上堆
phylum_order <- c(top_phyla, "Others")
data_long$Phylum_Display <- factor(data_long$Phylum_Display, levels = rev(phylum_order))

# 合并元数据
plot_df <- data_long %>%
  left_join(metadata_lab, by = "Sample") %>%
  mutate(Habitat = factor(Habitat, levels = c("Intertidal", "Subtidal")),
         Treatment = factor(Treatment, levels = c("Control", "Heat")))

# ==========================================
# 3. 颜色方案准备 (保持与野外图一致)
# ==========================================
# 使用同样的配色种子，确保跨图表颜色对应关系一致
my_cols <- pal_npg()(10)
if (top_n > 10) {
  my_cols <- colorRampPalette(my_cols)(top_n)
}
color_map <- setNames(c(my_cols[1:top_n], "#D3D3D3"), phylum_order)

# ==========================================
# 4. 绘图 (Science Advances 风格)
# ==========================================
p_exp_stack <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Phylum_Display)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.1) +
  # 双重分面：第一层按生境，第二层按处理
  facet_grid(. ~ Habitat + Treatment, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_map) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = "Experimental Samples", 
       y = "Relative Abundance", 
       fill = "Bacterial Phylum") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.y = element_text(size = 12),
    # 优化分面标签
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(size = 10, face = "bold", margin = margin(5, 5, 5, 5)),
    legend.text = element_text(face = "plain", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    panel.spacing = unit(0.2, "cm")
  )

# ==========================================
# 5. 展示与保存
# ==========================================
print(p_exp_stack)

# 建议在补充材料中标记为 Figure S3
ggsave("Figure_S3_Experimental_Phylum_Composition.pdf", p_exp_stack, width = 10, height = 6)

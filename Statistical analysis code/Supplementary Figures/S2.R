library(tidyverse)
library(ggsci)
library(patchwork)

# ==========================================
# 1. 载入数据
# ==========================================
data_rel <- read.csv("Phylum_Relative_Abundance_Field.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("Metadata_Field_Final.csv")

# ==========================================
# 2. 数据清洗与合并
# ==========================================
top_n <- 10

# 计算平均丰度
avg_abundance <- rowMeans(data_rel)

# 找出前 N 个门，同时显式排除 "Unclassified Bacteria"
# 这样确保 "Unclassified Bacteria" 一定会被归入 "Others"
top_phyla <- names(sort(avg_abundance, decreasing = TRUE))
top_phyla <- top_phyla[top_phyla != "Unclassified Bacteria"][1:top_n]

# 转换长表
data_long <- data_rel %>%
  rownames_to_column("Phylum") %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Abundance") %>%
  # 核心逻辑：不在 Top N 或者名字包含 Unclassified 的全部归为 Others
  mutate(Phylum_Display = ifelse(Phylum %in% top_phyla, Phylum, "Others"))

# 确定排序顺序：Top N 丰度降序排列，Others 放在最后
phylum_order <- c(top_phyla, "Others")
# 在 ggplot 堆叠图中，factor 的第一个 level 会在最下面
# 我们希望 Others 在最上面，所以反转顺序
data_long$Phylum_Display <- factor(data_long$Phylum_Display, levels = rev(phylum_order))

# 合并元数据
plot_df <- data_long %>%
  left_join(metadata, by = "Sample") %>%
  mutate(Habitat = factor(Habitat, levels = c("Intertidal", "Subtidal")))

# ==========================================
# 3. 颜色方案准备
# ==========================================
# 提取 NPG 配色并手动指定 Others 为灰色
# 使用 pal_npg 提取前 10 个颜色
library(scales)
my_cols <- pal_npg()(10)
# 如果 top_n 大于 10，则通过插值扩充颜色
if (top_n > 10) {
  my_cols <- colorRampPalette(my_cols)(top_n)
}

# 建立颜色映射表
# 这里的顺序要和 phylum_order 对应
color_map <- setNames(c(my_cols[1:top_n], "#D3D3D3"), phylum_order)

# ==========================================
# 4. 绘图 (Science Advances 风格)
# ==========================================
p_field_stack <- ggplot(plot_df, aes(x = Sample, y = Abundance, fill = Phylum_Display)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.1) + # width = 1 让柱子紧挨着，更显紧凑
  facet_grid(. ~ Habitat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_map) +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(x = "In-situ Samples", 
       y = "Relative Abundance", 
       fill = "Bacterial Phylum") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 5)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.text.y = element_text(size = 12),
    strip.background = element_rect(fill = "grey95", color = NA), # 去掉分面背景框
    strip.text = element_text(size = 12, face = "bold", margin = margin(5, 5, 5, 5)),
    legend.text = element_text(face = "plain", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = "right",
    panel.spacing = unit(0.2, "cm") # 缩小分面之间的距离
  )

# ==========================================
# 5. 展示与保存
# ==========================================
print(p_field_stack)

# 建议尺寸：宽10-12英寸，高6英寸
ggsave("Figure_S2_Field_Phylum_Composition.pdf", p_field_stack, width = 11, height = 6)


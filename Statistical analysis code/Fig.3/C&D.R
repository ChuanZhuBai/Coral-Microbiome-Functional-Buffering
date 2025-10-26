# 加载必要包
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(rstatix)

# 数据准备 ------------------------------------------------------------
species <- read.csv("relative-species.csv", header = TRUE, row.names = 1, check.names = FALSE)

species <- read.csv("CSS_relative_KEGG_bases_lab.csv", header =TRUE, row.names =1, check.names =FALSE)#功能组间的距离"Bray-Curtis Distance"
metadata <- read.csv("metadata.csv", header = TRUE)

# 确保样本顺序一致
species <- species[, metadata$Sample]

# 计算Bray-Curtis距离矩阵
bray_dist <- vegdist(t(species), method = "bray")

# 数据整理 ------------------------------------------------------------
dist_df <- as.data.frame(as.matrix(bray_dist))
dist_df$Sample <- rownames(dist_df)

dist_long <- dist_df %>%
  pivot_longer(cols = -Sample, names_to = "Sample2", values_to = "Distance") %>%
  left_join(metadata, by = "Sample") %>%
  left_join(metadata, by = c("Sample2" = "Sample"), suffix = c("", ".y")) %>%
  filter(Sample != Sample2) # 去除自身比较

# 筛选组内比较（Control vs. Heat）
stability_data <- dist_long %>%
  filter(
    (Group == "IT-C" & Group.y == "IT-H") | # 潮间带Control vs. Heat
      (Group == "ST-C" & Group.y == "ST-H")    # 潮下带Control vs. Heat
  ) %>%
  mutate(GroupPair = factor(ifelse(Group == "IT-C", "Intertidal", "Subtidal"),
                            levels = c("Intertidal", "Subtidal")))

# 统计分析 ------------------------------------------------------------
# 更严谨的统计检验：添加效应量计算
stat_test <- stability_data %>% 
  wilcox_test(Distance ~ GroupPair) %>%
  add_significance() %>%
  add_xy_position(x = "GroupPair")

# 计算组间差异的效应量
effect_size <- wilcox_effsize(stability_data, Distance ~ GroupPair)

# 可视化设置 ----------------------------------------------------------
fill_colors <- c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")
point_colors <- c("Intertidal" = "#E64B35B3", "Subtidal" = "#4DBBD5B3") # 添加透明度

# 科学可视化 ---------------------------------------------------------
ggplot(stability_data, aes(x = GroupPair, y = Distance)) +
  # 箱线图
  geom_boxplot(
    aes(fill = GroupPair),
    width = 0.6, 
    alpha = 0.8, 
    outlier.shape = NA,
    show.legend = FALSE
  ) +
  # 散点（颜色与组别一致，带透明度）
  geom_jitter(
    aes(color = GroupPair),
    width = 0.1, 
    size = 2,
    alpha = 0.6,
    show.legend = FALSE
  ) +
  # 颜色设置
  scale_fill_manual(values = fill_colors) +
  scale_color_manual(values = point_colors) +
  
  # 统计标注
  stat_pvalue_manual(
    stat_test, 
    label = "p = {p}{p.signif}",
    tip.length = 0.01,
    bracket.size = 0.6,
    label.size = 4.5
  ) +
  
  # 标题和坐标轴
  labs(
    x = "",
    y = "Bray-Curtis Dissimilarity\n(Control vs Heat)",
    title = "Microbial Community Stability Under Heat Stress",
    subtitle = paste0(
      "Wilcoxon rank-sum test: p = ", signif(stat_test$p, 3),
      ", Effect size (r) = ", signif(effect_size$effsize, 2)
    ),
    caption = "Data presented as median ± IQR\nPoints represent individual sample pairs"
  ) +
  
  # 主题设置
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    axis.text.x = element_text(size = 11, color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.margin = unit(c(10, 10, 10, 10), "pt"),
    plot.caption = element_text(size = 9, hjust = 0, color = "grey40")
  ) +
  
  # 保证y轴从0开始
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))





















# 加载必要包
library(vegan)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
dir()
# 读取数据
species <- read.csv("relative-species.csv", header =TRUE, row.names =1, check.names =FALSE)

species <- read.csv("CSS_relative_KEGG_bases_lab.csv", header =TRUE, row.names =1, check.names =FALSE)#功能组间的距离"Bray-Curtis Distance"

metadata <- read.csv("metadata.csv", header =TRUE)


#######################################
# 示例：使用vegan包计算功能冗余
library(vegan)
# 假设species为物种丰度表，kegg为功能丰度表
redundancy <- function(species, kegg) {
  spec_dist <- vegdist(t(species), "bray")
  func_dist <- vegdist(t(kegg), "bray")
  mantel(spec_dist, func_dist, method = "pearson")  # 检验群落与功能距离相关性
}
###########################################
# 确保样本顺序一致
species <- species[, metadata$Sample]

# 计算Bray-Curtis距离矩阵
bray_dist <- vegdist(t(species), method ="bray")

# 转换为数据框并合并分组信息
dist_df <- as.data.frame(as.matrix(bray_dist))
dist_df$Sample <- rownames(dist_df)
dist_long <- dist_df %>%
  pivot_longer(cols =-Sample, names_to ="Sample2", values_to ="Distance") %>%
  left_join(metadata, by ="Sample") %>%
  left_join(metadata, by =c("Sample2" ="Sample"), suffix =c("", ".y")) %>%
  filter(Sample != Sample2) # 去除自身比较

# 筛选组内比较（Control vs. Heat的同组内比较）
stability_data <- dist_long %>%
  filter(
    (Group =="IT-C" & Group.y =="IT-H") | # 潮间带Control vs. Heat
      (Group =="ST-C" & Group.y =="ST-H")    # 潮下带Control vs. Heat
  ) %>%
  mutate(GroupPair =ifelse(Group =="IT-C", "Intertidal", "Subtidal"))

# Wilcoxon秩和检验
wilcox_test <- wilcox.test(Distance ~ GroupPair, data =stability_data)
cat(paste0(
  "Wilcoxon p-value: ", round(wilcox_test$p.value, 10),
  ifelse(wilcox_test$p.value <0.05, " (显著)", " (不显著)")
))
#"IT-H" = "#E64B35", "ST-H" = "#4DBBD5"
# 输出组间距离中位数
stability_data %>%
  group_by(GroupPair) %>%
  summarise(Median_Distance =median(Distance))


write.csv(stability_data, "stability_data.csv")
ggplot(stability_data, aes(x =GroupPair, y =Distance, fill =GroupPair)) +
  geom_boxplot(width =0.6, alpha =0.8, outlier.shape =NA) +
  geom_jitter(width =0.1, size =2, alpha =0.5, color ="gray30") +
  scale_fill_manual(values =c("Intertidal" ="#E64B35", "Subtidal" ="#4DBBD5")) +
  labs(
    x ="",
    y ="Bray-Curtis Distance\n(Control vs. Heat)",
    title ="Stability of microbial community structure under heat stress",
    subtitle =paste0("Wilcoxon检验: p =", round(wilcox_test$p.value, 4))
  ) +
  theme_classic(base_size =14) +
  theme(
    legend.position ="none",
    plot.title =element_text(face ="bold", hjust =0.5),
    plot.subtitle =element_text(hjust =0.5)
  ) +
  stat_compare_means(
    method ="wilcox.test",
    comparisons =list(c("Intertidal", "Subtidal")),
    label ="p.signif",
    symnum.args =list(cutpoints =c(0, 0.001, 0.01, 0.05, 1), symbols =c("***", "**", "*", "ns"))
  )


# 计算各组内Beta多样性离散度（dispersion）
dispersion <- betadisper(bray_dist, group =metadata$Group)
anova(dispersion) # 检验组间离散度差异
boxplot(dispersion, col =c("#E69F00", "#56B4E9"), main ="Beta多样性离散度")





#评估潮间带(Intertidal)和潮下带(Subtidal)在热胁迫下的功能稳定性差异
#################################################
# 加载必要包
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)

# 1. 读取差异分析结果
it_res <- read.csv("KEGG_diff_results_IT-H_vs_IT-C.csv", stringsAsFactors = FALSE)
st_res <- read.csv("KEGG_diff_results_ST-H_vs_ST-C.csv", stringsAsFactors = FALSE)

# 2. 数据清洗与标注
it_res$Group <- "Intertidal"
st_res$Group <- "Subtidal"

# 合并数据并筛选显著差异KO (FDR < 0.05)
combined <- bind_rows(it_res, st_res) %>%
  filter(PValue < 0.05) %>%
  mutate(
    Direction = ifelse(logFC > 0, "Up", "Down"),
    Abs_logFC = abs(logFC)  # 计算变化幅度绝对值
  )

# 3. 统计检验
# 3.1 比较两组间logFC绝对值的差异
wilcox_test <- wilcox.test(Abs_logFC ~ Group, data = combined)
cat(sprintf(
  "变化幅度比较(中位数):\n潮间带 = %.2f\n潮下带 = %.2f\nWilcoxon p = %.3g",
  median(combined$Abs_logFC[combined$Group == "Intertidal"]),
  median(combined$Abs_logFC[combined$Group == "Subtidal"]),
  wilcox_test$p.value
))

# 3.2 按变化方向分组统计
change_stats <- combined %>%
  group_by(Group, Direction) %>%
  summarise(
    Count = n(),
    Median_logFC = median(logFC),
    Median_Abs_logFC = median(Abs_logFC),
    .groups = "drop"
  )

# 4. 可视化
# 4.1 变化幅度分布（箱线图）
p1 <- ggplot(combined, aes(x = Group, y = Abs_logFC, fill = Group)) +
  geom_boxplot(width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 1.5, alpha = 0.3) +
  scale_fill_manual(values = c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")) +
  labs(
    x = "",
    y = "|log2(Fold Change)|",
    title = "Functional Change Magnitude",
    subtitle = sprintf("Wilcoxon p = %.2g", wilcox_test$p.value)
  ) +
  stat_compare_means(
    comparisons = list(c("Intertidal", "Subtidal")),
    method = "wilcox.test",
    label = "p.signif"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")

print(p1)
# 4.2 差异KO数量与方向（条形图）#"Up" = "#619CFF", "Down" = "#F8766D"
p2 <- ggplot(change_stats, aes(x = Group, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Up" = "#F8766D", "Down" = "#619CFF")) +
  labs(
    x = "",
    y = "Number of Significant KOs",
    title = "Differential KO Counts",
    fill = "Regulation"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

print(p2)
# 4.3 按代谢通路分类展示变化（点图）
# 截断过长的通路名称
combined <- combined %>%
  mutate(Level3_short = ifelse(
    nchar(Level3) > 50,
    paste0(substr(Level3, 1, 47), "..."),
    Level3
  ))

p3_short <- ggplot(combined %>% filter(Level3 %in% top_pathways$Level3), 
                   aes(x = logFC, y = reorder(Level3_short, logFC))) +
  geom_point(aes(color = Group), alpha = 0.6) +
  scale_color_manual(values = c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")) +
  labs(y = "") +
  theme(axis.text.y = element_text(size = 10))

print(p3_short)
# 选择在至少一组中差异显著的top通路（按KO数量）
top_pathways <- combined %>%
  group_by(Level3) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count)) %>%
  head(20)  # 展示前20个最重要的通路

p3_filtered <- combined %>%
  filter(Level3 %in% top_pathways$Level3) %>%
  ggplot(aes(x = logFC, y = reorder(Level3, logFC, median), color = Group)) +
  geom_point(position = position_jitterdodge(), alpha = 0.7, size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")) +
  labs(
    x = "log2(Fold Change)", 
    y = "",
    title = "Top 20 Pathway Changes"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10),
    panel.grid.major.y = element_line(linetype = "dotted")
  )




print(p3_filtered)

p3_facet <- combined %>%
  filter(Level3 %in% top_pathways$Level3) %>%
  ggplot(aes(x = Group, y = logFC, color = Group)) +
  geom_boxplot(width = 0.5) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  facet_wrap(~ Level3, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c("Intertidal" = "#E64B35", "Subtidal" = "#4DBBD5")) +
  labs(x = "", y = "log2(Fold Change)") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    strip.text = element_text(size = 8)
  )

print(p3_facet)


library(plotly)
p3_interactive <- plot_ly(
  data = combined %>% filter(Level3 %in% top_pathways$Level3),
  x = ~logFC, 
  y = ~Level3, 
  color = ~Group,
  colors = c("#E64B35", "#4DBBD5"),
  type = 'scatter',
  mode = 'markers',
  hoverinfo = 'text',
  text = ~paste("KO:", KO, "<br>Pathway:", Level3)
) %>% 
  layout(
    yaxis = list(title = "", tickfont = list(size = 10)),
    xaxis = list(title = "log2(Fold Change)")
  )
htmlwidgets::saveWidget(p3_interactive, "pathway_changes.html")
# 5. 保存结果
# 5.1 保存统计结果
write.csv(change_stats, "functional_change_statistics.csv", row.names = FALSE)

# 5.2 保存图表
ggsave("change_magnitude_comparison.png", p1, width = 6, height = 5, dpi = 300)
ggsave("differential_KO_counts.png", p2, width = 6, height = 5, dpi = 300)
ggsave("pathway_changes.png", p3, width = 8, height = 6, dpi = 300)

# 6. 结果显示
gridExtra::grid.arrange(p1, p2, ncol = 2)
print(p3)









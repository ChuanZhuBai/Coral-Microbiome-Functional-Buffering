######Shannon指数
# 1. 加载必要包 
library(vegan) 
library(ggplot2) 
library(ggpubr) 
library(rstatix) 
library(ggbeeswarm) 
library(dplyr) 
dir() 

# 2. 数据预处理 --------------------------------------------------------------- 
species_data <- read.csv("absolute-species.csv", row.names = 1, check.names = FALSE) 
metadata <- read.csv("metadata.csv") 

# 将绝对丰度转换为相对丰度
rel_species_data <- t(apply(species_data, 2, function(x) x/sum(x)))

# 计算Shannon指数并合并数据 
alpha_data <- data.frame( 
  Sample = names(diversity(rel_species_data, index = "shannon")), 
  Shannon = diversity(rel_species_data, index = "shannon") 
) %>% 
  left_join(metadata, by = "Sample") %>% 
  mutate( 
    Group = factor(Group, levels = c("IT", "ST")),  # 先完成Group因子的转换 
    Group_label = ifelse(Group == "IT", "Intertidal", "Subtidal")  # 然后添加新列 
  ) 
# 查看前几行数据确认 
head(alpha_data) 

# 3. 异常值处理（IQR法）------------------------------------------------------ 
alpha_clean <- alpha_data %>% 
  group_by(Group) %>% 
  mutate( 
    Q1 = quantile(Shannon, 0.25), 
    Q3 = quantile(Shannon, 0.75), 
    IQR = Q3 - Q1, 
    lower_bound = Q1 - 1.5*IQR, 
    upper_bound = Q3 + 1.5*IQR 
  ) %>% 
  filter(Shannon >= lower_bound & Shannon <= upper_bound) %>% 
  ungroup() 

# 4. 统计检验 ---------------------------------------------------------------- 
stat_test <- alpha_clean %>% 
  wilcox_test(Shannon ~ Group) %>% 
  add_y_position(step.increase = 0.12) %>%  # 移除了错误的scale参数 
  mutate(label = case_when( 
    p < 0.001 ~ "***", 
    p < 0.01 ~ "**", 
    p < 0.05 ~ "*", 
    TRUE ~ paste0("p=", round(p, 2)) 
  )) 

# 5. 高级可视化（最终优化版）-------------------------------------------------- 
shannon_plot <- ggplot(alpha_clean, aes(x = Group, y = Shannon)) + 
  # 半透明箱线图 
  geom_boxplot( 
    aes(fill = Group), 
    width = 0.6, 
    alpha = 0.2, 
    linewidth = 0.5,  # 使用linewidth替代已弃用的size 
    outlier.shape = NA 
  ) + 
  # 小提琴图展示分布 
  geom_violin( 
    aes(fill = Group), 
    width = 0.7, 
    alpha = 0.3, 
    trim = FALSE, 
    scale = "width" 
  ) + 
  # 散点分布（带抖动） 
  geom_quasirandom( 
    aes(color = Group), 
    size = 2, 
    alpha = 0.8, 
    width = 0.15, 
    shape = 16, 
    method = "smiley" 
  ) + 
  # 统计标注 
  stat_pvalue_manual( 
    stat_test, 
    label = "label", 
    tip.length = 0.02, 
    bracket.size = 0.4, 
    vjust = 0.5 
  ) + 
  # 配色方案（使用指定颜色） 
  scale_fill_manual( 
    values = c("IT" = "#E64B35", "ST" = "#4DBBD5"), 
    guide = "none" 
  ) + 
  scale_color_manual( 
    values = c("IT" = "#E64B35", "ST" = "#4DBBD5"), 
    guide = "none" 
  ) + 
  # 坐标轴和标签 
  labs( 
    title = "Shannon Diversity Index Comparison", 
    subtitle = paste0("Wilcoxon rank-sum test: ", stat_test$label[1]), 
    x = NULL, 
    y = "Shannon Diversity Index" 
  ) + 
  # 主题设置（优化版） 
  theme_minimal(base_size = 12) +  # 增大基础字号 
  theme( 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5), 
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line = element_line(linewidth = 0.5, color = "black"),  # 加粗轴线 
    axis.ticks = element_line(linewidth = 0.5, color = "black"),  # 加粗刻度线 
    axis.text = element_text(color = "black", size = 11),  # 增大坐标轴文字 
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # 增大y轴标题 
    plot.margin = margin(10, 10, 10, 10)  # 调整边距 
  ) + 
  # Y轴扩展 
  scale_y_continuous( 
    expand = expansion(mult = c(0.05, 0.15)), 
    limits = c(NA, max(alpha_clean$Shannon) * 1.1) 
  ) 

# 6. 预览图形 
print(shannon_plot) 

# 7. 保存图形（高质量PDF） 
ggsave("Final_Shannon_Diversity_Optimized-new.pdf", 
       plot = shannon_plot, 
       width = 8.6,  # Nature单栏标准宽度(8.6 cm) 
       height = 6.5, # 黄金比例高度 
       units = "cm", 
       dpi = 600, 
       device = cairo_pdf) 


####################################################
####################################################
####################################################
#####Chao1指数
# ======================
# 修正后的Chao1指数分析代码
# ======================

# 1. 加载必要包
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(dplyr)
dir()
# 2. 数据预处理 ---------------------------------------------------------------
# 读取绝对丰度数据
abs_species_data <- read.csv("absolute-species.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv")

# 计算每个样本的Chao1指数
chao1_values <- apply(abs_species_data, 2, function(x) {
  est <- estimateR(x)
  est["S.chao1"]
})

# 创建数据框
alpha_data <- data.frame(
  Sample = names(chao1_values),
  Chao1 = unname(chao1_values),
  stringsAsFactors = FALSE
) %>% 
  left_join(metadata, by = "Sample") %>% 
  mutate(
    Group = factor(Group, levels = c("IT", "ST")),
    Group_label = ifelse(Group == "IT", "Intertidal", "Subtidal")
  )

# 3. 异常值处理（IQR法）------------------------------------------------------
alpha_clean <- alpha_data %>% 
  group_by(Group) %>% 
  mutate(
    Q1 = quantile(Chao1, 0.25),
    Q3 = quantile(Chao1, 0.75),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5*IQR,
    upper_bound = Q3 + 1.5*IQR
  ) %>% 
  filter(Chao1 >= lower_bound & Chao1 <= upper_bound) %>% 
  ungroup()

# 4. 统计检验 ----------------------------------------------------------------
stat_test <- alpha_clean %>% 
  wilcox_test(Chao1 ~ Group) %>% 
  add_y_position(step.increase = 0.12) %>% 
  mutate(label = case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ paste0("p=", round(p, 2))
  ))
  
  # 5. 高级可视化 --------------------------------------------------------------
  chao1_plot <- ggplot(alpha_clean, aes(x = Group, y = Chao1)) +
    geom_boxplot(aes(fill = Group), width = 0.6, alpha = 0.2, linewidth = 0.5, outlier.shape = NA) +
    geom_violin(aes(fill = Group), width = 0.7, alpha = 0.3, trim = FALSE, scale = "width") +
    geom_quasirandom(aes(color = Group), size = 2, alpha = 0.8, width = 0.15, shape = 16, method = "smiley") +
    stat_pvalue_manual(stat_test, label = "label", tip.length = 0.02, bracket.size = 0.4, vjust = 0.5) +
    scale_fill_manual(values = c("IT" = "#E64B35", "ST" = "#4DBBD5"), guide = "none") +
    scale_color_manual(values = c("IT" = "#E64B35", "ST" = "#4DBBD5"), guide = "none") +
    labs(
      title = "Chao1 Richness Index Comparison",
      subtitle = paste0("Wilcoxon rank-sum test: ", stat_test$label[1]),
      x = NULL,
      y = "Chao1 Richness Index"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(linewidth = 0.5, color = "black"),
      axis.ticks = element_line(linewidth = 0.5, color = "black"),
      axis.text = element_text(color = "black", size = 11),
      axis.title.y = element_text(size = 12, margin = margin(r = 10)),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  # 6. 保存图形
  ggsave("Chao1_Richness_Comparison.pdf", 
         plot = chao1_plot,
         width = 8.6, height = 6.5, units = "cm", dpi = 600, device = cairo_pdf)
  
  # 显示图形
  print(chao1_plot)














































###############绝对丰度
# 1. 加载必要包
library(vegan)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggbeeswarm)
library(dplyr)
dir()
# 2. 数据预处理 ---------------------------------------------------------------
species_data <- read.csv("absolute-species.csv", row.names = 1, check.names = FALSE)
metadata <- read.csv("metadata.csv")

# 计算Shannon指数并合并数据
# 计算Shannon指数并合并数据
alpha_data <- data.frame(
  Sample = names(diversity(t(species_data), index = "shannon")),
  Shannon = diversity(t(species_data), index = "shannon")
) %>% 
  left_join(metadata, by = "Sample") %>% 
  mutate(
    Group = factor(Group, levels = c("IT", "ST")),  # 先完成Group因子的转换
    Group_label = ifelse(Group == "IT", "Intertidal", "Subtidal")  # 然后添加新列
  )
# 查看前几行数据确认
head(alpha_data)
         
# 3. 异常值处理（IQR法）------------------------------------------------------
         alpha_clean <- alpha_data %>%
           group_by(Group) %>%
           mutate(
             Q1 = quantile(Shannon, 0.25),
             Q3 = quantile(Shannon, 0.75),
             IQR = Q3 - Q1,
             lower_bound = Q1 - 1.5*IQR,
             upper_bound = Q3 + 1.5*IQR
           ) %>%
           filter(Shannon >= lower_bound & Shannon <= upper_bound) %>%
           ungroup()
        
# 4. 统计检验 ----------------------------------------------------------------
stat_test <- alpha_clean %>%
  wilcox_test(Shannon ~ Group) %>%
  add_y_position(step.increase = 0.12) %>%  # 移除了错误的scale参数
  mutate(label = case_when(
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ paste0("p=", round(p, 2))
  ))
           
# 5. 高级可视化（最终优化版）--------------------------------------------------
shannon_plot <- ggplot(alpha_clean, aes(x = Group, y = Shannon)) +
  # 半透明箱线图
  geom_boxplot(
    aes(fill = Group),
    width = 0.6,
    alpha = 0.2,
    linewidth = 0.5,  # 使用linewidth替代已弃用的size
    outlier.shape = NA
  ) +
  # 小提琴图展示分布
  geom_violin(
    aes(fill = Group),
    width = 0.7,
    alpha = 0.3,
    trim = FALSE,
    scale = "width"
  ) +
  # 散点分布（带抖动）
  geom_quasirandom(
    aes(color = Group),
    size = 2,
    alpha = 0.8,
    width = 0.15,
    shape = 16,
    method = "smiley"
  ) +
  # 统计标注
  stat_pvalue_manual(
    stat_test,
    label = "label",
    tip.length = 0.02,
    bracket.size = 0.4,
    vjust = 0.5
  ) +
  # 配色方案（使用指定颜色）
  scale_fill_manual(
    values = c("IT" = "#E64B35", "ST" = "#4DBBD5"),
    guide = "none"
  ) +
  scale_color_manual(
    values = c("IT" = "#E64B35", "ST" = "#4DBBD5"),
    guide = "none"
  ) +
  # 坐标轴和标签
  labs(
    title = "Shannon Diversity Index Comparison",
    subtitle = paste0("Wilcoxon rank-sum test: ", stat_test$label[1]),
    x = NULL,
    y = "Shannon Diversity Index"
  ) +
  # 主题设置（优化版）
  theme_minimal(base_size = 12) +  # 增大基础字号
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, color = "gray40", hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, color = "black"),  # 加粗轴线
    axis.ticks = element_line(linewidth = 0.5, color = "black"),  # 加粗刻度线
    axis.text = element_text(color = "black", size = 11),  # 增大坐标轴文字
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),  # 增大y轴标题
    plot.margin = margin(10, 10, 10, 10)  # 调整边距
  ) +
  # Y轴扩展
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.15)),
    limits = c(NA, max(alpha_clean$Shannon) * 1.1)
  )

# 6. 预览图形
print(shannon_plot)

# 7. 保存图形（高质量PDF）
ggsave("Final_Shannon_Diversity_Optimized.pdf",
       plot = shannon_plot,
       width = 8.6,  # Nature单栏标准宽度(8.6 cm)
       height = 6.5, # 黄金比例高度
       units = "cm",
       dpi = 600,
       device = cairo_pdf)


























#################################################
##################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(rstatix)
dir()
# 假设已有数据框alpha_data，包含以下列：
# Sample, Group (IT/ST), Shannon, Chao1
library(tidyverse)
library(ggpubr)
library(rstatix)

# 1. 读取数据
alpha_data <- read.csv("diversity_indices_absolute.csv", row.names = 1, check.names = FALSE)

# 2. 数据清洗
clean_data <- alpha_data %>%
  group_by(Group) %>%
  mutate(
    Shannon = ifelse(Shannon > quantile(Shannon, 0.75) + 1.5*IQR(Shannon) | 
                       Shannon < quantile(Shannon, 0.25) - 1.5*IQR(Shannon), NA, Shannon),
    Chao1 = ifelse(Chao1 > quantile(Chao1, 0.75) + 1.5*IQR(Chao1) | 
                     Chao1 < quantile(Chao1, 0.25) - 1.5*IQR(Chao1), NA, Chao1)
  ) %>%
  na.omit()

# 3. 预处理
clean_data$Group <- factor(clean_data$Group, levels = c("IT", "ST"))
alpha_long <- clean_data %>%
  pivot_longer(cols = c(Shannon, Chao1), 
               names_to = "Index", 
               values_to = "Value") %>%
  mutate(Index = factor(Index, levels = c("Shannon", "Chao1")))
         
         # 4. 统计检验
         stat_test <- alpha_long %>%
           group_by(Index) %>%
           t_test(Value ~ Group) %>%
           adjust_pvalue(method = "BH") %>%
           add_xy_position(x = "Group", dodge = 0.8)
         
         # 5. 创建公共主题
         my_theme <- theme(
           plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
           axis.title = element_text(face = "bold"),
           axis.text = element_text(color = "black"),
           legend.position = "none",
           panel.background = element_blank(),
           panel.border = element_rect(fill = NA, color = "black"),
           strip.background = element_blank(),
           strip.text = element_text(face = "bold", size = 12)
         )
         
         # 6. 创建分面图（自动调整y轴范围）
         auto_plot <- ggplot(alpha_long, aes(x = Group, y = Value)) +
           geom_boxplot(aes(fill = Group), width = 0.6, alpha = 0.8, outlier.shape = NA) +
           geom_jitter(aes(fill = Group), width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.6) +
           facet_wrap(~ Index, scales = "free_y", nrow = 1) +
           stat_pvalue_manual(
             stat_test, 
             label = "p = {p.adj}", 
             tip.length = 0.01,
             size = 3.5
           ) +
           scale_fill_manual(values = c("IT" = "#E64B35", "ST" = "#00A087")) +
           labs(x = NULL, y = "Alpha Diversity Index", 
                title = "Microbial Diversity Comparison") +
           my_theme
         
         # 7. 创建手动控制y轴范围的组合图
         plot_shannon <- ggplot(filter(alpha_long, Index == "Shannon"), 
                                aes(x = Group, y = Value)) +
           geom_boxplot(aes(fill = Group), width = 0.6, alpha = 0.8, outlier.shape = NA) +
           geom_jitter(aes(fill = Group), width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.6) +
           stat_pvalue_manual(
             filter(stat_test, Index == "Shannon"),
             label = "p = {p.adj}", 
             tip.length = 0.01,
             size = 3.5
           ) +
           scale_fill_manual(values = c("IT" = "#E64B35", "ST" = "#00A087")) +
           labs(x = NULL, y = "Shannon Index") +
           ylim(1, 4) +
           my_theme
         
         plot_chao1 <- ggplot(filter(alpha_long, Index == "Chao1"), 
                              aes(x = Group, y = Value)) +
           geom_boxplot(aes(fill = Group), width = 0.6, alpha = 0.8, outlier.shape = NA) +
           geom_jitter(aes(fill = Group), width = 0.15, size = 2, shape = 21, color = "black", alpha = 0.6) +
           stat_pvalue_manual(
             filter(stat_test, Index == "Chao1"),
             label = "p = {p.adj}", 
             tip.length = 0.01,
             size = 3.5
           ) +
           scale_fill_manual(values = c("IT" = "#E64B35", "ST" = "#00A087")) +
           labs(x = NULL, y = "Chao1 Index") +
           ylim(100, 400) +
           my_theme
         
         # 组合两个图
         manual_plot <- ggarrange(plot_shannon, plot_chao1, 
                                  ncol = 2, common.legend = TRUE) %>%
           annotate_figure(top = text_grob("Microbial Diversity Comparison", 
                                           face = "bold", size = 14))
         
         # 8. 保存图形
         ggsave("alpha_diversity_auto.png", auto_plot, width = 10, height = 5, dpi = 300)
         ggsave("alpha_diversity_manual.png", manual_plot, width = 10, height = 5, dpi = 300)
         
         # 同时显示两个版本
         auto_plot
         manual_plot
# 保存为PDF（投稿质量）
ggsave(
  "alpha_diversity_auto.pdf",auto_plot,
  width = 10, height = 5,
  units = "cm",
  device = cairo_pdf,
  dpi = 600
)
ggsave(
  "alpha_diversity_manual.pdf",manual_plot,
  width = 10, height = 5,
  units = "cm",
  device = cairo_pdf,
  dpi = 600
)




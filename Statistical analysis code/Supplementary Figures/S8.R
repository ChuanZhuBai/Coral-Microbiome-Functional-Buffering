# 加载必要包
library(tidyverse)
library(vegan)
library(ggpubr)
library(rstatix)
library(multcompView) 
library(patchwork)

# ==========================================
# 1. 定义带字母标记的 Alpha 绘图函数 (修复版)
# ==========================================
plot_shannon_with_letters_v2 <- function(data_path, metadata, title, y_label) {
  
  # A. 读取数据
  df <- read.csv(data_path, row.names = 1, check.names = FALSE)
  df <- df[, metadata$Sample]
  
  # B. 计算 Shannon 指数
  shannon_val <- diversity(t(df), index = "shannon")
  
  # C. 整合数据
  plot_df <- data.frame(Sample = names(shannon_val), value = as.numeric(shannon_val)) %>%
    left_join(metadata, by = "Sample") %>%
    mutate(Group = factor(Group, levels = c("IT-C", "IT-H", "ST-C", "ST-H")))
  
  # D. 统计分析
  # 1. 整体检验 (Kruskal-Wallis)
  kw_test <- plot_df %>% kruskal_test(value ~ Group)
  
  # 2. 两两比较 (Wilcoxon + BH校正)
  pwc <- plot_df %>% 
    wilcox_test(value ~ Group, p.adjust.method = "BH")
  
  # E. 生成字母标记 (核心修复部分)
  # 为了防止 multcompLetters 报错，临时替换组名中的横杠为下划线
  p_val_vec <- pwc$p.adj
  # 构造 IT_C-IT_H 这种格式
  names(p_val_vec) <- paste0(gsub("-", "_", pwc$group1), "-", gsub("-", "_", pwc$group2))
  
  # 计算字母 (只有当整体有差异时字母才有意义，但通常全画)
  cld_res <- multcompLetters(p_val_vec)
  cld <- cld_res$Letters
  
  # 将下划线还原回横杠，匹配回原始 Group
  cld_df <- data.frame(
    Group = gsub("_", "-", names(cld)), 
    Letters = as.character(cld),
    stringsAsFactors = FALSE
  )
  
  # 计算每组的最大值，用于确定字母位置
  plot_summary <- plot_df %>%
    group_by(Group) %>%
    summarise(max_val = max(value), .groups = "drop") %>%
    left_join(cld_df, by = "Group")
  
  # F. 配色
  color_pal <- c("IT-C" = "#FC8D62", "IT-H" = "#E41A1C", "ST-C" = "#4DBBD5", "ST-H" = "#377EB8")
  
  # G. 绘图
  p <- ggplot(plot_df, aes(x = Group, y = value, fill = Group)) +
    geom_violin(alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.3, outlier.shape = NA, alpha = 0.8, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.5) +
    # 动态调整字母高度：在最大值基础上增加 5% 的高度，防止压线
    geom_text(data = plot_summary, aes(x = Group, y = max_val + (max(plot_df$value)*0.05), label = Letters), 
              vjust = 0, fontface = "bold", size = 5) +
    scale_fill_manual(values = color_pal) +
    labs(title = title, x = "", y = y_label, 
         subtitle = paste0("Kruskal-Wallis, p = ", signif(kw_test$p, 3))) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "none"
    )
  
  return(list(plot = p, stats = pwc, data = plot_df))
}

# ==========================================
# 2. 重新运行分析
# ==========================================
metadata <- read.csv("metadata.csv")
rownames(metadata) <- metadata$Sample

# 这里的 metadata 和文件名请确保正确
res_gen <- plot_shannon_with_letters_v2(
  "Experimental_Genus_RelAbund_Final.csv", 
  metadata, 
  "Taxonomic Alpha Diversity", 
  "Shannon Index (Genus level)"
)

res_ko <- plot_shannon_with_letters_v2(
  "Experimental_KO_CSS_Normalized_Final.csv", 
  metadata, 
  "Functional Alpha Diversity", 
  "Shannon Index (KO level)"
)

# 拼图展示
library(patchwork)
final_alpha_plot <- res_gen$plot + res_ko$plot + plot_annotation(tag_levels = 'A')
print(final_alpha_plot)

# ==========================================
# 3. 结果导出
# ==========================================
ggsave("Figure_S_Experimental_Shannon_Letters.pdf", final_alpha_plot, width = 10, height = 6)
write.csv(res_gen$stats, "Stats_Wilcox_Shannon_Genus.csv")
write.csv(res_ko$stats, "Stats_Wilcox_Shannon_KO.csv")

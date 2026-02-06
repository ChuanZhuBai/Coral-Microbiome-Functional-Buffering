# 加载必要包
library(vegan)
library(ggplot2)
library(tidyverse)
library(patchwork)

# ==========================================
# 1. 核心处理函数定义
# ==========================================

# 数据过滤与预处理：保留整数计数并过滤极低丰度噪音
preprocess_data <- function(file_path, is_ko = FALSE) {
  data <- read.csv(file_path, header = TRUE, row.names = 1, check.names = FALSE)
  if(is_ko) {
    data <- data[!rownames(data) %in% c("Unclassified", "unclassified", "Unassigned"), ]
  }
  # 四舍五入为整数（rarefy要求）
  data <- round(data)
  
  # 过滤掉在所有样本中相对丰度都低于 0.001% 的项，提高计算效率
  rel_abund <- sweep(data, 2, colSums(data), "/")
  keep <- rowSums(rel_abund > 0.00001) >= 1
  return(data[keep, ])
}

# 快速步进式计算稀释曲线数据
calculate_rare_table <- function(data, steps = 30) {
  data_t <- t(data)
  results <- list()
  
  for(i in 1:nrow(data_t)) {
    total_reads <- sum(data_t[i, ])
    # 生成步进序列
    seqs <- seq(100, total_reads, length.out = steps)
    richness <- sapply(seqs, function(x) rarefy(data_t[i, ], x))
    
    results[[rownames(data_t)[i]]] <- data.frame(
      reads = seqs,
      richness = richness,
      sample = rownames(data_t)[i],
      group = ifelse(grepl("IT", rownames(data_t)[i]), "Intertidal", "Subtidal")
    )
  }
  return(do.call(rbind, results))
}

# 计算 Good's Coverage (用于客观评价饱和度)
get_coverage_stats <- function(data) {
  data_t <- t(data)
  coverage <- sapply(1:nrow(data_t), function(i) {
    n <- sum(data_t[i, ])
    f1 <- sum(data_t[i, ] == 1) # 单例数
    1 - (f1 / n)
  })
  return(data.frame(Sample = rownames(data_t), Coverage = coverage))
}

# ==========================================
# 2. 执行计算
# ==========================================

# 读取并过滤数据
cat("正在读取并过滤数据...\n")
genus_filtered <- preprocess_data("genus_abund.csv")
ko_filtered <- preprocess_data("ko_bases.csv", is_ko = TRUE)

# 计算稀释曲线
cat("正在计算稀释曲线数据 (81 样本)...\n")
df_genus_rare <- calculate_rare_table(genus_filtered)
df_ko_rare <- calculate_rare_table(ko_filtered)

# 计算覆盖率
cov_genus <- get_coverage_stats(genus_filtered)
cov_ko <- get_coverage_stats(ko_filtered)

# ==========================================
# 3. 顶级期刊风格绘图
# ==========================================

plot_rare_publication <- function(df_rare, title, y_lab) {
  # 计算组平均线
  group_avg <- df_rare %>%
    group_by(group, reads = round(reads, -3)) %>%
    summarise(mean_rich = mean(richness), .groups = "drop")
  
  ggplot() +
    # 绘制 81 条背景细线 (半透明)
    geom_line(data = df_rare, aes(x = reads, y = richness, group = sample, color = group), 
              alpha = 0.15, size = 0.2) +
    # 绘制组平均趋势线 (加粗)
    geom_smooth(data = df_rare, aes(x = reads, y = richness, color = group), 
                method = "gam", size = 1.2, se = FALSE) +
    # 颜色设置 (SA常用配色)
    scale_color_manual(values = c("Intertidal" = "#E64B35FF", "Subtidal" = "#4DBBD5FF")) +
    scale_x_continuous(labels = scales::scientific) +
    labs(title = title, x = "Sequencing Depth (Reads)", y = y_lab, color = "Habitat") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "bold", size = 10),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )
}

p1 <- plot_rare_publication(df_genus_rare, "Bacterial Genus Saturation", "Number of Genera")
p2 <- plot_rare_publication(df_ko_rare, "Functional KO Saturation", "Number of KOs")

# 合并展示
final_plot <- (p1 / p2) + plot_layout(guides = 'collect') & theme(legend.position = "bottom")
print(final_plot)

# ==========================================
# 4. 生成自动评估报告 (用于 Discussion 写作)
# ==========================================

cat("\n=== 数据饱和度报告 ===\n")
cat("Genus 平均覆盖率:", round(mean(cov_genus$Coverage)*100, 2), "%\n")
cat("KO 功能平均覆盖率:", round(mean(cov_ko$Coverage)*100, 2), "%\n")

Mif(mean(cov_genus$Coverage) > 0.98) {
  cat("结论：数据已高度饱和，当前 10G 深度已足以代表群落多样性。\n")
}

# 保存图片
ggsave("Figure_S1_Rarefaction_Professional.pdf", final_plot, width = 7, height = 9)

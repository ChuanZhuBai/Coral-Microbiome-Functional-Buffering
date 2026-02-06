#IT和ST两组的PCoA分析
# 1. 加载必要包
library(vegan)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(Cairo)
dir()
# 2. 数据准备
# 读取物种绝对丰度表（物种×样本）
species <- as.matrix(read.csv("genus.csv", row.names = 1, check.names = FALSE))
species_mat <- species[rowSums(species) > 0.0001, ]
# 读取元数据 (更新为使用SampleID列)
metadata <- read.csv("metadata.csv")
rownames(metadata) <- metadata$SampleID  # 注意这里改为SampleID
# 3. 计算Bray-Curtis距离
bray_dist <- vegdist(t(species_mat), method = "bray")
bray_dist_matrix <- as.matrix(bray_dist)
# 4. PCoA分析
pcoa_result <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_scores <- as.data.frame(pcoa_result$points)
colnames(pcoa_scores) <- c("PCoA1", "PCoA2")
# 5. 添加分组信息并去除离群点
pcoa_scores$Group <- metadata[rownames(pcoa_scores), "Group"]
pcoa_scores$Habitat <- metadata[rownames(pcoa_scores), "Habitat"]
# 离群点检测（基于马氏距离）
mah_dist <- mahalanobis(pcoa_scores[,1:2], 
                        colMeans(pcoa_scores[,1:2]), 
                        cov(pcoa_scores[,1:2]))
pcoa_scores$Outlier <- mah_dist > qchisq(0.975, df=2)
pcoa_clean <- subset(pcoa_scores, !Outlier)
# 添加分组信息
pcoa_scores$Group <- metadata[rownames(pcoa_scores), "Group"]
pcoa_scores$Habitat <- metadata[rownames(pcoa_scores), "Habitat"]
# 5. PERMANOVA检验
permanova <- adonis2(bray_dist ~ Group + Habitat, data = metadata, permutations = 9999)
# 6. 绘图优化
# 6. 绘图优化
p <- ggplot(pcoa_scores, aes(x = PCoA1, y = PCoA2, color = Group, shape = Habitat)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.2) +
  geom_point(size = 3) +
  labs(x = sprintf("PCoA1 (%.1f%%)", pcoa_result$eig[1]/sum(pcoa_result$eig)*100),
       y = sprintf("PCoA2 (%.1f%%)", pcoa_result$eig[2]/sum(pcoa_result$eig)*100)) +
  scale_color_manual(values = c("IT" = "#E64B35", 
                                "ST" = "#4DBBD5")) +
  theme_classic()

# 提取 PERMANOVA 结果（直接使用 permanova 对象）
permanova_r2 <- permanova$R2[1]  # 第一个因子（Group）的 R²
permanova_p <- permanova$`Pr(>F)`[1]  # 第一个因子（Group）的 P 值

# 格式化 P 值显示（科学计数法或 <0.001）
format_pval <- function(p) {
  if (p < 0.001) {
    return("P < 0.001")
  } else {
    return(sprintf("P = %.3f", p))
  }
}

# 添加 PERMANOVA 结果注释
p <- p + 
  annotate("text", 
           x = min(pcoa_scores$PCoA1), 
           y = max(pcoa_scores$PCoA2), 
           label = sprintf("PERMANOVA:\nR² = %.2f\n%s", 
                           permanova_r2, 
                           format_pval(permanova_p)),
           hjust = 0, 
           vjust = 1, 
           size = 3)

print(p)
# 7. 保存图形
ggsave("PCoA_Plot-new.pdf", plot = p, width = 8.6, height = 5, units = "cm", device = cairo_pdf)

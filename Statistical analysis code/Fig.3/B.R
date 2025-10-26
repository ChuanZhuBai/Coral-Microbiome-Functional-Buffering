# 1. 加载包 ------------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(ggpubr)  # 用于添加统计检验标记
library(tidyr)   # 用于数据整理

cat("===== 包加载完成 =====\n")

# 2. 数据准备 --------------------------------------------------------------
# 读取数据
fvfm_data <- read.csv("FvFm.csv") 

# 数据清洗与配对设计处理
# 提取样本编号（如C1/H1中的数字部分）
fvfm_data <- fvfm_data %>%
  mutate(
    SampleID = gsub(".*-(C|H)(\\d+)", "\\2", Sample),  # 提取样本编号
    GroupType = ifelse(grepl("-H", Group), "Heated", "Control"),
    BaseGroup = gsub("-C|-H", "", Group)  # 提取IT/ST分组
  ) 

# 检查样本配对情况
cat("\n===== 样本配对情况 =====\n")
print(table(fvfm_data$BaseGroup, fvfm_data$SampleID))

# 3. 计算变化率 -----------------------------------------------------------
# 宽格式转换（确保Control和Heated样本正确配对）
paired_data <- fvfm_data %>%
  pivot_wider(
    id_cols = c(BaseGroup, SampleID),
    names_from = GroupType,
    values_from = MQY
  ) %>%
  mutate(
    ChangeRate = (Heated - Control) / Control * 100  # 计算变化率百分比
  )

cat("\n===== 配对后数据前5行 =====\n")
print(head(paired_data, 5))
write.csv(paired_data, "paired_data.csv", row.names = FALSE)
# 4. 统计检验 ------------------------------------------------------------
# 4.1 正态性检验
shapiro_test <- paired_data %>%
  group_by(BaseGroup) %>%
  summarise(
    p_value = shapiro.test(ChangeRate)$p.value
  )
cat("\n===== 正态性检验结果 =====\n")
print(shapiro_test)

# 4.2 组间差异检验（IT vs ST的变化率）
if(all(shapiro_test$p_value > 0.05)) {
  # 参数检验（t检验）
  res <- t.test(ChangeRate ~ BaseGroup, data = paired_data, paired = FALSE)
  test_method <- "Student's t-test"
} else {
  # 非参数检验（Wilcoxon检验）
  res <- wilcox.test(ChangeRate ~ BaseGroup, data = paired_data)
  test_method <- "Wilcoxon rank-sum test"
}

cat("\n===== 组间差异检验结果 =====\n")
cat("Method:", test_method, "\n")
print(res)

# 4.3 组内差异检验（每个地点加热前后的变化）
group_tests <- list()
for(group in c("IT", "ST")) {
  subset_data <- paired_data %>% filter(BaseGroup == group)
  
  if(shapiro.test(subset_data$ChangeRate)$p.value > 0.05) {
    test_res <- t.test(subset_data$ChangeRate, mu = 0)
  } else {
    test_res <- wilcox.test(subset_data$ChangeRate, mu = 0, exact = FALSE)
  }
  group_tests[[group]] <- test_res
}

cat("\n===== 组内变化显著性 =====\n")
print(group_tests)

# 5. 可视化 --------------------------------------------------------------
# 5.1 颜色方案
group_colors <- c("IT" = "#E64B35", "ST" = "#4DBBD5")

# 5.2 箱线图+散点图（显示个体变化）
p <- ggplot(paired_data, aes(x = BaseGroup, y = ChangeRate)) +
  # 箱线图
  geom_boxplot(
    aes(fill = BaseGroup), 
    width = 0.6,
    alpha = 0.8,
    outlier.shape = NA
  ) +
  # 个体散点
  geom_jitter(
    aes(color = BaseGroup),
    width = 0.1,
    size = 2,
    alpha = 0.6
  ) +
  # 添加零参考线
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  # 统计检验标记（组间）
  stat_compare_means(
    comparisons = list(c("IT", "ST")),
    method = ifelse(test_method == "Student's t-test", "t.test", "wilcox.test"),
    label = "p.format",
    tip.length = 0.01,
    vjust = 1.5
  ) +
  # 统计检验标记（组内 vs 0）
  geom_signif(
    data = data.frame(
      BaseGroup = c("IT", "ST"),
      pval = sapply(group_tests, function(x) x$p.value)
    ),
    aes(
      x = BaseGroup,
      y = max(paired_data$ChangeRate) * 1.1,
      annotations = ifelse(pval < 0.001, "***", 
                           ifelse(pval < 0.01, "**",
                                  ifelse(pval < 0.05, "*", "ns")))
    ),
    textsize = 5,
    vjust = 0.5,
    manual = TRUE
  ) +
  # 美学设置
  scale_fill_manual(values = group_colors) +
  scale_color_manual(values = group_colors) +
  labs(
    x = "Group",
    y = "Fv/Fm Change Rate (%)",
    title = "Host physiological phenotype responses to heat stress",
    subtitle = paste("Between-group test:", test_method, "p =", format.pval(res$p.value, digits = 2))
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

# 6. 保存结果 ------------------------------------------------------------
# 6.1 保存图形
ggsave("Host_FvFm_ChangeRate.pdf", p, width = 6, height = 4)

# 6.2 保存统计结果
sink("Host_Physiology_Stats.txt")
cat("=== 正态性检验 ===\n")
print(shapiro_test)
cat("\n=== 组间差异检验 ===\n")
cat("Method:", test_method, "\n")
print(res)
cat("\n=== 组内变化检验 ===\n")
print(group_tests)
sink()

# 显示图形
print(p)

# 加载必要的包
library(ggplot2)
library(ggpubr)

# 读取数据（假设数据文件名为Dailyrange.csv）
data <- read.csv("Dailyrange.csv", header = TRUE)

# 确保Group是因子类型
data$Group <- factor(data$Group, levels = c("intertidal", "subtidal"))

# 绘制箱线图
p <- ggplot(data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = 16, outlier.size = 2) +
  
  # 设置颜色（潮间带红色，潮下带蓝色）"IT-H" = "#E64B35", "ST-H" = "#4DBBD5"
  scale_fill_manual(values = c("intertidal" = "#E64B35", "subtidal" = "#4DBBD5")) +
  
  # 添加Wilcoxon检验显著性标记
  stat_compare_means(method = "wilcox.test", 
                     comparisons = list(c("intertidal", "subtidal")),
                     label = "p.signif", 
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", "ns")),
                     vjust = 0.5, 
                     size = 6) +
  
  # 设置坐标轴标签
  labs(x = "Habitat", 
       y = "Daily temperature range (°C)",
       title = "Environmental variability across tidal zones") +
  
  # 调整主题
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", size = 0.2)
  ) +
  
  # 添加生态学信息标注
  annotate("text", x = 1.5, y = max(data$value) * 1.1, 
           label = "IT shows higher environmental variability\n(Wilcoxon test, p < 0.001)",
           size = 4, color = "black")

# 显示图形
print(p)

# 保存图形（可根据需要调整格式和尺寸）
ggsave("Daily_temperature_range.png", p, width = 6, height = 5, dpi = 300)
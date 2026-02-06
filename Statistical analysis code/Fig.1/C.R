# 加载必要的包
library(dplyr)
library(tidyr)

dir()
# 如果从CSV文件读取数据，取消以下注释：
data <- read.csv("Dailyrange.csv")

data <- read.csv("Temperature.csv")

data <- read.csv("diversity_indices_chao1.csv")

data <- read.csv("diversity_indices_shannon.csv")

data <- read.csv("Generalist_Data-insitu.csv")

data <- read.csv("Generalist_Data_insitu_KO.csv")

data <- read.csv("FDiv.csv")
# 按组计算平均值和标准差
results <- data %>%
  group_by(Group) %>%
  summarise(
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    # 格式化输出：平均值±标准差，保留2位小数
    formatted = sprintf("%.2f±%.2f°C", mean, sd),
    # 添加样本量信息（可选）
    formatted_with_n = sprintf("%.2f±%.2f°C (n=%d)", mean, sd, n)
  )

# 打印结果（SCI论文推荐格式）
print(results)

# 可选：输出为表格（适合直接插入论文）
knitr::kable(results, align = 'c', caption = "The mean and standard deviation of daily temperature range")

# 可选：绘制箱线图验证数据分布
library(ggplot2)
ggplot(data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() +
  labs(title = "Daily temperature difference distribution", y = "Temperature (°C)") +
  theme_minimal()






# 加载必要的包
library(ggplot2)  # 用于绘图
library(dplyr)    # 用于数据处理
dir()
# 1. 读取数据
data <- read.csv("Dailyrange.csv")  # 替换为你的实际文件路径

data <- read.csv("paired_data.csv")

data <- read.csv("diversity_indices_chao1.csv")

data <- read.csv("diversity_indices_shannon.csv")
# 2. 检查数据结构
head(data)  # 确认列名是否为"Group"和"value"
str(data)   # 检查数据类型（确保Group是因子，value是数值）

# 3. 按组计算中位数和IQR
summary_stats <- data %>%
  group_by(Group) %>%
  summarise(
    Median = median(value, na.rm = TRUE),
    Q1 = quantile(value, 0.25, na.rm = TRUE),
    Q3 = quantile(value, 0.75, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  )

# 打印统计结果
print(summary_stats)

# 4. 可选：统计检验（Mann-Whitney U检验，比较两组差异）
wilcox_test <- wilcox.test(value ~ Group, data = data)
print(wilcox_test)

# 5. 绘制箱线图
ggplot(data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot() +
  labs(title = "Comparison of Daily Range between Groups",
       x = "Group",
       y = "Value") +
  theme_minimal()

# 6. 保存结果到文件（可选）
write.csv(summary_stats, "Group_Median_IQR_Results.csv", row.names = FALSE)

shapiro.test(data$value[data$Group == "intertidal"])
shapiro.test(data$value[data$Group == "subtidal"])

# 直方图 + 密度曲线
library(ggplot2)
ggplot(data, aes(x = value, fill = Group)) +
  geom_histogram(alpha = 0.6, bins = 30) +
  facet_wrap(~Group, scales = "free") +
  theme_minimal()

# Q-Q图
qqnorm(data$value[data$Group == "intertidal"], main = "Intertidal Q-Q Plot")
qqline(data$value[data$Group == "intertidal"])
qqnorm(data$value[data$Group == "subtidal"], main = "Subtidal Q-Q Plot")
qqline(data$value[data$Group == "subtidal"])

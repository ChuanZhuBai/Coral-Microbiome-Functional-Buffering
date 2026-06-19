# =============================================================================
# Diurnal (daily) temperature range: Intertidal vs Subtidal
# =============================================================================
library(ggplot2)
library(ggpubr)
library(dplyr)

habitat_colors <- c("intertidal" = "#E64B35FF", "subtidal" = "#4DBBD5FF")

data <- read.csv("Dailyrange.csv", header = TRUE)
data$Group <- factor(data$Group, levels = c("intertidal", "subtidal"))
stopifnot(!anyNA(data$Group))

data$value <- as.numeric(data$value)
if (anyNA(data$value)) stop("Non-numeric or missing values in 'value' column.")
message(sprintf("Rows: %d (intertidal=%d, subtidal=%d)",
                nrow(data), sum(data$Group=="intertidal"), sum(data$Group=="subtidal")))

# =============================================================================
# Summary table: per-habitat diurnal temperature range statistics + Wilcoxon
#   (mean +/- SD, median, IQR, min-max range, n) 
# =============================================================================
summary_tab <- data %>%
  group_by(Group) %>%
  summarise(
    n        = n(),
    mean     = mean(value),
    sd       = sd(value),
    median   = median(value),
    Q1       = quantile(value, 0.25),
    Q3       = quantile(value, 0.75),
    min      = min(value),
    max      = max(value),
    .groups  = "drop"
  ) %>%
  mutate(
    `mean +/- SD` = sprintf("%.2f +/- %.2f", mean, sd),
    `range (min-max)` = sprintf("%.1f - %.1f", min, max),
    Habitat = ifelse(Group == "intertidal", "Intertidal", "Subtidal")
  ) %>%
  select(Habitat, n, mean, sd, `mean +/- SD`, median, Q1, Q3, min, max, `range (min-max)`)


wt <- wilcox.test(value ~ Group, data = data)
eff_n <- nrow(data)

z <- qnorm(wt$p.value / 2, lower.tail = FALSE)
r_eff <- z / sqrt(eff_n)

cat("\n===== Diurnal temperature range — per-habitat summary =====\n")
print(as.data.frame(summary_tab), row.names = FALSE, digits = 4)
cat(sprintf("\nWilcoxon rank-sum: W=%.1f, p=%.3g | effect size r=%.3f\n",
            wt$statistic, wt$p.value, r_eff))
cat(sprintf("Intertidal: %.2f +/- %.2f C (range %.1f-%.1f)\n",
            summary_tab$mean[summary_tab$Habitat=="Intertidal"],
            summary_tab$sd[summary_tab$Habitat=="Intertidal"],
            summary_tab$min[summary_tab$Habitat=="Intertidal"],
            summary_tab$max[summary_tab$Habitat=="Intertidal"]))
cat(sprintf("Subtidal:   %.2f +/- %.2f C (range %.1f-%.1f)\n",
            summary_tab$mean[summary_tab$Habitat=="Subtidal"],
            summary_tab$sd[summary_tab$Habitat=="Subtidal"],
            summary_tab$min[summary_tab$Habitat=="Subtidal"],
            summary_tab$max[summary_tab$Habitat=="Subtidal"]))


summary_export <- summary_tab %>%
  mutate(Wilcoxon_p = signif(wt$p.value, 3), effect_size_r = round(r_eff, 3))
write.csv(summary_export, "TableS_diurnal_temperature_range_summary.csv", row.names = FALSE)
message("Wrote TableS_diurnal_temperature_range_summary.csv")


p <- ggplot(data, aes(x = Group, y = value, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = 16, outlier.size = 2) +
  scale_fill_manual(values = habitat_colors) +
  scale_x_discrete(labels = c("intertidal" = "Intertidal",
                              "subtidal"   = "Subtidal")) +
  
  stat_compare_means(
    method      = "wilcox.test",
    comparisons = list(c("intertidal", "subtidal")),
    label       = "p.signif",
    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                       symbols   = c("***", "**", "*", "ns")),
    size = 6
  ) +
  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Habitat", y = "Diurnal temperature range (\u00B0C)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    axis.text          = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2)
  )

print(p)
ggsave("Figure_1B_daily_temperature_range.png", p, width = 4, height = 5, dpi = 600)
ggsave("Figure_1B_daily_temperature_range.pdf", p, width = 4, height = 5)
# =============================================================================
# Fig. 1E — Heat-induced change in Fv/Fm relative to control, Intertidal vs Subtidal.
# DESIGN: split-colony. Each colony was split into one Control (29 C) and one
#   Heat (33 C) fragment.
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

habitat_colors <- c("IT" = "#E64B35FF", "ST" = "#4DBBD5FF")
habitat_labels <- c("IT" = "Intertidal", "ST" = "Subtidal")


sep <- if (grepl("\t", readLines("FvFm.csv", n = 1))) "\t" else ","
d <- read.csv("FvFm.csv", header = TRUE, sep = sep, check.names = FALSE)

d <- d %>%
  mutate(
    BaseGroup = gsub("-C$|-H$", "", Group),            # IT / ST
    Treatment = ifelse(grepl("-H$", Group), "Heated", "Control"),
    ColonyID  = sub("^.*[CH]", "", Sample)             
  )
stopifnot(all(d$BaseGroup %in% c("IT", "ST")),
          all(d$Treatment %in% c("Control", "Heated")))


wide <- d %>%
  pivot_wider(id_cols = c(BaseGroup, ColonyID),
              names_from = Treatment, values_from = MQY)
if (anyNA(wide$Control) || anyNA(wide$Heated))
  stop("Some colonies lack a matching Control/Heat fragment (pairing not 1:1). ",
       "Check the sample IDs in FvFm.csv.")

wide <- wide %>% mutate(ChangeRate = (Heated - Control) / Control * 100)
wide$BaseGroup <- factor(wide$BaseGroup, levels = c("IT", "ST"))


between_res <- wilcox.test(ChangeRate ~ BaseGroup, data = wide)


within_res <- lapply(c("IT", "ST"), function(g) {
  sub <- wide[wide$BaseGroup == g, ]
  wilcox.test(sub$Heated, sub$Control, paired = TRUE, exact = FALSE)
})
names(within_res) <- c("IT", "ST")

med <- wide %>% group_by(BaseGroup) %>%
  summarise(n = n(), median_change_pct = median(ChangeRate), .groups = "drop")

sink("Figure_1E_stats.txt")
cat("Design: split-colony (paired within habitat)\n\nMedian % change by habitat:\n")
print(as.data.frame(med))
cat("\nBetween-habitat Wilcoxon rank-sum (change rates):\n"); print(between_res)
cat("\nWithin-habitat paired Wilcoxon signed-rank (heat vs control):\n")
print(within_res)
sink()
print(as.data.frame(med)); print(between_res)
write.csv(wide, "Figure_1E_change_rates.csv", row.names = FALSE)


p <- ggplot(wide, aes(x = BaseGroup, y = ChangeRate)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_boxplot(aes(fill = BaseGroup), width = 0.5, alpha = 0.85,
               outlier.shape = NA) +
  geom_jitter(aes(color = BaseGroup), width = 0.12, size = 2, alpha = 0.6) +
  stat_compare_means(method = "wilcox.test",
                     comparisons = list(c("IT", "ST")),
                     label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols   = c("***", "**", "*", "ns")),
                     size = 6) +
  scale_fill_manual(values = habitat_colors) +
  scale_color_manual(values = habitat_colors) +
  scale_x_discrete(labels = habitat_labels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  labs(x = "Habitat", y = expression(Delta*"Fv/Fm relative to control (%)")) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        axis.text = element_text(color = "black"))

print(p)
ggsave("Figure_1E_FvFm_change.png", p, width = 4, height = 5, dpi = 600)
ggsave("Figure_1E_FvFm_change.pdf", p, width = 4, height = 5)
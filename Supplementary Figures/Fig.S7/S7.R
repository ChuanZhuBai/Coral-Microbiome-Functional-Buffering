# =============================================================================
# Fig. S5 — robustness of the generalist-richness vs FRI relationship to the
# =============================================================================
library(tidyverse); library(patchwork)

hab_cols <- c(Intertidal = "#E64B35FF", Subtidal = "#4DBBD5FF")

read_auto <- function(file) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- load -------------------------------------------------------------------
niche   <- read_auto("Figure_2D_niche_breadth_values.csv")          # genus, B, Group, Class
genus_abun <- read_auto("genus.csv"); rownames(genus_abun) <- genus_abun[[1]]; genus_abun[[1]] <- NULL
genus_abun <- as.matrix(genus_abun); storage.mode(genus_abun) <- "numeric"
micro <- read_auto("SourceData_2HJ_per_sample.csv") %>%             # Sample, Group, Habitat, Global_FRI
  select(Sample, Group, Habitat, Global_FRI)

# detected genera (relative abundance > 0) per sample
present_long <- as.data.frame(sweep(genus_abun, 2, colSums(genus_abun), "/")) %>%
  rownames_to_column("genus") %>%
  pivot_longer(-genus, names_to = "Sample", values_to = "ab") %>%
  filter(ab > 0) %>% select(genus, Sample)

# ---- generalist richness at a given pooled-B percentile ---------------------
gen_rich_at <- function(pct) {
  thr <- as.numeric(quantile(niche$B, pct, na.rm = TRUE))
  gens <- niche %>% filter(B >= thr) %>% pull(genus) %>% unique()
  rich <- present_long %>% filter(genus %in% gens) %>%
    group_by(Sample) %>% summarise(Richness = n_distinct(genus), .groups = "drop")
  list(threshold = thr, n_gen = length(gens), rich = rich)
}

# Top 10% = 90th percentile ; Top 20% = 80th percentile
defs <- list("Top 10%" = 0.90, "Top 20%" = 0.80)

run_one <- function(label, pct) {
  g <- gen_rich_at(pct)
  d <- micro %>% inner_join(g$rich, by = "Sample") %>%
    mutate(Richness = replace_na(Richness, 0),
           Habitat = factor(Habitat, levels = c("Intertidal","Subtidal")))
  fit <- summary(lm(Global_FRI ~ Richness, data = d))
  message(sprintf("[%s] threshold(B)=%.3f, %d generalist genera | pooled R2=%.3f p=%.3g",
                  label, g$threshold, g$n_gen, fit$r.squared, fit$coefficients[2,4]))
  for (grp in c("IT","ST")) {
    dd <- d[d$Group == grp, ]
    ct <- cor.test(dd$Richness, dd$Global_FRI, method = "pearson")
    message(sprintf("   within %s: r=%.3f, p=%.3g (n=%d)", grp, ct$estimate, ct$p.value, nrow(dd)))
  }
  list(data = d, fit = fit, label = label)
}

res <- imap(defs, ~ run_one(.y, .x))

# ---- plot: one panel per threshold ------------------------------------------
make_panel <- function(r) {
  d <- r$data
  lab <- sprintf("Pooled: R² = %.2f, p = %.3g", r$fit$r.squared, r$fit$coefficients[2,4])
  ggplot(d, aes(Richness, Global_FRI)) +
    geom_smooth(method = "lm", color = "black", fill = "grey85", alpha = 0.5, linewidth = 1) +
    geom_smooth(aes(color = Habitat), method = "lm", se = FALSE,
                linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(fill = Habitat), shape = 21, size = 3.2, stroke = 0.4, alpha = 0.85) +
    annotate("text", x = min(d$Richness), y = max(d$Global_FRI),
             label = lab, hjust = 0, vjust = 1, size = 3.6, fontface = "bold") +
    scale_fill_manual(values = hab_cols) + scale_color_manual(values = hab_cols) +
    labs(title = sprintf("Generalists: %s niche breadth", r$label),
         x = "Richness of generalist genera (per sample)",
         y = "Mean Functional Redundancy Index") +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
          axis.text = element_text(color = "black"), legend.position = "top",
          legend.title = element_blank())
}

p_s5 <- make_panel(res[[1]]) + make_panel(res[[2]]) + plot_layout(guides = "collect") &
  theme(legend.position = "top")

# ---- save -------------------------------------------------------------------
src <- imap_dfr(res, ~ .x$data %>% mutate(Threshold = .y) %>%
                  select(Sample, Group, Habitat, Threshold, Richness, Global_FRI))
write.csv(src, "SourceData_S5_generalist_robustness.csv", row.names = FALSE)
ggsave("FigureS5_generalist_robustness.pdf", p_s5, width = 9, height = 4.4)
print(p_s5)
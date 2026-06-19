# =============================================================================
# regressions of coral thermotolerance on microbiome properties
# =============================================================================
library(tidyverse)

hab_cols <- c(Intertidal = "#E64B35FF", Subtidal = "#4DBBD5FF")

read_auto <- function(file) {
  sep <- if (grepl("\t", readLines(file, n = 1))) "\t" else ","
  read.csv(file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE)
}

# ---- Fv/Fm: build Sample name and retention rate ----------------------------
fvfm <- read_auto("Figure_1E_change_rates.csv") %>%
  mutate(Sample    = paste0(BaseGroup, ColonyID),          # IT + 1 -> "IT1"
         Retention = 100 + as.numeric(ChangeRate)) %>%      # = Heated/Control*100
  select(Sample, Retention)

# ---- microbiome metrics (FRI count + effective, generalist richness) --------
micro <- read_auto("SourceData_2HJ_per_sample.csv") %>%
  select(Sample, Group, Habitat, Generalist_Richness, Global_FRI, Global_FRI_eff)

# ---- merge: keep only samples present in BOTH (inner join) ------------------
dat <- inner_join(micro, fvfm, by = "Sample")
dat$Habitat <- factor(dat$Habitat, levels = c("Intertidal","Subtidal"))
message(sprintf("Matched samples (metagenome AND Fv/Fm): %d (IT=%d, ST=%d)",
                nrow(dat), sum(dat$Group=="IT"), sum(dat$Group=="ST")))

# ---- DIAGNOSTIC: compare FRI_count vs FRI_eff as predictor of retention -----
# (decides whether switching the x metric changes the within-habitat picture)
message("\n--- predictor comparison for retention (pooled R2 / within-habitat r,p) ---")
for (xv in c("Global_FRI", "Global_FRI_eff", "Generalist_Richness")) {
  fit <- summary(lm(reformulate(xv, "Retention"), data = dat))
  rIT <- with(dat[dat$Group=="IT",], cor.test(get(xv), Retention))
  rST <- with(dat[dat$Group=="ST",], cor.test(get(xv), Retention))
  message(sprintf("%-20s pooled R2=%.3f p=%.3g | IT r=%.3f p=%.3g | ST r=%.3f p=%.3g",
                  xv, fit$r.squared, fit$coefficients[2,4],
                  rIT$estimate, rIT$p.value, rST$estimate, rST$p.value))
}
message("---------------------------------------------------------------------\n")

# ---- generic regression reporter: pooled + within-habitat -------------------
report_reg <- function(d, xvar, yvar, tag) {
  f  <- as.formula(sprintf("%s ~ %s", yvar, xvar))
  fit <- summary(lm(f, data = d))
  message(sprintf("[%s] POOLED: R2=%.3f, p=%.3g, slope=%.3f",
                  tag, fit$r.squared, fit$coefficients[2,4], fit$coefficients[2,1]))
  for (g in c("IT","ST")) {
    dd <- d[d$Group==g,]
    ct <- cor.test(dd[[xvar]], dd[[yvar]], method = "pearson")
    message(sprintf("   within %s: r=%.3f, p=%.3g (n=%d)", g, ct$estimate, ct$p.value, nrow(dd)))
  }
  invisible(fit)
}

# ---- plotting helper: pooled line + habitat-coloured points + within R2 -----
make_plot <- function(d, xvar, yvar, xlab, ylab) {
  fit <- summary(lm(as.formula(sprintf("%s ~ %s", yvar, xvar)), data = d))
  lab_pool <- sprintf("Pooled: R² = %.2f, p = %.3g", fit$r.squared, fit$coefficients[2,4])
  ggplot(d, aes(.data[[xvar]], .data[[yvar]])) +
    geom_smooth(method = "lm", color = "black", fill = "grey85",
                alpha = 0.5, linewidth = 1) +                     # pooled
    geom_smooth(aes(color = Habitat), method = "lm", se = FALSE,
                linewidth = 0.7, linetype = "dashed") +           # within-habitat
    geom_point(aes(fill = Habitat), shape = 21, size = 3.5, stroke = 0.4, alpha = 0.85) +
    annotate("text", x = min(d[[xvar]]), y = max(d[[yvar]]),
             label = lab_pool, hjust = 0, vjust = 1, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = hab_cols) +
    scale_color_manual(values = hab_cols) +
    labs(x = xlab, y = ylab) +
    theme_classic(base_size = 12) +
    theme(axis.text = element_text(color = "black"), legend.position = "top",
          legend.title = element_blank())
}

# ---- 3G: retention ~ FRI ----------------------------------------------------
report_reg(dat, "Global_FRI", "Retention", "3G FRI")
p_3g <- make_plot(dat, "Global_FRI", "Retention",
                  "Whole-community FRI", "Fv/Fm retention rate (%)")

# ---- 3H: retention ~ generalist richness ------------------------------------
report_reg(dat, "Generalist_Richness", "Retention", "3H generalist")
p_3h <- make_plot(dat, "Generalist_Richness", "Retention",
                  "Richness of generalist genera (per sample)", "Fv/Fm retention rate (%)")

# ---- save -------------------------------------------------------------------
write.csv(dat, "SourceData_3G_3H_regression.csv", row.names = FALSE)
ggsave("Figure_3G_FRI_vs_retention.pdf",        p_3g, width = 6, height = 4.5)
ggsave("Figure_3H_generalist_vs_retention.pdf", p_3h, width = 6, height = 4.5)
print(p_3g); print(p_3h)
# Fig. S1. Rarefaction curves

library(vegan)
library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
if (length(script_file) == 1) {
  setwd(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
}

habitat_colors <- c(
  Intertidal = "#E64B35FF",
  Subtidal = "#4DBBD5FF"
)

read_count_table <- function(file) {
  counts <- read.csv(file, row.names = 1, check.names = FALSE)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "numeric"

  if (anyNA(counts)) {
    stop("Non-numeric or missing values detected in ", file)
  }
  if (any(counts < 0)) {
    stop("Negative counts detected in ", file)
  }
  if (any(abs(counts - round(counts)) > sqrt(.Machine$double.eps))) {
    stop("Rarefaction requires integer count data: ", file)
  }
  if (any(colSums(counts) == 0)) {
    stop("Samples with zero total counts detected in ", file)
  }

  counts <- round(counts)
  counts[rowSums(counts) > 0, , drop = FALSE]
}

summarise_counts <- function(counts, dataset) {
  data.frame(
    dataset = dataset,
    features = nrow(counts),
    samples = ncol(counts),
    min_reads = min(colSums(counts)),
    median_reads = median(colSums(counts)),
    max_reads = max(colSums(counts)),
    row.names = NULL
  )
}

get_habitat <- function(sample) {
  case_when(
    grepl("^IT", sample) ~ "Intertidal",
    grepl("^ST", sample) ~ "Subtidal",
    TRUE ~ NA_character_
  )
}

make_rarefaction_data <- function(counts, dataset, step = 1000) {
  rare_df <- rarecurve(t(counts), step = step, label = FALSE, tidy = TRUE)

  rare_df %>%
    transmute(
      dataset = dataset,
      sample = as.character(Site),
      reads = Sample,
      richness = Species,
      habitat = factor(get_habitat(as.character(Site)), levels = names(habitat_colors))
    )
}

get_coverage <- function(counts, dataset) {
  sample_table <- t(counts)
  depth <- rowSums(sample_table)
  singletons <- rowSums(sample_table == 1)

  data.frame(
    dataset = dataset,
    sample = rownames(sample_table),
    habitat = get_habitat(rownames(sample_table)),
    reads = depth,
    goods_coverage = 1 - singletons / depth
  )
}

plot_rarefaction <- function(data, title, y_label) {
  ggplot(data, aes(reads, richness, group = sample, color = habitat)) +
    geom_line(linewidth = 0.3, alpha = 0.65) +
    scale_color_manual(values = habitat_colors, drop = FALSE) +
    scale_x_continuous(labels = label_comma(accuracy = 1)) +
    labs(
      title = title,
      x = "Subsampled reads",
      y = y_label,
      color = "Habitat"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
}

genus_counts <- read_count_table("genus_abund.csv")
ko_counts <- read_count_table("ko_abund.csv")

input_summary <- bind_rows(
  summarise_counts(genus_counts, "Genus"),
  summarise_counts(ko_counts, "KO")
)

rarefaction_data <- bind_rows(
  make_rarefaction_data(genus_counts, "Genus"),
  make_rarefaction_data(ko_counts, "KO")
)

coverage_data <- bind_rows(
  get_coverage(genus_counts, "Genus"),
  get_coverage(ko_counts, "KO")
)

coverage_summary <- coverage_data %>%
  group_by(dataset) %>%
  summarise(
    samples = n(),
    mean_goods_coverage = mean(goods_coverage),
    min_goods_coverage = min(goods_coverage),
    .groups = "drop"
  )

write.csv(input_summary, "Figure_S1_input_summary.csv", row.names = FALSE)
write.csv(coverage_data, "Figure_S1_sample_depth_coverage.csv", row.names = FALSE)
write.csv(coverage_summary, "Figure_S1_coverage_summary.csv", row.names = FALSE)
write.csv(rarefaction_data, "Figure_S1_rarefaction_points.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Figure_S1_sessionInfo.txt")

print(input_summary)
print(coverage_summary)

p_genus <- plot_rarefaction(
  filter(rarefaction_data, dataset == "Genus"),
  "Genus-level rarefaction",
  "Observed genera"
)

p_ko <- plot_rarefaction(
  filter(rarefaction_data, dataset == "KO"),
  "KO-level rarefaction",
  "Observed KOs"
)

fig_s1 <- (p_genus / p_ko) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(fig_s1)

ggsave("Figure_S1_Rarefaction.png", fig_s1, width = 7, height = 8.5, dpi = 300)
ggsave("Figure_S1_Rarefaction.pdf", fig_s1, width = 7, height = 8.5)

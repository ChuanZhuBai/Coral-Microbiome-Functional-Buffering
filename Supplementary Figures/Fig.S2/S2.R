# Fig. S2. Field phylum-level taxonomic composition

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(scales)
library(ggsci)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
if (length(script_file) == 1) {
  setwd(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
}

top_n <- 10
field_groups <- c("Intertidal", "Subtidal")

phylum_counts <- read.csv("Phylum_abund.csv", row.names = 1, check.names = FALSE)
phylum_counts <- as.matrix(phylum_counts)
storage.mode(phylum_counts) <- "numeric"

metadata <- read_excel("Metadata_SampleKey.xlsx") %>%
  as.data.frame() %>%
  filter(Group_Description %in% field_groups) %>%
  mutate(
    Habitat = factor(Group_Description, levels = field_groups),
    Sample = Paper_Sample_ID
  )

if (!all(metadata$Raw_Sample_ID_in_Figshare %in% colnames(phylum_counts))) {
  stop("Some field samples in Metadata_SampleKey.xlsx are absent from Phylum_abund.csv.")
}

field_counts <- phylum_counts[, metadata$Raw_Sample_ID_in_Figshare, drop = FALSE]
if (any(colSums(field_counts) == 0)) {
  stop("Samples with zero phylum-level reads detected.")
}

sample_summary <- metadata %>%
  transmute(
    Raw_sample_ID = Raw_Sample_ID_in_Figshare,
    Sample = Sample,
    Habitat = Habitat,
    Phylum_assigned_reads = colSums(field_counts)[Raw_Sample_ID_in_Figshare]
  )

field_rel <- sweep(field_counts, 2, colSums(field_counts), "/")
colnames(field_rel) <- metadata$Sample
field_rel <- field_rel[rowSums(field_rel) > 0, , drop = FALSE]

unclassified_phyla <- rownames(field_rel)[grepl("unclassified", rownames(field_rel), ignore.case = TRUE)]
mean_abundance <- sort(rowMeans(field_rel), decreasing = TRUE)
top_phyla <- names(mean_abundance)[!names(mean_abundance) %in% unclassified_phyla][1:top_n]

input_summary <- data.frame(
  input_phyla = nrow(phylum_counts),
  input_samples = ncol(phylum_counts),
  field_samples = ncol(field_counts),
  intertidal_samples = sum(metadata$Habitat == "Intertidal"),
  subtidal_samples = sum(metadata$Habitat == "Subtidal"),
  field_phyla_with_nonzero_counts = nrow(field_rel),
  field_all_zero_phyla = nrow(phylum_counts) - nrow(field_rel),
  row.names = NULL
)

phylum_summary <- data.frame(
  Phylum = names(mean_abundance),
  Mean_relative_abundance = as.numeric(mean_abundance),
  Rank = seq_along(mean_abundance),
  Display = ifelse(names(mean_abundance) %in% top_phyla, names(mean_abundance), "Others"),
  row.names = NULL
)

plot_data <- field_rel %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Phylum") %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Relative_abundance") %>%
  left_join(
    metadata %>% select(Sample, Habitat),
    by = "Sample"
  ) %>%
  mutate(Phylum_plot = if_else(Phylum %in% top_phyla, Phylum, "Others")) %>%
  group_by(Sample, Habitat, Phylum_plot) %>%
  summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop") %>%
  mutate(
    Sample = factor(Sample, levels = metadata$Sample),
    Phylum_plot = factor(Phylum_plot, levels = c(top_phyla, "Others"))
  )

if (any(abs(tapply(plot_data$Relative_abundance, plot_data$Sample, sum) - 1) > 1e-8)) {
  stop("Relative abundance does not sum to 1 for all samples.")
}

write.csv(input_summary, "Figure_S2_input_summary.csv", row.names = FALSE)
write.csv(sample_summary, "Figure_S2_field_sample_summary.csv", row.names = FALSE)
write.csv(field_rel, "Figure_S2_field_phylum_relative_abundance.csv")
write.csv(plot_data, "Figure_S2_plot_data_top10_plus_others.csv", row.names = FALSE)
write.csv(phylum_summary, "Figure_S2_phylum_mean_relative_abundance.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Figure_S2_sessionInfo.txt")

phylum_colors <- pal_npg("nrc")(length(top_phyla))
color_map <- setNames(c(phylum_colors, "#BDBDBD"), c(top_phyla, "Others"))

fig_s2 <- ggplot(plot_data, aes(Sample, Relative_abundance, fill = Phylum_plot)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.08) +
  facet_grid(. ~ Habitat, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_map, drop = FALSE) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = "Field samples",
    y = "Relative abundance",
    fill = "Phylum"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.title.y = element_text(margin = margin(r = 6)),
    strip.background = element_rect(fill = "grey95", color = NA),
    strip.text = element_text(margin = margin(4, 4, 4, 4)),
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.35, "cm"),
    panel.spacing.x = unit(0.18, "cm")
  )

print(input_summary)
print(head(phylum_summary, 15))

print(fig_s2)

ggsave("Figure_S2_Field_Phylum_Composition.png", fig_s2, width = 11, height = 6, dpi = 300)
ggsave("Figure_S2_Field_Phylum_Composition.pdf", fig_s2, width = 11, height = 6)

# =============================================================================
# In-lab experiment phylum-level taxonomic composition
# =============================================================================
library(readxl); library(ggplot2); library(dplyr); library(tidyr)
library(tibble); library(scales); library(ggsci)

script_args <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_args[grepl("^--file=", script_args)])
if (length(script_file) == 1) {
  setwd(dirname(normalizePath(script_file, winslash = "/", mustWork = TRUE)))
}

top_n <- 10

exp_groups <- c("Intertidal Control", "Intertidal Heat",
                "Subtidal Control",   "Subtidal Heat")
group_short <- c("Intertidal Control" = "IT-C", "Intertidal Heat" = "IT-H",
                 "Subtidal Control"   = "ST-C", "Subtidal Heat"   = "ST-H")

phylum_counts <- read.csv("Phylum_abund.csv", row.names = 1, check.names = FALSE)
phylum_counts <- as.matrix(phylum_counts); storage.mode(phylum_counts) <- "numeric"


metadata <- read_excel("Metadata_SampleKey.xlsx") %>%
  as.data.frame() %>%
  filter(Group_Description %in% exp_groups) %>%
  mutate(Group   = factor(group_short[Group_Description],
                          levels = c("IT-C","IT-H","ST-C","ST-H")),
         Habitat = factor(Group_Description, levels = exp_groups),
         Sample  = Paper_Sample_ID)

if (nrow(metadata) != 24)
  message("NOTE: expected 24 experimental samples, found ", nrow(metadata))
if (!all(metadata$Raw_Sample_ID_in_Figshare %in% colnames(phylum_counts)))
  stop("Some experimental samples in SampleKey are absent from Phylum_abund.csv: ",
       paste(setdiff(metadata$Raw_Sample_ID_in_Figshare, colnames(phylum_counts)), collapse=", "))


exp_counts <- phylum_counts[, metadata$Raw_Sample_ID_in_Figshare, drop = FALSE]
if (any(colSums(exp_counts) == 0)) stop("Samples with zero phylum-level reads detected.")

sample_summary <- metadata %>%
  transmute(Raw_sample_ID = Raw_Sample_ID_in_Figshare, Sample = Sample,
            Group = Group, Habitat = Habitat,
            Phylum_assigned_reads = colSums(exp_counts)[Raw_Sample_ID_in_Figshare])

exp_rel <- sweep(exp_counts, 2, colSums(exp_counts), "/")
colnames(exp_rel) <- metadata$Sample
exp_rel <- exp_rel[rowSums(exp_rel) > 0, , drop = FALSE]


unclassified_phyla <- rownames(exp_rel)[grepl("unclassified", rownames(exp_rel), ignore.case = TRUE)]
mean_abundance <- sort(rowMeans(exp_rel), decreasing = TRUE)
top_phyla <- names(mean_abundance)[!names(mean_abundance) %in% unclassified_phyla][1:top_n]
top_phyla <- top_phyla[!is.na(top_phyla)]


input_summary <- data.frame(
  input_phyla = nrow(phylum_counts), input_samples = ncol(phylum_counts),
  experimental_samples = ncol(exp_counts),
  IT_C = sum(metadata$Group=="IT-C"), IT_H = sum(metadata$Group=="IT-H"),
  ST_C = sum(metadata$Group=="ST-C"), ST_H = sum(metadata$Group=="ST-H"),
  phyla_with_nonzero_counts = nrow(exp_rel),
  all_zero_phyla = nrow(phylum_counts) - nrow(exp_rel), row.names = NULL)

phylum_summary <- data.frame(
  Phylum = names(mean_abundance), Mean_relative_abundance = as.numeric(mean_abundance),
  Rank = seq_along(mean_abundance),
  Display = ifelse(names(mean_abundance) %in% top_phyla, names(mean_abundance), "Others"),
  row.names = NULL)


plot_data <- exp_rel %>% as.data.frame(check.names = FALSE) %>%
  rownames_to_column("Phylum") %>%
  pivot_longer(-Phylum, names_to = "Sample", values_to = "Relative_abundance") %>%
  left_join(metadata %>% select(Sample, Group, Habitat), by = "Sample") %>%
  mutate(Phylum_plot = if_else(Phylum %in% top_phyla, Phylum, "Others")) %>%
  group_by(Sample, Group, Habitat, Phylum_plot) %>%
  summarise(Relative_abundance = sum(Relative_abundance), .groups = "drop") %>%
  mutate(Sample = factor(Sample, levels = metadata$Sample),
         Phylum_plot = factor(Phylum_plot, levels = c(top_phyla, "Others")))

if (any(abs(tapply(plot_data$Relative_abundance, plot_data$Sample, sum) - 1) > 1e-8))
  stop("Relative abundance does not sum to 1 for all samples.")


write.csv(input_summary, "Figure_input_summary.csv", row.names = FALSE)
write.csv(sample_summary, "Figure_experimental_sample_summary.csv", row.names = FALSE)
write.csv(exp_rel, "Figure_experimental_phylum_relative_abundance.csv")
write.csv(plot_data, "Figure_plot_data_top10_plus_others.csv", row.names = FALSE)
write.csv(phylum_summary, "Figure_phylum_mean_relative_abundance.csv", row.names = FALSE)
writeLines(capture.output(sessionInfo()), "Figure_sessionInfo.txt")


phylum_colors <- pal_npg("nrc")(length(top_phyla))
color_map <- setNames(c(phylum_colors, "#BDBDBD"), c(top_phyla, "Others"))

fig_s6 <- ggplot(plot_data, aes(Sample, Relative_abundance, fill = Phylum_plot)) +
  geom_col(width = 0.95, color = "white", linewidth = 0.08) +
  facet_grid(. ~ Group, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = color_map, drop = FALSE) +
  scale_y_continuous(labels = label_percent(accuracy = 1), limits = c(0,1), expand = c(0,0)) +
  labs(x = "Experimental samples", y = "Relative abundance", fill = "Phylum") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_text(margin = margin(t = 6)),
        axis.title.y = element_text(margin = margin(r = 6)),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(margin = margin(4,4,4,4), face = "bold"),
        legend.position = "right", legend.title = element_text(size = 9),
        legend.text = element_text(size = 8), legend.key.size = unit(0.35, "cm"),
        panel.spacing.x = unit(0.18, "cm"))

print(input_summary); print(head(phylum_summary, 15)); print(fig_s6)
ggsave("Figure_Experiment_Phylum_Composition.png", fig_s6, width = 11, height = 6, dpi = 300)
ggsave("Figure_Experiment_Phylum_Composition.pdf", fig_s6, width = 11, height = 6)
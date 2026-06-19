# Functional buffering by bacterial generalists associated with coral holobiont thermotolerance

## Overview
This repository contains the custom R scripts and metadata necessary to reproduce the analysis presented in our paper. 

## Directory Structure
The repository is organized by figure numbers corresponding to the manuscript to ensure transparency and reproducibility:

*   `Statistical analysis code/` : Contains all the customized R scripts used for statistical analyses and figure generation.
    *   `Fig.1/` : Scripts for physiological data analysis (Figure 1).
    *   `Fig.2/` : Scripts for baseline microbiome architectures, niche breadth, and Functional Redundancy Index (FRI) (Figure 2).
    *   `Fig.3/` : Scripts for functional reorganization, Functional Retention Rate (FRR) (Figure 3).
    *   `Supplementary Figures/` : Scripts for supplementary analyses.

*   `Metadata_SampleKey.xlsx` : A comprehensive sample key explicitly mapping all raw bioinformatic output IDs (e.g., restI1A) to the standardized Sample IDs used in the manuscript (e.g., IT-C1).

## Data Availability
*   **Raw Sequence Data**: All raw metagenomic FASTQ files have been deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession number **PRJNA1474294**.
*   **Processed Data & Matrices**: Due to GitHub file size limits, all processed taxonomic and functional abundance matrices, ORF annotations, and mapping statistics are openly available in Figshare at: **[https://doi.org/10.6084/m9.figshare.30511487]**.

## Prerequisites
The scripts were developed and tested in R version 4.3.3. 
Required packages include: `tidyverse`, `vegan`, `ggpubr`, `metagenomeSeq`, `patchwork`, `igraph`.
## Repository structure

```
Statistical analysis R scripts/
├── Fig.1/
│   ├── 1.B/          # Diurnal temperature range (field)
│   └── 1.E/          # Heat-induced ΔFv/Fm (split-colony experiment)
├── Fig.2/
│   ├── 2A/           # Genus-level PCoA (field)
│   ├── 2B/           # Functional (KO) PCoA (field)
│   ├── 2C/           # Genus Shannon alpha diversity (field)
│   ├── 2D/           # Niche breadth (Levins B) + specialist/generalist classification
│   └── 2E&F&G/       # Functional Redundancy Index (FRI) and drivers
├── Fig.3/
│   ├── 3A&B/         # CPCoA of experimental heat-stress communities
│   ├── 3C&D/         # Microbiome stability (Bray-Curtis distance, control vs heat)
│   ├── 3E/           # Functional retention rate (FRR) under heat stress
│   └── 3F&G/         # Regression: coral thermotolerance ~ FRI / generalist richness
├── Supplementary Figures/
│   ├── Fig.S1/       # Rarefaction curves (field)
│   ├── Fig.S2/       # Field phylum composition
│   ├── Fig.S3/       # Genus co-occurrence networks (field)
│   ├── Fig.S4/       # KO functional diversity (experimental, CSS-normalised)
│   ├── Fig.S5/       # Differential abundance volcano plot (edgeR)
│   ├── Fig.S6/       # KEGG pathway enrichment (clusterProfiler ORA)
│   ├── Fig.S7/       # Generalist threshold robustness sensitivity
│   ├── Fig.S8/       # Experimental phylum composition
│   ├── Fig.S9/       # Experimental Shannon alpha diversity
│   └── Fig.S10/      # KO co-occurrence networks (experimental heat groups)
├── Bacterial metagenomic data/
│   └── rename.R      # File-renaming utility for raw metagenomic outputs

## Script ↔ figure mapping

> **Note on figure numbering:** The filenames written inside scripts reflect an earlier draft numbering. The *Figure label in code* column shows what the script actually writes to disk; the *FINAL manuscript figure* column gives the published numbering. Where the FINAL column is blank, please confirm against the final submitted manuscript.

| Script | What it generates | Input data | Output files | Figure label in code | FINAL manuscript figure (TO CONFIRM) |
|---|---|---|---|---|---|
| `Fig.1/1.B/1.B.R` | Diurnal temperature range boxplot (IT vs ST) + Wilcoxon summary table | `Dailyrange.csv` | `Figure_1B_daily_temperature_range.png/.pdf`; `TableS_diurnal_temperature_range_summary.csv` | Figure_1B | Fig. 1B |
| `Fig.1/1.E/1.E.R` | Heat-induced ΔFv/Fm boxplot (split-colony design, IT vs ST) | `FvFm.csv` | `Figure_1E_FvFm_change.png/.pdf`; `Figure_1E_change_rates.csv`; `Figure_1E_stats.txt` | Figure_1E | Fig. 1E |
| `Fig.2/2A/2A.R` | Genus-level Bray-Curtis PCoA (field communities, IT vs ST) + PERMANOVA | `genus.csv`; `metadata.csv` | `Figure_2A_genus_PCoA.pdf/.png` | Figure_2A | Fig. 2A |
| `Fig.2/2B/2B.R` | Functional (KO) Bray-Curtis PCoA (field, IT vs ST) + PERMANOVA | `kegg_abund.csv`; `metadata.csv` | `Figure_2B_KO_PCoA.pdf/.png` | Figure_2B | Fig. 2B |
| `Fig.2/2C/2C.R` | Genus-level Shannon alpha diversity, rarefied (IT vs ST) | `genus_abund.csv`; `metadata.csv` | `Figure_2C_genus_Shannon.pdf/.png`; `Figure_2C_shannon_genus.csv` | Figure_2C | Fig. 2C |
| `Fig.2/2D/2D.R` | Levins niche-breadth density + ECDF, plus the specialist-vs-generalist boxplot shown as the embedded inset within Fig. 2D; niche classification used downstream | genus.csv; genus_abund.csv; metadata.csv | Figure_2D_niche_breadth.pdf/.png; Figure_2E_specialist_generalist.pdf/.png; Figure_2D_niche_breadth_values.csv; threshold- and prevalence-sensitivity CSVs | Figure_2D, Figure_2E | Fig. 2D (the Figure_2E_specialist_generalist output is the embedded inset of Fig. 2D; Fig. 2E itself is produced by 2E&F&G.R) |
| `Fig.2/2E&F&G/2E&F&G.R` | Global FRI (2H), stress-pathway FRI bar chart (2I), and FRI ~ generalist-richness regression (2J); outputs FRI source data for 3F&G | `genus.csv`; `extracted_taxon_kegg_bacteria.csv`; `KEGG_paths.csv`; `metadata.csv`; `Figure_2D_niche_breadth_values.csv` | `Figure_2H_global_FRI.pdf`; `Figure_2I_stress_FRI.pdf`; `Figure_2J_regression.pdf`; `SourceData_2HJ_per_sample.csv`; `TableS5_FRI_stress_stats.csv`; `TableSx_all_pathways_IT_vs_ST.csv` | Figure_2H, 2I, 2J | Fig. 2E, 2F, 2G (TO CONFIRM) |
| `Fig.3/3A&B/3A&B.R` | CPCoA (db-RDA) of experimental heat-stress communities at genus (3A) and KO (3B) levels | `genus_abund.csv`; `ko_abund.csv`; `metadata.csv` | `Figure_3A_3B_CPCoA.pdf`; `SourceData_3AB_CPCoA_scores.csv`; `SourceData_3AB_global_summary.csv`; pairwise PERMANOVA CSVs | Figure_3A, Figure_3B | Fig. 3A, 3B |
| `Fig.3/3C&D/3C&D.R` | Microbiome compositional stability (mean BC distance, control vs heat) at genus (3C) and KO (3D) levels | `genus_abund.csv`; `ko_abund.csv`; `metadata.csv` | `Figure_3C_3D_stability.pdf`; `SourceData_3C_stability_genus.csv`; `SourceData_3D_stability_KO.csv`; `Stats_3C_3D_stability.csv` | Figure_3C, Figure_3D | Fig. 3C, 3D |
| `Fig.3/3E/3E.R` | Functional retention rate (FRR) under heat per KEGG Level-2 pathway. The published Fig. 3E is the between-habitat FRR-difference bar plot (Figure_3E_SameSide_BarPlot.pdf); the dumbbell, per-pathway, and overall-retention plots are additional/exploratory outputs | KO abundance table; metadata.csv; KEGG pathway annotation | Figure_3E_SameSide_BarPlot.pdf (published Fig. 3E); Figure_3E_direction_dumbbell.pdf; Figure_3E_overall_retention.pdf; Figure_3E_perpathway_retention.pdf; FRR source-data CSVs | Figure_3E | Fig. 3E |
| `Fig.3/3F&G/3F&G.R` | Regressions of coral thermotolerance (ΔFv/Fm) on (3G) FRI and (3H) generalist richness | `Figure_1E_change_rates.csv`; `SourceData_2HJ_per_sample.csv`; `Figure_2D_niche_breadth_values.csv` | `Figure_3G_FRI_vs_retention.pdf`; `Figure_3H_generalist_vs_retention.pdf`; `SourceData_3G_3H_regression.csv` | Figure_3G, Figure_3H | Fig. 3F, 3G (TO CONFIRM) |
| `Supplementary Figures/Fig.S1/S1.R` | Rarefaction curves for field samples (genus level) | Genus count CSVs (one per sample or combined) | `Figure_S1_Rarefaction.png/.pdf`; coverage and rarefaction-points CSVs | Figure_S1 | Fig. S1 |
| `Supplementary Figures/Fig.S2/S2.R` | Field phylum-level relative-abundance stacked bar chart (top 10 + others) | `Phylum_abund.csv`; `metadata.csv` | `Figure_S2_Field_Phylum_Composition.png/.pdf`; phylum summary CSVs | Figure_S2 | Fig. S2 |
| `Supplementary Figures/Fig.S3/S3.R` | Genus co-occurrence networks (field IT vs ST) and network topology table | `IT.csv`; `ST.csv` (genus relative-abundance tables) | `IT_network.gml`; `ST_network.gml`; `network_topology_metrics.csv` | Figure_S3 | Fig. S3 |
| `Supplementary Figures/Fig.S4/S4.R` | KO-level functional diversity (Shannon + Levins B, CSS-normalised) in the experiment | `genus_abund.csv`; `ko_abund.csv`; `metadata.csv` | `Figure_2F_KO_diversity.pdf/.png`; `Figure_2F_KO_shannon.csv`; `Figure_2F_KO_levinsB_css.csv` | Figure_2F | Fig. S4 (TO CONFIRM) |
| `Supplementary Figures/Fig.S5/S5.R` | Differential-abundance volcano plot (edgeR TMM, IT vs ST, KO level) | KO/genus count matrix + metadata | `Figure_S4_Volcano.png/.pdf`; `Figure_S4_DA_results.csv`; `Figure_S4_summary.csv` | Figure_S4 | Fig. S5 (TO CONFIRM) |
| `Supplementary Figures/Fig.S6/S6.R` | KEGG pathway enrichment facet plot (clusterProfiler ORA on DA genes) | DA results from S5; KEGG pathway annotations | `Figure_2G_KEGG_enrichment_FacetStyle.pdf/.png`; `Figure_2G_KEGG_enrichment_full.csv` | Figure_2G | Fig. S6 (TO CONFIRM) |
| `Supplementary Figures/Fig.S7/S7.R` | Sensitivity of generalist classification across niche-breadth thresholds | `genus.csv`; `metadata.csv`; niche breadth data | `FigureS5_generalist_robustness.pdf`; `SourceData_S5_generalist_robustness.csv` | FigureS5 | Fig. S7 (TO CONFIRM) |
| `Supplementary Figures/Fig.S8/S8.R` | Experimental phylum-level relative-abundance stacked bar chart | `Phylum_abund.csv` (experimental); `metadata.csv` | `Figure_Experiment_Phylum_Composition.png/.pdf`; phylum summary CSVs | (unlabelled in code) | Fig. S8 (TO CONFIRM) |
| `Supplementary Figures/Fig.S9/S9.R` | Experimental Shannon alpha diversity at genus and KO levels (4 groups: IT-C/H, ST-C/H) | `genus_abund.csv`; `ko_abund.csv`; `metadata.csv` | `Figure_S_experimental_Shannon.pdf`; `Stats_KW_experimental_Shannon.csv` | Figure_S | Fig. S9 (TO CONFIRM) |
| `Supplementary Figures/Fig.S10/S10.R` | KO co-occurrence networks for heat-treated groups (IT-H vs ST-H) with KEGGREST annotation | `ko_abund.csv`; KEGGREST API (internet required) | `Figure_3F_network_{IT-H,ST-H}_preview.pdf`; `Figure_3F_hubs_{IT-H,ST-H}.csv`; `Figure_3F_network_topology.csv` | Figure_3F | Fig. S10 (TO CONFIRM) |
| `Bacterial metagenomic data/rename.R` | Utility: renames/reformats raw metagenomic output files | Excel files (readxl) | (renamed files in-place) | — | — |

---

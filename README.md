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

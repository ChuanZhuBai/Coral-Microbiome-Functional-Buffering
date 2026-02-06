# Functional buffering by redundant bacterial generalists underpins coral holobiont thermotolerance

## Overview
This repository contains the custom R scripts and processed data necessary to reproduce the analysis presented in the paper:  
**"Functional buffering by redundant bacterial generalists underpins coral holobiont thermotolerance in fluctuating environments"**

## Directory Structure
- `Data/`: Processed taxonomic and functional abundance matrices (due to size limits, raw FASTQ/BAM files are deposited in NCBI SRA PRJNAXXXXXX).
- `Scripts/`: R scripts for data normalization, ecological indices calculation (FRI, Niche breadth), and figure generation.
- `Results/`: Output figures and statistical reports.

## Prerequisites
The scripts were developed and tested in R version 4.2.0+. Required packages:
- `tidyverse`, `vegan`, `ggpubr`, `metagenomeSeq`, `patchwork`, `rstatix`, `multcompView`

## Workflow
1. Run `Scripts/01_Data_Processing.R` to perform CSS normalization and ID mapping.
2. Run `Scripts/02_Field_Analysis.R` to generate baseline patterns (Figure 2).
3. Run `Scripts/03_Experiment_Analysis.R` to analyze heat stress responses (Figure 3).

## Contact
Your Name (your_email@institution.edu)

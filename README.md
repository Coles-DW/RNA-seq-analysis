RNA-seq Analysis

This repository contains R scripts and workflows for the analysis of RNA-seq datasets, including differential gene expression, visualization of up- and down-regulated genes, functional enrichment analysis, and network-level correlation analysis.

Repository Structure
RNA-seq-analysis/
├── scripts/
│   ├── DEG_barplot_timecourse.R   # Plot bar charts of differentially expressed genes over time/conditions  
│   ├── DEseq2.R                   # Differential expression analysis using DESeq2  
│   ├── Filtering_counts.R         # Preprocessing and filtering of raw count data  
│   ├── GO_enrichment_barplot.R    # Visualization of GO enrichment results (e.g. g:Profiler, topGO)  
│   ├── Merge.R                    # Combine and organize outputs from different steps  
│   ├── PLS-DA.R                   # Partial Least Squares Discriminant Analysis for expression profiles  
│   ├── TopGO.R                    # Functional enrichment using the topGO package  
│   └── utils/                     # (Optional) helper functions  
├── data/                          # Input data (not included for confidentiality)  
├── results/                       # Example outputs (plots, tables)  
└── README.md                      # Documentation  

Requirements

The following R packages are used throughout the scripts:

tidyverse

ggplot2

DESeq2

topGO

car

emmeans

multcompView

multcomp

visreg

ggeffects

gridExtra

viridis

Install missing packages with:

install.packages(c("tidyverse", "ggplot2", "car", "emmeans", 
                   "multcompView", "multcomp", "visreg", 
                   "ggeffects", "gridExtra", "viridis"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("DESeq2", "topGO"))

Example Workflows

Data Preprocessing

Use Filtering_counts.R to filter raw RNA-seq counts for quality control.

Differential Expression

Run DEseq2.R for differential gene expression analysis.

Use DEG_barplot_timecourse.R to visualize DEGs across time or treatment.

Functional Enrichment

Perform GO analysis with TopGO.R.

Plot results with GO_enrichment_barplot.R.

Advanced Analyses

Run PLS-DA.R for multivariate separation of experimental groups.

Combine multiple outputs with Merge.R.

Usage

Clone the repository and run scripts in R:

git clone https://github.com/Coles-DW/RNA-seq-analysis.git
cd RNA-seq-analysis/scripts
Rscript DEseq2.R


Or open the scripts in RStudio for interactive analysis.

Notes

Input data is not included for confidentiality.

Scripts are written to be modular and customizable.

Example figures and outputs can be stored in the results/ directory.

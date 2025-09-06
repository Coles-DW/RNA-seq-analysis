RNA-seq Analysis

This repository contains R scripts and workflows for the analysis of RNA-seq datasets, including differential gene expression, visualization of up- and down-regulated genes, functional enrichment analysis, and network-level correlation analysis.

Repository Structure
RNA-seq-analysis/
├── scripts/
│   ├── DEG_Analysis.R                 # Differential gene expression analysis  
│   ├── DEG_Proportions_Plot.R         # Bar plot of up- and down-regulated DEGs  
│   ├── GO_Enrichment_Visualization.R  # Visualization of GO enrichment (g:Profiler results)  
│   ├── Spearman_Correlation_Network.R # Correlation analysis for network visualization  
│   └── utils/                         # Helper functions (if needed)  
├── data/                              # Input data (not included in repo for confidentiality)  
├── results/                           # Example output figures and tables  
└── README.md                          # Documentation  

Requirements

The scripts rely on the following R packages:

tidyverse

ggplot2

car

emmeans

multcompView

multcomp

visreg

ggeffects

lme4

gridExtra

viridis

Install missing packages with:

install.packages(c("tidyverse", "ggplot2", "car", "emmeans", 
                   "multcompView", "multcomp", "visreg", 
                   "ggeffects", "lme4", "gridExtra", "viridis"))

Example Workflows

Differential Expression Analysis

Run statistical models (linear or mixed effects).

Perform pairwise contrasts.

Generate ANOVA tables and effect plots.

DEG Proportion Plotting

Visualize numbers of up- and down-regulated genes across timepoints or conditions.

GO Term Enrichment

Import results from g:Profiler or other tools.

Generate dot plots of enriched terms by category.

Correlation Network Analysis

Calculate Spearman correlations across samples.

Export network edge lists for Cytoscape or other visualization platforms.

Usage

Clone the repository and run scripts in R:

git clone https://github.com/Coles-DW/RNA-seq-analysis.git
cd RNA-seq-analysis/scripts
Rscript DEG_Analysis.R


Alternatively, open scripts in RStudio for interactive exploration.

Notes

Input data files are not provided for confidentiality.

Scripts are written to be modular and can be adapted to other datasets.

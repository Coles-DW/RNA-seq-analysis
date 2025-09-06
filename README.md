RNA-seq Analysis Pipeline

A modular RNA-seq analysis workflow for preprocessing, differential expression, and functional annotation.

Features

Filter and preprocess RNA-seq count data

Differential expression analysis using DESeq2

Gene Ontology (GO) enrichment and visualization

Partial Least Squares Discriminant Analysis (PLS-DA)

Reproducible and adaptable R scripts

Repository Structure
RNA-seq-analysis/
├── scripts/          # R scripts for analysis
├── data/             # Input counts and annotation files
├── plots/            # Generated figures
└── README.md

Quick Start

Clone the repo

git clone https://github.com/Coles-DW/RNA-seq-analysis.git
cd RNA-seq-analysis


Install dependencies

install.packages(c("ggplot2", "dplyr", "tidyr"))
BiocManager::install(c("DESeq2", "topGO", "gProfileR"))


Run scripts

source("scripts/Filtering_counts.R")
source("scripts/Merge.R")
source("scripts/DEseq2.R")
source("scripts/TopGO.R")
source("scripts/goanalysis_gprofilier_plot.R")
source("scripts/PLS-DA.R")


View results in plots/.

Notes

Scripts assume correctly formatted input files in data/

Modular workflow: adapt scripts to your dataset as needed

License

MIT License © Donovin Coles

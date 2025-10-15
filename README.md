# scRNA-seq Roots - Single-cell RNA-seq analysis for root tissues

This repository contains a reproducible R workflow to download and analyze public single-cell RNA-seq data (example: Arabidopsis / maize root datasets).
Files:
- data_download.R : downloads required public data (GEO/10x)
- analysis.R : full Seurat pipeline (QC, normalization, PCA, UMAP, clustering, markers)
- renv_install.R : initialize renv environment
- README.md : this file

Run:
1. Open R in project folder
2. source("renv_install.R")  # optional
3. source("data_download.R")
4. source("analysis.R")


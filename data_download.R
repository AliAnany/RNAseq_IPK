# data_download.R
# Downloads example scRNA-seq data for roots.
# This script attempts to download from GEO (GSE accession) or 10x samples.
# Edit paths/URLs as needed.

if (!requireNamespace("GEOquery", quietly=TRUE)) install.packages("GEOquery", repos="https://cloud.r-project.org")
library(GEOquery)

outdir <- "data"
dir.create(outdir, showWarnings = FALSE)

# Example GEO accession (Arabidopsis root single-cell): E-MTAB-9355 is on ArrayExpress,
# we'll try GEO GSE if available; user can change accession below.
# Attempt to download GSE155928 (maize scRNA-seq) as an example:
accession <- "GSE155928"

cat("Attempting to download", accession, "via GEOquery ...\n")
try({
  gset <- getGEO(accession, GSEMatrix = FALSE, destdir = outdir)
  saveRDS(gset, file = file.path(outdir, paste0(accession, "_raw.rds")))
  cat("Saved raw GEO object to", outdir, "\n")
}, silent = TRUE)

cat("If GEO download fails, please download 10x-formatted matrix and place in data/\n")

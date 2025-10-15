# analysis.R - scRNA-seq analysis using Seurat
# Full reproducible script. Assumes data/ contains either a 10X directory or a saved GEO object.

# Install and load packages
packages <- c("Seurat","patchwork","dplyr","ggplot2","Matrix")
installed <- rownames(installed.packages())
for (p in packages) if (!p %in% installed) install.packages(p, repos="https://cloud.r-project.org")

library(Seurat); library(patchwork); library(dplyr); library(ggplot2)

data_dir <- "data"
# Try to load 10X
tenx_dirs <- list.files(data_dir, pattern = "filtered_feature_bc_matrix|10x|matrix.mtx", recursive = TRUE, full.names = TRUE)
seu <- NULL

if (length(tenx_dirs) > 0) {
  cat("Found possible 10x files. Attempting to read 10X data...\n")
  # If there's a directory named features/ or matrix.mtx
  try({
    seu <- Read10X(data.dir = data_dir)
    seu <- CreateSeuratObject(counts = seu, project = "Root_scRNA")
  }, silent = TRUE)
}

# If not, try GEO raw rds
if (is.null(seu)) {
  gse_files <- list.files(data_dir, pattern = "GSE|gse|_raw.rds", full.names = TRUE)
  if (length(gse_files) > 0) {
    cat("Loading GEO rds file:", gse_files[1], "\n")
    gset <- readRDS(gse_files[1])
    # Attempt extraction - depends on the GEO object structure
    # Here we assume expression matrix available as a table in gset
    # User may need to adapt this block for the specific GEO object
    if (length(gset) >= 1) {
      eset <- try(exprs(gset[[1]]), silent = TRUE)
      if (!inherits(eset, "try-error")) {
        seu <- CreateSeuratObject(counts = as.matrix(eset), project = "Root_scRNA")
      }
    }
  }
}

if (is.null(seu)) stop("No input data found in data/. Place a 10X folder or GEO rds file there.")

# Basic QC
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-|^mt-")
VlnPlot(seu, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3)

# Filter
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 20)

# Normalization and feature selection
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
ElbowPlot(seu)

# UMAP + clustering
seu <- RunUMAP(seu, dims = 1:20)
seu <- FindNeighbors(seu, dims = 1:20)
seu <- FindClusters(seu, resolution = 0.5)
DimPlot(seu, reduction = "umap", label = TRUE) + ggtitle("scRNA-seq UMAP")

# Marker identification
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, file = "output/top10_markers_per_cluster.csv", row.names = FALSE)

# Save Seurat object
saveRDS(seu, file = "output/seurat_root_scRNA.rds")
cat("Analysis complete. Outputs in output/\n")

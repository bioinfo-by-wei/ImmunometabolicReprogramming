
rm(list=ls())
# Load necessary libraries
library(Seurat)
library(qs)

# Read in the 10X Genomics data
data <- Read10X("datadir")

# Create a Seurat object with metadata
dt <- CreateSeuratObject(counts = data, meta.data = meta, project = 'data')

# Calculate mitochondrial gene percentage
dt <- PercentageFeatureSet(dt, pattern = "^MT-", col.name = 'percent.mt')

# Filter cells based on quality control metrics
dt <- subset(dt, subset = percent.mt < 20 & nFeature_RNA > 200 & nFeature_RNA < 6000)

# Normalize the data
dt <- NormalizeData(dt)

# Scale the data
dt <- ScaleData(dt, features = rownames(dt))

# Identify variable features
dt <- FindVariableFeatures(dt)

# Perform PCA
dt <- RunPCA(dt)

# Plot elbow plot to determine the number of principal components
ElbowPlot(dt, ndims = 50)

# UMAP dimensionality reduction and clustering
dims <- 1:30  # Choose appropriate dimensions based on ElbowPlot
dt <- RunUMAP(dt, dims = dims)
dt <- FindNeighbors(dt, dims = dims)
dt <- FindClusters(dt, resolution = 0.1)


# Define major cell types for clusters
major_type <- data.frame(ClusterID = 0:32, major_type = NA)

# Assign major types to specific clusters
major_type[major_type$ClusterID %in% c(2, 8, 12, 3, 18), "major_type"] <- 'B/plasma'
major_type[major_type$ClusterID %in% c(0, 1, 5, 7, 30, 28), "major_type"] <- 'T/NK'
major_type[major_type$ClusterID %in% c(4, 10, 20, 27), "major_type"] <- 'myeloids'
major_type[major_type$ClusterID %in% c(25), "major_type"] <- 'mast'
major_type[major_type$ClusterID %in% c(11, 31, 32), "major_type"] <- 'endothelial'
major_type[major_type$ClusterID %in% c(6, 14, 16, 26, 29, 13), "major_type"] <- 'fibroblast'
major_type[major_type$ClusterID %in% c(23), "major_type"] <- 'entericglial'
major_type[major_type$ClusterID %in% c(9, 19, 21, 22, 15, 17, 24), "major_type"] <- 'epithelial'

# Add major type annotation to Seurat object's metadata
dt@meta.data$major_type <- "NA"
for (i in seq_len(nrow(major_type))) {
  dt@meta.data[which(dt@meta.data$RNA_snn_res.0.1 == major_type$ClusterID[i]), 'major_type'] <- major_type$major_type[i]
}

# Save the updated Seurat object
qsave(dt, "crc_major.qs", nthreads = 10)

# Subset Seurat object by major types
subsets <- list(
  epi = dt[, dt@meta.data$major_type %in% c("epithelial")],
  mye = dt[, dt@meta.data$major_type %in% c("myeloids")],
  tnk = dt[, dt@meta.data$major_type %in% c("T/NK")],
  endo = dt[, dt@meta.data$major_type %in% c("endothelial")],
  fib = dt[, dt@meta.data$major_type %in% c("fibroblast")],
  mast = dt[, dt@meta.data$major_type %in% c("mast")],
  entericglial = dt[, dt@meta.data$major_type %in% c("entericglial")],
  b = dt[, dt@meta.data$major_type %in% c("B/plasma")]
)

# Function to process each subset
process_subset <- function(subset, name) {
  subset <- qread(paste0(name, ".qs"), nthreads = 10)
  
  # Re-normalize, find variable features, scale, and run PCA
  subset <- NormalizeData(subset)
  subset <- FindVariableFeatures(subset, selection.method = "vst", nfeatures = 2000)
  subset <- ScaleData(subset, features = rownames(subset))
  subset <- RunPCA(subset, features = VariableFeatures(object = subset))
  ElbowPlot(subset)
  
  # Batch correction using Harmony if needed
  if (!is.null(subset@meta.data$orig.ident)) {
    subset <- RunHarmony(subset, reduction = "pca", group.by.vars = "orig.ident", reduction.save = "harmony")
  }
  
  # Neighbors, clusters, and UMAP
  dims <- 1:20  # Adjust based on ElbowPlot results
  subset <- FindNeighbors(subset, dims = dims, reduction = ifelse(!is.null(subset@reductions$harmony), "harmony", "pca"))
  subset <- FindClusters(subset, resolution = 0.15)
  subset <- RunUMAP(subset, dims = dims, reduction = ifelse(!is.null(subset@reductions$harmony), "harmony", "pca"), reduction.name = "umap")
  
  # Dimension plots
  DimPlot(subset, reduction = "umap", label = TRUE)
  DimPlot(subset, reduction = 'umap', group.by = "sample.ID", label = TRUE, pt.size = 0.5) + NoLegend()
  
  # Find markers
  markers <- FindAllMarkers(subset, logfc.threshold = 0.25, only.pos = TRUE)
  write.csv(markers, file = paste0(name, "_markers_0.15.csv"))
  
  # Save processed subset
  qsave(subset, paste0(name, "_cluster.qs"), nthreads = 10)
}

# Apply the function to each subset
lapply(names(subsets), function(name) process_subset(subsets[[name]], name))

# For mast and entericglial, rename subtypes after clustering
for (subset_name in c("mast", "entericglial")) {
  subset <- qread(paste0(subset_name, "_cluster.qs"), nthreads = 10)
  
  # Define subtypes
  sub_type <- data.frame(ClusterID = unique(subset@meta.data$RNA_snn_res.0.1), sub_type = NA)
  if (subset_name == "mast") {
    sub_type[sub_type$ClusterID %in% c(0, 1), "sub_type"] <- 'mast'
    sub_type[sub_type$ClusterID %in% c(2, 3, 4, 5), "sub_type"] <- 'del'
  } else if (subset_name == "entericglial") {
    sub_type[sub_type$ClusterID %in% c(0, 1, 2, 3), "sub_type"] <- 'entericglial'
    sub_type[sub_type$ClusterID %in% c(4, 5), "sub_type"] <- 'del'
  }
  
  # Update metadata
  subset@meta.data$sub_type <- "NA"
  for (i in seq_len(nrow(sub_type))) {
    subset@meta.data[which(subset@meta.data$RNA_snn_res.0.1 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
  }
  
  # Save updated subset
  qsave(subset, paste0(subset_name, "_subtype.qs"), nthreads = 10)
}


# Define subtypes for myeloid cells based on cluster IDs
sub_type <- data.frame(ClusterID = 0:10, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0), "sub_type"] <- 'APOE+SPP1+TAM'
sub_type[sub_type$ClusterID %in% c(1), "sub_type"] <- 'Monocyte'
sub_type[sub_type$ClusterID %in% c(2), "sub_type"] <- 'Neutrophils'
sub_type[sub_type$ClusterID %in% c(3), "sub_type"] <- 'FOLR2+TAM'
sub_type[sub_type$ClusterID %in% c(4), "sub_type"] <- 'SELENOP+TAM'
sub_type[sub_type$ClusterID %in% c(6), "sub_type"] <- 'cDC1'
sub_type[sub_type$ClusterID %in% c(8), "sub_type"] <- 'mature cDC'
sub_type[sub_type$ClusterID %in% c(9), "sub_type"] <- 'pDC1'
sub_type[sub_type$ClusterID %in% c(5,7,10), "sub_type"] <- 'del'

# Update metadata with subtype annotations
mye@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  mye@meta.data[which(mye@meta.data$RNA_snn_res.0.15 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(mye, "mye_subtype.qs", nthreads = 10)



# Define subtypes for B cells/plasma cells based on cluster IDs
sub_type <- data.frame(ClusterID = 0:9, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0,2,4,6,7,8,9), "sub_type"] <- 'Plasma'
sub_type[sub_type$ClusterID %in% c(1,3), "sub_type"] <- 'B'
sub_type[sub_type$ClusterID %in% c(5), "sub_type"] <- 'del'

# Update metadata with subtype annotations
b@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  b@meta.data[which(b@meta.data$RNA_snn_res.0.1 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(b, "b_subtype.qs", nthreads = 10)


# Define subtypes for T/NK cells based on cluster IDs
sub_type <- data.frame(ClusterID = 0:8, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0), "sub_type"] <- 'Cytotoxic CD8+T'
sub_type[sub_type$ClusterID %in% c(1), "sub_type"] <- 'Central memory CD4+T cell'
sub_type[sub_type$ClusterID %in% c(2), "sub_type"] <- 'Naïve CD4+T'
sub_type[sub_type$ClusterID %in% c(3), "sub_type"] <- 'Treg'
sub_type[sub_type$ClusterID %in% c(4), "sub_type"] <- 'Effector memory CD8+T'
sub_type[sub_type$ClusterID %in% c(5), "sub_type"] <- 'NK'
sub_type[sub_type$ClusterID %in% c(6), "sub_type"] <- 'Proliferating'
sub_type[sub_type$ClusterID %in% c(7,8), "sub_type"] <- 'del'

# Update metadata with subtype annotations
tnk@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  tnk@meta.data[which(tnk@meta.data$RNA_snn_res.0.2 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(tnk, "tnk_subtype.qs", nthreads = 10)


# Define subtypes for fibroblasts based on cluster IDs
sub_type <- data.frame(ClusterID = 0:10, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0), "sub_type"] <- 'INHBA+Fib'
sub_type[sub_type$ClusterID %in% c(1), "sub_type"] <- 'Pericytes'
sub_type[sub_type$ClusterID %in% c(2), "sub_type"] <- 'C7+Fib'
sub_type[sub_type$ClusterID %in% c(3), "sub_type"] <- 'ABCA8+Fib'
sub_type[sub_type$ClusterID %in% c(4), "sub_type"] <- 'ATP5F1D+Fib'
sub_type[sub_type$ClusterID %in% c(5), "sub_type"] <- 'ADH1B+Fib'
sub_type[sub_type$ClusterID %in% c(6), "sub_type"] <- 'VSTM2A+Fib'
sub_type[sub_type$ClusterID %in% c(7), "sub_type"] <- 'Myofibroblasts'
sub_type[sub_type$ClusterID %in% c(8), "sub_type"] <- 'del'
sub_type[sub_type$ClusterID %in% c(9), "sub_type"] <- 'Smooth muscle'
sub_type[sub_type$ClusterID %in% c(10), "sub_type"] <- 'MEST+Fib'

# Update metadata with subtype annotations
fib@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  fib@meta.data[which(fib@meta.data$RNA_snn_res.0.1 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(fib, "fib_subtype.qs", nthreads = 10)


# Define subtypes for endothelial cells based on cluster IDs
sub_type <- data.frame(ClusterID = 0:6, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0), "sub_type"] <- 'Venous'
sub_type[sub_type$ClusterID %in% c(1), "sub_type"] <- 'Angiogenic'
sub_type[sub_type$ClusterID %in% c(2), "sub_type"] <- 'Capillary'
sub_type[sub_type$ClusterID %in% c(3), "sub_type"] <- 'Lymphatic'
sub_type[sub_type$ClusterID %in% c(4), "sub_type"] <- 'Pericytes like'
sub_type[sub_type$ClusterID %in% c(5), "sub_type"] <- 'del'
sub_type[sub_type$ClusterID %in% c(6), "sub_type"] <- 'Proliferating'

# Update metadata with subtype annotations
endo@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  endo@meta.data[which(endo@meta.data$RNA_snn_res.0.1 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(endo, "endo_subtype.qs", nthreads = 10)


# Define subtypes for epithelial cells based on cluster IDs
sub_type <- data.frame(ClusterID = 0:9, sub_type = NA)
sub_type[sub_type$ClusterID %in% c(0), "sub_type"] <- 'KRT23+epi'
sub_type[sub_type$ClusterID %in% c(1), "sub_type"] <- 'REG4+epi'
sub_type[sub_type$ClusterID %in% c(2), "sub_type"] <- 'Colonocyte'
sub_type[sub_type$ClusterID %in% c(3), "sub_type"] <- 'Goblet cell'
sub_type[sub_type$ClusterID %in% c(4), "sub_type"] <- 'Paneth'
sub_type[sub_type$ClusterID %in% c(5), "sub_type"] <- 'SPON2+epi'
sub_type[sub_type$ClusterID %in% c(6), "sub_type"] <- 'SPINK4+epi'
sub_type[sub_type$ClusterID %in% c(7), "sub_type"] <- 'del'
sub_type[sub_type$ClusterID %in% c(8), "sub_type"] <- 'PGF+epi'
sub_type[sub_type$ClusterID %in% c(9), "sub_type"] <- 'Enteroendocrine'

# Update metadata with subtype annotations
epi@meta.data$sub_type <- "NA"
for (i in seq_len(nrow(sub_type))) {
  epi@meta.data[which(epi@meta.data$RNA_snn_res.0.15 == sub_type$ClusterID[i]), 'sub_type'] <- sub_type$sub_type[i]
}

# Save the updated Seurat object
qsave(epi, "epi_subtype.qs", nthreads = 10)

# Load the subtype-annotated Seurat objects
epi <- qread("epi_subtype.qs", nthreads = 10)
endo <- qread("endo_subtype.qs", nthreads = 10)
fib <- qread("fib_subtype.qs", nthreads = 10)
b <- qread("b_subtype.qs", nthreads = 10)
tnk <- qread("tnk_subtype.qs", nthreads = 10)
mye <- qread("mye_subtype.qs", nthreads = 10)
ent <- qread("ent_subtype.qs", nthreads = 10)
mast <- qread("mast_subtype.qs", nthreads = 10)
crc <- qread("crc_tnk_epi_end_fib.qs", nthreads = 10)

# Merge the datasets into one Seurat object
merged <- Reduce(function(x, y) merge(x, y = y, add.cell.ids = TRUE, project = "CRC"), list(epi, endo, fib, b, tnk, mye, ent, mast))

# Ensure unique cell barcodes
merged <- merged[, !duplicated(colnames(merged))]

# Join layers to consolidate information
crc_1 <- JoinLayers(merged)

# Save the final combined Seurat object
qsave(crc_1, "crc_subtype_all.qs", nthreads = 48)





































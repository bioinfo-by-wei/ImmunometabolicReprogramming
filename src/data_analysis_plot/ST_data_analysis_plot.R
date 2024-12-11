# Clear workspace and load necessary libraries
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(jsonlite)
library(png)
library(GSVA)
library(qs)
library(ggpubr)
library(rstatix)
library(future)

### Spatial Transcriptomics Data Analysis

#### Sample 2
# Load spatial transcriptomics data
st1 <- readRDS("/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015.rds")
spe1 <- st1

# Normalize and transform the data using SCTransform
spe1 <- SCTransform(spe1, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# Plot gene expression
p3 <- SpatialFeaturePlot(spe1, features = c("APOE"), pt.size = 1.6, alpha = c(0.1, 1))

# Dimensionality reduction, clustering, and visualization
spe1 <- RunPCA(spe1, assay = "SCT", verbose = FALSE)
spe1 <- FindNeighbors(spe1, reduction = "pca", dims = 1:30)
spe1 <- FindClusters(spe1, verbose = FALSE)
spe1 <- RunUMAP(spe1, reduction = "pca", dims = 1:30)

# UMAP plot with cluster labels
p4 <- DimPlot(spe1, reduction = "umap", label = TRUE)

# Subset interesting clusters
spe2 <- subset(spe1, idents = c(1, 2, 3, 5, 6))

# Spatial dimension plot with cluster labels
p5 <- SpatialDimPlot(spe1, label = TRUE, label.size = 3)

# Highlight specific clusters
p6 <- SpatialDimPlot(spe1, cells.highlight = CellsByIdentities(object = spe1, idents = c(2, 0, 3)), 
                     facet.highlight = TRUE, ncol = 3, pt.size.factor = 3)

# Differential gene expression analysis
de_markers <- FindMarkers(spe1, ident.1 = 0, ident.2 = 3)
p7 <- SpatialFeaturePlot(object = spe1, features = rownames(de_markers)[1:3], 
                         alpha = c(0.1, 1), pt.size.factor = 3, ncol = 3)

# Perform differential expression for all clusters
de_markers2 <- FindAllMarkers(spe1)
write.csv(de_markers2, "/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015_marker.csv")

# Save processed object
saveRDS(spe1, "/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015_run.rds")

#### Visualization
# Reload processed object
spe1 <- readRDS("/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015_run.rds")

# Spatial dimension plot with cluster labels
p1 <- SpatialDimPlot(spe1, label = TRUE, label.size = 4)
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/ST_plot/ST_000015_st_cluster.pdf", 
       plot = p1, device = "pdf", width = 4, height = 4, units = "in")

# Gene expression plot
p2 <- SpatialFeaturePlot(spe1, features = c("CYP27A1"), pt.size = 1.6, alpha = c(0.1, 1))
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/ST_plot/ST_000015_st_CYP27A1.pdf", 
       plot = p2, device = "pdf", width = 4, height = 4, units = "in")

#### Gene Set Scoring
# Load gene sets and clean up empty entries
geneset <- read.csv("/public2022/weijinfen/flux/tam_new_analysis/st_data/st_geneset.csv")
cleaned_geneset <- lapply(geneset, function(x) x[x != ""])

# Ensure genes list is not empty before adding module score
genes_list <- as.character(cleaned_geneset[["Treg_1"]])
if (length(genes_list) > 0) {
  spe1 <- AddModuleScore(object = spe1, features = list(genes_list), name = "Treg_1")
  score_column <- "Treg_1_Score"
  spe1@meta.data[[score_column]] <- spe1[[score_column]]
  print(colnames(spe1@meta.data))
} else {
  warning("No genes found for gene set: Treg_1")
}

# Save scored object
qsave(spe1, "/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015_run_score.qs")

#### Plotting Gene Set Scores
score_column <- "Phagocytosis1"
p3 <- SpatialFeaturePlot(spe1, features = score_column, alpha = c(0.1, 1), min.cutoff = -0.3, max.cutoff = 0.8) +
  theme_bw() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank()) +
  ggtitle(" ")
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/ST_plot/ST_000015_st_Phagocytosis1.pdf", 
       plot = p3, device = "pdf", width = 3, height = 3, units = "in")

#### Correlation Scatter Plot
# Load scored object
spe1 <- qread("/public2022/weijinfen/flux/tam_new_analysis/st_data/SCAR_ST_000015_run_score.qs")

# Extract relevant columns for correlation analysis
df <- as.matrix(spe1@meta.data[, 8:16])
df1 <- as.data.frame(df)

# Define ranges for axes
xlim_range <- range(df1$Treg_11, na.rm = TRUE)
ylim_range <- range(df1$Lipid_loaded_TAM1, na.rm = TRUE)

# Create scatter plot with regression line and correlation coefficient
p4 <- ggscatter(df1, x = "Treg_11", y = "Lipid_loaded_TAM1",
                size = 1, color = "gray", 
                xlim = xlim_range, ylim = ylim_range,
                add = "reg.line", add.params = list(color = "#FC4E07", fill = "#FC4E07", size = 1),
                conf.int = TRUE) +
  stat_cor(method = "spearman", label.x = min(df1$Treg_11, na.rm = TRUE), label.y = max(df1$Lipid_loaded_TAM1, na.rm = TRUE), size = 6) +
  xlab("Treg score") + ylab("Lipid loaded TAM score") +
  theme(axis.title = element_text(size = 16), text = element_text(size = 16))

# Save plot
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/ST_plot/cor_scatter/score_cor_8.pdf", 
       plot = p4, device = "pdf", width = 4, height = 4, units = "in")

#### Violin Plots for Score Comparisons Across Clusters
# Extract metadata and select columns of interest
meta_data <- spe1@meta.data
selected_columns <- c("Anti_inflammatory1", "Lipid_loaded_TAM1", "Treg1", "CD8Tex1", "CD8Teff1", "Immunosuppression1", "Phagocytosis1", "Proliferation1", "SCT_snn_res.0.2")
df <- meta_data[, selected_columns]

# Convert data to long format
df_long <- df %>%
  pivot_longer(cols = starts_with(c("Anti_inflammatory1", "Lipid_loaded_TAM1", "Treg1", "CD8Tex1", "CD8Teff1", "Immunosuppression1", "Phagocytosis1", "Proliferation1")),
               names_to = "Score_Type", values_to = "Score")

# Define custom colors
custom_colors <- c("#f6bcfd", "#ff7d00", "#ffc512", "#ffa300", "#8dd3c6", "#ff6581", "#f8d90d", "#a5da6b", "#e578d6", "#ffd2d8", "#90e4cd", "#f06152")[1:5]

# Loop to create violin plots for each score type
for (score_type in unique(df_long$Score_Type)) {
  p <- ggplot(df_long[df_long$Score_Type == score_type, ], aes(x = factor(SCT_snn_res.0.2), y = Score, fill = factor(SCT_snn_res.0.2))) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
    labs(title = paste(score_type, "Scores in Cluster"),
         x = "Cluster", y = "Score") +
    theme_minimal() +
    theme(legend.position = "right", axis.title = element_text(size = 14), axis.text = element_text(size = 12),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(color = "black", size = 0.5)) +
    scale_fill_manual(values = custom_colors)

  # Save plot
  ggsave(filename = paste0("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/ST_plot/", score_type, "_violin_plot.pdf"), 
         plot = p, width = 10, height = 5)
}


# Clear workspace and load necessary libraries
rm(list = ls())
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(qs)

### Load New Spatial Transcriptomics Data (GSE226997)

#### Conda Environment: singlecell
# Read 10x count data
ct <- Read10X("/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/Ajou_Visium_P1")

# Read image data
img <- Read10X_Image(image.dir = "/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/Ajou_Visium_P1/Spatial", filter.matrix = TRUE)

# Create Seurat object
sc1 <- CreateSeuratObject(counts = ct, assay = "Spatial")
img <- img[Cells(x = sc1)]
DefaultAssay(sc1) <- "Spatial"
sc1[['slice1']] <- img

# Normalize and transform the data using SCTransform
sc1 <- SCTransform(sc1, assay = "Spatial", verbose = FALSE)

# Save processed object
qsave(sc1, "/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/sc1_v1.qs")

#### Conda Environment: flux
# Reload processed object
sc1 <- qread("/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/sc1_v1.qs")

# Visualization of raw data
p1 <- VlnPlot(sc1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
p2 <- SpatialFeaturePlot(sc1, features = "nCount_Spatial")
p3 <- SpatialFeaturePlot(sc1, features = "APOE", pt.size.factor = 2)
p4 <- SpatialFeaturePlot(sc1, features = "APOE", pt.size.factor = 2, alpha = c(0.5, 1))

# Dimensionality reduction and clustering
sc1 <- RunPCA(sc1, assay = "SCT", verbose = FALSE)
sc1 <- FindNeighbors(sc1, reduction = "pca", dims = 1:30)
sc1 <- FindClusters(sc1, verbose = FALSE)
sc1 <- RunUMAP(sc1, reduction = "pca", dims = 1:30)

# UMAP plot with cluster labels
p5 <- DimPlot(sc1, reduction = "umap", label = TRUE)

# Spatial dimension plot with cluster labels
p6 <- SpatialDimPlot(sc1, label = TRUE, label.size = 3)

# Differential gene expression analysis
markers1 <- FindAllMarkers(sc1)
write.csv(markers1, "/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/p1_marker.csv")

# Save processed object
qsave(sc1, "/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/sc1_v2.qs")

#### Lollipop Plot for Gene Correlation Analysis
# Load correlation data
df <- qread("/public2022/weijinfen/flux/tam_new_analysis/st_data/GSE226997/st_p3_gene_cor.csv")
colnames(df) <- c("X", "Y", "Correlation", "P_Value")

# Calculate absolute correlation and categorize
df$Abs_Correlation <- abs(df$Correlation)
df$cor1 <- cut(df$Abs_Correlation,
               breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, Inf),
               labels = c("< 0.1", "0.1 - 0.2", "0.2 - 0.3", "0.3 - 0.4", "0.4 - 0.5", "> 0.5"),
               right = FALSE)

df$pvalue1 <- cut(df$P_Value,
                  breaks = c(0, 0.001, 0.01, 0.05, Inf),
                  labels = c("< 0.001", "< 0.01", "< 0.05", "> 0.05"),
                  right = FALSE)

# Order by correlation
df <- df[order(df$Correlation, decreasing = TRUE), ]
df$Y_Order <- factor(df$Y, levels = df$Y)

# Create lollipop plot
p <- ggplot(df, aes(y = Correlation, x = Y_Order, color = Y_Order)) +
  scale_color_manual(name = "Variable",
                     values = c(
                       "Immunosuppression" = "#003366",
                       "Phagocytosis" = "#A020F0",
                       "Treg" = "#FFC0CB",
                       "Anti_inflammatory" = "#FF69B4",
                       "CD3D" = "#8B008B",
                       "CD68" = "#FFD700",
                       "Proliferation" = "#F0E68C",
                       "CEACAM5" = "#00CED1",
                       "EPCAM" = "#ADD8E6"
                     )) +
  geom_segment(aes(x = Y_Order, xend = Y_Order, y = 0, yend = Correlation), size = 1) +
  geom_point(aes(size = cor1)) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 14)) +
  labs(y = "Correlation coefficient with SPI1 expression", x = "") +
  coord_flip()

# Add p-value labels
df$pvalue2 <- cut(df$P_Value,
                  breaks = c(0, 0.05, Inf),
                  labels = c("< 0.05", "> 0.05"),
                  right = FALSE)

p1 <- ggplot(df, aes(y = 0, x = Y_Order, label = ifelse(P_Value < 0.05, sprintf("%.3e", P_Value), ""))) +
  geom_text(aes(color = pvalue2), vjust = -0.2, size = 4) +
  scale_color_manual(name = "", values = c("< 0.05" = "red", "> 0.05" = "black")) +
  theme_void() +
  guides(color = F) +
  coord_flip()

# Combine plots
final_plot <- plot_grid(p, p1, nrow = 1)

# Save combined plot
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/SC_plot_new/cor/p3_cor_spi_genes.pdf", 
       plot = final_plot, width = 7, height = 8)


# Load necessary libraries
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(qs)

### Load and Visualize Spatial Transcriptomics Data

#### Load Spatial Transcriptomics Data
sc1 <- readRDS("/wangzhenzenlab/slb/stereoseq/STT0000036/CRCP56_T/CRCP56_T_bin1.rds")
spe1 <- sc1

#### Cluster Visualization
p1 <- SpatialDimPlot(spe1, label = TRUE, label.size = 3, pt.size.factor = 2) 
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/SC_plot_new/sc_cluster/p4_st_cluster.pdf", plot = p1, width = 4, height = 4, units = "in")

#### Multiple Gene Expression Visualization
genes <- c("SIRPA", "SOAT1", "LIPA", "PSAP", "CYP27A1", "CH25H", "ABCA1", "ABCG1", "SIGLEC10", "SCARB2", "MSR1", "CD36", "TREM2", "OLR1", "CD163", "APOE", "SPP1", "STAT5A", "SPI1", "CD68", "NPC1", "NPC2", "FOXP3", "CD3D", "CD8A", "CD8B", "EPCAM", "CEACAM5")

plots <- list()
for (gene in genes) {
  p <- SpatialFeaturePlot(spe1, features = gene, pt.size.factor = 2, alpha = c(0.5, 1))
  plots[[gene]] <- p
}

combined_plot <- wrap_plots(plots, ncol = 6)
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/SC_plot_new/sc_cluster/p4_multiple_genes_expression.pdf", plot = combined_plot, width = 10, height = 10, units = "in")

#### Key Gene Expression Visualization
p2 <- SpatialFeaturePlot(spe1, features = c("SPP1"), pt.size.factor = 2, alpha = c(0.5, 1), max.cutoff = 2)
ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/SC_plot_new/sc_cluster/p4_st_SPP1.pdf", plot = p2, width = 4, height = 4, units = "in")

### Calculate Gene Set Scores Across Multiple Samples

#### Define Sample File Paths and Read Geneset
file_paths <- c(
  "/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/qs_data/CRCP112_bin1_cluster.qs",
  "/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/qs_data/CRCP59_bin1_cluster.qs",
  "/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/qs_data/CRCP95_bin1_cluster.qs",
  "/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/qs_data/CRCP56_bin1_cluster.qs",
  "/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/qs_data/CRCP67_bin1_cluster.qs"
)

geneset <- read.csv("/public2022/weijinfen/flux/tam_new_analysis/st_data/st_geneset.csv")
cleaned_geneset <- lapply(geneset, function(x) x[x != ""])
gene_set_names <- names(cleaned_geneset)

#### Process Each Sample
for (file_path in file_paths) {
  sc1 <- qread(file_path)
  sample_name <- sub("^.*/(CRC[^_]+).*", "\\1", file_path)
  
  spe1 <- sc1
  
  for (gene_set_name in gene_set_names) {
    genes <- cleaned_geneset[[gene_set_name]]
    
    if (length(genes) > 0) {
      spe1 <- AddModuleScore(object = spe1, features = list(genes), name = gene_set_name)
      
      score_column <- gene_set_name
      
      if (!score_column %in% colnames(spe1@meta.data)) {
        score_columns <- grep(paste0("^", gene_set_name, "\\d+$"), colnames(spe1@meta.data), value = TRUE)
        
        if (length(score_columns) > 0) {
          score_column <- score_columns[1]
        } else {
          stop(paste0("Score column for ", gene_set_name, " not found in meta.data"))
        }
      }
      
      spe1@meta.data[[gene_set_name]] <- spe1@meta.data[[score_column]]
    }
  }
  
  save_path <- file.path("/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds", paste0(sample_name, "_bin1_cluster_genescore.qs"))
  qsave(spe1, save_path)
  cat(paste("Processed and saved:", sample_name, "\n"))
}

### Plotting

#### Load Processed Data
sc1 <- qread("/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/CRCP112_bin1_cluster_genescore.qs")
spe1 <- sc1
plot_df <- as.data.frame(spe1@meta.data)
spring_palette <- c("#f6bcfd", "#ff7d00", "#ffc512", "#ffa300", "#8dd3c6", "#ff6581", "#f8d90d", "#a5da6b", "#e578d6", "#ffd2d8", "#90e4cd", "#f06152", "#7fb8ff", "#d6eaf8","#ff00ff")

#### Cluster UMAP Plot
p3 <- ggplot(plot_df, aes(x = x, y = y, color = seurat_clusters)) +
  geom_point(size = 0.2, alpha = 1) + 
  scale_color_manual(values = spring_palette) + 
  theme_minimal() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

pdf("/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/p112_plot/p112_major_umap_new1.pdf", height = 6, width = 7.5)
print(p3)
dev.off()

#### Spatial Score Plots
score_columns <- c("Anti_inflammatory", "Lipid_loaded_TAM", "Treg", "CD8Tex", "CD8Teff", "Immunosuppression", "Phagocytosis", "Proliferation", "SPI1")

for (score_column in score_columns) {
  p3 <- ggplot(plot_df, aes(x = x, y = y, color = plot_df[[score_column]])) +
    geom_point(size = 0.2, alpha = 1) +
    scale_color_gradient2(low = "#0000FF", mid = "#ADD8E6", high = "#FF0000", midpoint = median(plot_df[[score_column]]), limits = range(plot_df[[score_column]])) +
    theme_minimal() +
    theme(
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  pdf_file_path <- file.path("/public2022/weijinfen/flux/tam_new_analysis/st_data/crc_stereoseq/rds/p112_plot", paste0("p112_", score_column, "_umap_new.pdf"))
  pdf(pdf_file_path, height = 6, width = 7.5)
  print(p3)
  dev.off()
}

#### Violin Plots for Gene Set Scores
df <- plot_df[, c(score_columns, "SCT_snn_res.0.8")]
df_long <- df %>%
  pivot_longer(cols = all_of(score_columns), names_to = "Score_Type", values_to = "Score")

custom_colors <- spring_palette

for (score_type in unique(df_long$Score_Type)) {
  p <- ggplot(df_long[df_long$Score_Type == score_type, ], aes(x = factor(SCT_snn_res.0.8), y = Score, fill = factor(SCT_snn_res.0.8))) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", color = "black", outlier.shape = NA) +
    labs(title = paste(score_type, "Scores in Cluster"), x = "Cluster", y = "Score") +
    theme_minimal() +
    theme(
      legend.position = "right",
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", size = 0.5)
    ) +
    scale_fill_manual(values = custom_colors)
  
  ggsave(filename = paste0(score_type, "_violin_plot.pdf"), plot = p, width = 10, height = 5)
}










































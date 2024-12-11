
library(Seurat)
library(qs)
library(biovizBase)
library(EnsDb.Mmusculus.v79)
library(GenomicRanges)
library(ArchR)
addArchRGenome("hg38")

# Define paths and list files
pathFragments <- "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac_pbs/"
inputFiles <- list.files(pathFragments, pattern = "fragments.tsv.gz", full.names = TRUE)
fileNames <- list.files(pathFragments, pattern = "fragments.tsv.gz$", full.names = FALSE)
cleanFileNames <- gsub(".fragments.tsv.gz", "", fileNames)

sampleNames <- sapply(cleanFileNames, function(x) {
  parts <- unlist(strsplit(x, "_"))
  if (length(parts) >= 2) {
    return(parts[2])
  } else {
    return(x)
  }
})

# Create Arrow files from fragment files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 4, # Minimum TSS fragments per cell
  minFrags = 1000, # Minimum total fragments per cell
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# Initialize ArchR project
proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  geneAnnotation = geneAnnoHg38,
  genomeAnnotation = genomeAnnoHg38,
  outputDirectory = "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac_pbs/"
)

# Load existing project and get marker features
projmouse <- loadRDS("/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/project/projmouse_1.rds")
markersGS <- getMarkerFeatures(
  ArchRProj = projmouse, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS)
saveRDS(markerList, "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/project/marker.rds")

# Combine marker data into a single dataframe and write to CSV
markerDataFrames <- lapply(names(markerList), function(cluster) {
  data <- as.data.frame(markerList[[cluster]])
  data$Cluster <- cluster
  data
})
combinedDataFrame <- do.call(rbind, markerDataFrames)
write.csv(combinedDataFrame, file = "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/project/marker.csv", row.names = FALSE)

# Define marker genes for heatmap
markerGenes <- c(
  "CD34", "GATA1", "PAX5", "MS4A1", "EBF1", "MME", "CD14", "CEBPB", "MPO", "IRF8", "CD3D", "CD8A", "TBX21", "IL7R",
  "EPCAM", "KRT19", "CDH1", "KRT18", "CD68", "MARCO", "FCGR3A", "LYZ", "NCAM1", "NKG7", "GNLY", "KLRD1", "KIT", "MS4A2", "GATA2",
  "DCN", "COL1A1", "COL1A2", "THY1", "PECAM1", "CLDN5", "FLT1", "RAMP2"
)
markerGenes <- tolower(markerGenes)

# Generate and draw heatmap
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# Capitalize first letter of each word in marker genes
capitalizeFirstLetter <- function(x) {
  paste0(toupper(substring(x, 1, 1)), substring(x, 2))
}
capitalizedMarkerGenes <- sapply(markerGenes, capitalizeFirstLetter, USE.NAMES = FALSE)
print(capitalizedMarkerGenes)

# Plot UMAP embedding colored by GeneScoreMatrix for Col1a1
p <- plotEmbedding(
  ArchRProj = projmouse,
  colorBy = "GeneScoreMatrix",
  name = "Col1a1",
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

# Add Impute Weights for MAGIC imputation
projmouse <- addImputeWeights(projmouse)

# Save plots as PDFs
plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Col1a1-Imputation.pdf", 
        ArchRProj = projmouse, 
        addDOC = FALSE, width = 5, height = 5)

# Plot browser track for Epcam gene
p <- plotBrowserTrack(
  ArchRProj = projmouse, 
  groupBy = "Clusters", 
  geneSymbol = "Epcam", 
  upstream = 50000,
  downstream = 50000
)
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projmouse, 
        addDOC = FALSE, width = 5, height = 5)

# Launch ArchR Browser and save updated project
ArchRBrowser(projmouse)
saveRDS(projmouse, "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/project/projmouse_2.rds")


library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

# Load previously saved ATAC project
atac2 <- addGroupCoverages(ArchRProj = atac1, groupBy = "Clusters2")

# Ensure MACS2 is installed and available
pathToMacs2 <- findMacs2()

# 10.2 Calling Peaks with MACS2
atac2 <- addReproduciblePeakSet(
    ArchRProj = atac2, 
    groupBy = "Clusters2", 
    pathToMacs2 = pathToMacs2
)

# 10.3 Calling Peaks with TileMatrix
atac2 <- addReproduciblePeakSet(
    ArchRProj = atac2, 
    groupBy = "Clusters2",
    peakMethod = "Tiles",
    method = "p"
)
getPeakSet(atac2)

# 10.4 Add Peak Matrix
atac2 <- addPeakMatrix(atac2)
getAvailableMatrices(atac2)

# 11.1 Identifying Marker Peaks
table(atac2$Clusters2)
markersPeaks <- getMarkerFeatures(
    ArchRProj = atac2, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clusters2",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05 & Log2FC >= 1", returnGR = TRUE)

# 11.2 Plotting Marker Peaks
heatmapPeaks <- markerHeatmap(
    seMarker = markersPeaks, 
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
    transpose = TRUE
)
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap.pdf", width = 8, height = 6, ArchRProj = atac2, addDOC = FALSE)

pma <- markerPlot(seMarker = markersPeaks, name = "apoe+tam", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markersPeaks, name = "apoe+tam", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "apoe-Markers-MA-Volcano.pdf", width = 5, height = 5, ArchRProj = atac2, addDOC = FALSE)

# Browser track visualization
p <- plotBrowserTrack(
    ArchRProj = atac2, 
    groupBy = "Clusters2", 
    geneSymbol = "Apoe",
    features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)[["apoe+tam"]],
    upstream = 50000,
    downstream = 50000
)
plotPDF(p, name = "Plot-Tracks-With-Features_Apoe.pdf", width = 5, height = 5, ArchRProj = atac2, addDOC = FALSE)

# 11.3 Pairwise Testing Between Groups
markerTest <- getMarkerFeatures(
    ArchRProj = atac2, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters2",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "apoe+tam",
    bgdGroups = "cycling tam"
)

pma <- markerPlot(seMarker = markerTest, name = "apoe+tam", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- markerPlot(seMarker = markerTest, name = "apoe+tam", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "apoe-vs-cycling-tam-Markers-MA-Volcano.pdf", width = 5, height = 5, ArchRProj = atac2, addDOC = FALSE)

# 12.1 Motif Enrichment in Differential Peaks
atac2 <- addMotifAnnotations(ArchRProj = atac2, motifSet = "cisbp", name = "Motif")

# apoe+tam vs cycling tam
markerTest <- getMarkerFeatures(
    ArchRProj = atac2, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters2",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "apoe+tam",
    bgdGroups = "cycling tam"
)

# apoe+tam vs Monocyte
markerTestMonocyte <- getMarkerFeatures(
    ArchRProj = atac2, 
    useMatrix = "PeakMatrix",
    groupBy = "Clusters2",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = "apoe+tam",
    bgdGroups = "Monocyte"
)

# Motif Enrichment Analysis
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = atac2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.4"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 3,
        nudge_x = 2,
        color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(P-adj) Motif Enrichment") + 
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 18)
    )

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = atac2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
    geom_point(size = 1) +
    ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 3,
        nudge_x = 2,
        color = "black"
    ) + theme_ArchR() + 
    ylab("-log10(FDR) Motif Enrichment") +
    xlab("Rank Sorted TFs Enriched") +
    scale_color_gradientn(colors = paletteContinuous(set = "comet")) +
    theme(
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        plot.title = element_text(size = 20),
        plot.subtitle = element_text(size = 18)
    )
plotPDF(ggUp, ggDo, name = "apoe-vs-Monocyte-Markers-Motifs-Enriched.pdf", width = 8, height = 8, ArchRProj = atac2, addDOC = FALSE)

# 12.2 Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = atac2,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

# 13.1 Motif Deviations
saveRDS(atac2, "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/project/atac3.rds")

if(!"Motif" %in% names(atac2@peakAnnotation)){
    atac2 <- addMotifAnnotations(ArchRProj = atac2, motifSet = "cisbp", name = "Motif")
}

atac2 <- addBgdPeaks(atac2)
atac2 <- addDeviationsMatrix(
    ArchRProj = atac2, 
    peakAnnotation = "Motif",
    force = TRUE
)    

plotVarDev <- getVarDeviations(atac2, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores.pdf", width = 10, height = 10, ArchRProj = atac2, addDOC = FALSE)

motifs <- c("Sfpi1", "Spic", "Bcl11b", "Elf5", "Spib", "Elf1","Stat2", "Irf1","Ehf", "Elf3", "Irf1", "cl11a", "Stat1", "Stat3", "Stat5a", "Stat2")

markerMotifs <- getFeatures(atac2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

p <- plotGroups(ArchRProj = atac2, 
    groupBy = "Clusters2", 
    colorBy = "MotifMatrix", 
    name = markerMotifs,
    imputeWeights = getImputeWeights(atac2)
)
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation.pdf", width = 6, height = 6, ArchRProj = atac2, addDOC = FALSE)

p2 <- lapply(seq_along(p), function(x){
    if(x != 1){
        p[[x]] + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank()) + ylab("")
    }else{
        p[[x]] + guides(color = FALSE, fill = FALSE) + 
        theme_ArchR(baseSize = 6) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        theme(axis.ticks.y=element_blank(), axis.title.y=element_blank()) + ylab("")
    }
})
combined_plot <- do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))
ggsave(filename = "Plot-Groups-Deviations-w-Imputation_p2.pdf", plot = combined_plot, width = 15, height = 15, units = "in")

# UMAP-TF score
p <- plotEmbedding(
    ArchRProj = atac2, 
    colorBy = "MotifMatrix", 
    name = sort(markerMotifs), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(atac2)
)
plotPDF(p, name = "umap-Groups-Deviations-w-Imputation.pdf", width = 6, height = 6, ArchRProj = atac2, addDOC = FALSE)

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
})
combined_plot <- do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
ggsave(filename = "umap1-TF-score-Groups-Deviations-w-Imputation_p2.pdf", plot = combined_plot, width = 12, height = 15, units = "in")

# UMAP-TF ATAC expression
markerRNA <- getFeatures(atac2, select = paste(motifs, collapse="|"), useMatrix = "GeneScoreMatrix")
markerRNA <- c("Stat1","Stat2", "Stat5a", "Stat3", "Spi1")
p <- plotEmbedding(
    ArchRProj = atac2, 
    colorBy = "GeneScoreMatrix", 
    name = sort(markerRNA), 
    embedding = "UMAP",
    imputeWeights = getImputeWeights(atac2)
)
plotPDF(p, name = "umap-gene-Groups-Deviations-w-Imputation.pdf", width = 6, height = 6, ArchRProj = atac2, addDOC = FALSE)


library(ArchR)
library(GenomicRanges)

# 14.1 Motif Footprinting

# Retrieve motif positions from the ATAC project
motifPositions <- getPositions(atac2)

# Define motifs of interest
motifs <- c("Sfpi1", "Spic", "Bcl11b", "Elf5", "Spib", "Elf1", "Stat2", "Irf1", "Ehf", "Elf3", "Irf1", "cl11a")

# Extract marker motifs by matching motif names to available motif positions
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))

# Get footprints for selected motifs grouped by clusters
seFoot <- getFootprints(
    ArchRProj = atac2, 
    positions = motifPositions[markerMotifs], 
    groupBy = "Clusters2"
)

# 14.2 Subtracting the Tn5 Bias

# Plot footprints with subtraction normalization to remove Tn5 bias
plotFootprints(
    seFoot = seFoot,
    ArchRProj = atac2, 
    normMethod = "Subtract",
    plotName = "Footprints-Subtract-Bias.pdf",
    addDOC = FALSE,
    smoothWindow = 5
)

# Plot footprints with division normalization to adjust for Tn5 bias
plotFootprints(
    seFoot = seFoot,
    ArchRProj = atac2, 
    normMethod = "Divide",
    plotName = "Footprints-Divide-Bias.pdf",
    addDOC = FALSE,
    smoothWindow = 5
)

# Plot footprints without any normalization
plotFootprints(
    seFoot = seFoot,
    ArchRProj = atac2, 
    normMethod = "None",
    plotName = "Footprints-No-Normalization.pdf",
    addDOC = FALSE,
    smoothWindow = 5
)

# 14.3 Feature Footprinting Around Transcription Start Sites (TSS)

# Add group coverages if not already added
if (!"GroupCoverages" %in% names(atac2)) {
    atac2 <- addGroupCoverages(ArchRProj = atac2, groupBy = "Clusters2")
}

# Get footprints around TSSs with specified flanking regions
seTSS <- getFootprints(
    ArchRProj = atac2, 
    positions = GRangesList(TSS = getTSS(atac2)), 
    groupBy = "Clusters2",
    flank = 2000
)

# Plot TSS footprints without normalization but with flank normalization
plotFootprints(
    seFoot = seTSS,
    ArchRProj = atac2, 
    normMethod = "None",
    plotName = "TSS-No-Normalization.pdf",
    addDOC = FALSE,
    flank = 2000,
    flankNorm = 100
)


library(ArchR)
library(GenomicRanges)
library(ggplot2)
library(ggrepel)

# 15.2 Co-accessibility Analysis with ArchR

# Add co-accessibility to the ATAC project using IterativeLSI dimensions
atac2 <- addCoAccessibility(
    ArchRProj = atac2,
    reducedDims = "IterativeLSI"
)

# Get co-accessibility loops with correlation cutoff and specified resolution
cA <- getCoAccessibility(
    ArchRProj = atac2,
    corCutOff = 0.5,
    resolution = 10000,
    returnLoops = TRUE
)

# Define marker genes for plotting browser tracks
markerGenes <- c("Mmp12", "Apoe", "Arg1", "Fabp5", "Pf4", "Lyz2", "Hmox1", "Wfdc17", "Ctsb", "Ccl6", "Fth1", "Lgals3", "Ctsd", "Ctss")

# Plot browser tracks including co-accessibility loops
p <- plotBrowserTrack(
    ArchRProj = atac2, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = cA
)

# Save the plot as PDF
plotPDF(plotList = p, 
        name = "Plot-1-Tracks-Marker-Genes-with-CoAccessibility_2.pdf", 
        ArchRProj = atac2, 
        addDOC = FALSE, width = 5, height = 5)

# 15.3 Peak-to-Gene Linkage with ArchR

# Add peak-to-gene links using IterativeLSI dimensions
atac2 <- addPeak2GeneLinks(
    ArchRProj = atac2,
    reducedDims = "IterativeLSI"
)

# Get peak-to-gene links with correlation cutoff and specified resolution
p2g <- getPeak2GeneLinks(
    ArchRProj = atac2,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = TRUE
)

# Define marker genes for plotting browser tracks
markerGenes <- c("Nr1h3", "Npc1", "Npc2", "Soat1", "Lipa", "Psap", "Cyp27a1", "Ch25h", "Abca1", "Abcg1", "Sirpa", "Olr1", "Trem2", "Scarb2", "Msr1", "Cd36")

# Plot browser tracks including peak-to-gene links
p <- plotBrowserTrack(
    ArchRProj = atac2, 
    groupBy = "Clusters2", 
    geneSymbol = markerGenes, 
    upstream = 50000,
    downstream = 50000,
    loops = p2g
)

# Save the plot as PDF
plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
        ArchRProj = atac2, 
        addDOC = FALSE, width = 5, height = 5)

# Plot peak-to-gene heatmap
p <- plotPeak2GeneHeatmap(ArchRProj = atac2, groupBy = "Clusters2")
plotPDF(p, 
        name = "plot-Peak2-Gene-Heatmap.pdf", 
        ArchRProj = atac2, 
        addDOC = FALSE, width = 10, height = 20)

# Replot heatmap with k=10
p <- plotPeak2GeneHeatmap(ArchRProj = atac2, groupBy = "Clusters2", k=10)
plotPDF(p, 
        name = "plot-Peak2-Gene-Heatmap_k10.pdf", 
        ArchRProj = atac2, 
        addDOC = FALSE, width = 15, height = 20)

# 15.4 Identification of Positive TF-Regulators

# Get grouped sparse encoding (SE) from MotifMatrix
seGroupMotif <- getGroupSE(ArchRProj = atac2, useMatrix = "MotifMatrix", groupBy = "Clusters2")

# Filter SE for 'z' seqnames and calculate maxDelta
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- Reduce("cbind", lapply(seq_len(ncol(seZ)), function(x) rowMaxs(assay(seZ) - assay(seZ)[,x]))) %>% rowMaxs

# Correlate GeneScoreMatrix with MotifMatrix
corGSM_MM <- correlateMatrices(
    ArchRProj = atac2,
    useMatrix1 = "GeneScoreMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

# Correlate GeneIntegrationMatrix with MotifMatrix
corGIM_MM <- correlateMatrices(
    ArchRProj = atac2,
    useMatrix1 = "GeneIntegrationMatrix",
    useMatrix2 = "MotifMatrix",
    reducedDims = "IterativeLSI"
)

# Annotate maxDelta for correlated matrices
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

# Identify significant TF regulators for GeneScoreMatrix correlation
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.2 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"

# Plot volcano plot for GeneScoreMatrix correlation
p <- ggplot(data = data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(corGSM_MM$maxDelta)*1.05)) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(TFRegulator == "YES", GeneScoreMatrix_name, "")),
    size = 3,
    nudge_x = 0.1,
    nudge_y = 0.1,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = "grey50",
    inherit.aes = TRUE
  ) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/Plots/volcano_plot_with_labels.pdf", plot = p, width = 13, height = 13, units = "in")

# Identify significant TF regulators for GeneIntegrationMatrix correlation
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.3 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"

# Plot volcano plot for GeneIntegrationMatrix correlation
p <- ggplot(data = data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(expand = c(0,0), limits = c(0, max(corGSM_MM$maxDelta)*1.05)) +
  ggrepel::geom_text_repel(
    aes(label = ifelse(TFRegulator == "YES" | (maxDelta > 5 & padj < 0.05 & cor > 0.3), GeneScoreMatrix_name, "")),
    size = 3,
    nudge_x = 0.1,
    nudge_y = 0.1,
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.5, "lines"),
    segment.color = "grey50",
    inherit.aes = TRUE
  ) +
  theme(axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 14), axis.text.y = element_text(size = 14))

ggsave(filename = "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/atac/Plots/volcano_plot_with_labels-apoe.pdf", plot = p, width = 13, height = 13, units = "in")



library(Seurat)
library(dplyr)
library(ggplot2)

# 1. Preprocessing and Clustering scRNA-seq Data
# Normalize data
sce <- NormalizeData(sce, 
                     normalization.method = "LogNormalize", 
                     scale.factor = 1e4)

# Identify variable features
sce <- FindVariableFeatures(sce)

# Scale data
sce <- ScaleData(sce)

# Perform PCA
sce <- RunPCA(sce, features = VariableFeatures(object = sce))

# Perform UMAP reduction
sce <- RunUMAP(sce, dims = 1:20, reduction = "pca")

# Plot UMAP with cluster labels
p <- DimPlot(sce, reduction = "umap", label = TRUE)

# Find neighbors and clusters
sce <- FindNeighbors(sce, dims = 1:20)
sce <- FindClusters(sce, resolution = 0.2)

# Identify all markers using Wilcoxon test
markers <- FindAllMarkers(object = sce, test.use = "wilcox", only.pos = TRUE, logfc.threshold = 0.25)

# 2. Annotating Clusters Based on Marker Genes

# Load single-cell data
sc <- qread("/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/sc_data_1.qs")
sce <- sc

# Plot initial UMAP with Seurat clusters labeled
DimPlot(sce, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)

# Define cell types based on clusters
marker <- data.frame(
  cluster = 0:7,
  cell = as.character(0:7)
)

marker$cell <- case_when(
  marker$cluster %in% c(0, 2, 7) ~ 'Fibroblast',
  marker$cluster %in% c(1, 4) ~ 'Epithelial',
  marker$cluster %in% c(5) ~ 'Endothelial',
  marker$cluster %in% c(3) ~ 'Myeloid',
  marker$cluster %in% c(6) ~ 'T NK cell',
  TRUE ~ as.character(marker$cluster)
)

# Add major cell type annotations to metadata
sce@meta.data$major <- sapply(sce@meta.data$seurat_clusters, function(x) {
  marker[marker$cluster == x, 'cell']
})

# Extract Myeloid cells for further analysis
mye <- sce[, sce$major == "Myeloid"]

# Re-run preprocessing and clustering for Myeloid cells
mye <- NormalizeData(mye, normalization.method = "LogNormalize", scale.factor = 1e4)
mye <- FindVariableFeatures(mye)
mye <- ScaleData(mye)
mye <- RunPCA(mye, features = VariableFeatures(object = mye))
mye <- RunUMAP(mye, dims = 1:20, reduction = "pca")
p <- DimPlot(mye, reduction = "umap", label = TRUE)
mye <- FindNeighbors(mye, dims = 1:20)
mye <- FindClusters(mye, resolution = 0.1)
markers_mye <- FindAllMarkers(object = mye, test.use = "wilcox")

# Save markers to CSV
write.csv(markers_mye, "/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/sc_mye_marker.csv")

# Define subtypes for Myeloid cells
marker_subtype <- data.frame(
  cluster = 0:4,
  cell = as.character(0:4)
)

marker_subtype$cell <- case_when(
  marker_subtype$cluster %in% c(0) ~ 'apoe+tam',
  marker_subtype$cluster %in% c(1) ~ 'Fibroblast',
  marker_subtype$cluster %in% c(2) ~ 'cycling tam',
  marker_subtype$cluster %in% c(3) ~ 'B cell',
  marker_subtype$cluster %in% c(4) ~ 'monocyte',
  TRUE ~ as.character(marker_subtype$cluster)
)

# Add subtype annotations to metadata
mye@meta.data$subtype <- sapply(mye@meta.data$RNA_snn_res.0.1, function(x) {
  marker_subtype[marker_subtype$cluster == x, 'cell']
})

# Update sce object with subtype information
sce@meta.data$subtype <- sce@meta.data$major
sce@meta.data$subtype[sce$major == "Myeloid"] <- mye@meta.data$subtype

# Add BioClassification column for integration with ATAC data
sce@meta.data$BioClassification <- sce@meta.data$subtype

# 3. Visualization of Cluster Annotations

# Load annotated single-cell data
sc <- qread("/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE220635/sc_data_celltype.qs")

# Define colors for plotting
color1 <- c("#55b7e6", "#56ba77", "#9faa3f", "#bd79b6", "#930E14", "#F8D5F4", "#96D1C6", "#BAD65D", "#EE7D6B", "#C8C1DE", "#2681B6","#ed6e69")

# Plot UMAP colored by subtype
p1 <- DimPlot(sc, reduction = "umap", group.by = "subtype", cols = color1) +
      NoAxes() + 
      theme(panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))

pdf("/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE22065/atac/Plots/umap_scdata_celltype.pdf", height=4, width=5)
print(p1)
dev.off()

# Plot feature expression for genes of interest
genes_of_interest <- c("Spi1")
color <- c("#f4f2ef", "blue")

p <- FeaturePlot(sc, features = genes_of_interest, cols = color, pt.size = 0.1, max = 2) +
     theme(panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"))

pdf("/public2022/weijinfen/flux/tam_new_analysis/atac_seq/GSE22065/atac/Plots/umap_scdata_Spi1.pdf", height=4, width=4)
print(p)
dev.off()






























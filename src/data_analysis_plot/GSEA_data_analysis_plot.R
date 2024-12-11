










# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(GSVA)
library(Seurat)
library(qs)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(future)
library(pheatmap)
library(ggsignif)
library(cowplot)
library(paletteer)
library(gplots)
library(stringr)
library(tibble)
library(dplyr)
library(tidyr)
library(data.table)


# Read the DEG data
gene <- read.csv("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/flt_tam.csv")

# Convert gene symbols to Entrez IDs
df_id <- bitr(gene$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
df_all <- merge(gene, df_id, by = "SYMBOL", all = FALSE)

# Write the combined data to CSV
write.csv(df_all, "/public2022/weijinfen/flux/cesc/tcga_result/edger_tumor_normal_all_DEG_enid.csv")

# Filter significant DEGs based on PValue
df_all_sort <- df_all[df_all$PValue < 0.05, ]

# Sort by logFC in descending order
df_all_sort <- df_all_sort[order(-df_all_sort$logFC), ]

# Create a named vector of logFC values using Entrez IDs
gene_fc <- setNames(df_all_sort$logFC, df_all_sort$ENTREZID)


# Perform KEGG pathway enrichment analysis
KEGG <- gseKEGG(gene_fc, organism = "hsa")

# Perform GO enrichment analysis (biological process)
GO_BP <- gseGO(geneList = gene_fc,
               ont = "BP",
               maxGSSize = 500,
               keyType = "ENTREZID",
               OrgDb = org.Hs.eg.db)

# Save sorted GO results
sortGO_BP <- GO_BP@result[order(GO_BP@result$NES, decreasing = TRUE), ]
write.csv(sortGO_BP, "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/flt_tam_GO_gsea.csv")


# Select specific GO terms for visualization
paths <- c("GO:0006954")  # Example GO term

# Plot GO enrichment results
p1 <- gseaplot2(GO_BP, paths, pvalue_table = FALSE)

p2 <- gseaplot2(
  GO_BP,
  paths,
  color = 'red',
  rel_heights = c(1.5, 0.5, 1),
  subplots = 1:2,
  pvalue_table = FALSE,
  title = 'inflammatory response',
  ES_geom = 'line'
)

p3 <- gseaplot(
  GO_BP,
  paths,
  color.line = 'red',
  by = "runningScore",
  title = 'lipid localization'
)
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/apoe_gesa_lipid_localization.pdf", plot = p3, width = 8, height = 8, units = "in")

# Combine multiple pathways for visualization
top_go <- c("GO:0010876", "GO:0006869", "GO:0006629", "GO:0008202", "GO:0030301", "GO:0006955", "GO:0015918", "GO:0006897")
top_go_gsea <- GO_BP
top_go_gsea@result <- top_go_gsea@result[top_go, ]
top_go_gsea@result$Description <- paste0('NES=', round(top_go_gsea@result$NES, 2), " ", top_go_gsea@result$Description)

pdf(file = "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/apoe_gesa_lipid_localization_2.pdf", width = 8, height = 9)
p4 <- gseaplot2(
  top_go_gsea,
  top_go,
  color = c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#00ff00', '#ff00ff')
)
print(p4)
dev.off()


# Read sorted GO terms data
dt <- read.csv("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/cycle_tam_GO_gsea_sort.csv")
dt$Description <- factor(dt$Description, levels = rev(dt$Description))

# Define custom theme
mytheme <- theme(
  legend.position = 'none',
  plot.title = element_text(size = 14, face = 'bold'),
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.ticks.y = element_blank()
)

# Create bar plot for GO terms
p <- ggplot() +
  geom_bar(data = dt,
           aes(x = NES, y = Description),
           width = 0.7,
           fill = "#E41A1C",
           alpha = 0.7,
           stat = 'identity') +
  labs(title = 'Biological processes of Cycling+TAM') +
  theme_classic() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(axis.text.y = element_blank()) +
  geom_text(data = dt,
            aes(x = 0.1,
                y = Description,
                label = Description),
            size = 7,
            hjust = 0) +
  mytheme

# Save the plot
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/tam_gene_gsea_pdf/cycle_tam_bar.pdf", plot = p, width = 8, height = 6, units = "in")


# Load necessary library
library(dplyr)

# Read the gene set data from CSV
gene_set <- read.csv("/public2022/weijinfen/flux/cesc/gsea/glycosylation.csv", stringsAsFactors = FALSE)

# Split genes by their set names
gene_sets <- split(gene_set$gene, gene_set$name)

# Define write.gmt function to write gene sets to a GMT file
write.gmt <- function(gs, file) {
  # Open connection to file
  con <- file(file, open = "wt")
  
  # Iterate over each gene set and write to file
  for (setName in names(gs)) {
    genes <- gs[[setName]]
    if (length(genes) > 0) {
      line <- paste(setName, "tmp", paste(genes, collapse = "\t"), sep = "\t")
      writeLines(line, con)
    }
  }
  
  # Close the file connection
  close(con)
}

# Specify the output GMT file path
file <- "/public2022/weijinfen/flux/cesc/gsea/glycosylation_1.gmt"

# Write the gene sets to the GMT file
write.gmt(gene_sets, file)

# Confirm the operation
cat("GMT file has been created successfully at:", file, "\n")














































































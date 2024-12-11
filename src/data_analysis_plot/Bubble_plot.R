# Load necessary libraries
library(Seurat)
library(qs)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggh4x)
library(rstatix)
library(future)
library(pheatmap)

# Read the data
df <- qread("/public2022/weijinfen/flux/singlecell_data/gene_exp/crc8_all_analysis_gene_exp_inallcells.qs")

# Filter data for tumor samples
df_tumor <- df[df$sampletype_rev2 == "Tumor", ]
modules <- list(
  biosynthesis = c("SQLE", "LSS", "LBR", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "ERG28", "SC5D", "EBP", "DHCR24", "DHCR7", "HMGCS1", "HMGCS2", "GGPS1", "FDFT1", "FDPS", "IDI1", "MVK", "PMVK", "MVD", "ACAT1", "ACAT2", "HMGCR", "ACLY", "ACSS1", "ACSS2"),
  uptake_and_regulation = c("OLR1", "TREM2", "MSR1", "CD36", "LDLRAP1", "SREBF1", "SREBF2", "NPC1L1", "NUMB", "LIMA1", "SCARB1", "SCARB2", "LDLR", "LRP1"),
  efflux = c("ABCG1", "ABCA1"),
  esterification = c("SOAT1", "SOAT2", "LIPA"),
  transport = c("NPC1", "NPC2"),
  regulation = c("PCSK9", "NR1H3", "NR1H2"),
  oxidization = c("CYP27A1", "CH25H"),
  negative_regulation = c("SIRPA", "SIGLEC10"),
  prosaposin = c("PSAP")
)

# Combine all module data into one dataframe
df_all_modules <- do.call(rbind, lapply(names(modules), function(module_name) {
  module_genes <- modules[[module_name]]
  df_filtered <- df_tumor[, c("major", module_genes)]
  df_long <- reshape2::melt(df_filtered, id.vars = "major")
  df_long$gene <- df_long$variable
  df_long$module <- module_name
  return(df_long)
}))

# Summarize mean expression and expression percentage
df_summary <- df_all_modules %>%
  group_by(major, gene, module) %>%
  summarise(
    mean_expression = mean(value, na.rm = TRUE),
    percentage = sum(value > 0) / n()
  ) %>%
  ungroup()

# Create a new variable containing module and gene information
df_summary$module_gene <- paste(df_summary$module, df_summary$gene, sep = ": ")

# Calculate overall mean expression for normalization
mean_expression_all <- tapply(df_summary$mean_expression, df_summary$gene, mean)

# Normalize expression
df_summary$normalized_expression <- df_summary$mean_expression / mean_expression_all[df_summary$gene]

# Preprocess percentage column
df_summary$percentage <- pmin(pmax(df_summary$percentage, 0.1), 0.6)

# Define module colors
module_colours <- c("#e8736b", "#3a7c38", "#459093", "#e8c2aa", "#7cb9e8", "#d73027", "#fdae61", "#66c2a5", "#fc8d62")

# Create bubble plot
bubble_plot <- ggplot(df_summary, aes(x = interaction(gene, module), y = major)) + 
  guides(x = "axis_nested") +
  geom_point(aes(fill = normalized_expression, size = percentage, colour = module), shape = 21) +
  scale_colour_manual(values = module_colours) + # Add module color scheme
  theme_bw(base_size = 14) +
  xlab('Gene') + ylab('Cell Type') +  
  scale_fill_gradient2(
    low = 'white', 
    mid = '#EB1D36', 
    high = '#990000',
    midpoint = (min(df_summary$normalized_expression, na.rm = TRUE) + max(df_summary$normalized_expression, na.rm = TRUE)) / 2,
    limits = c(0, max(df_summary$normalized_expression, na.rm = TRUE)),
    breaks = c(0, max(df_summary$normalized_expression, na.rm = TRUE)),
    labels = sprintf("%.2f", c(0, max(df_summary$normalized_expression, na.rm = TRUE))),
    name = 'Normalized Expression'
  ) +  
  scale_size(
    range = c(0.5, 3),  # Reduce range to make bubbles smaller
    limits = c(0.1, 0.6),
    breaks = seq(0.1, 0.6, 0.1),
    labels = seq(10, 60, 10)
  ) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(color = 'black'),
    aspect.ratio = 0.5,
    plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = 'cm'),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = 'italic', size = 6),
    axis.text.y = element_text(size = 10),
    ggh4x.axis.nestline.x = element_line(size = 0.6),
    ggh4x.axis.nesttext.x = element_text(colour = module_colours, size = 10)
  ) +
  coord_cartesian(clip = 'off') +
  guides(
    fill = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      barwidth = unit(5, "cm")
    ),
    size = guide_legend(
      title = "Percent of Cells",
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      override.aes = list(color = "black", fill = "grey")
    ),
    colour = guide_legend(title = "Module") # Add module legend
  )

# Save the plot
pdf("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_diff_allcell_pdf/all_progression_diff_major_bubble_new.pdf", width = 10, height = 6)
print(bubble_plot)
dev.off()

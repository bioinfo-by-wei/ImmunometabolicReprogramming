library(METAFlux)
library(qs)
library(Seurat)
library(dplyr)
library(tidyr)
library(reshape2)
data("human_blood")
data("human_gem")
data("nutrient_lookup_files")

# Add cell type names to the flux file
# Single file processing
df_flux <- qread("/public2022/weijinfen/flux/scflux/crc8/C109_sc_flux_raw.qs", nthreads = 48)
df_flux <- as.data.frame(df_flux)
df_flux$celltype <- sub("^([^ ]+ [^ ]+).*", "\\1", rownames(df_flux))
df_flux$HMR <- sub("^[^ ]+ [^ ]+ ", "", rownames(df_flux))
# Move HMR and celltype to the first two columns
df_flux <- select(df_flux, celltype, HMR, everything())
# Score file
df_score <- qread("/public2022/weijinfen/flux/scores/crc8/C109_scores.qs", nthreads = 48)
# Map cell types from df_score$cell_type_fractions to df_flux$celltype
celltype_names <- names(df_score$cell_type_fractions)
celltype_map <- setNames(celltype_names, paste0("celltype", 1:length(celltype_names)))
df_flux <- df_flux %>%
  mutate(celltype_num = as.integer(gsub("celltype", "", celltype))) %>%
  mutate(celltype = celltype_map[paste0("celltype", celltype_num)]) %>%
  select(-celltype_num)
qsave(df_flux, "C108_flux_cell.qs", nthreads = 48)

# Batch process multiple files
# Define file paths
flux_path <- "/public2022/weijinfen/flux/tam_new_analysis/scflux_raw/crc8/"
score_path <- "/public2022/weijinfen/flux/tam_new_analysis/scores/"
output_path <- "/public2022/weijinfen/flux/tam_new_analysis/scflux_cell/"
flux_files <- list.files(flux_path, pattern = "_sc_flux_raw.qs$", full.names = TRUE)
score_files <- list.files(score_path, pattern = "_scores.qs$", full.names = TRUE)
# Extract sample names
sample_names <- sub("_sc_flux_raw.qs$", "", basename(flux_files))
# Define processing function
process_sample <- function(sample_name) {
  # Build file paths
  flux_file <- file.path(flux_path, paste0(sample_name, "_sc_flux_raw.qs"))
  score_file <- file.path(score_path, paste0(sample_name, "_scores.qs"))
  output_file <- file.path(output_path, paste0(sample_name, "_flux_cell.qs"))
  # Read flux file
  df_flux <- qread(flux_file, nthreads = 48)
  df_flux <- as.data.frame(df_flux)
  # Add two new columns
  df_flux$celltype <- sub("^([^ ]+ [^ ]+).*", "\\1", rownames(df_flux))
  df_flux$HMR <- sub("^[^ ]+ [^ ]+ ", "", rownames(df_flux))
  # Move HMR and celltype to the first two columns
  df_flux <- select(df_flux, celltype, HMR, everything())
  # Read score file
  df_score <- qread(score_file, nthreads = 48)
  # Map cell types from df_score$cell_type_fractions to df_flux$celltype
  celltype_names <- names(df_score$cell_type_fractions)
  celltype_map <- setNames(celltype_names, paste0("celltype", 1:length(celltype_names)))
  df_flux <- df_flux %>%
    mutate(celltype_num = as.integer(gsub("celltype", "", celltype))) %>%
    mutate(celltype = celltype_map[paste0("celltype", celltype_num)]) %>%
    select(-celltype_num)
  # Save processed file
  qsave(df_flux, output_file, nthreads = 48)
}
lapply(sample_names, process_sample)

# Integrate major cell type flux information - original data is calculated for subtypes
# Define a function to extract sample ID and the following character (normal/tumor marker)
extract_sample_and_status <- function(batch_id) {
  underscore_pos <- regexpr("_", batch_id)
  if (underscore_pos != -1) {
    sample_info <- substr(batch_id, 1, underscore_pos)
    sample_and_status <- paste0(sample_info, substr(batch_id, underscore_pos + 1, underscore_pos + 1))
    return(sample_and_status)
  } else {
    stop("Invalid batch ID format.")
  }
}
# Apply function to batchID vector
simplified_batchID <- sapply(crc@meta.data$batchID, extract_sample_and_status)
# Create a data frame with original and simplified batchID
batchID_mapping <- data.frame(
  original_batchID = names(simplified_batchID),
  simplified_batchID = simplified_batchID
)
# Split cells by original batchID
cells_by_batchID <- split(seq_len(ncol(crc)), crc@meta.data$batchID)
# Initialize a vector for simplified batchID per cell
simplified_batchID_per_cell <- character(ncol(crc))
# Iterate over each batchID group
for (batchID in unique(crc@meta.data$batchID)) {
  simplified_ID <- batchID_mapping[batchID_mapping$original_batchID == batchID, "simplified_batchID"]
  simplified_batchID_per_cell[cells_by_batchID[[batchID]]] <- simplified_ID
}
# Add simplified batchID to Seurat object metadata
crc$sampleBatchID <- simplified_batchID_per_cell

# Major and revised subtype cell relationship
# Create a data frame if not already created
cell_classification <- data.frame(
  major = crc$major,
  revised_subtype = crc$revised_subtype
)
# Group by major and revised_subtype, count cell numbers for each revised_subtype
unique_cell_relations <- cell_classification %>%
  group_by(revised_subtype, major) %>%
  summarise(cell_count = n()) %>%
  distinct(revised_subtype, .keep_all = TRUE)
# Write results to CSV
write.csv(unique_cell_relations, "/public2022/weijinfen/flux/cellflux/crc8/crc8_major_subtype_relation.csv")

# Calculate cell proportions for crc single-cell data
# Cell proportions
metadata <- as.data.frame(crc@meta.data)
metadata$sampleBatchID <- as.character(metadata$sampleBatchID)
metadata$major <- as.factor(metadata$major)
# Count each major cell type in each sample
cell_counts_by_sample <- metadata %>%
  group_by(sampleBatchID, major) %>%
  summarise(count = n())
# Total cell count per sample
total_cells_per_sample <- metadata %>%
  group_by(sampleBatchID) %>%
  summarise(total_count = n())
# Join two data frames to calculate proportions
cell_proportions <- cell_counts_by_sample %>%
  left_join(total_cells_per_sample, by = "sampleBatchID") %>%
  mutate(proportion = count / total_count)
# Convert result to desired dataframe format
proportions_df <- cell_proportions %>%
  select(sampleBatchID, major, proportion) %>%
  spread(key = major, value = proportion)
# View results
print(proportions_df)
# Optionally save results to CSV
write.csv(proportions_df, file = "/public2022/weijinfen/flux/cellflux/crc8/crc8_majorcell_proportions_major.csv")



# Comparison of metabolite accumulation/consumption across cell subtypes within a single sample —— Boxplot with significance test —— 10.25
# Differences of a specific metabolite among cell subpopulations across multiple samples
#### Create a random color mapping for each cell type
random_colors <- runif(n_distinct(long_data$Celltype), min = 0, max = 1)
color_mapping <- setNames(hcl(h = random_colors * 360, c = 100, l = 65, alpha = 1), unique(long_data$Celltype))

#### Read in data
df_all_normalized <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_subtype_tumor/phospholipids extracellular pool[s]_uptake_score_tumor.qs")
library(tidyr)
# Convert data to long format
long_data <- df_all_normalized %>%
  pivot_longer(cols = -Sample, names_to = "Celltype", values_to = "Value")
# Calculate global Kruskal-Wallis test p-value
kruskal_result <- kruskal.test(Value ~ Celltype, data = long_data)
p_value <- signif(kruskal_result$p.value, digits = 3)
# Create boxplot
p <- ggplot(long_data, aes(x = Celltype, y = Value, colour = Celltype)) + 
  geom_boxplot(aes(colour = Celltype), width = 0.5, fill = NA) +  # Set boxes hollow
  geom_jitter(width = 0.2, height = 0, size = 0.8, alpha = 0.4) + 
  scale_color_manual(values = color_mapping) + 
  guides(color = FALSE) +  # Remove legend
  # Remove gray grid background
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 16),  # Adjust X-axis tick label font size
    axis.text.y = element_text(size = 13),  # Adjust Y-axis tick label font size
    axis.title.x = element_text(size = 16),  # Adjust X-axis title font size
    axis.title.y = element_text(size = 16),  # Adjust Y-axis title font size
    panel.background = element_rect(fill = "white", colour = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Add reference lines for X and Y axes
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(colour = "black") 
  ) +  
  labs(title = "Phospholipids Extracellular Pool Uptake Score by Cell Type",
       x = "",
       y = "Value",
       color = "Celltype")
# Manually add p-value to the plot, ensuring it is fully visible
y_max <- max(long_data$Value)
y_offset <- diff(range(long_data$Value)) * 0.15
p <- p + 
  geom_text(aes(x = 3, y = y_max + y_offset,
                label = paste0("p-value = ", p_value)),
            vjust = 2.5, hjust = 0.5, size = 5, color = "black")  # Set p-value label font size to 12
# Save the plot
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_diff_pdf/phospholipids extracellular pool[s]_tumor_diff1.pdf", plot = p, width = 18, height = 8, dpi = 300)



# Batch processing for multiple metabolites, saving the score data of each metabolite in each cell type and plotting boxplots.
setwd("/public2022/weijinfen/flux/tam_new_analysis/scflux_mean/tumor_flux_mean")
# Load built-in dataset
data("human_gem")
# Read all samples into a list
sample_files <- list.files(pattern = "*.qs")
# Get all cell types across all samples
samples <- lapply(sample_files, function(file) qread(file))
all_cells <- unique(unlist(lapply(samples, colnames)))

# Metabolites from Exchange/demand reactions
exchange_reactions <- human_gem[human_gem$SUBSYSTEM == "Exchange/demand reactions", ]
exchange_reactions_equation <- exchange_reactions$EQUATION
# Extract metabolite names
metabolites <- gsub("\$s\$ <=>$|\$s\$$", "", exchange_reactions_equation)
# Create mapping from reaction ID to metabolite name
reaction_to_metabolite_map <- setNames(metabolites, exchange_reactions$ID)

# Define a function to read a single sample file and extract specific metabolite flux values
extract_metabolite_flux <- function(filename, metabolite) {
  data <- qread(filename)
  new_row_names <- sub("^internal_medium ", "", rownames(data))
  rownames(data) <- new_row_names
  common_ids <- intersect(rownames(data), exchange_reactions$ID)
  data <- data[common_ids, ]
  data$metabolite <- reaction_to_metabolite_map[rownames(data)]
  rownames(data) <- data$metabolite
  data <- data[, -which(colnames(data) == "metabolite")]
  metabolite_data <- data[rownames(data) == metabolite, ]
  # Ensure the dataframe contains all possible cell types
  for (cell_type in all_cells) {
    if (!(cell_type %in% colnames(metabolite_data))) {
      metabolite_data[cell_type] <- NA
    }
  }  
  # Add sample name column
  metabolite_data$Sample <- gsub("_subcell_mean_flux.qs", "", basename(filename))
  return(metabolite_data)
}

# Initialize a list to store data for each metabolite
all_metabolite_data <- list()
# List of metabolites of interest
metabolites_of_interest <- metabolites

for (metabolite in metabolites_of_interest) {
  # Extract data for the current metabolite from each sample file
  metabolite_data_list <- lapply(sample_files, function(file) extract_metabolite_flux(file, metabolite)) 
  # Merge all sample data
  all_metabolite_data[[metabolite]] <- do.call(rbind, metabolite_data_list)
  # Data preprocessing
  all_metabolite_data[[metabolite]] <- all_metabolite_data[[metabolite]] %>%
    replace(is.na(.), 0) %>%
    mutate(across(where(is.numeric), ~sign(.) * abs(.)^(1/3)))  # Apply cube root transformation
  # Set 'Sample' column as row names
  rownames(all_metabolite_data[[metabolite]]) <- all_metabolite_data[[metabolite]]$Sample
  output_dir <- "/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_incells/"
  filename <- paste0(output_dir, metabolite, "_uptake_score_tumor.qs")
  qsave(all_metabolite_data[[metabolite]], filename)
     
  # Convert data to long format
  long_data <- all_metabolite_data[[metabolite]] %>% pivot_longer(cols = -Sample, names_to = "Celltype", values_to = "Value")
  # Create random color mapping for each sample
  random_colors <- runif(n_distinct(long_data$Sample), min = 0, max = 1)
  color_mapping <- setNames(hcl(h=random_colors * 360, c=100, l=65, alpha=1), unique(long_data$Sample))
  # Calculate global Kruskal-Wallis test p-value
  kruskal_result <- kruskal.test(Value ~ Celltype, data = long_data)
  p_value <- signif(kruskal_result$p.value, digits = 3)
  # Create boxplot
  p <- ggplot(long_data, aes(x = Celltype, y = Value, colour = Sample)) +
    geom_boxplot(aes(colour = NULL), width = 0.5) +
    geom_jitter(width = 0.2, height = 0, size = 0.8, alpha = 0.4) +
    scale_color_manual(values = color_mapping) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(colour = "black")
    ) +
    labs(title = paste0(metabolite, " Uptake Score by Cell Type"),
         x = "Celltype",
         y = "Value",
         color = "Sample") 
  # Manually add p-value to the plot
  y_max <- max(long_data$Value)
  y_offset <- diff(range(long_data$Value)) * 0.15
  p <- p + 
    geom_text(aes(x = 3, y = y_max + y_offset,
                  label = paste0("p-value = ", p_value)),
              vjust = 2.5, hjust = 0.5, size = 4, color = rgb(0, 0, 0, 0.8))
  # Save the plot
  ggsave(paste0("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_diff/", metabolite, "subcell_allsample_tumor_box.png"), 
          plot = p, width = 20, height = 18, dpi = 300) 
}


# Calculate correlations for multiple metabolites between two cell types —— SPP1+APOE+TAM Microfold
# Load necessary libraries
setwd("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_incells/")
library(ggplot2)
library(data.table)
library(stringr)

# Get all .qs files in the current directory
sample_files <- list.files(pattern = "*.qs")

# Initialize two empty data.tables
spp1_apoe_tam_df <- data.table(Sample = character(), stringsAsFactors = FALSE)
microfold_df <- data.table(Sample = character(), stringsAsFactors = FALSE)

# Define a function to extract metabolite name from filename
extract_metabolite_name <- function(filename) {
  # Use str_extract from stringr package to extract content before [s]
  metabolite_name <- str_extract(filename, "^(.*)\$s\$")
  if (is.na(metabolite_name)) {
    stop(paste0("Cannot extract metabolite name from filename ", filename))
  }
  return(metabolite_name)
}

# Iterate over all files
for (file in sample_files) {
  # Read file
  metabolite_data <- qread(file) 
  # Extract metabolite name
  metabolite_name <- extract_metabolite_name(file)  
  # Extract contents of SPP1+APOE+TAM cell type
  spp1_apoe_tam_content <- metabolite_data[, c("SPP1+APOE+TAM", "Sample")]
  setnames(spp1_apoe_tam_content, c("SPP1+APOE+TAM", "Sample"), c(metabolite_name, "Sample")) 
  # Extract contents of Microfold cell type
  microfold_content <- metabolite_data[, c("Microfold", "Sample")]
  setnames(microfold_content, c("Microfold", "Sample"), c(metabolite_name, "Sample")) 
  # If first iteration, initialize data.tables
  if (nrow(spp1_apoe_tam_df) == 0) {
    spp1_apoe_tam_df <- spp1_apoe_tam_content
  } else {
    # Merge data.tables, add new column
    spp1_apoe_tam_df <- merge(spp1_apoe_tam_df, spp1_apoe_tam_content, by = "Sample", all.x = TRUE)
  }  
  if (nrow(microfold_df) == 0) {
    microfold_df <- microfold_content
  } else {
    # Merge data.tables, add new column
    microfold_df <- merge(microfold_df, microfold_content, by = "Sample", all.x = TRUE)
  }
}

# Save results
qsave(spp1_apoe_tam_df,"spp1_tam_metabolism.qs")
qsave(microfold_df, "microfold_metabolism.qs")

# Calculate correlations
correlation_df <- data.frame()
for (i in 2:ncol(spp1_apoe_tam_df)) {
  for (j in 2:ncol(microfold_df)) {
    # Calculate correlation coefficient and its p-value
    cor_test <- cor.test(spp1_apoe_tam_df[, i], microfold_df[, j], method = "spearman")   
    # Extract correlation coefficient and p-value
    correlation <- cor_test$estimate
    p_value <- cor_test$p.value    
    correlation_df <- rbind(correlation_df, data.frame(
      Metabolite_SPP1_APOE_TAM = names(spp1_apoe_tam_df)[i],
      Metabolite_Microfold = names(microfold_df)[j],
      Correlation = correlation,
      P_Value = p_value
    ))
  }
}

# View correlation results
write.csv(correlation_df,"/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_cor/spp1_microfolf_metabolism_cor.csv")

# Correlation between two metabolites across multiple samples —— Scatter plot
lactate <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_incells/HMR_9135_lactate_uptake_score_tumor.qs")
cholesterol <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_incells/HMR_9285_cholesterol_uptake_score_tumor.qs")
library(ggplot2)

# Extract lactate content in Microfold cells
microfold_lactate <- lactate[, c("Microfold","Sample")]
# Extract cholesterol content in SPP1+APOE+TAM cells
spp1_apoe_tam_cholesterol <- cholesterol[, c("SPP1+APOE+TAM","Sample")]

# Create a new data.table for plotting
dt_cor <- merge(microfold_lactate, spp1_apoe_tam_cholesterol, by = "Sample")
correlation <- cor(dt_cor$Microfold, dt_cor$SPP1+APOE+TAM, method = "spearman")
names(dt_cor) <- c("Sample", "Microfold_Lactate", "SPP1_APOE_TAM_Cholesterol")

ggscatter(dt_cor, x = "Microfold_Lactate", y = "SPP1_APOE_TAM_Cholesterol",
          size = 1,
          add = "reg.line",  # Add regression line
          add.params = list(color = "#77C034", fill = "#C5E99B", size = 1),  # Customize regression line color
          conf.int = TRUE  # Add confidence interval
) +
  stat_cor(method = "spearman", label.x = 0.2, label.y = 6, label.sep = "\n") +
  xlab("Lactate content in Microfold cells") +
  ylab("Cholesterol content in SPP1+APOE+TAM cells")



# Compare metabolite consumption or accumulation in a single cell type —— Bar plot
# Plot the top 15 accumulated or consumed metabolites for a single cell type
# Merge multiple samples by median, and process multiple samples in a loop

# Load all samples into a list
library(data.table)
sample_files <- list.files(pattern = "*.qs")
samples <- lapply(sample_files, function(file) {
  df <- qread(file)
  cbrt <- function(x) sign(x) * abs(x)^(1/3)
  cbrt_df <- apply(df, 2, cbrt)  # Apply cube root transformation
  new_row_names <- sub("^internal_medium ", "", rownames(cbrt_df))
  rownames(cbrt_df) <- new_row_names
  common_ids <- intersect(row.names(cbrt_df), nutrient_lookup_files$ID)
  # Keep only rows with common IDs
  df_common <- cbrt_df[common_ids, ]
  # Add EQUATION column
  df_common$EQUATION <- nutrient_lookup_files$EQUATION[match(common_ids, nutrient_lookup_files$ID)]
  # Remove extra characters
  df_common$EQUATION <- gsub("\$s\$ <=>", "", df_common$EQUATION)
  rownames(df_common) <- df_common$EQUATION
  df_common
})

# Get all cell types across all samples
all_cells <- unique(unlist(lapply(samples, colnames)))

# Create a data.frame to store all pathway values for all cell types and samples, initialized as NA
merged_data <- data.frame(matrix(NA, nrow=length(all_cells), ncol=length(samples)))
setnames(merged_data, paste0("Sample_", seq_along(samples)))
rownames(merged_data) <- all_cells

# Iterate over each sample, adding its values to merged_data
for (i in seq_along(samples)) {
  sample <- samples[[i]]
  common_cells <- intersect(colnames(sample), all_cells)
  for (cell in common_cells) {
    merged_data[cell, i] <- sample[cell, ]
  }
}

# Create a data.frame for storing final results, columns are all cell types, rows are all pathways
final_result <- data.frame(matrix(NA, nrow=length(rownames(samples)), ncol=length(all_cells)))
setnames(final_result, all_cells)
rownames(final_result) <- rownames(samples)

# Iterate over each pathway
for (metabolite in rownames(samples)) {
  # Iterate over each cell type
  for (cell_type in all_cells) {
    # Get values for the current cell type from all samples
    cell_values <- unlist(lapply(samples, function(sample) {
      if (cell_type %in% colnames(sample)) {
        return(sample[metabolite, ][[cell_type]])
      } else {
        return(NA)
      }
    }))   
    # Remove NA values
    cell_values <- cell_values[!is.na(cell_values)]    
    # If there are at least two values, calculate the median
    if (length(cell_values) >= 2) {
      final_result[metabolite, ][[cell_type]] <- median(cell_values)
    } else if (length(cell_values) == 1) {
      # If only one value, keep that value
      final_result[metabolite, ][[cell_type]] <- cell_values
    }
  }
}
final_result$EQUATION <- rownames(final_result)

# Output results
write.csv(final_result, file="/public2022/weijinfen/flux/cellflux/test/cell_diff_flux/major/final_merged_major_flux_scores_tumor.csv", row.names=TRUE)

### Bar plot
library(ggplot2)
library(purrr)

plot_cell_type <- function(cell_type) {
  data_full <- final_result[, c("EQUATION", cell_type)]
  max_15 <- data_full[order(data_full[, cell_type], decreasing = TRUE), ][1:15, ]
  min_15 <- data_full[order(data_full[, cell_type]), ][1:15, ]
  top_15 <- rbind(max_15, min_15)
  
  # Attempt to convert the value column to numeric, ignoring non-numeric entries
  data_long <- data.frame(
    value = as.numeric(na.omit(top_15[, cell_type])),
    rowname = top_15[, "EQUATION"][!is.na(top_15[, cell_type])]
  )
  
  # Ensure all values are numeric before ordering by absolute value of value
  data_long <- data_long[order(abs(data_long$value)), ]
  
  # Create bar plot
  p <- ggplot(data_long) +
    geom_col(aes(x = value, y = reorder(rowname, value), fill = value > 0),
             size = 0.25, color = "white") +
    geom_text(aes(x = ifelse(value > 0, -.005, .005), y = rowname, 
                  label = rowname,
                  hjust = ifelse(value > 0, 1, 0)), size = 5) +
    geom_vline(xintercept = 0, size = 1, color = "grey40") +
    scale_x_continuous(expand = expansion(add = c(0,0)),
                       breaks = seq(-.4,.2, by = .2)) + 
    scale_y_discrete(expand = c(.025,.025)) +
    scale_fill_manual(values = c("TRUE" = "#2376b7","FALSE" = "#45b787")) +
    coord_cartesian(clip = "off",
                    xlim = c(min(data_long$value) * 1.3, max(data_long$value) * 1.3),
                    ylim = c(1, nrow(data_long))) +  
    theme_minimal() + 
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          axis.text.x = element_text(face = "bold", size = rel(1), color = "black"),
          panel.border = element_rect(color = "black", linewidth = 1, linetype = "solid", fill = 'transparent'))
  p
}

purrr::map(all_cells[-length(all_cells)], function(cell_type) {
  p <- plot_cell_type(cell_type)
  ggsave(paste0("/public2022/weijinfen/flux/cellflux/test/cell_diff_flux/major/", cell_type, "_bar.pdf"), plot = p, width = 10, height = 10, units = "in", device = "pdf")
})


###
### Correlation between metabolites and cell proportions (Table + Scatter Plot)
### Calculate the correlation between metabolite levels and cell proportions
proportions_df <- read.csv("/public2022/weijinfen/flux/cellflux/test/crc8_majorcell_proportions_major.csv") # Proportions of various cells in all samples, rows are samples
df_result <- read.csv("/public2022/weijinfen/flux/cellflux/test/flux_diff_cell/major_allsample_HMR_9063.csv") # Values of a specific metabolite in all cell types across all samples, rows are samples

# Extract common samples
common_samples <- intersect(rownames(df_result), rownames(proportions_df))

# Extract Epithelial cell related data
df_result_B <- df_result[common_samples, "Epithelial"]
proportions_B <- proportions_df[common_samples, "Epithelial"]

# Calculate correlation
correlation_test <- cor.test(df_result_B, proportions_B)

# Output correlation coefficient and p-value
cat("Correlation coefficient (r): ", correlation_test$estimate, "\n")
cat("P-value: ", correlation_test$p.value, "\n")

# Plot scatter plot
library(ggpubr)
dt1 <- data.frame(df_result_B, proportions_B)
p1 <- ggscatter(dt1, x = "df_result_B", y = "proportions_B", 
                color = "#eba0b3", size=1, fill = "lightgray",
                add = "reg.line", conf.int = TRUE, 
                add.params = list(color = "black", fill = "lightgray"),
                cor.coef = TRUE,
                cor.method = "pearson")    
ggsave("/public2022/weijinfen/flux/cellflux/test/cor_figure/cor_epipro_HMR_9063_plot.pdf", plot = p1, width = 3, height = 3, units = "in") 


###
### Compare differences in metabolite levels in a specific cell type between two groups of samples
### df_result contains values for a single metabolite HMR_9063 in each sample's major cell types
### Compare differences between normal and cancer samples
# Divide samples
df_result <- read.csv("/public2022/weijinfen/flux/cellflux/test/flux_diff_cell/major_allsample_HMR_9063.csv") # Values of a specific metabolite in all cell types across all samples, rows are samples
rownames(df_result) <- df_result[,1]
df_result <- df_result[,-1]

pre_C123_samples <- subset(df_result, rownames(df_result) <= "C123")
post_C123_samples <- subset(df_result, rownames(df_result) > "C123")

# Calculate statistics
pre_stats <- apply(pre_C123_samples, 2, summary)
post_stats <- apply(post_C123_samples, 2, summary)

# Perform statistical tests
results <- list()
for (cell_type in colnames(df_result)) {
  pre_values <- pre_C123_samples[, cell_type]
  post_values <- post_C123_samples[, cell_type] 
  # Use Wilcoxon rank sum test (Mann-Whitney U test)
  test_result <- wilcox.test(pre_values, post_values)
  results[[cell_type]] <- test_result
}

# Print p-values for each cell type
for (cell_type in colnames(df_result)) {
  cat(cell_type, ": p-value = ", results[[cell_type]]$p.value, "\n")
}

# Visualize results
library(ggplot2)
library(ggpubr)

# Convert datasets to long format
pre_data <- melt(pre_C123_samples)
pre_data$Group <- "Pre-C123"
post_data <- melt(post_C123_samples)
post_data$Group <- "Post-C123"

# Combine datasets
combined_data <- rbind(pre_data, post_data)

p <- ggboxplot(
  combined_data,
  x = "Group",
  y = "value",
  color = "Group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter", # Optionally add points to show raw data distribution
  facet.by = "variable",
  ylab = "Metabolite A Value",
  xlab = "Sample Group"
) +
  stat_compare_means(method = "wilcox") + # Add p-values
  theme_minimal() +
  theme( # Customize theme
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Set background to white
    strip.background = element_rect(fill = "white", colour = NA), # Set facet label background to white
    axis.text.x = element_text(size = rel(1.2), angle = 45, hjust = 1), # Rotate X-axis text and increase its size
    axis.text.y = element_text(size = rel(1.2)), # Increase Y-axis text size
    strip.text = element_text(size = rel(1.2)), # Increase facet label text size
    axis.line = element_line(colour = "black", size = 1) # Add black axis lines with width of 1
  ) +
  facet_wrap(~variable, ncol = 6) # Adjust facet layout so that 6 comparisons are placed in one row
ggsave("/public2022/weijinfen/flux/cellflux/test/flux_diff_cell/major_box/test_majorcell_group_HMR_9063_plot.pdf", plot = p, width = 10, height = 4, units = "in")


###
### Compare differences in metabolite levels in a specific cell type between tumor and normal samples
### df_result contains values for a single metabolite HMR_9063 in each sample's major cell types
### Compare differences between normal and cancer samples
### SPP1+TAM

library(ggplot2)
library(ggpubr)
library(dplyr)

# Read data
df_normal <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_subtype_normal/fatty acid pool[s] uptake_score_normal.qs")
df_tumor <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_subtype_tumor/fatty_acid pool[s]_uptake_score_tumor.qs")

# Extract SPP1+APOE+TAM data
SPP1_APOE_TAM_col <- "SPP1+APOE+TAM"
normal_values <- df_normal[, SPP1_APOE_TAM_col, drop = FALSE]
tumor_values <- df_tumor[, SPP1_APOE_TAM_col, drop = FALSE]

# Convert data to long format
normal_data <- normal_values %>%
  as.data.frame() %>%
  rename(value = !!sym(SPP1_APOE_TAM_col)) %>%
  mutate(Group = "Normal", Sample = rownames(.))
tumor_data <- tumor_values %>%
  as.data.frame() %>%
  rename(value = !!sym(SPP1_APOE_TAM_col)) %>%
  mutate(Group = "Tumor", Sample = rownames(.))

# Combine datasets
combined_data <- bind_rows(normal_data, tumor_data) 

# Plot boxplot
p <- ggboxplot(
  combined_data,
  x = "Group",
  y = "value",
  color = "Group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter", # Optionally add points to show raw data distribution
  ylab = "fatty acid level \n in SPP1+APOE+TAM",
  xlab = NULL
) +
  stat_compare_means(method = "wilcox") + # Add p-values
  theme_minimal() +
  theme( 
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Set background to white
    strip.background = element_rect(fill = "white", colour = NA), # Set facet label background to white
    axis.text.x = element_text(size = rel(1.2), angle = 45, hjust = 1), # Rotate X-axis text and increase its size
    axis.text.y = element_text(size = rel(1.2)), # Increase Y-axis text size
    strip.text = element_text(size = rel(1.2)), # Increase facet label text size
    axis.line = element_line(colour = "black", size = 1) # Add black axis lines with width of 1
  )
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_diff_tumor_normal_pdf/fatty acid-pool_spp1_tumorvsnormal_plot.pdf", plot = p, width = 5, height = 4, units = "in")


##
## Epithelial cells
library(ggplot2)
library(ggpubr)
library(dplyr)

# Read data
df_normal <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_subtype_normal/26-hydroxycholesterol[s] uptake_score_normal.qs")
df_tumor <- qread("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_score_subtype/26-hydroxycholesterol[s] uptake_score_tumor.qs")

# Extract Microfold data
SPP1_APOE_TAM_col <- "Microfold"
normal_values <- df_normal[, SPP1_APOE_TAM_col, drop = FALSE]
tumor_values <- df_tumor[, SPP1_APOE_TAM_col, drop = FALSE]

# Convert data to long format
normal_data <- normal_values %>%
  as.data.frame() %>%
  rename(value = !!sym(SPP1_APOE_TAM_col)) %>%
  mutate(Group = "Normal", Sample = rownames(.))
tumor_data <- tumor_values %>%
  as.data.frame() %>%
  rename(value = !!sym(SPP1_APOE_TAM_col)) %>%
  mutate(Group = "Tumor", Sample = rownames(.))

# Combine datasets
combined_data <- bind_rows(normal_data, tumor_data) 

# Plot boxplot
p <- ggboxplot(
  combined_data,
  x = "Group",
  y = "value",
  color = "Group",
  palette = c("#00AFBB", "#E7B800"),
  add = "jitter", # Optionally add points to show raw data distribution
  ylab = "26-hydroxycholesterol Value",
  xlab = NULL
) +
  stat_compare_means(method = "wilcox") + # Add p-values
  theme_minimal() +
  theme( 
    panel.grid.major = element_blank(), # Remove major grid lines
    panel.grid.minor = element_blank(), # Remove minor grid lines
    panel.background = element_rect(fill = "white", colour = NA), # Set background to white
    strip.background = element_rect(fill = "white", colour = NA), # Set facet label background to white
    axis.text.x = element_text(size = rel(1.2), angle = 45, hjust = 1), # Rotate X-axis text and increase its size
    axis.text.y = element_text(size = rel(1.2)), # Increase Y-axis text size
    strip.text = element_text(size = rel(1.2)), # Increase facet label text size
    axis.line = element_line(colour = "black", size = 1) # Add black axis lines with width of 1
  )
ggsave("/public2022/weijinfen/flux/tam_new_analysis/scflux_result/metabolism_diff_tumor_normal/cholesterol_microfold_tumorvsnormal_plot.pdf", plot = p, width = 4, height = 4, units = "in")































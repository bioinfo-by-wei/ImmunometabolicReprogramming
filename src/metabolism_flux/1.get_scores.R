
library(Seurat)
library(METAFlux)
library(parallel)
options(mc.cores=10)
source('./my_functions.R')
data("human_blood")
library(qs)

# Create output directory
dir.create('scores/crc', recursive = TRUE)

# Define a function to calculate scores for parallel processing. 
# The function takes samplename and cancertype as input, calculates the corresponding scores,
# and saves the scores along with other information into a qs file.
calc_scores <- function(sample_name, cancertype){
  sc <- obj.list[[sample_name]]
  # Calculate scores using the rewritten functions
  mean_exp <- parallel_calculate_avg_exp(myseurat = sc, myident = 'revised_subtype', n_bootstrap = 100, seed = 1)
  scores <- parallel_calculate_reaction_score(data = mean_exp)
  cell_type_fractions <- table(sc$revised_subtype) / nrow(sc@meta.data)
  all_cell_types <- unique(sc$revised_subtype)
  num_cell_types <- length(all_cell_types)
  out_obj <- list()
  out_obj[['scores']] <- scores
  out_obj[['cell_type_fractions']] <- cell_type_fractions
  out_obj[['num_cell_types']] <- num_cell_types
  qsave(out_obj, file = paste0('scores/', cancertype, '/', sample_name, "_scores.qs"), nthreads = 48)
}

file <- "/public2022/weijinfen/flux/script/crc_new/crc_res_barcode.qs"
message(file, ' reading')
data <- qread(file, nthreads = 48)
message("reading finished, splitting file by sample")
obj.list <- SplitObject(data, split.by = "patientbarcode")
message('calculating scores')
patients <- names(obj.list)
cancertype <- "crc"
mclapply(patients, calc_scores, cancertype)


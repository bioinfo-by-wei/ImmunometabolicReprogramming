library(METAFlux)
library(qs)
library(Seurat)
source('./my_functions.R')
data("human_blood")

# Create output directory for storing results separately by sample due to long computation time
dir.create('scflux', recursive = TRUE)

# Calculate flux for crc8
cancer <- 'crc8'
for (file in dir(paste0('scores/', cancer))) {
    message(cancer, ':', file, '   processing')
    sample_name <- strsplit(file, '_scores.qs')[1]
    # Read data
    tmpdata <- qread(paste0('scores/', cancer, '/', file), nthreads = 48)
    num_cell_types <- tmpdata[['num_cell_types']]
    cell_type_fractions <- tmpdata[['cell_type_fractions']]
    scores <- tmpdata[['scores']]
    # Compute sc_flux
    flux <- fast_compute_sc_flux(num_cell = num_cell_types, fraction = cell_type_fractions, fluxscore = scores, medium = human_blood)
    # Save result file for each sample separately
    qsave(flux, paste0('scflux/', cancer, '/', sample_name, '_sc_flux_raw.qs'))
}


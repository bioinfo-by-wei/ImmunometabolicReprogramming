#!/bin/bash
# Activate the Conda environment
conda activate pyscenic

# Run change.py to convert CSV data into a Loom file
python change.py

# Set the directory for cisTarget databases (Change to your own directory)
dir=/data/index_genome/cisTarget_databases/

# Define paths to necessary database files
tfs=\$dir/hs_hgnc_tfs.txt
feather=\$dir/hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather
tbl=\$dir/motifs-v9-nr.hgnc-m0.001-o0.0.tbl

# Ensure the database files are complete and correct
input_loom=./sample.loom
ls \$tfs \$feather \$tbl

# 2.1 Gene Regulatory Network (GRN) Inference
pyscenic grn \
--num_workers 10 \
--output adj.sample.tsv \
--method grnboost2 \
\$input_loom \
\$tfs  # Transcription factor file, a list of 1839 gene names

# 2.2 cisTarget Analysis
pyscenic ctx \
adj.sample.tsv \$feather \
--annotations_fname \$tbl \
--expression_mtx_fname \$input_loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 20 \
--mask_dropouts

# 2.3 AUCell Analysis
pyscenic aucell \
\$input_loom \
reg.csv \
--output out_SCENIC.loom \
--num_workers 10
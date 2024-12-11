#!/bin/bash

data_names=(
    'CRCP71_T' 'CRCP95_T' 'CRCP99_T' 'CRCP100_T' 'CRCP104_T' 'CRCP107_T' 'CRCP112_T'
)

r_script_path="/wangzhenzenlab/slb/stereoseq/h5ad2rds.R"
input_dir="/wangzhenzenlab/slb/stereoseq/STT0000036"
output_dir="/wangzhenzenlab/slb/stereoseq/STT0000036"

for data_name in "${data_names[@]}"; do
    input_file="${input_dir}/${data_name}/${data_name}_bin1.h5ad"
    output_file="${output_dir}/${data_name}/${data_name}_bin1.rds"

    echo "Running R script for $data_name..."
    Rscript "$r_script_path" --infile "$input_file" --outfile "$output_file"
done
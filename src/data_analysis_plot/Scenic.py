import os
import loompy as lp
import numpy as np
import scanpy as sc

# Get the current working directory
current_dir = os.getcwd()
print(f"Current working directory: {current_dir}")

# List files in the current directory
files = os.listdir(current_dir)
print("Files in current directory:")
print(files)

# Read CSV file into an AnnData object
csv_file = "for.scenic.data.csv"
adata = sc.read_csv(csv_file)

# Prepare row attributes and column attributes for Loom file creation
row_attrs = {"Gene": np.array(adata.var_names)}
col_attrs = {"CellID": np.array(adata.obs_names)}

# Create Loom file from the AnnData object
loom_file = "sample.loom"
lp.create(loom_file, adata.X.transpose(), row_attrs, col_attrs)

print(f"Loom file '{loom_file}' created successfully.")
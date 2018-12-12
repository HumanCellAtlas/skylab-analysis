# npz2rds 
Converts Optimus output to an output suitable for reading in with R (rds file). Uses two scripts to do the conversion via txt files.


# Example Usage

npz2rds.sh -c input/sparse_counts_col_index.npy -r input/sparse_counts_row_index.npy -d input/sparse_counts.npz -o counts.rds

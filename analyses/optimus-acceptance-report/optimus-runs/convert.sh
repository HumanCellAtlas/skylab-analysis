#!/bin/bash

source ~/default_venv/bin/activate

npz2rds.sh -c ../MergeCountFiles/sparse_counts_col_index.npy \
	   -r ../MergeCountFiles/sparse_counts_row_index.npy \
	   -d ../MergeCountFiles/sparse_counts.npz \
	   -o counts.rds

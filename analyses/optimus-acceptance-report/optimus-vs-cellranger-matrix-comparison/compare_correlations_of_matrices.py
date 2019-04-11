#!/usr/bin/env python
import argparse
import os
import re
import sys
import csv
import gzip
import numpy as np
from scipy import sparse
import scipy.io
from scipy.io import mmread
import scipy.io as sp_io
import pandas as pd


def main():
    description = """The script compares numerical equivalency of the count matrices from Cell Ranger 
                  and Optimus by subsampling a prescribed number of barcodes and computing correlation.
                  The correlations are computing for pairs of row vector (one from Cell Ranger count matrix
                  and the other from the Optimus counts matrix) of gene counts across the  genes 
                  with or without random permutations (columnwise) and reordering (across barcodes) to check
                  for robustness of the correlations (which is done by downsteam codes written in R)"""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        '--optimus-output',
        dest='optimus_dir',
        required=True,
        help='Optimus output dir containing files sparse_counts.npz, sparse_counts_(col/row)_index.npy')

    parser.add_argument(
        '--cellranger-output',
        dest='cellranger_dir',
        required=True,
        help='Cell Ranger output dir containing files genes.tsv(features.tsv), barcodes.tsv and matrix.mtx')

    parser.add_argument(
        '--normalize',
        dest='normalize',
        action="store_true",
        help='normalize by gene counts across a columns by dividing  by the sum total of counts in that column')

    parser.add_argument(
        '--shuffle',
        dest='shuffle',
        action="store_true",
        help='randomly shuffle gene counts for on of the row vectors for  individual pairs')

    parser.add_argument(
        '--sub-sample-size',
        dest='sub_sample_size',
        type=int,
        default=5000,
        help='number random barcodes picked for computing the correlations')

    parser.add_argument(
        '--cross-correlation',
        dest='cross_correlation',
        action="store_true",
        help='random barcode pairs correlation')

    args = parser.parse_args()

    # read the optimus count matrix  (in sparse format files), gene names and barcodes
    opt_genes = np.load(args.optimus_dir + "/sparse_counts_col_index.npy")
    opt_barcodes = np.load(args.optimus_dir + "/sparse_counts_row_index.npy")
    opt_counts = np.load(args.optimus_dir + "/sparse_counts.npz")

    # reads the gene names
    opt_genes_dict = {}
    for i, gene in enumerate(opt_genes):
        opt_genes_dict[gene] = i

    # reads the barcode names for Optimus
    opt_barcodes_dict = {}
    for i, barcode in enumerate(opt_barcodes):
        opt_barcodes_dict[barcode] = i

    # create the csr matrix  for Optimus
    opt_csr_matrix = sparse.csr_matrix((opt_counts['data'], opt_counts['indices'],
                                        opt_counts['indptr']), shape=opt_counts['shape'])

    # print('opt matrix shape', opt_csr_matrix.shape)

    # read the cell ranger files in tsv and mtx format
    cell_genes, cell_barcodes, cell_mat = load_cellranger_mtx(
        args.cellranger_dir)

    # reads the barcode names for Cell Ranger
    cell_genes_dict = {}
    for i, gene in enumerate(cell_genes):
        cell_genes_dict[gene] = i

    # create the csr matrix  for Cell Ranger
    cell_barcodes_dict = {}
    for i, barcode in enumerate(cell_barcodes):
        cell_barcodes_dict[barcode] = i

    cell_csr_matrix = cell_mat.tocsr().transpose()
    print('cell matrix shape', cell_csr_matrix.shape)

    # compute the set of genes in Cell Ranger and Optimus
    cell_genes_set = set(cell_genes)
    opt_genes_set = set(opt_genes)

    # compute the list of common genes
    com_genes_set = cell_genes_set.intersection(opt_genes_set)
    print('common genes set   :', len(com_genes_set))

    # cmpute the set of common barcodes
    cell_barcodes_set = set(cell_barcodes)
    opt_barcodes_set = set(opt_barcodes)
    com_barcodes_set = cell_barcodes_set.intersection(opt_barcodes_set)
    print('common barcodes set:', len(com_barcodes_set))

    """ compte the submatrix from the countr matrix from Cell Ranger with the 
    set of common gene and barcodes comupted above"""
    row_list = [cell_barcodes_dict[x] for x in sorted(list(com_barcodes_set))]
    col_list = [cell_genes_dict[x] for x in sorted(list(com_genes_set))]
    cell_csr_submatrix = cell_csr_matrix[row_list, :][:, col_list]
    print('cell csr sub-matrix', cell_csr_submatrix.shape)

    """ compte the submatrix from the countr matrix from Optimus with the 
    set of common gene and barcodes comupted above"""
    row_list = [opt_barcodes_dict[x] for x in sorted(list(com_barcodes_set))]
    col_list = [opt_genes_dict[x] for x in sorted(list(com_genes_set))]
    opt_csr_submatrix = opt_csr_matrix[row_list, :][:, col_list]

    print('opt csr sub-matrix', opt_csr_submatrix.shape)
    # cell_csr_submatrix_arr = cell_csr_submatrix.toarray()
    # opt_csr_submatrix_arr = opt_csr_submatrix.toarray()

    tot = len(row_list)

    # these are normalizing divisors for the Optimus and Cell Ranger count matrix
    opt_norm = opt_csr_submatrix.sum(axis=0)
    cell_norm = cell_csr_submatrix.sum(axis=0)

    # now compute the array of inverses to accomplish the division for sparse matrix operations
    opt_norm_inv = []
    cell_norm_inv = []
    for i in range(opt_norm.shape[1]):
        if opt_norm.item(i) > 0.9:
            opt_norm_inv.append(1.0 / opt_norm.item(i))
        else:
            opt_norm_inv.append(1.0)

        if cell_norm.item(i) > 0.9:
            cell_norm_inv.append(1.0 / cell_norm.item(i))
        else:
            cell_norm_inv.append(1.0)

    # this is and random permutation across the cell barcodes or rows for sampling
    row_indices = np.random.permutation(range(len(row_list)))

    k = 0

    print("correlation:{}\t{}\t{}".format("correl", "cell_ranger", "optimus"))
    # process one row at a time 
    for i in row_indices:
        # sample size requirement has been met
        if k > args.sub_sample_size:
            break

        _row1 = np.asarray(cell_csr_submatrix[i, :].todense()).reshape(-1)
        # pick a random 2nd barcode if cross correlation
        if args.cross_correlation:
            j = np.random.choice(row_indices)
        else:
            j = i
        _row2 = np.asarray(opt_csr_submatrix[j, :].todense()).reshape(-1)

        if args.normalize:
            # normalize by the total column count or gene count
            row1 = np.multiply(
                np.asarray(_row1, dtype=np.float32),
                np.asarray(cell_norm_inv, dtype=np.float32))
            row2 = np.multiply(
                np.asarray(_row2, dtype=np.float32),
                np.asarray(opt_norm_inv, dtype=np.float32))
        else:
            # do not normalize at all 
            row1 = _row1
            row2 = _row2

        # shuffle the gene counts for the 2nd row
        if args.shuffle:
            np.random.shuffle(row2)

        try:
            cell_sum = np.sum(
                np.asarray(cell_csr_submatrix[i, :].todense()).reshape(-1))
            opt_sum = np.sum(
                np.asarray(opt_csr_submatrix[j, :].todense()).reshape(-1))
            """ if the cells are empty 0's then ignore them, these are anyway removed in the 
            R script for barcodes or random barcode pairs with less than 10 total counts"""
            if cell_sum == 0 or opt_sum == 0:
                continue

            # compute the correlations
            correl = np.corrcoef(row1, row2)[0][1]
            k += 1
            print("correlation:{}\t{}\t{}".format(correl, cell_sum, opt_sum))
        except:
            pass


def print_values(com_barcodes_set, com_genes_set, cell_barcodes_dict, cell_genes_dict, cell_csr_matrix,
                 opt_barcodes_dict, opt_genes_dict, opt_csr_matrix):
    fout = sys.stdout
    b = 0
    for i, barcode in enumerate(com_barcodes_set):
        # print(i)
        rowcount = cell_csr_matrix[
                   cell_barcodes_dict[barcode], :].count_nonzero()
        if rowcount <= 52:
            continue
        # print('rowcount 10', rowcount)

        # skip zero read barcodes
        fout.write("{}\t{}\t".format(i, barcode))

        # collect a set of genes with non-zero counts
        genes_disp = []
        k = 0
        for j, gene in enumerate(com_genes_set):
            if cell_csr_matrix[cell_barcodes_dict[barcode],
                               cell_genes_dict[gene]] > 0:
                genes_disp.append(gene)
                k += 1
            if k > 35:
                break
        fout.write('\t'.join(genes_disp) + "\n")

        fout.write("\t\t\t")
        for gene in genes_disp:
            fout.write("\t{}".format(
                cell_csr_matrix[cell_barcodes_dict[barcode],
                                cell_genes_dict[gene]]))
        fout.write("\n")

        fout.write("\t\t\t")
        for gene in genes_disp:
            fout.write("\t{}".format(opt_csr_matrix[opt_barcodes_dict[barcode],
                                                    opt_genes_dict[gene]]))
        fout.write("\n")

        b += 1
        if b > 2:
            break


def load_cellranger_mtx(matrix_dir):
    mat = scipy.io.mmread(os.path.join(matrix_dir, "matrix.mtx"))

    features_path = os.path.join(matrix_dir, "genes.tsv")
    feature_ids = [
        row[0] for row in csv.reader(open(features_path), delimiter="\t")
    ]
    gene_names = [
        row[1] for row in csv.reader(open(features_path), delimiter="\t")
    ]

    barcodes_path = os.path.join(matrix_dir, "barcodes.tsv")
    barcodes = [
        re.sub('-\d', '', row[0])
        for row in csv.reader(open(barcodes_path), delimiter="\t")
    ]
    return gene_names, barcodes, mat


if __name__ == "__main__":
    main()

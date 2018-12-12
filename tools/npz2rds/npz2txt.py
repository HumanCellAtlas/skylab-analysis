#!/usr/bin/env python

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Description: A python script that will read the output count matrix from
##   the Optimus pipeline and save it as discrete text files in the selected output
##   directory

from optparse import OptionParser

parser = OptionParser()
parser.add_option('--col-index',dest='colindex', help='column index npy file')
parser.add_option('--row-index',dest='rowindex', help='row index npy file')
parser.add_option('--counts',dest='counts', help='counts npz file')
parser.add_option('--output-dir',dest='output',help='output directory')

(options, args) = parser.parse_args()

## Check that the options are valid
import numpy as np

## Column indices
col_indices = np.load(options.colindex)
np.savetxt(options.output + '/sparse_counts_col_index.txt', col_indices, fmt='%s')

## Row indices
row_indices = np.load(options.rowindex)
np.savetxt(options.output + '/sparse_counts_row_index.txt', row_indices, fmt='%s')

## The count data
data = np.load(options.counts)

## Save each component as individual txt file
np.savetxt(options.output + '/sparse_counts_indices.txt',data['indices'],fmt='%i')
np.savetxt(options.output + '/sparse_counts_indptr.txt',data['indptr'],fmt='%i')
np.savetxt(options.output + '/sparse_counts_shape.txt',data['shape'],fmt='%i')
np.savetxt(options.output + '/sparse_counts_data.txt',data['data'],fmt='%i')

### End of Script

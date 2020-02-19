#!/usr/bin/env python3

import os
import xlrd
import csv
import argparse

parser = argparse.ArgumentParser(description="Convert an cumulus differential expression Excel sheet into a csv files")
parser.add_argument('--input-xlsx-path', dest='input_xslx_path', help="input xlsx files")
parser.add_argument('--output-prefix', dest='output_prefix', help="output prefix")
args = parser.parse_args()

input_xlsx_path = args.input_xslx_path
output_prefix = args.output_prefix

print("Opening workbook " + input_xlsx_path)
open_book = xlrd.open_workbook(filename=input_xlsx_path, on_demand=True)

## Loop over all sheets
for cur_sheet_name in open_book.sheet_names():
    print("Loading sheet " + cur_sheet_name)
    curr_sheet = open_book.sheet_by_name(cur_sheet_name)
    output_file_name = output_prefix + cur_sheet_name + '.csv'
    print("Opening file " + output_file_name)
    f = open(output_file_name, 'w')
    c = csv.writer(f)
    for r in range(curr_sheet.nrows):
        c.writerow(curr_sheet.row_values(r))
    f.close()

print("Done")

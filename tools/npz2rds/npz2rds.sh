#!/bin/bash

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Date: December 12, 2018
## Description: A wrapper script that converts optimus
##   output to RDS R files


show_help() {
    echo "Usage: $0 [arguments]"
    echo ""
    echo "Arguments:"
    echo "  -c\t\tcolumn name file (.npy)"
    echo "  -r\t\trow name file (.npy)"
    echo "  -d data file (.npz)"
    echo "  -o output file (.rds)"
    echo "  -h print this helpful message"
    echo ""
}


## Process command line arguments
OPTIND=1

## Init variables
colindexfile=""
rowindexfile=""
countsfile=""
outputfile=""

while getopts "hc:r:d:o:" opt; do
    case "$opt" in
	h)
	    show_help
	    exit 0
	    ;;
	c)
	    colindexfile=$OPTARG
	    ;;
	r)
	    rowindexfile=$OPTARG
	    ;;
	d)
	    countsfile=$OPTARG
	    ;;
	o)
	    outputfile=$OPTARG
	    ;;
    esac
done

shift $((OPTIND-1))

## DEBUG check args
# echo colindex $colindexfile
# echo rowindex $rowindexfile
# echo counts $countsfile
# echo outputfile $outputfile


## Make tmp directory
tmpdir=`mktemp -d`

## Convert the npz to text
./npz2txt.py --col-index $colindexfile \
	     --row-index $rowindexfile \
	     --counts $countsfile \
	     --output-dir $tmpdir

## Convert the text to rds
./sparseTxt2Rds.R --input-dir $tmpdir --output-file $outputfile


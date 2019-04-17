import pysam
import sys
import argparse

def read_whitelist(whitelist_file):
    """This function reads the white list file and returns the list"""

    # create an empty dictionary to store the white lists
    whitelistdict = {}
    with open(whitelist_file, 'r') as fp:
        for line in fp:
            whitelistdict[line.strip()] = True

    # return the list of whitelist
    return whitelistdict.keys()


def barcode_counting():
    description = """This script computes statistics for cell barcode correction, against a whitelist, 
                   for the reads in a BAM file from Optimus after the Attach10Barcodes task. 
                   The script also takes a cutoff parameter to count barcodes with 
                   this minimum number of reads in the sample"""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--bam',
                        dest="input_bam",
                        required=True,
                        help='any of the BAM files following the Attach10XBarcodes task in an Optimums run')

    parser.add_argument('--whitelist',
                        dest="whitelist",
                        required=True,
                        help='the whitelist containing a set of barcodes from Cell Ranger')

    parser.add_argument('--cutoff',
                        dest="cutoff",
                        required=True,
                        help='minium number of reads to consider any corrected barcode to be not spurious')

    args = parser.parse_args()

    whitelist = read_whitelist(args.whitelist)

    barcodes = {}
    corrected_barcodes = {}

    samfile = pysam.AlignmentFile(args.input_bam, "rb")

    tot = 0  # total reads
    i = 0  # reads with 0 error
    j = 0  # reads with 1 error
    k = 0  # reads with 2 or more
    l = 0

    # retrieve reads one by one
    for read in samfile:
        tot += 1
        # if there is a cell barcode
        if read.has_tag('CB'):
            cbtag = read.get_tag('CB')

            if not cbtag in barcodes:
                barcodes[cbtag] = 0

            barcodes[cbtag] += 1

            crtag = read.get_tag('CR')
            # if raw cell barcode is already in the white list
            if crtag in whitelist:
                i = i + 1
            else:
                if not crtag in corrected_barcodes:
                    corrected_barcodes[crtag] = [0, cbtag]
                corrected_barcodes[crtag][0] += 1
                j = j + 1
        else:
            # print('d', mindistance(read.get_tag('CR')))
            # ed = mindistance(read.get_tag('CR'))
            # errors[ed] +=1
            k = k + 1

        if tot % 3000000 == 0:
            print(tot, i, j, k)
            print_statistics(barcodes, corrected_barcodes)

    print("{}\t{}\t{}\t{}".format(tot, i, j, k))

    print_statistics(barcodes, corrected_barcodes, args.cutoff)


def print_statistics(barcodes, corrected_barcodes, cutoff):
    """ This script computes and reports the following: 
       (a) total number of barcodes, 
       (b) number of correcte barcodes,
       (c) number of barcodes with at least 'cutoff' reads ,
       (d) number of reads in corrected barcodes (with min cutoff)
    """
    barcode_counts = sorted(barcodes.values())
    num_barcodes = len(barcode_counts)

    corr_to_bar_with_min_count = 0
    reads_in_min_count_corrected_barcodes = 0
    for barcode in corrected_barcodes:
        # if corrected barcode has been seen at least 'cutoff times'
        if barcodes[corrected_barcodes[barcode][1]] >= cutoff:
            corr_to_bar_with_min_count += 1
            reads_in_min_count_corrected_barcodes += corrected_barcodes[barcode][0]

    print("Cutoff: {}".format(cutoff))
    print(
        "Total barcodes: {}    Corrected barcodes : {}  Corrected barcodes: {} with min #reads {}  Reads in corrected barcodes: {}".format(
            num_barcodes, len(corrected_barcodes), corr_to_bar_with_min_count, cutoff,
            reads_in_min_count_corrected_barcodes))


if __name__ == "__main__":
    barcode_counting()

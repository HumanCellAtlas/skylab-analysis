#!/usr/bin/env bash
# This script:
# Takes an Optimus multi-aligned bam file as input
# Extracts aligned reads from chr1, multialigned reads from chr2, and unaligned reads from the input bam
# Reverts the extracted data to fastq and exits.
# These inputs are then run through Optimus and Cell Ranger as part of a pending scientific testing framework.

set -eo pipefail

# grab the picard jar
curl -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.18.4/picard.jar

# assumes you've grabbed samtools with one of these commands
# suo apt-get install samtools
# sudo yum install samtools

INPUT_BAM=$1  # bam file input from Optimus, HCA uses the pbmc8k dataset found here:
# http://s3-us-west-2.amazonaws.com/10x.files/samples/cell-exp/2.1.0/pbmc8k/pbmc8k_fastqs.tar
N_UNIQUE=${2:-1000000}
N_UNALIGNED=${3:-100000}
N_MULTIALIGNED=${4:-100000}

# sort and index the bam file
samtools sort -o sorted.bam ${INPUT_BAM}
samtools index sorted.bam

# get some unique-aligned reads
samtools view -H sorted.bam > header.txt
samtools view sorted.bam 1 | head -n ${N_UNIQUE} > unique.sam

# get some unaligned reads
samtools view -f 4 sorted.bam | head -n "${N_UNALIGNED}" > unaligned.sam

# get some multi-aligned reads
# the awk command here retrieves records that were aligned between 2-9 times, but then extracts
# only the primary records (HI:i:1) because we're converting this back to fastq, so we don't
# need the alternative alignment locations.
samtools view sorted.bam 2 |\
  awk '/NH:i:[2-9]/&&/HI:i:1/' |\
  head -n "${N_MULTIALIGNED}" > multialigned.sam

# put things back together
cat header.txt unique.sam unaligned.sam multialigned.sam | samtools view -b > subset.bam
samtools sort -o subset_sorted.bam subset.bam
samtools index subset_sorted.bam

# revert to fastq
java -Xmx4g -jar /usr/picard/picard.jar SamToFastqWithTags \
  I="${INPUT_BAM}" \
  FASTQ="r2.fastq.gz" \
  SEQUENCE_TAG_GROUP=CR,UR \
  QUALITY_TAG_GROUP=CY,UY \
  COMPRESS_OUTPUTS_PER_TAG_GROUP=true \
  INCLUDE_NON_PRIMARY_ALIGNMENTS=false
mv CR_UR.fastq.gz r1.fastq.gz

# all done
echo "finished subsetting $2 unique, $3 unaligned and $4 multi-aligned reads from $1"


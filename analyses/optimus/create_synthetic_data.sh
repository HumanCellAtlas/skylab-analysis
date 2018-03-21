#!/usr/bin/env bash

# this script takes as input the PBMC8k dataset, multi-aligned with Optimus.
# it generates a bam file containing unique, unaligned, and multialigned reads from different
# chromosomes
# it then regenerates fastq data

# grab the picard jar
curl -o picard.jar https://github.com/broadinstitute/picard/releases/download/2.18.4/picard.jar

# assumes you've grabbed samtools with one of these commands
# suo apt-get install samtools
# sudo yum install samtools

INPUT_BAM=$1  # pbmc8k resequenced output
N_UNIQUE=$2  # 1,000,000
N_UNALIGNED=$3  # 100,000
N_MULTIALIGNED=$4  # 100,000

# sort and index the bam file
samtools sort -o sorted.bam ${INPUT_BAM}
samtools index sorted.bam

# get some unique-aligned reads
samtools view -H sorted.bam > header.txt
samtools view sorted.bam 1 | head -n ${N_UNIQUE} > unique.sam

# get some unaligned reads
samtools view -f 4 sorted.bam | head -n "${N_UNALIGNED}" > unaligned.sam

# get some multi-aligned reads
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
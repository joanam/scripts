#!/bin/bash

# Map the reads against the PhiX genome creating:
# - a sam file (<prefix>_phix.sam) with all PhiX reads aligned
# - a fastq file (<prefix>.noPhiXreads.fastq) with all non-PhiX reads
# written by Joana Meier (Nov 2016)

# Note, this line needs to be adjusted to show the correct path and prefix of the bowtie2 index of the PhiX genome
phix_reference="/cluster/project/gdc/shared/p129/ref-genome/phix174"

# This script requires that bowtie2 is installed, exectable the command 'bowtie2'. 
# Else, the script must be modified with the correct path and/or name of bowtie 2.

# Check if two arguments were given
if [[ $1 = "-h" ]]; then
  echo "Usage: PhiX_removal.sh <raw reads fastq file> <outfile prefix>"
  exit 0
fi

if (( "$#" != 2 )); then
  echo "Please provide 2 arguments: <raw reads fastq file> <outfile prefix>"
  echo "Usage: PhiX_removal.sh  <raw reads fastq file> <outfile prefix>"
  exit 0
fi


file=$1
prefix=$2
pu=`awk -F: '{print $1":"$2}' $file | head -1 | sed 's/@//'`

# Run bowtie2 to write out the PhiX reads into a sam file and generate a fastq file without PhiX reads
bowtie2 -x ${phix_reference} \
 -q $file --phred33 --end-to-end \
 -p 4 -N 1 --no-unal \
 --un $prefix".noPhiXreads.fastq" \
 -S $prefix"_phix.sam" \
 --rg-id $prefix"_phix" \
 --rg "ID:"$prefix"_phix" \
 --rg "PG:Bowtie2_end-to-end_default_settings" \
 --rg "PL:ILLUMINA" \
 --rg "PU:"$pu \
 --rg "SM:"$prefix"_phix"


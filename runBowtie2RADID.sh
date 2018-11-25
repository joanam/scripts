#! /bin/bash
# Usage: bsub -n 4 "runBowtie2RADID.sh <path_to_bowtie2_refgenome_index>"

# This script runs bowtie2 on all files ending on .fastq in the current directory
# This script assumes the usage of the standard naming scheme of the aquatic ecology and evolution group at EAWAG. 
# Users utilizing a different naming scheme must adjust the script

# The alignment settings are default end-to-end alignment and
# the .bam file contains only aligned sequences (option --no-unal)
# Check bowtie2 manual for the explanation of the arguments and for parameter optimization
# Check bowtie.log files after running the script
# created by David Marques and Joana Meier, May 2014
# Note runBowtie2RADID.sh now directly generates sorted and indexed bam files, samToBamSortIndex.sh is deprecated!

# Load the modules (adjust depending on the cluster/server setup)
module load bowtie2 samtools

# This clause checks if there was a path to the bowtie2 index file given
if [ $# -lt 1 ]
then
	echo "Error: Please provide path to bowtie2 index like /cluster/project/gdc/ref-genome/someGenome (no suffix)"
	exit 1
fi

# This clause double-checks if the bowtie2-index is available where specified
if [ -f $1.1.bt2 ]
then
	echo ""
else
	echo "Error: bowtie2-index not found, please add the correct file/path combination (or remove the file suffixes if given before)"
	exit 1
fi

# Loops over all files in the current directory with the ending GQI*.fastq
for i in *GQI*.fastq
	do
	echo -e "\nbowtie2 aligning: "$i

	# These variables capture different informations from the RADID string
	# example: "12345.gasacuCHL1.GQI14.fastq"
	# f captures "12345.gasacuCHL1" and l captures "GQI14"
	f=$(echo $i | cut -f 1-2 -d ".")
	l=$(echo $i | cut -f 3 -d ".")
	ind=${i%.fastq}

	# The variable m captures the name of the sequencing machine
	# and the run-number (important for base-quality score recalibration
	m=$(head -1 $i | cut -f 1-2 -d ":" | cut -b 2-)

	bowtie2 \
	  -x $1 \
	  -q $i \
	  --phred33 \
	  --end-to-end \
	  --no-unal \
	  -p 4 \
	  -N 1 \
	  --rg-id $f \
	  --rg "ID:"$f \
	  --rg "LB:"$l \
	  --rg "PG:bowtie2" \
	  --rg "PL:ILLUMINA" \
	  --rg "PU:"$m \
	  --rg "SM:"$f \
	  2> $ind.bowtie2.log | samtools view -bS - > $ind.bowtie2.bam
	  
	  # sorting and indexing
	  samtools sort -o $ind.bowtie2.bam -O bam -T $ind.tmp $ind.bowtie2.bam
	  samtools index $ind.bowtie2.bam
done

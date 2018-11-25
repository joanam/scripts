#!/bin/bash

# This bash script allows to demultiplex RAD data with barcodes of different lengths
# It uses process_radtags from Stacks, starting with the longest barcodes
# It creates for each barcode length in the output directory:
# - a log file that is called e.g. process_radtags_6bp.log
# - the demultiplexed fastq files (e.g. sample_AGTGC.fq)
# - an output files with reads that did not contain a barcode in the provided list (*.discards)
# arguments: readsFastqFile(can be gzipped) barcodesFile outputFolder trimmingLength
# The barcodes file should be in the following format, one barcode per line, tab-delimited: <individual name> "tab" <barcode>

# created by Joana Meier, July 2014


# Load stacks (needs to be adjusted according to the cluster/server setup):
module load gcc/4.8.2 gdc perl/5.18.4 stacks/1.40

# This script assumes, that process_radtags is executable from the current directory. If that is not the case, specify the path here:
path=""


#*******************************************
#   Check correct usage of the bash Script
#*******************************************

# Check if help is requested
if [[ $1 = "-h" ]]; then
  echo -e "\nUsage: demultiplexing.sh readsFastqFile barcodesFile outputFolder trimmingLength"
  echo -e "\nThe fastq file can be gzipped or not"
  echo -e "\nNote: the barcodes file should contain per line a fishecRADid than a tab, and a barcode of 5-8 bp length"
  echo -e "\nNote: If you set the trimming value to a value exceeding the read length, it will discard all reads as low quality" 
  exit 0
fi

# Check if the right number of arguments is provided
if (( "$#" != 4 )); then
  echo -e "\nPlease provide 4 arguments: reads.fastq barcodes.file path/to/outputFiles trimmingLength"
  exit 0
fi

# Check if the fastq file exists and contains data
if [[ ! -s $1 ]]
then
        echo -e "\n"$1" not found or empty. Please provide file with reads to demultiplex"
        exit 1
fi

# Check if the barcodes file exists and contains data
if [[ ! -s $2 ]]
then
	echo -e "\n"$2" not found or empty. Please provide file with barcodes"
	exit 1
fi

# Check if the given output path exists
if [[ ! -d $3 ]]
then
	echo -e "\noutput path "$3" does not exist! Please provide a valid path"
	exit 1
fi

# Check that the barcodes file contains barcodes in the second column
if [[ $(awk '/[[:graph:]]/ {if(!match($2,/[ATGC][ATGC][ATGC][ATGC]+/)) print $2}' $2 | wc -l )  -gt 0 ]]
then
	echo -e "\n$(awk '/[[:graph:]]/ {if(!match($2,/[ATGC][ATGC][ATGC][ATGC]+/)) print $2}' $2 | head -1) is not a barcode!\n" 
	echo "Please provide a barcode file in the following format: fishecRADid <tab> barcode" 
	exit 1
fi


#************************************
# Check if the input file is zipped
#************************************

if [[ $1 == *"gz" ]]
then
	format=gzfastq
else
	format=fastq
fi


#**************************
# Read in arguments
#**************************

readsFile=$1
barcodes=$2
out=$3
# if present, remove slash at the end of the path name
out=`echo $out | sed -e 's-/$--'`
readsPath="."
trim=$4


#*******************************************************
# in case that there are DOS line endings, remove them
#*******************************************************

dos2unix -q -k $barcodes



#************************************************************
# demultiplex the reads starting with the longest barcodes
#************************************************************

for i in 8 7 6 5
do
 # print all barcodes of length $i into a temporary barcodes file
 awk -v len=$i '{if(length($2)==len) print $2}' $barcodes > $barcodes"_"$i".tmp"

 # check if the temporary barcodes file is empty
 if [[ ! -s $barcodes"_"$i".tmp" ]]
   then 
	rm $barcodes"_"$i".tmp"
	echo -e "\nno barcodes of "$i"bp length"

   # if the temporary barcodes file contains barcodes: extract reads with these barcodes
   else 
	echo -e "\ndemultiplexing reads for "$i"bp long barcodes"
	
	# demultiplex the reads with the temporary barcodes file
${path}process_radtags \
         -f $readsPath"/"$readsFile \
         -o $out \
         -t $trim \
         -b $barcodes"_"$i".tmp" \
         -i $format -r -e "sbfI" -E 'phred33' -D
	
	# rename the process_radtags log file to prevent overwriting	
	mv $out/process_radtags.log $out/process_radtags_$i"bp.log"

	# remove temporary barcodes file
	rm $barcodes"_"$i".tmp"

	# change the readsFileName and path
	readsFile=$readsFile".discards"
	readsPath=$out
	format="fastq"

 fi
done


#************************************************************
# Suggest next steps to the user
#********************************************************

echo -e "\n\n finished!"
echo "Check the process_radtags log files for additional barcodes"
echo "If everything looks ok, delete the .discards file"
echo "To rename the files and recover the header information, use: header_recoverfromstacks.sh on the Github page of David Marques"

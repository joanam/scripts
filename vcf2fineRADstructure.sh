#!/bin/bash

# Shell script to convert the vcf file to an input format for fineRADstructure
# Written by Joana Meier, 2018
# Sorry, this is a very slow script as it extracts each RAD locus from the vcf file with vcftools

# usage: vcf2fineRADstructure.sh <vcf file> <path to bam files>

# This script requires the createRADmappingReport.sh script which requires the folder with the bam files
# It gets the RADloci present in at least 10 individuals with at least 10 reads each
# It runs files of 1000 lines in parallel (increase number of lines if it uses too many CPU).


# get the vcf file name (without suffix)
file=$1
file=${file%.gz}
file=${file%.vcf}
bamfilesFolder=$2
bamfilesFolder=${bamfilesFolder%/}

# Extract the RADtags present in at least 10 inds (positions of RAD sites from mapping report)
currentDir=`pwd`

# Go to the directory containing the bam files, run createMappingReport; come back
cd $bamfilesFolder; createRADmappingReport.sh; cd $currentDir

# Extract RADloci with at least 10 individuals from seq_depth_min10.txt file
cut -d" " -f 3,4 $bamfilesFolder/seq_depth_min10.txt | sort | uniq -c > RADpos.c
awk '{if($1>10) print $2,$3}' RADpos.c > RADpos
sort -V RADpos > RADpos.sorted
awk '{print $1" "$2-100" "$2"\n"$1" "$2" "$2+100}' RADpos.sorted | \
  grep -v "locus" > RADloci

# Get the number of individuals
if [ -s $file.vcf.gz ]
then
  nind=`zgrep ^#CH ${file}.vcf.gz | awk '{print NF-9}'`
  suff=".gz" # suffix
  vcfgz="gz" # for vcftools
else if [ -s $file.vcf ]
  nind=`grep ^#CH ${file}.vcf | awk '{print NF-9}'`
  suff="";vcfgz=""
else
  echo -e "Error: file $file.vcf[.gz] not found!\nexiting..."
  exit 1
fi

# Split the RADloci file into files of 10000 lines for parallelisation
split -l 10000 RADloci RADloci.

# Function to generate the $RADlocifile.file
function generateFile {
 a=0
 while read i
 do
   # increment the index
   a=$((a+1))

   # get the SNPs of the RADlocus
   vcftools --${vcfgz}vcf $file.vcf$suff --chr `echo $i | cut -d" " -f1` --plink-tped \
          --from-bp `echo $i | cut -d" " -f2` --to-bp `echo $i | cut -d" " -f3` \
          --out ${RADlocifile}.$a

   # if the RADlocus contains SNPs: append the info in the correct format to the ${RADlocifile}.file
   if [[ -s ${RADlocifile}.$a.tped ]]
   then
   echo -e "locus"$a"\t"`head -1 ${RADlocifile}.$a.tped | cut -f 1,2`"\t" | tr -d '\n' >> ${RADlocifile}.file
   for((c=5;c<($nind*2+5);c+=2))
   do
    echo -e `cut -f $c ${RADlocifile}.$a.tped | \
      tr -d '\n'`"/"`cut -f $((c+1)) ${RADlocifile}.$a.tped | \
      tr -d '\n'`"\t" |  tr -d '\n' >> ${RADlocifile}.file
   done
   echo "" >> ${RADlocifile}.file
   fi

   # remove temporary files
   rm ${RADlocifile}.$a.log ${RADlocifile}.$a.tfam ${RADlocifile}.$a.tped
 
 # Read in the RADlocifile given as argument when running the function
 done < $1

}


# run each file separately (in parallel) with the function above
for RADlocifile in RADloci.*
do
  generateFile $RADlocifile &
done


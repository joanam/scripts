#!/bin/bash

# Script to convert vcf to eigenstrat format for ADMIXTOOLS
# Written by Joana Meier
# It takes a single argument: the vcf file (can be gzipped)

if [ $# -ne 1 ]
 then
 echo "Please provide the vcf file to parse"
 exit 1
fi

file=$1
file=${file%.gz}


# if the vcf file is gzipped:

if [ -s $file.gz ]
then

 file=${1%.vcf}

 # Get a .map and .ped file
 vcftools15 --gzvcf $file".vcf.gz" \
         --plink --mac 1.0 \
         --out $file

else
 vcftools15 --vcf $file".vcf" \
         --plink --mac 1.0 \
         --out $file

fi


# Change the .map file to match the requirements of ADMIXTOOLS
awk -F"\t" 'BEGIN{scaff="";add=0}{
        split($2,newScaff,":")
        if(!match(newScaff[1],scaff)){
                scaff=newScaff[1]
                add=lastPos
        }
        count+=0.0001
        pos=add+$4
        print "1\t"$2"\t"count"\t"pos
        lastPos=pos
}' ${file}.map > better.map
mv better.map ${file}.map

# Change the .ped file to match the ADMIXTOOLS requirements
awk 'BEGIN{ind=1}{printf ind"\t"$2"\t0\t0\t0\t1\t"; for(i=7;i<=NF;++i) printf $i"\t";ind++;printf "\n"}' ${file}.ped  > tmp.ped
mv tmp.ped ${file}.ped


# create an inputfile for convertf
echo "genotypename:    ${file}.ped" > par.PED.EIGENSTRAT.${file}
echo "snpname:         ${file}.map" >> par.PED.EIGENSTRAT.${file}
echo "indivname:       ${file}.ped" >> par.PED.EIGENSTRAT.${file}
echo "outputformat:    EIGENSTRAT" >> par.PED.EIGENSTRAT.${file}
echo "genotypeoutname: ${file}.eigenstratgeno" >> par.PED.EIGENSTRAT.${file}
echo "snpoutname:      ${file}.snp" >> par.PED.EIGENSTRAT.${file}
echo "indivoutname:    ${file}.ind" >> par.PED.EIGENSTRAT.${file}
echo "familynames:     NO" >> par.PED.EIGENSTRAT.${file}


# Use CONVERTF to parse PED to eigenstrat
convertf -p par.PED.EIGENSTRAT.${file}

# change the snp file for ADMIXTOOLS:
awk 'BEGIN{i=0}{i=i+1; print $1"\t"$2"\t"$3"\t"i"\t"$5"\t"$6}' $file.snp > $file.snp.tmp
mv $file.snp.tmp $file.snp

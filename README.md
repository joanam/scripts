# scripts
Scripts to handle NGS data and other biological data

# HaploABBABABA_multithreaded # Attention: I am working on a bug fix, please do not currently use this software
Java program to calculate ABBA-BABA (D-statistics) to infer gene flow and the five-population test to infer the direction of gene flow. 

For more info, use: java -jar HaploABBABABA_multithreaded.jar -h

Usage: java -jar HaploABBABABA_multithreaded.jar \
         -i <input.vcf> \
         -p <populations file> -c <file with combinations to test> \
         -o <output name> -t <n threads> -b <n bootstraps>


# ldPruning.sh
Bash script to prune SNPs in high linkage disequilibrium from a vcf file. This is important for all down-stream analyses that assume no linkage among SNPs, i.e. that require independent SNPs. 

Usage: ldPruning.sh <vcffile[.gz]> [optional: <LD threshold (R^2), default 0.1> <outputformat vcf/plink> <if 0/1 recoding requested, write 01>]

Requires vcftools and plink


# vcf2phylip.py
Python2 script to convert a vcf file to phylip format, e.g. for RAxML

Usage: vcf2phylip.py -i <input.vcf> -o <output.py> [optional: -r -f -e]
        if -r is specified, the reference sequence will be included in the phylip file
        if -f is specified, all sites not in the vcf file will be printed as missing (N)
        if -e is specified, indels are not printed (else replaced by N)

If -f is specified, there will be no frameshifts if a site in the vcf file was filtered out due to low quality and also indels will be handled so that the positions and the length of the sequence is the same as the reference sequence length. This is very useful for data partitioning (e.g. specifying first, second and third codons) without having to adjust the positions because of sites lost during filtering steps or due to indels.



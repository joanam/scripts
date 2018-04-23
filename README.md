# scripts
Java, R, Python and Bash scripts to handle NGS data and other biological data


## ldPruning.sh
Bash script to prune SNPs in high linkage disequilibrium from a vcf file. This is important for all down-stream analyses that assume no linkage among SNPs, i.e. that require independent SNPs, e.g. STRUCTURE, PCA, treemix...

Requires vcftools and plink

```
Usage: ldPruning.sh <vcffile[.gz]> [optional: <LD threshold (R^2), default 0.1> <outputformat vcf/plink> <if 0/1 recoding requested, write 01>]
```


## vcf2phylip.py
Python2 script to convert a vcf(.gz) file to phylip format, e.g. for RAxML

```
Usage: vcf2phylip.py -i <input.vcf(.gz)> -o <output.py> [optional: -r -f -e -m]
        if -r is specified, the reference sequence will be included in the phylip file
        if -f is specified, all sites not in the vcf file will be printed as missing (N)
        if -e is specified, indels are not printed (else replaced by N)
        if -m is specified, haploid genotypes as e.g. for mitochondrial DNA are expected 
           (use e.g. when genotype calling was performed with a GATK tool and the option --sample_ploidy 1

If -f is specified, there will be no frameshifts if a site in the vcf file was filtered out due to low 
quality and also indels will be handled so that the positions and the length of the sequence is the 
same as the reference sequence length. This is very useful for data partitioning (e.g. specifying first, 
second and third codons) without having to adjust the positions because of sites lost during filtering 
steps or due to indels.
```


## convertVCFtoEigenstrat.sh
Bash script which converts a vcf file to eigenstrat format for ADMIXTOOLS.

Requires vcftools and converf

```
Usage: convertVCFtoEigenstrat.sh <vcffile>
```


## allelicBalance.py
Python script that handles genotypes with strong allelic disbalance, i.e. heterozygote genotypes with one allele covered by most reads and another allele covered by much fewer reads. This is usually an indication of contamination but it can also be present if some reads are PCR duplicates. The user can specify if disbalanced heterozygote genotypes should be set as missing data or turned into homozygote genotypes of the allele with more reads.

```
usage: allelicBalance.py [-h] -i I -o O (-hom | -excl | -two) [-p P] [-p2 P2]
                         [-r R]

optional arguments:
           -h, --help           show this help message and exit
           -i I, --input I      input file in vcf format [required]
           -o O, --output O     output file [required]
           -hom, --homozygote   set failing genotypes as homozygous
           -excl, --exclude     set failing genotypes as missing
           -two, --twoSteps     set failing genotypes as missing if pvalue<p and as homozygous if pvalue<p2
           -p P, --pvalue P     p-value threshold for binomial test [default: 0.01]
           -p2 P2, --pHomoz P2  second p-value threshold for binomial test if -two is specified, 
                                threshold to make failing genotypes homozygote [default=0.005]
           -r R, --ratio R      hard cutoff for allelic ratio [default=0.2]
```
 

## addRecombRates.r

R script to generate a map in plink format, e.g. for phasing with BEAGLE. It requires a vcf file and a linkage map as input.
```
Usage:  addRecombRates.r <vcffile> <linkage map>
```


## vcf2fineSTR.lsf

Script that can be submitted to the Euler cluster (ETHZ) generating the fineSTRUCTURE input files from a vcf file and a linkage map. This script is not intended to be useful for people without Euler access but rather as an example of how one could generate a fineSTRUCTURE file from a vcf file using vcftools and the makeRecombMap4fineSTRUCTURE.sh and convertIndsFineSTR.sh scripts.

```
Usage:  bsub < vcf2fineSTR.lsf
```

 
## HaploABBABABA_multithreaded # Attention: I am working on a bug fix, please do not currently use this software
Java program to calculate ABBA-BABA (D-statistics) to infer gene flow and the five-population test to infer the direction of gene flow. 

For more info, use: java -jar HaploABBABABA_multithreaded.jar -h
```
Usage: java -jar HaploABBABABA_multithreaded.jar \
         -i <input.vcf> \
         -p <populations file> -c <file with combinations to test> \
         -o <output name> -t <n threads> -b <n bootstraps>
```


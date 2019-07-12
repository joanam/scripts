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

## PhiX_removal.sh
Bash script to remove PhiX (virus) reads from NGS data. This script generates a sam file with the aligned PhiX reads and a fastq file with all non-PhiX reads. This script may be useful to remove PhiX reads from RAD libraries.

```
Usage:  PhiX_removal < raw reads file > < output files prefix>
```

## demultiplexing.sh
Bash script to demultiplex RAD reads with barcodes of different lengths. It uses process_radtags from Stacks. It also trimms the reads if the trimming length is shorter than the read length. This may make sense to remove read ends of bad sequencing quality and 
to ensure that samples with different barcode lengths end up with the same demultiplexed read length. Note, the barcodes file needs to have one line per barcode giving the sample name and then tab-delimited the barcode.

```
Usage:  demultiplexing.sh <readsFastqFile> <barcodesFile> <outputFolder> <trimmingLength>
```

## runBowtie2RADID.sh
This script is meant for RAD data mostly for the fish ecology and evolution group at EAWAG as it requires a naming scheme used by that group. However, it may also be useful to other users as a guide to setup a script that automatically aligns all fastq files in a directory.


 
## HaploABBABABA_multithreaded # Attention: I am working on a bug fix, please do not currently use this software
Java program to calculate ABBA-BABA (D-statistics) to infer gene flow and the five-population test to infer the direction of gene flow. 

For more info, use: java -jar HaploABBABABA_multithreaded.jar -h
```
Usage: java -jar HaploABBABABA_multithreaded.jar \
         -i <input.vcf> \
         -p <populations file> -c <file with combinations to test> \
         -o <output name> -t <n threads> -b <n bootstraps>
```

## createRADmappingReport.sh
Bash script that generates a mapping report of RAD data in bam files including for each individual the number of RAD loci sequenced, the mean sequencing depth, the RAD loci covered by at least 10 reads, and the mean sequencing depth of these RAD loci. It also generates a file with one line per RAD locus and individual giving the sequencing depth for each. The script needs to be run in the directory containing the bam files. The report can be used to test how good a RAD library is. If samples with few reads have few loci sequenced at high depth instead of many loci at low depth, this is strong indication for high duplication levels.

```
Usage: cd <directory of bam files>; createRADmappingReport.sh
```

## vcf2fineRADstructure.sh
Bash script that converts a vcf file into the input format required for fineRADstructure (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html). It is rather slow as it uses vcftools to extract each RAD locus from the vcf file and then converts it to the correct format. It uses createRADmappingReport.sh to define the RAD loci. This script should be run with a phased vcf file as fineRADstructure needs haplotypes. For RAD data, I would use e.g. GATK ReadBackedPhasing (https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_phasing_ReadBackedPhasing.php).

```
Usage: vcf2fineRADstructure.sh <vcf file> <path to bam files>
```

#! /usr/bin/env python

# Author: Joana

from sys import *
import os, time, argparse, re, scipy, numpy
from collections import defaultdict
from scipy import stats

parser = argparse.ArgumentParser(description='Filter out genotypes with strong allic disbalance in the given VCF file')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-hom', '--homozygote', action='store_true', help="set failing genotypes as homozygous",default=False)
group.add_argument('-excl', '--exclude', action='store_true', help="set failing genotypes as missing", default=False) 
group.add_argument('-two', '--twoSteps', action='store_true', help="set failing genotypes as missing if pvalue<p and as homozygous if pvalue<p2", default=False)
parser.add_argument('-p', '--pvalue', type=float, dest='p', help="p-value threshold for binomial test [default: 0.01]", required=False, default=0.01)
parser.add_argument('-p2', '--pHomoz', type=float, dest='p2', help="second p-value threshold for binomial test if -two is specified, threshold to make failing genotypes homozygote [default=0.005]", default=0.005)
parser.add_argument('-r', '--ratio', type=float, dest='r', help="hard cutoff for allelic ratio [default=0.2]", default=0.2)

args = parser.parse_args()

readbytes=0
threshold=args.p
threshold2=args.p2
ratio=args.r

inputF=open(args.i,'r')
outputF=open(args.o, 'w')


for Line in inputF:
    # DATA SECTION: clause checks if the header section is over
    if re.match('^#',Line) is None:

        # Get the columns of that line
        columns=Line.strip("\n").split("\t")

        # Only check SNPs (monomorphic sites are not modified)
        if columns[5]!=".":
            # Add the info to the site
            result=columns[0:9]

            # Get the genotypes
            genotypecolumns=range(9,len(columns))

            # Check each individual if it is a heterozygote
            for ind in genotypecolumns:
                genotype=columns[ind]
                genotype=genotype.split(":")

                if "/" in genotype[0]:
                    alleles=genotype[0].split("/") 
                elif "|" in genotype[0]:
                    alleles=genotype[0].split("|") 
                else:
                    result.append(":".join(genotype))       

                # If the genotype is heterozygous check the allelic balance
                if alleles[0]!=alleles[1]:
                    reads=genotype[1].split(",")

                    # Calculate the probability for the observed allele distribution with a binomial test
                    pval=scipy.stats.binom_test(reads[0], n=int(reads[0])+int(reads[1]), p=0.5)

                    # if one of the alleles has no reads (weirdly this happens)
                    if int(reads[0])==0 or int(reads[1])==0:
                        result.append("./.")

                    # If significant allelic disbalance
                    elif pval<threshold:

                        # If -excl is specified, set failing genotypes as missing
                        if args.exclude:
                            result.append("./.")

                        # if -hom is specified, set failing genotypes as homozygous
                        elif args.homozygote:
                        
                            # Replace the genotype by the more common allele homozygote
                            if reads[0]>reads[1]:
                                genotype[0]=alleles[0]+"/"+alleles[0]
                            else:
                                genotype[0]=alleles[1]+"/"+alleles[1]

                            result.append(":".join(genotype))

                        # if -two is specified, set genotypes failing threshold as missing, and those failing threshold2 as homozygous
                        elif args.twoSteps:
                            if pval>threshold2:
                                result.append("./.")
                            else:
                                if reads[0]>reads[1]:
                                        genotype[0]=alleles[0]+"/"+alleles[0]
                                else:
                                        genotype[0]=alleles[1]+"/"+alleles[1]
                        else: 
                            exit("either -hom or -excl or -two is required, please specify how failing genotypes should be handled")

                    # If the binomial test is not significant but the allic ratio very small (at low depth possible) 
                    elif float(reads[0])/float(reads[1])<ratio or float(reads[1])/float(reads[0])<ratio:
                        result.append("./.")


                    # If the genotype makes the test
                    else:
                        result.append(":".join(genotype))

                # If the genotype is homozygous, just append as is
                else:
                    result.append(":".join(genotype))

            outputF.write('\t'.join(result)+"\n")

        # If it is not a SNP (monomorphic) just write the line out
        else:
            outputF.write(Line)

    # If it is a header line, just write it out
    else:
        outputF.write(Line)



inputF.close()
outputF.close()

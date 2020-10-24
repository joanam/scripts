#! /usr/bin/env python

# Authors: Samuel Wittwer, Joana Meier, David Marques, Irene Keller

from sys import *
import os, time, argparse, re, numpy
from collections import defaultdict

parser = argparse.ArgumentParser(description='Filter all sites in the given VCF file')

parser.add_argument('-i', '--input', dest='i', help="input file in vcf format [required]", required=True)
parser.add_argument('-o', '--output', dest='o', help="output file [required]", required=True)
parser.add_argument('-q', '--QUAL', dest='qual', type=int, help="minimum quality of the SNP (QUAL field of vcf) [default 30]", default=30)
parser.add_argument('-p', '--minGQ', dest='gq', type=int, help="minimum quality for each single genotype (GQ of genotype) [default 30]", default=30)
parser.add_argument('-d', '--readDepth', dest='depth', type=int, help="minimal read depth per single genotype (DP of genotype) [default 1]", default=1)
parser.add_argument('-n', '--minInd', dest='minInd', type=int, help="minimal number of individuals with gentoype per site [default 1]", default=1)
parser.add_argument('-m', '--minAllele', dest='minAllele', type=int, help="minimum number of times an allele must be observed [default 1]", default=1)
parser.add_argument('-v', '--verbose', action='store_true', help="if -v is specified, progress info will be written to stdout (slows the program down)", default=True)
parser.add_argument('-rindel', '--removeIndels', action='store_true', help="if -rindel is specified, indels will be filtered", default=False)
parser.add_argument('-risnp', '--removeIndelSNPs', action='store_true', help="if -risnp is specified, SNPs around indels (distance specified by -idist) will be filtered", default=False)
parser.add_argument('-idist', '--indelDist', type=int, help="SNPs at the given distance or smaller from indels will be filtered if -risnp given as additional argument [default 10]", default=10)
parser.add_argument('-rmono', '--removeMono', action='store_true', help="if -rmono is specified, monomorphic sites (except SNPs to the reference only) will be filtered", default=False)
parser.add_argument('-rmsnp', '--removeMultiallelicSNPs', action='store_true', help="if -rmsnp is specified, multiallelic (>2 alleles) SNPs will be filtered", default=False)
parser.add_argument('-rrsnp', '--removeRefSNPs', action='store_true', help="if -rrsnp is specified, SNPs to the reference only (=monomorphic within sample) will be filtered", default=False)
parser.add_argument('-np', '--minIndPop', dest='minIndPop', type=int, help="minimum number of individuals per population-code", default=0)
parser.add_argument('-e', '--exclude', dest='e', help="exclude popcodes containing the following string(s) from minIndPop-filtering", required=False, default=["!!!!!!!!!"], nargs='*')

args = parser.parse_args()

if args.verbose:
    print "\n\tfiltering the VCF file with QUAL>=%d, GQ>=%d, DP>=%d, minInd>=%d, minAllele>=%d" % (args.qual,args.gq,args.depth,args.minInd,args.minAllele)
    if args.minIndPop>0 and args.e==["!!!!!!!!!"]:
        print"\t + sites with <%d genotypes per population-code" %(args.minIndPop)
    if args.minIndPop>0 and args.e!=["!!!!!!!!!"]:
        print"\t + sites with <%d genotypes per population-code (except for populations containing the string %s) will be filtered" %(args.minIndPop,args.e)
    if args.removeIndels:
        print"\t + all indels will be filtered"
    if args.removeIndelSNPs:
        print"\t + SNPs %d bp near indels will be filtered" %(args.indelDist)
    if args.removeMono:
        print"\t + all monomorphic sites (except SNPs to the reference only) will be filtered"
    if args.removeRefSNPs:
        print"\t + all SNPs to the reference only (monomorphic within sample) will be filtered"
    if args.removeMultiallelicSNPs:
        print"\t + all multiallelic (>2 alleles) SNPs will be filtered"

filesize=os.stat(args.i)[6]
readbytes=0

inputF=open(args.i,'r')
outputF=open(args.o, 'w')

LineNumber=0 
relLineNumber=0
prevScaffold=""
pos=0
prevPos=0
sites=[]
positions=[]
indel=False
snp=False
exclu=0
if args.minIndPop>0 and args.e!=["!!!!!!!!!"]:
    exclu=1
SNPtomonoFail=False
snpCounter=0
refsnpCounter=0
multiallelicCounter=0
indelCounter=0
monoCounter=0
baseCounter=0
printedBases=0
indelPos=1000000000000
indelDist=100000000000

#failed sites counters:
qualFail=0
SNPtomonoCount=0
minIndFail=0
minIndPopFail=0
minAlleleFail=0
nearIndelSNPFail=0
refsnpfiltered=0
multiallelicfiltered=0
monofiltered=0
snpsPassed=0
refsnpsPassed=0
multiallelicPassed=0
indelsretained=0

for Line in inputF:
    # DATA SECTION: clause checks if the header section is over
    if re.match('^#',Line) is None:
        # REPORTING TO STDOUT: writes out every millionth line number and the percentage processed if -v selected
        relLineNumber+=1
        readbytes+=len(Line)
        done=str(int((float(readbytes)/filesize)*100))
        if args.verbose and relLineNumber>9999:
            LineNumber+=10000
            stdout.write("Line %d, %s%% done %s"%(LineNumber,done,"\r"))
            stdout.flush()
            time.sleep(1)
            relLineNumber=0

        columns=Line.strip("\n").split("\t")
        scaffold=columns[0]
        pos=int(columns[1])

        if scaffold == prevScaffold:
            indelDist+=(pos-prevPos)
        else:
            indelDist=10000000000000    # if it's a new scaffold prevent effects of indels in the previous scaffold

        # QUAL filtering - checks if overall quality of the site is good (-q)
        # If filtered, a site is not printed.
        if columns[5]=="." or float(columns[5])>=args.qual:
            result=columns[0:9]
            genotypecolumns=range(9,len(columns))
            indCounter = 0
            alleleCounter = 0
            countOne=0
            countZero=0
            countTwo=0
            countThree=0
            refsnp = False
            snp = False
            indel = False
            SNPtomonoFail=False
            multiallelic=False
            minIndAlleleCheck = False

            # GQ/DP filtering - checks if genotype quality (-p) and depth (-d) are sufficient for each single genotype. 
            # If filtered out, genotype is set to missing "./."
            for ind in genotypecolumns:
                genotype=columns[ind]
                genotype=genotype.split(":")

                # Case 1: Genotype missing "./." or "./.:.:4" or "0/0" without further information -> no filtering necessary
                if ('./.' in genotype) or (genotype==["0/0"]):
                    result.append("./.")

                # Case 2: Non-variant site equal to reference, has only 2 fields (GT:DP) -> only depth filtering possible
                elif len(genotype)==2:
                    if float(genotype[1])>=args.depth:
                        result.append(columns[ind])
                        indCounter+=1
                        countZero+=1
                    else:
                        result.append("./.")

                # Case 3: SNP/indel against on reference, but non-variant within samples, has only 3 fields (GT:AD:DP) -> only depth filtering possible
                elif len(genotype)==3:
                    if float(genotype[2])>=args.depth:
                        result.append(columns[ind])
                        indCounter+=1
                        countOne+=1
                    else:
                        result.append("./.")

                # Case 4: SNP/indel with genotype qualities, > 3 fields present (GT:AD:DP:GQ:PL) -> genotype quality and depth filtering
                elif len(genotype)>3:
                    if genotype[2]!=".":
                        if float(genotype[2])>=args.depth: 
                            if float(genotype[3])>=args.gq:
                                result.append(columns[ind])
                                indCounter+=1
                                if "/" in genotype[0]:
                                   alleles=genotype[0].split("/") # counts the occurrences of each allele (up to 3 alt alleles -> LIMITATION)
                                elif "|" in genotype[0]:
                                   alleles=genotype[0].split("|") # counts the occurrences of each allele (up to 3 alt alleles -> LIMITATION)
                                if int(alleles[0])==0:
                                    countZero+=1
                                elif int(alleles[0])==1:
                                    countOne+=1
                                elif int(alleles[0])==2:
                                    countTwo+=1
                                elif int(alleles[0])==3:
                                    countThree+=1
                                if int(alleles[1])==0:
                                    countZero+=1
                                elif int(alleles[1])==1:
                                    countOne+=1
                                elif int(alleles[1])==2:
                                    countTwo+=1
                                elif int(alleles[1])==3:
                                    countThree+=1
                            else:
                                result.append("./.")
                        else:
                            result.append("./.")
                    else:
                        result.append("./.")
            # allele counts after genotype-based filtering
            countArray=[countZero, countOne, countTwo, countThree]

            # Identification of indels, SNPs and refSNPs
            # also counts bases, indels, SNPs and refSNPs
            baseCounter+=1
            if len(columns[3])>1: # in case of a deletion (reference allele is more than one base) -> indel
                indel = True
                indelPos = pos
            elif len(columns[4])>1: # could be an indel or >1 alternative alleles
                alleles = columns[4].split(",")  # if there are more than 1 alternative alleles, they are separated by commas
                for alt in alleles:
                    if len(alt)>1: # in case of an insertion (more than one base in alt allele) -> indel
                        indel = True
                        indelPos = pos
                        snp = False
                        break
            if indel:
                indelCounter+=1
            elif (columns[4]!="."):
                snp = True
                snpCounter+=1
                if (len(columns[3])>1 or len(columns[4])>1):
                    multiallelic = True
                    multiallelicCounter+=1
                # In case variable genotypes got filtered (SNP->monomorphic) or SNP is a refSNP (SNP against reference, but monomorphic within sample)
                if (sum(i>0 for i in countArray)<2):
                    if countArray[0]>0:
                        snp = False
                        SNPtomonoFail = True
                    else:
                        refsnp = True
                        refsnpCounter+=1
            else:
                monoCounter+=1

            # minIND and minAllele filtering (-n & -m), only for non-INDELs
            # Checks that there are enough individuals that passed and that there are enough alternative allele counts
            if indel:
                minIndAlleleCheck = True
            elif indCounter >= args.minInd:
                if args.minIndPop>0:
                    for m in range(exclu,len(NbrPerSpecies)):
                        present=0
                        absent=0
                        for n in range(0,len(result)-9):
                            if levels[n]==m:
                                if './.' in result[n+9]:
                                    absent+=1
                                else:
                                    present +=1
                        if (present>=args.minIndPop) or (absent==0):
                            continue
                        else:
                            break
                    if m==(len(NbrPerSpecies)-1) and (present>=args.minIndPop or absent==0):	# checks if the loop was exited early. Check of present/absent is necessary in case the last factor level is the only one with insufficiant number of genotypes
                        if (not snp):
                            minIndAlleleCheck = True
                        elif (snp and (sum(i>=args.minAllele for i in countArray) >= sum(i>0 for i in countArray))):
                            minIndAlleleCheck = True
                            snpsPassed+=1
                            if refsnp:
                                refsnpsPassed+=1
                            if multiallelic:
                                multiallelicPassed+=1
                        else:
                            minAlleleFail+=1
                    else:
                        minIndAlleleCheck = False
                        minIndPopFail+=1
                else:
                    if (not snp):
                        minIndAlleleCheck = True
                    elif (snp and (sum(i>=args.minAllele for i in countArray) >= sum(i>0 for i in countArray))):
                        minIndAlleleCheck = True
                        snpsPassed+=1
                        if refsnp:
                            refsnpsPassed+=1
                        if multiallelic:
                            multiallelicPassed+=1

                    else:
                        minAlleleFail+=1
            else:
                minIndFail+=1

            # Counts the number of SNPs that were converted to monomorphic loci and that passed quality filtering.
            if SNPtomonoFail and minIndAlleleCheck:
                SNPtomonoCount+=1

            # Direct OUTPUT without INDEL-proximity SNP filter: Writes locus to file, if -risnp option is not specified
            if(not args.removeIndelSNPs) and minIndAlleleCheck:
                if indel:
                    if (not args.removeIndels):
                        outputF.write('\t'.join(result)+"\n")
                        printedBases+=1
                        indelsretained+=1
                elif (not snp):
                    if args.removeMono:
                        monofiltered+=1
                    else:
                        outputF.write('\t'.join(result)+"\n")
                        printedBases+=1
                elif snp:
                    if refsnp and args.removeRefSNPs:
                        refsnpfiltered+=1
                    elif multiallelic and args.removeMultiallelicSNPs:
                        multiallelicfiltered+=1
                    else:
                        outputF.write('\t'.join(result)+"\n")
                        printedBases+=1

            # Looped OUTPUT with INDEL-proximity SNP filter: If option -ri is specified, INDEL-proximity SNPs are filtered
            elif args.removeIndelSNPs and minIndAlleleCheck:
                # start over if scaffold is new -> write out former sites & start memorizing new sites
                if scaffold!=prevScaffold:      # if it's a new scaffold, write out the previous lines and empty the lists
                    if len(sites)>0:            # to not create an empty line if nothing contained in sites
                        outputF.write('\n'.join(sites)+"\n")
                        printedBases+=len(sites)
                        del sites[:]
                        del positions[:]
                    if indel:
                        indelDist=0
                        if (not args.removeIndels):
                            outputF.write('\t'.join(result)+"\n") #print '\t'.join(result)+"\n" 
                            printedBases+=1
                            indelsretained+=1
                    else:  # if it's not an indel
                        if (not args.removeMono) and (not args.removeRefSNPs) and (not args.removeMultiallelicSNPs):
                            sites.append('\t'.join(result))
                            if snp:
                                positions.append(pos)
                            else:
                                positions.append(-1000) # to force printing of monomorphic sites even near indels
                        elif (not snp):
                            if args.removeMono:
                                monofiltered+=1
                            else:
                                sites.append('\t'.join(result))
                                positions.append(-1000)
                        elif snp:
                            if refsnp and args.removeRefSNPs:
                                refsnpfiltered+=1
                            elif multiallelic and args.removeMultiallelicSNPs:
                                multiallelicfiltered+=1
                            else:
                                sites.append('\t'.join(result))
                                positions.append(pos)

                # If current site is an INDEL, writes out non-SNPs and SNPs outside the proximity limit
                elif indel:    # if the current site is an indel: write out previous lines that are distant enough and empty the lists
                    indelDist=0
                    for index in range(len(positions)):
                        if positions[index]<(indelPos-args.indelDist):
                            outputF.write(sites[index]+"\n") #print sites[index]+"\n" # 
                            printedBases+=1
                        else:
                            nearIndelSNPFail+=1
                    if (not args.removeIndels):
                        outputF.write('\t'.join(result)+"\n") #print '\t'.join(result)+"\n" 
                        printedBases+=1
                        indelsretained+=1
                    del sites[:]
                    del positions[:]

                # For current site outside proximity range of INDEL, appends to "sites" and "positions"
                elif indelDist>args.indelDist:      # if the distance to the nearest indel is large enough
                    if len(sites)>args.indelDist:   # if there are enough sites in the list, write out the first site and remove it
                        outputF.write(sites[0]+"\n") # print sites[0]+"\n" #
                        printedBases+=1    
                        sites.pop(0)
                        positions.pop(0)
                    if (not args.removeMono) and (not args.removeRefSNPs) and (not args.removeMultiallelicSNPs):
                        sites.append('\t'.join(result))
                        if snp:
                            positions.append(pos)
                        else:
                            positions.append(-1000) # to force printing of monomorphic sites even near indels
                    elif (not snp):
                        if args.removeMono:
                            monofiltered+=1
                        else:
                            sites.append('\t'.join(result))
                            positions.append(-1000)
                    elif snp:
                        if refsnp and args.removeRefSNPs:
                            refsnpfiltered+=1
                        elif multiallelic and args.removeMultiallelicSNPs:
                            multiallelicfiltered+=1
                        else:
                            sites.append('\t'.join(result))
                            positions.append(pos)

                # right after the indel, prints only non-SNPs to the sites variable
                elif (not snp):
                    if args.removeMono:
                        monofiltered+=1
                    else:
                        sites.append('\t'.join(result))
                        positions.append(-1000)

                # Counts filtered SNPs because of INDEL proximity
                else:
                    if refsnp and args.removeRefSNPs:
                        refsnpfiltered+=1
                    elif multiallelic and args.removeMultiallelicSNPs:
                        multiallelicfiltered+=1
                    else:
                        nearIndelSNPFail+=1
        else:
            qualFail+=1

        prevScaffold=scaffold
        prevPos=pos

    
    # HEADER
    # check if it is a header line and if yes, write it to the output file
    else:  # writes header information into outfile
        outputF.write(Line)
        if re.match('^##',Line) is None:
            if args.minIndPop>0:
                header=Line.strip("\n").split("\t")
                header=header[9:len(header)]  # header now contains all the individual IDs
                species=[]
                excluded=[]
                for i in range(0,len(header)):
                    ID=header[i]
                    sub=ID.split('.')[1]	# extracts popcode from our standard RAD naming convention, i.e. "FishecNumber"."popcode"."library"
                    regex=re.compile("|".join(args.e))	# search pattern for excluded popcode (-e option)
                    if regex.search(sub) is not None:
                        species.append("!")	# gives "!" value to excluded popuations -> this gives excluded popcodes the first factor position in the function categorical
                        excluded.append(sub)
                    else:
                        species.append(sub)
                a=numpy.array(species)
                poplist=list(set(a))
                levels=[dict([list(reversed(i)) for i in list(enumerate(sorted(list(set(a)))))])[j] for j in a]
                NbrPerSpecies=defaultdict(int)
                for k in levels:
                    NbrPerSpecies[k]+=1		# NbrPerSpecies contains counts of how often each factor level is observed. There must be a simpler way to do this. (counting nbr of columns in design matrix produced by categorical should do the trick)
                if len(NbrPerSpecies)==1:	# Aborts with warning message if all samples were excluded from population-based filtering by -e option
                    outputF.close()
                    os.remove(args.o)
                    exit("\n\tWARNING: all samples were excluded from filtering for minimum number of individuals per population-code [-e option] \n\t\tor the sample names are not in standard format [fishecnumber.popcode.library]\n\tFiltering aborted! Disable population-code based filtering (-np 0) or choose another -e value!")
                del header
                del ID
                del sub
                del regex
                del a
        
        
# Outputs the site variable after the last line in the input, if the Indel-associated SNP removal option was chosen
if args.removeIndelSNPs:
    for index in range(len(positions)):
        outputF.write(sites[index]+"\n") #print sites[index]+"\n" # 
        printedBases+=1

if args.verbose:
    print "\n\tInput file contains: \n\t %d\ttotal sites with QUAL>=%d, consisting of:" %(baseCounter,args.qual)
    print "\t %d\tmonomorphic sites" %(monoCounter)
    print "\t %d\tSNPs (%d SNPs to the reference only and %d SNPs polymorphic in sample, incl. %d multiallelic SNPs)" %(snpCounter,refsnpCounter,snpCounter-refsnpCounter,multiallelicCounter)
    print "\t %d\tindels" %(indelCounter)

    print "\n\tFiltering report:"
    print "\t %d\tsites filtered out due to low quality (QUAL)" %(qualFail)
    print "\t %d\tsites filtered out because <%d individuals covered (-minInd)" %(minIndFail,args.minInd)
    if args.minIndPop>0:
        print "\t %d\tsites filtered out because <%d individuals covered per population-code (-minIndPop) " %(minIndPopFail,args.minIndPop)
    if args.e!=["!!!!!!!!!"]:
        print "\t \tfollowing population-codes not filtered out for min. number of individuals (-excluded):"
        print "\t \t "+", ".join(list(set(excluded)))
    print "\t %d\tSNPs filtered out due to <%d minor allele counts (-minAllele)" %(minAlleleFail,args.minAllele)
    print "\t %d\tSNPs converted to monomorphic sites due to filtered genotypes (GQ and/or DP)" %(SNPtomonoCount)
    if args.removeIndels:
        print "\t %d\tindels were filtered out" %(indelCounter)
    if args.removeIndelSNPs:
        print "\t %d\tvalid SNPs were filtered out due to within %d bp proximity to indels" %(nearIndelSNPFail,args.indelDist)
    if args.removeMono:
        print "\t %d\tvalid monomorphic sites were filtered out" %(monofiltered)
    if args.removeRefSNPs:
        print "\t %d\tvalid SNPs to the reference only (monomorphic within sample) were filtered out" %(refsnpfiltered)
    if args.removeMultiallelicSNPs:
        print "\t %d\tvalid multiallelic (>2 alleles) SNPs were filtered out" %(multiallelicfiltered)

    print "\n\tOutput file %s contains: \n\t %d\ttotal sites, consisting of:" %(args.o,printedBases)
    print "\t %d\tmonomorphic sites" %(printedBases-indelsretained-(snpsPassed-nearIndelSNPFail-refsnpfiltered-multiallelicfiltered))
    print "\t %d\tSNPs (%d SNPs to the reference only and %d SNPs polymorphic in sample, incl. %d multiallelic SNPs)" %(snpsPassed-nearIndelSNPFail-refsnpfiltered-multiallelicfiltered,refsnpsPassed-refsnpfiltered,snpsPassed-refsnpsPassed-nearIndelSNPFail-multiallelicfiltered,multiallelicPassed-multiallelicfiltered)
    print "\t %d\tindels\n" %(indelsretained)

inputF.close()
outputF.close()

#!/cluster/apps/r/3.1.2_openblas/x86_64/bin/Rscript --vanilla

# Written by Joana Meier, May-2017
# usage: addRecombRates.r <vcffile> <linkage map>
# Generates a map file in plink format for the positions provided in the vcf file

# Load the required packages
require(data.table)
require(gtools)

# Read in the vcf file name and linkage map name (including path) from the command line argument
args<-commandArgs(TRUE)
file<-args[1]
map<-args[2]
prefix<-sub('\\.vcf.*', '', file) 

# If the file is gzipped, add info to gunzip
if(grepl(file,pattern=".gz")) file=paste("gunzip -c ",file,sep="")

# Read in the data (just the first 5 columns)
data<-fread(cmd=file,header=T,skip="#CHR",select=c(1:5),data.table=F)

# Extract the first five columns with position information and add a column for recombination distances
data<-cbind(data,"Morgan"=vector(length = length(data[,1])))

# Read in the table of physical vs recombination distance
recomb<-read.table(map,header=T,sep="\t")

# For each chromosome separately, add the recombination distance
# For unmapped scaffolds, use the empirical mean of 2cM/Mb

chrom<-mixedsort(levels(as.factor(data[,1])))

for(i in 1:length(chrom)){

    # get the name of the ith chromosome/scaffold
    chr=chrom[i]

    # get the lines of the dataset of this specific chromosome
    dataChr<-data[data[,1]==chr,]

    # Get physical and recombination distances for that chromosome:
    recChr<-recomb[recomb$CHROM==chr,]
    if(length(recChr$CHROM)>0){
     d<-smooth.spline(recChr$POS,recChr$cM,spar=0.7)
     recPos<-predict(d,dataChr[,2])$y
     recPos[recPos<0]<-0

     # Make the recombination distance monotonously increasing:
     for(i in 2:length(recPos)){
         if(recPos[i]<recPos[i-1]) recPos[i]=recPos[i-1]+0.0000001
     }
    }
    else{
     recPos<-dataChr$POS/0.5e6
    }
    data[data[,1]==chr,"Morgan"]<-recPos/100
}

toPrint<-cbind(data[,1],paste(data[,1],data[,2],sep=":"),data[,"Morgan"],data[,2])

write.table(toPrint,paste(prefix,".plink.map",sep=""),row.names=F,quote=F,sep=" ",col.names=F)



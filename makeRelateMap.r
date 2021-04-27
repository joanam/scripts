#!Rscript --vanilla

# Written by Joana Meier, Oct 2019
# usage: makeRelateMap.r <vcffile> <linkage map>
# Generates a map file with the position, recomb rate, cM position as required by Relate

# Load the required packages
require(data.table)
require(gtools)

# Read in the vcf file name and linkage map name (including path) from the command line argument
args<-commandArgs(TRUE)
file<-args[1]
map<-args[2]
prefix<-sub('\\.vcf.*', '', file) 

# avoid scientific notation
options(scipen=999)

# If the file is gzipped, add info to gunzip
if(grepl(file,pattern=".gz")) file=paste("gunzip -c ",file,sep="")

# Read in the data (just the first 5 columns)
data<-fread(file,header=T,skip="#CHR",select=c(1:5),data.table=F)

# Extract the first five columns with position information and add a column for recombination distances
data<-cbind(data,"rec"=vector(length = length(data[,1])),"cM"=vector(length = length(data[,1])))

# Read in the table of physical vs recombination distance
recomb<-read.table(map,header=T,sep="\t")

# For each chromosome separately, add the recombination distance
# For unmapped scaffolds, use the empirical mean of 2cM/Mb

chrom<-mixedsort(levels(as.factor(data[,1])))

for(i in length(chrom)){

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
         if(recPos[i]<=recPos[i-1]) recPos[i]=recPos[i-1]+0.000000001
     }
     # predict the recombination rate for each position
     recRate<-stats:::predict.smooth.spline(d,as.integer(dataChr[,2]),deriv=1)$y

     # Set 0 / negative recombination rates to tiny value
     recRate[recRate<=0]<-0.000000000001
    }
    else{
     recPos<-dataChr$POS/0.5e6
     recRate<-rep(2,times=length(dataChr[,1]))
    }
    data[data[,1]==chr,length(names(data))-1]<-recPos
    data[data[,1]==chr,length(names(data))]<-recRate
}

toPrint<-cbind(data[,2],data[,length(names(data))],data[,length(names(data))-1])

write.table(toPrint,paste(prefix,".relate.map",sep=""),row.names=F,quote=F,sep=" ",col.names=F)



#!/cluster/apps/r/3.1.2_openblas/x86_64/bin/Rscript --vanilla

# Written by Joana Meier, May-2017
# usage: makeRecombMap4fineSTRUCTURE.sh <vcffile> <linkage map>
# Generates a recomb file for fineSTRUCTURE with recombination rates in M/bp


# Check if there were two arguments given
if [ $# -lt 2 ]
then
        echo -e "ERROR: Not enough arguments provided!\nUsage:  makeRecombMap4fineSTRUCTURE.sh <vcffile> <linkage map>"
        exit 1
fi


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

# Read in the data (just the column with the positions)
data<-fread(file,header=T,skip="#CHR",select=c(1,2),data.table=F)

# add a column for recombination distances in Morgan as required by fineSTRUCTURE
data<-cbind("chrom"=data[,1],"start.pos"=data[,2],"recom.rate.perbp"=rep(0.00,length(data[,1])))
data<-as.data.frame(data)
data$start.pos<-as.integer(as.character(data$start.pos))
data$recom.rate.perbp<-as.double(data$recom.rate.perbp)

# Read in the linkage map, i.e. table of physical vs recombination distance, teb-delimited
recomb<-read.table(map,header=T,sep="\t")

# For each chromosome separately, add the recombination distance
# For unmapped scaffolds, use the empirical mean of 2M/Mb = 0.0002

chrom<-mixedsort(levels(as.factor(data[,1])))
add=0

for(i in 1:length(chrom)){

    # get the name of the ith chromosome/scaffold
    chr=chrom[i]

    # get the lines of the dataset of this specific chromosome
    dataChr<-data[data[,1]==chr,]
    n<-length(dataChr)

    # Get physical position and the local recombination rate, if no info in linkage map use mean of 2 cM/Mb
    recChr<-recomb[recomb$CHROM==chr,]
    if(length(recChr$CHROM)>0){
     d<-smooth.spline(recChr$POS,recChr$cM,spar=0.7)
     recRate<-stats:::predict.smooth.spline(d,as.integer(dataChr[,2]),deriv=1)$y*100
     recRate[recRate<0]<-0.0000001
    }
    else{
     recRate<-rep(0.00025,times=length(dataChr[,1]))
    }
    recRate[length(recRate)]<-(-9.00)
    data[data[,1]==chr,3]<-recRate
    data[data[,1]==chr,2]<-data[data[,1]==chr,2]+add

    add=add+dataChr[length(dataChr[,1]),2]
    print(paste(chr," finished, starting at ",add,sep=""))
}

write.table(data[,2:3],paste(prefix,".recomb",sep=""),row.names=F,quote=F,sep=" ",col.names=T)

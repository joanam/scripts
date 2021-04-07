#!~/bin/Rscript --vanilla

# Usage: removeNonsignFbranch.r -i <fbranch file with Zscores> -z <zscore threshold> -o <output filename>

# Load libaries
library(optparse)

# Read input arguments
option_list = list(
  make_option(c("-z", "--zthreshold"), type="double", default=3.0,
              help="z score threshold", metavar="double"),
  make_option(c("-i", "--infile"), type="character", default=NULL,
              help="file with fbranch scores with z scores produced by running Dsuite Fbranch with -Z True", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="output file name", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

threshold<-opt$zthreshold

# Get the line number where z scores start
test<-read.table(file=opt$infile,sep=";",as.is=T,comment.char = ":")
lines<-which(grepl(test[,1],pattern="#"))-2

# Read in the fbranch values and z scores
fbranch<-read.table(opt$infile,header=T,as.is=T,comment.char = "#",nrows=lines)
zscore<-read.table(opt$infile,header=T,as.is=T,comment.char = "#",skip=lines+1)

# Extract values entries
fMatrix<-as.matrix(fbranch[-1,-c(1,2,3)])
zMatrix<-as.matrix(zscore[-1,-c(1,2,3)])

# Set fbranch values to zero if z is below threshold
fMatrix[apply(zMatrix,MARGIN=1:2,FUN = function(x) x<threshold)]<-0

# New matrix with non-significant f4 ratio values set to 0
fbranchNew<-cbind(fbranch[-1,c(1,2,3)],fMatrix)
fbranchNew<-rbind(fbranch[1,],fbranchNew)

# Write out the modified fbranch matrix
write.table(fbranchNew,file=opt$outfile,row.names = F,quote=F,na = "nan",sep = "\t")

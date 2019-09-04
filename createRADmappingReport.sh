#!/bin/bash

# calculates the number of reads, number of loci, mean sequencing depth and number of loci with min 10 reads and their mean depth  

# create the report.txt file with a line for each individual (each bam file) with read counts
echo -e "sample\tsampleLib\tmappedReads" > report.txt
for i in *${1}*.bam
do
echo "counting reads in "$i
totalReads=`samtools view $i | wc -l`
echo -e ${i%.GQI*}"\t"${i%.bam}"\t"$totalReads >> report.txt
done

# creates a file with a line for each locus with mapped reads in each contig for each individual

echo "sample marker scaffold locus orientation depth" > seq_depth.txt
echo "the following samples will be included: "
echo *${1}*.bam

for i in *${1}*.bam
do
  	echo "working on "$i

        samtools view $i | awk \
        'BEGIN{
         markerold=""
	 marker=""
	 scaffold=""
         locus=""
	 orientation=""
	 depth=""
	 cigar=""
	 first="false"
	 sample=""	
        }
	{
	 if($5>30){		# mapping quality MAPQ>30 to filter out badly aligned reads

	  # for the first scaffold: 
	  if(scaffold==""){

	    # get the sample name	
            for(i=1;i<=NF;i++) if($i ~ /RG:Z:.*/){sample=$i} 

	    # if reverse check for indels and use the last site - 4
	    if($2==16){
	     b=0
             cigar=$6
	     gsub(/M/,"+",cigar);gsub(/D/,"+",cigar);gsub(/\d*.I/,"",cigar);sub(/$/,"0",cigar)
             split(cigar,indexx,"+");for (i in indexx) b+=indexx[i]
	     locus=$4+b-4
	    }
	    # if forward, set locus to the first base position
	    else{
	     locus=$4
	    }

	    # set the variables marker, scaffold, orientation and depth
	    marker=$3"_"locus"_"$2
	    scaffold=$3
	    orientation=$2
  	    depth=1
	  }
	
 	  # if the locus is from the same scaffold as the last one:
	  else if(scaffold==$3){
	    markerold=marker
	    locusold=locus
	    
  	    # for reverse reads
	    if($2==16){
             b=0
             cigar=$6
             gsub(/M/,"+",cigar);gsub(/D/,"+",cigar);gsub(/\d*.I/,"",cigar);sub(/$/,"0",cigar)
             split(cigar,indexx,"+");for (i in indexx) b+=indexx[i]
             locus=$4+b-4
            }
	    # for forward reads
            else{
             locus=$4
            }
	    marker=$3"_"locus"_"$2

	    # if it is the same locus increase the variable depth by 1
	    if(marker==markerold){
	     depth+=1
	    }

 	    # if it is a new locus:
	    else{
	     # print the information about the previous locus:
	     print sample,markerold,scaffold,locusold,orientation,depth
	     
	     scaffold=$3
             if($2==16){
              b=0
              cigar=$6
              gsub(/M/,"+",cigar);gsub(/D/,"+",cigar);gsub(/\d*.I/,"",cigar);sub(/$/,"0",cigar)
              split(cigar,indexx,"+");for (i in indexx) b+=indexx[i]
              locus=$4+b-4
             }
             else{
              locus=$4
             }
             marker=$3"_"locus"_"$2
             orientation=$2
             depth=1
	    }
	  }
	  
	  # if it is a new scaffold:
	  else{
	    print sample,marker,scaffold,locus,orientation,depth
	    scaffold=$3
	    if($2==16){
             b=0
             cigar=$6
             gsub(/M/,"+",cigar);gsub(/D/,"+",cigar);gsub(/\d*.I/,"",cigar);sub(/$/,"0",cigar)
             split(cigar,indexx,"+");for (i in indexx) b+=indexx[i]
             locus=$4+b-4
            }
            else{
             locus=$4
            }
            marker=$3"_"locus"_"$2
	    orientation=$2
            depth=1
	   }
	  }
	}END{
	  # print the information of the last locus
          print sample,marker,scaffold,locus,orientation,depth
        }' | sed 's/RG:Z://g' >> seq_depth.txt
done

awk '{if($6>9) print $0}' seq_depth.txt > seq_depth_min10.txt


# get per locus information
 
echo -e "locus\ttotalInds\ttotalReads" > contigs.txt
cat seq_depth.txt | \
awk '{
 counter[$2]++
 totReads[$2]+=$6
}
END{
 for (i in counter) print i"\t"counter[i]"\t"totReads[i]
}' >> contigs.txt 

# get per locus information (count only individuals with at least 10 reads)

echo -e "locus\ttotalInds\ttotalReads" > lociMin10reads.txt
cat seq_depth_min10.txt | \
awk '{
 counter[$2]++
 totReads[$2]+=$6
}
END{
 for (i in counter) print i"\t"counter[i]"\t"totReads[i]
}' >> lociMin10reads.txt



# create a file with lociN and meanDepth per ind
echo -e "sample lociN meanDepth" > indInfo.txt
awk '{
counter[$1]++
depth[$1]+=$6
}
END{
for(i in counter) print i,counter[i],depth[i]/counter[i]
}' seq_depth.txt >> indInfo.txt


# create a file with lociN and meanDepth per ind counting loci with min 10 reads
echo -e "sample lociN meanDepth" > indInfoMin10reads.txt
awk '{
counter[$1]++
depth[$1]+=$6
}
END{
for(i in counter) print i,counter[i],depth[i]/counter[i]
}' seq_depth_min10.txt >> indInfoMin10reads.txt



# sort the stats files skipping the header line (careful: "sample" cannot be included in any sample name)
sort indInfoMin10reads.txt | awk '!/sample/' > sorted_indInfoMin10reads.txt
sort indInfo.txt | awk '!/sample/' > indInfoSorted
sort report.txt | awk '!/sample/' > sortedReport


# create a file with Nreads, lociN, meanDepth... for each individual
join -a 1 -a 2 sortedReport indInfoSorted -1 1 -2 1 > indInfoAdv.txt
sort indInfoAdv.txt > indInfoAdvSorted
echo "sample sampleLib mappedReads lociN meanDepth lociNmin10reads meanDepthMin10reads" > mappingReport.txt
join -a 1 -a 2 indInfoAdvSorted sorted_indInfoMin10reads.txt -1 1 -2 1 >> mappingReport.txt


# Remove files not needed anymore:
rm indInfoAdvSorted indInfoAdv.txt sortedReport indInfoSorted sorted_indInfoMin10reads.txt contigs.txt \
# rm	lociMin10reads.txt indInfoMin10reads.txt indInfo.txt report.txt  

# rm seq_depth*

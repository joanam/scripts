#! /usr/bin/env python

# python version 2.7.2+
# by Joana, script to convert vcf to phylip


import sys, getopt, getpass, RAD

def main(argv):
	input = sys.stdin
	output = sys.stdout
	writeref = False
	fill = False	
	noIndels = False
	prev=100000000000000000000000000000000000
	
	try:
		opts, args = getopt.getopt(argv, "i:o:", ["input","output"])
		for opt, arg in opts:
			if opt in ("-i", "--input"):
				input = open(arg,'r')
			if opt in ("-o", "--output"):
				output = open(arg,'w')
			if opt in ("-f", "--fill"):
				fill=True
                        if opt in ("-r", "--reference"):
                                writeref = True
			if opt in ("-e", "--excludeIndels"):
				noIndels = True

	except getopt.GetoptError:
		print "\nUsage: vcf2phylip.py -i <input.vcf> -o <output.py> [optional: -r -f -e]"
		print "\tif -r is specified, the reference sequence will be included in the phylip file"
		print "\tif -f is specified, all sites not in the vcf file will be printed as missing (N)"
		print "\tif -e is specified, indels are not printed (else replaced by N)"
		sys.exit(2)

	if noIndels and fill:
		print "-f and -e are incompatible! Please decide if you want all sites or not"
		sys.exit(2)

	headerinfo = RAD.vcf_ExtractHeaderInfo(input)	#skips header, retains IDs in headerinfo

	IDs = []	#IDs holds all sample names without equal spacing (used in nexus)
	IDs.append("reference ")
	resultsequences = []	#resultsequences holds all complete sequences to be written to phylip or nexus file
	resultsequences.append([])		#resultsequences[0] holds reference
	for entry in headerinfo[1]:
		IDs.append(entry.replace("-",".")+" ")
		resultsequences.append([])
	
	samplenames = RAD.FillUp(IDs)
	linecounter = 0
	print "\ngenerating phylip file with ",len(samplenames)-1," individuals"
	for line in input:
                site = line.strip('\n').split('\t')
                indel = False
		pos = int(site[1])
	
		# If missing positions should be filled up with Ns (-f specified, e.g. sites of low quality that were filtered out)
		if pos > (prev + 1) and fill:
			addLine=pos-(prev+1)
	                individualcounter = 1
                        for individual in site[9:]:
                                resultsequences[individualcounter]+= "N" * addLine
                                individualcounter += 1
			linecounter += addLine
	                resultsequences[0] += "N" * addLine

		# site contains a deletion, replace by missing data
		if len(site[3])>1:
			indel=True

		else:
			alleles = site[4].split(",")  # if there are more than 1 alternative alleles, they are separated by commas
        	        for alt in alleles:
	                	if len(alt)>1 or '*' in alt:  # in case of an insertion
 					indel=True
					break

		alternativeslist = []
		alternativeslist.append(site[3])
		if site[4] != ".":
			for entry in site[4].split(","):
				alternativeslist.append(entry)
                individualcounter = 1


		if not indel:
			for individual in site[9:]:
				if './.' in individual:
					resultsequences[individualcounter]+="N"
				else:
					resultsequences[individualcounter] += GetGenotype(individual[:3].split("/"), alternativeslist)
                                individualcounter += 1
	                linecounter += 1
	                resultsequences[0] += site[3]


		# If the site is an indel and noIndels is not specified, print as missing data (else not printed)
		elif not noIndels:
                        for individual in site[9:]:
                                resultsequences[individualcounter]+= "N"
                                individualcounter += 1
			linecounter += 1
			resultsequences[0] += site[3][:1]

		prev=int(site[1])

	input.close()
	
	RAD.phy_WritePhylipSequences2(samplenames, resultsequences, output, False)
	output.write("\n")
	output.close()


# vcf_NumberToBase
# Takes GT field [list, sep="/"]and alternativeslist [REF, Alt1, .... , AltN] and returns string from Ambiguity Matrix according to standard ambiguity codes
# Ambiguity matrix applied is:				A	C	G	T	.
# See CoordinatesDictionary.			A	A	M	R	W	N
#						C	M	C	S	Y	N
#						G	R	S	G	K	N
#						T	W	Y	K	T	N
#						.	N	N	N	N	N


AmbiguityMatrix = [["A","M","R","W","N"],["M","C","S","Y","N"],["R","S","G","K","N"],["W","Y","K","T","N"],["N","N","N","N","N","N"]]
CoordinatesDictionary = { "A":0 , "C":1, "G":2, "T":3, ".":4 }
GetGenotype = lambda individual, altlist: AmbiguityMatrix[CoordinatesDictionary[altlist[int(individual[0])]]][CoordinatesDictionary[altlist[int(individual[1])]]]

def vcf_NumberToBase(GTlist, alternativeslist):
	if GTlist[0] == ".":
		return "N"
	else:
		return str(GetGenotype(GTlist, alternativeslist))


if __name__ == '__main__':
	main(sys.argv[1:])


import argparse
from argparse import RawTextHelpFormatter
import sys
import re
import os

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Changes Caroline's vcflike format to a vcf format readable by vcftools to calculate genotype correlations (LD).

""")
parser.add_argument('--vcflike', required=True, dest='vcflike', type=str, help="vcflike input file (Caroline's format) to be changed to vcftools readable format (usable for calculating genotype correlations).")
parser.add_argument('--sample_names', required=False, dest='snames', type=str, help="Comma separated list of samples names in the vcflike input file that should be present in the vcf output file.", default="")

args = parser.parse_args()


#-------------------
# function
#-------------------


#def function_name(x,y):


####################
# main
####################

f=open(args.vcflike,"r")

# print vcf header
print "##fileformat=VCFv4.0"

for l in f:
	a=l.split()
		
	if re.match("CHR",l):
		asamp=a[8:len(a)]
# 		print asamp
		if (args.snames!=""):
			snames=args.snames.split(",") # sample names
			# index sample names, has the length of all samples in the input file plus one for all the header information combined
			isnames=[i+1 for i, j in enumerate(asamp) if j in snames] 
			isnames.append(0)
			isnames=sorted(isnames)
		else:
			isnames=range(0,(len(a)-7)) # number of samples + 1
# 		print isnames
		h=["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"]
		h.extend(asamp)
		h=list(h[i] for i in isnames)
		print "\t".join(h)
			
	else:
		CHROM,POS,REF,ALT= list( a[i] for i in [0,1,2,3])		
		out=["\t".join([CHROM,POS,".",REF,ALT,".","PASS",".","GT:DP:GQ"])]
		
		samp=[]
		asamp=a[8:len(a)]
		for s in range(0,len(asamp)):
			ss=asamp[s].split(":")
			sss=list( ss[i] for i in [0,4,5] )
			samp.append( ":".join(sss) )
		out.extend(samp)
		
		out=list(out[i] for i in isnames)
		print "\t".join(out)
	

#CHR	POS.CHR	REF	ALT	NUM.ALTS	NUM.POLY	NUM.MISS	FORMAT
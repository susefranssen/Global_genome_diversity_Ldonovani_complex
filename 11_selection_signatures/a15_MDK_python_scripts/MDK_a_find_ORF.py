
import argparse
from argparse import RawTextHelpFormatter
import sys
import re
import os

from Bio import Seq
from Bio import SeqIO
import Bio
import regex 



parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Note: This script is specific for the Ldon analysis as file paths and group names
for this project are hard coded in the file. :(

Takes fasta file with gene sequences. It searches the longest ORF on the 5'-3' strand 
(assumes data is stranded) and returns it as fasta to stout.

Format of input fasta -CDS:
Each sequence header contains three elements separates by '_'
1. gene name
2. chr name
3. start position of gene
4. end position of the gene
5. strand

The nt sequence is in one line per gene.

Format of the output fasta in sdout:
Each sequence header contains four elements separates by '_'
1. gene name
2. chr name
3. new start position of gene starting at longest ORF
4. new end position of gene ending at the longest ORF

The created log file contains some information for gene gene with respect to old and new 
start and end positions and the length of the longest found ORF.
The number of genes with and without ORF is stated in the end. (Those are ommitted in the fasta
sdout.)
""")
parser.add_argument('--CDS', required=True, dest='CDS', type=str, help="fasta file containing cds sequences")
parser.add_argument('--log', required=True, dest='log', type=str, help="file name where intermediate stats are written to")

args = parser.parse_args()


#-------------------
# function
#-------------------


#def function_name(x,y):


####################
# main
####################

# print("biopython version", Bio.__version__)

f=open(args.CDS,"r")
f_log=open(args.log,"w")

no_ORF_samples=0
ORF_samples=0

for l in f:
    if re.match('>',l):
        # headers are expected to contain the gene name, the chr name and the start pos
        gene,chr,spos,epos,strand=l.rstrip('\n').lstrip('>').split("_")
#        print(gene)
    else:
        nuc = l.replace('\n','')
        startP = regex.compile('ATG')
        longest = (0,0,0,0,0,0)
        for m in startP.finditer(nuc, overlapped=True):
#            print("")
#            print(m)
#            print(len(Seq.Seq(nuc)[m.start():]) )
#            print(len(Seq.Seq(nuc)[m.start():]) % 3)
#            print(Seq.Seq(nuc)[m.start():])
            if len(Seq.Seq(nuc)[m.start():]) % 3 ==0:
                nuc_mod=nuc
            elif len(Seq.Seq(nuc)[m.start():]) % 3 ==1: # not a multiple of 3, remove last nt
#                nuc_mod=nuc+"NN"
                nuc_mod=nuc[:-1]
            elif len(Seq.Seq(nuc)[m.start():]) % 3 ==2: # not a multiple of 3, remove last 2 nts
#                nuc_mod=nuc+"N"
                nuc_mod=nuc[:-2]
            pro = Seq.Seq(nuc_mod)[m.start():].translate(to_stop=True)
            if len(pro) > longest[0]:
                # protein_length nt_length, spos_ORF, epos_ORF, nt_ORF, aa_seq
#                print(nuc_mod)
#                print(len(Seq.Seq(nuc)[m.start():]) % 3)
                longest = (len(pro), len(nuc), m.start(), m.start()+len(pro)*3+3, nuc_mod[m.start():m.start()+len(pro)*3+3], len(nuc_mod[m.start():m.start()+len(pro)*3+3]))
# 		print("")

        if (strand=="+"):
            new_spos=str(int(spos)+longest[2]) # old spos of the gene with offset to the ORF
            new_epos=str(int(spos)+longest[3]-1)
        else:
#            print(longest)
            new_spos=str(int(spos)+(longest[1]-longest[3])) # old spos plus difference in the ending position (as dealing with revcompl)
            new_epos=str(int(spos)+(longest[1]-longest[2]-1)) # old spos plus ORF length minus offset (as dealing with revcompl)
#            print(new_spos)
#            print(new_epos)
        if (longest[0]>0):
            print("".join([ ">",gene,"_",chr,"_",new_spos,"_",new_epos,"_",strand ]) )
            print(longest[4])
            ORF_samples=ORF_samples+1
        else:
            no_ORF_samples=no_ORF_samples+1
		
        f_log.write("".join([ ">",gene,"_",chr,"_",new_spos,"_",new_epos,"_",strand, " " ]) )
        f_log.write("".join([ "spos ",str(spos), " new_spos ",new_spos, " epos ",str(int(spos)+longest[1]-1), " new_epos ",new_epos, " len_ORF ", str(longest[5]), "\n" ]) )
f_log.write("".join(["\nNumber of genes with no ORF ",str(no_ORF_samples),"\n"]) )
f_log.write("".join(["Number of genes with ORF ",str(ORF_samples),"\n"]) )

f.close()
f_log.close()






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



""")
parser.add_argument('--fasta', required=True, dest='fasta', type=str, help="fasta file with alignmed sequences, e.g. with muscle")
args = parser.parse_args()


#-------------------
# function
#-------------------




####################
# main
####################

# print("biopython version", Bio.__version__)


f_in=open(args.fasta,"r")

header=""
nt=""
aln={}
for l in f_in:
    if re.match('>',l):
        if header!="":
            aln[header]=nt
        header=l.rstrip('\n').split("\t")[0].split(" ")[0]
        nt=""
    else:
        nt = nt+l.replace('\n','')
f_in.close()        
aln[header]=nt
#print(aln) 
ids=list(aln.keys())
ids.sort()
#print(ids)


seq1spos="" # starting position of seq1
seq2spos="" # starting position of seq2
seq1epos="" # end position of seq1
seq2epos="" # end position of seq2
#
# find out start and endpos of each sequence
for i in range(0,len(aln[ids[0]])): #go through all positions in the alignment
#    print(aln[ids[0]][i] +" "+ aln[ids[1]][i])
    if aln[ids[0]][i]!="-" and seq1spos=="": # start of respective seq (first time an nt is observed)
        seq1spos=i+1
    if aln[ids[0]][i]!="-": # any time another nt is observed the respective epos is updated
        seq1epos=i+1
        
    if aln[ids[1]][i]!="-" and seq2spos=="": # start of respective seq (first time an nt is observed)
        seq2spos=i+1
    if aln[ids[1]][i]!="-": # any time another nt is observed the respective epos is updated
        seq2epos=i+1 


overlap=0 # length where both sequences overlap
match=0 # number of of matches between both sequences over the entire alignment
mismatch=0 # number of of mismatches between both sequences over the entire alignment
seq1del=0
seq2del=0
seq1longer=0 # number of bases sequence 1 is longer that the other one
seq2longer=0 # number of bases sequence 2 is longer that the other one
#
# find out number of match and mismatches within the overlapping regions of both sequences 
for i in range(0,len(aln[ids[0]])): #go through all positions in the alignment
    inseq1=False
    inseq2=False
    
    if seq1spos <= (i+1) and (i+1) <= seq1epos: # position within seq1
        inseq1=True        
    if seq2spos <= (i+1) and (i+1) <= seq2epos: # position within seq2
        inseq2=True        
        
    if inseq1 and inseq2:
        overlap=overlap+1
        if aln[ids[0]][i] == aln[ids[1]][i]:
            match=match+1
        elif aln[ids[0]][i] != aln[ids[1]][i] and aln[ids[0]][i]=="-":
            seq1del=seq1del+1
        elif aln[ids[0]][i] != aln[ids[1]][i] and aln[ids[1]][i]=="-":
            seq2del=seq2del+1
        else:
            mismatch=mismatch+1
    elif inseq1 and not inseq2:
        seq1longer=seq1longer+1
    elif not inseq1 and inseq2:
        seq2longer=seq2longer+1
            
seq1length=seq1epos-seq1spos+1        
seq2length=seq2epos-seq2spos+1  

seq1name=ids[0].split("_")[0][1:]
seq2name=ids[1].split("_")[0][1:]
#print(seq1name)
#print(seq2name)

print("".join(map(str,["seq1\tseq2\tseq1length\tseq2length\tseq1spos\tseq2spos\tseq1epos\tseq2epos\toverlap\tmatch\tmismatch\tseq1del\tseq2del\tseq1longer\tseq2longer\tmatchpc"])))
print("\t".join(map(str,[ids[0][1:],ids[1][1:],seq1length,seq2length,seq1spos,seq2spos,seq1epos,seq2epos, overlap,match,mismatch,seq1del,seq2del, seq1longer, seq2longer, match/overlap*100])))
         


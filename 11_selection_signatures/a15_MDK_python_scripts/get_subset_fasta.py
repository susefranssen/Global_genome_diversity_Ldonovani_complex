
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

Extract fasta sequences from a fasta file that contain tags in their header.
The tags are provided in an additional input file.

""")
parser.add_argument('--fasta', required=True, dest='fasta', type=str, help="original fasta file")
parser.add_argument('--tag', required=True, dest='tag', type=str, help="file containing tags separated by rows")
parser.add_argument('--out', required=True, dest='out', type=str, help="name of the output fasta")
parser.add_argument('--sep', required=False, dest='sep', type=str, default="_", help="separator in the header, tag will only be searched in the part of the header before the first separator occurs (default '_')")
args = parser.parse_args()


#-------------------
# function
#-------------------


def print_seq(header):
    print_seq=False
    for tt in tags:
        if re.search(tt, header.split(args.sep)[0]):
            print_seq=True
#            print(header +" "+tt)

    if print_seq:
        f_out.write(header+"\n")
        f_out.write(nuc+"\n")


####################
# main
####################

# print("biopython version", Bio.__version__)


f_in=open(args.fasta,"r")
f_out=open(args.out,"w")
f_tags=open(args.tag,"r")

tags=[]

for l in f_tags:
    tags.append(l.replace('\n',''))

print(tags)
print("number of tags "+str(len(tags)))

header=""
nuc=""

for l in f_in:
    if re.match('>',l):
        if header!="":
            print_seq(header)
        header=l.rstrip('\n')
        nuc=""
    else:
        nuc = nuc+l.replace('\n','')
        
print_seq(header)

f_in.close()
f_out.close()
f_tags.close()



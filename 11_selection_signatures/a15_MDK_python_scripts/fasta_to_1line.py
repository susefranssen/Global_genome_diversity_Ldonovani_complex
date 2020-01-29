
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

Changes nt seq of each fasta entry to one line.

""")
parser.add_argument('--fasta', required=True, dest='fasta', type=str, help="original fasta file")
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
for l in f_in:
    if re.match('>',l):
        if header!="":
            print(header)
            print(nt)
        header=l.rstrip('\n')
        nt=""
    else:
        nt = nt+l.replace('\n','')
        
print(header)
print(nt)

f_in.close()



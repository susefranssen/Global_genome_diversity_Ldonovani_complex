
import argparse
from argparse import RawTextHelpFormatter
import re
import os
import math

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC




parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Adds  summary figures and results parameters to the plain results from the MK.pl script.


""")
parser.add_argument('--MKres', required=True, dest='MKres', type=str, help="Output file from the MK.pl script")

args = parser.parse_args()


#-------------------
# classes
#-------------------

class Gene:
    def __init__(self, file,Pn,Dn,Ps,Ds,MKcodons,avgsample,FETpval):
        self.file=file
        self.Pn=int(Pn)
        self.Dn=int(Dn)
        self.Ps=int(Ps)
        self.Ds=int(Ds)
        self.MKcodons=MKcodons
        self.avgsample=avgsample
        self.FETpval=FETpval
        self.NI=None
        self.set_NI() # neutraliy index
        self.neglogNI= math.log10(self.NI)*-1
        self.alpha=1-self.NI
        self.mean_alpha_counter_contr=None
        self.mean_alpha_denom_contr=None
        self.set_mean_alpha_contr()
        
    # neutraliy index
    def set_NI(self):
        self.NI=((self.Pn+1)/(self.Dn+1))/((self.Ps+1)/(self.Ds+1))
        
    # set contribtion of a gene for calc mean alpha across genes
    def set_mean_alpha_contr(self):
        # counter (Zaehler) contribution for mean alpha
        self.mean_alpha_counter_contr= (self.Ds+1) * (self.Pn+1) / ((self.Ps+1) + (self.Ds+1))
        # denominator (Nenner) contribution for mean alpha
        self.mean_alpha_denom_contr= (self.Ps+1) * (self.Dn+1) / ((self.Ps+1) + (self.Ds+1))


        
    def __str__(self):
        return "\t".join(map(str,[self.file.split("_")[0], self.Pn,self.Dn,self.Ps,self.Ds,self.MKcodons,self.avgsample, self.FETpval, self.NI, self.neglogNI, self.alpha, self.mean_alpha_counter_contr, self.mean_alpha_denom_contr]))
        

    def str_summary(self):
        return "\t".join(map(str,[self.file.split("_")[0], self.Pn,self.Dn,self.Ps,self.Ds,self.MKcodons,self.avgsample, self.FETpval, self.NI, self.neglogNI, self.alpha, self.mean_alpha_counter_contr, self.mean_alpha_denom_contr]))




#-------------------
# variables
#-------------------



#-------------------
# function
#-------------------


#def function_name(x,y):


####################
# main
####################

# print("biopython version", Bio.__version__)

f=open(args.MKres,"r")

tag=args.MKres.rsplit(".txt")[0] # tag for out files
f_table=open(tag+"_table.txt","w")
f_summary=open(tag+"_summary.txt","w")

# contained the summed up contribution of each gene to the mean alpha calculation method (Stoletzki and Eyre-Walker 2011)
mean_alpha_counter_sum=0
mean_alpha_denom_sum=0

# read in
for l in f: # for each gene in the input file
    if re.match('file',l): # header line
#        print("\t".join(["gene","Pn","Dn","Ps","Ds","MKcodons","avgsample","FETpval","NI","neglogNI","alpha","mean_alpha_counter_contr","mean_alpha_denom_contr"]))    
        f_table.write("\t".join(["gene","Pn","Dn","Ps","Ds","MKcodons","avgsample","FETpval","NI","neglogNI","alpha","mean_alpha_counter_contr","mean_alpha_denom_contr"])+"\n")
    else:
        file,Pn,Dn,Ps,Ds,MKcodons,avgsample,FETpval= l.replace('\n','').split("\t")
        gene=Gene(file,Pn,Dn,Ps,Ds,MKcodons,avgsample,FETpval)
#        print(gene)
        f_table.write(gene.str_summary()+"\n")
        mean_alpha_counter_sum=mean_alpha_counter_sum+gene.mean_alpha_counter_contr
        mean_alpha_denom_sum=mean_alpha_denom_sum+gene.mean_alpha_denom_contr

mean_alpha=1 - (mean_alpha_counter_sum/mean_alpha_denom_sum)
#print("Mean alpha: "+ str(mean_alpha))
f_summary.write("Mean alpha: "+ str(mean_alpha)+"\n")

f.close()        
f_table.close()
f_summary.close()

   

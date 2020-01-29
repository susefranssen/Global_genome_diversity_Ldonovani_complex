
import argparse
from argparse import RawTextHelpFormatter
import re
import os

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC




parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter,description="""

Note: This script is specific for the Ldon analysis as file paths and group names
for this project are hard coded in the file. :(

This script generates the input files for the McDonald-Kreitman script 
(~/software/McDonald_Kreitman_script/MK.pl).
It requires a fasta file of ORFs as input 
(obtained by ~/00scripts/general_scripts/python3/MDK-test_Ldon_A15/MDK_a_find_ORF.py)
and group information to be used as in- and out-groups for the MK test 
(several groupings of interest can be given).
In the output all sequences from the ingroup will obtain the same name tag. 
All sequences in the outgroup will be reduced to two 'unphased' haplotypes 
and will obtain two seprate name tags to be compatible with the MK.pl script.
(note: all positions polymorphic in the outgroup are ignorged in the MK script.)

This script takes a fasta files with open reading frames of genes of interest (--ORF)
in the right orientation of the gene.
The header needs to contain the following information:

Each sequence header contains four elements separates by '_'
1. gene name
2. chr name
3. start position on the ref genome of gene starting at longest ORF (in orientation of the chr in the ref fasta)
4. end position on the ref genome of gene ending at the longest ORF (in orientation of the chr in the ref fasta)
5. orientation of the gene onthe reference genome, i.e. - or +

For each gene a folder will be created for each gene with the respective gene name.
each folder will contain the folloing information:
1. vcf files for each population with genotypes for all SNPs across the populations,
    i.e. if a SNP is monomorphic in one population but polymorphic in the entire
    dataset this SNP position will be present in all population vcf files to not
    ignore SNPs that are differentially fixed between populations.
2. fasta file containing the edited ORF sequences for each sample of each population;
    For each sample two sequences are present to represent both SNPs if the sample
    is heterozygous. SNPs are not correctly phased.
    The header contains the group and the sample name + tag for 
    (unphased haplotype 'A' and 'B')
3. fasta output file for each in-/out-group pairing


""")
parser.add_argument('--ORF', required=True, dest='ORF', type=str, help="fasta file containing ORF sequences")
parser.add_argument('--ingroup', required=True, dest='ingroup', type=str, help="comma separated list of all ingroup populations; if several ingroup and outgroup combinations are given they should be separated by '_'")
parser.add_argument('--outgroup', required=True, dest='outgroup', type=str, help="comma separated list of all outgroup populations; if several ingroup and outgroup combinations are given they should be separated by '_'")

args = parser.parse_args()


#-------------------
# classes
#-------------------

class Gene:
    def __init__(self, name, chr, spos, epos, nt, strand):
        self.name=name
        self.chr=chr
        self.spos=int(spos)
        self.epos=int(epos)
        self.nt=nt
        self.strand=strand
        self.revcompl=""
        
    def __str__(self):
        return "\t".join(map(str,[self.name, self.chr, self.spos, self.epos, self.nt]))
        
    def set_revcompl(self, revcompl):
        self.revcompl=revcompl
        

class SNP:
    def __init__(self, chr, pos, ref, alt, samples, values):
        self.chr=chr
        self.pos=int(pos)
        self.ref=ref
        self.alt=alt
        self.samples=samples
        self.values=values

    def __str__(self):
        return "\t".join(map(str,[self.chr, self.pos, self.ref, self.alt, self.samples, self.values]))




#-------------------
# variables
#-------------------

# vcf files
# "/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.",chr,"_",pop,"_diploid.vcf.gz"

# includes all pops we want to get the SNP info from, are part of the file name
pops=["Ldon1_noBPK157A1_44","Ldon5_7","Ldon4_noGILANI_18","Ldon3_noLRCL53_7","Linf_noISSextr_noInf152_43","Ldon4_4","CH","TurkeyH_11"]
pops_transl={"Ldon1_noBPK157A1_44":"Ldon1", "Ldon5_7":"Ldon2", "Ldon4_noGILANI_18":"Ldon3","Ldon4_4":"Ldon4", "Ldon3_noLRCL53_7":"Ldon5", "Linf_noISSextr_noInf152_43":"Linf1", "CH":"CHLinf", "TurkeyH_11":"CUKLinf"}

ingroups=args.ingroup.split("_")
outgroups=args.outgroup.split("_")

#-------------------
# function
#-------------------


#def function_name(x,y):


####################
# main
####################

# print("biopython version", Bio.__version__)

f=open(args.ORF,"r")

# generate all "unphased haplotypes"
for l in f: # gor each gene in the input file
    if re.match('>',l):
        # headers are expected to contain the gene name, the chr name and the start pos
        gname,chr,spos,epos,strand=l.rstrip('\n').lstrip('>').split("_")
    else:
        nuc = l.replace('\n','')
        gene = Gene(gname,chr,spos,epos,nuc,strand)

#        print(gene)
        folder=gene.name#+"_tmp"
        os.system("mkdir "+folder+"_tmp/")
        
        # fasta output with edited SNPs for all samples
        f_gene=open(folder+"_tmp/"+gene.name+".fa",'w')
        
        # for all pops, extracts the part of the vcf file corresponding to the current gene
        for pop in pops:
            count=1 # count to be added to group name
           
            ORF_edits={}
            new_pop=True # going through a new population
#            print(pop)
            
            vcf_in="/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex."+gene.chr+"_"+pop+".vcf.gz"
            vcf_out=folder+"_tmp/gene_"+pop+".vcf"            
            os.system("zcat "+vcf_in+" | grep '^#CHROM' > "+vcf_out)            
            cmd="tabix "+vcf_in+" "+gene.chr+":"+str(gene.spos)+"-"+str(gene.epos)
            os.system(cmd+" >> "+vcf_out)

            vcf=open(vcf_out, mode='r')
            # go through every SNP in this gene
            SNPpres=False # indicates if a SNP is present in the respective population
            for r in vcf:
#                print(r)
                if re.match('#',r):
                    samples=r.split()[9:] # all sample names
                else:
                    SNPpres=True
                    chr,pos=r.split()[0:2]
                    ref,alt=r.split()[3:5]
                    values=r.split()[9:] # values of the respective samples
                    snp=SNP(chr,pos,ref,alt,samples,values)
                    if snp.alt.count(",") == 0: # bialleleic
                        alleles={"0":snp.ref, "1":snp.alt, ".":"N"}
                    elif snp.alt.count(",") == 1: # trialleleic
                        alt1,alt2=snp.alt.split(",")
                        alleles={"0":snp.ref, "1":alt1, "2":alt2, ".":"N"}
                    elif snp.alt.count(",") == 2: # tetra-alleleic
                        alt1,alt2,alt3=snp.alt.split(",")
                        alleles={"0":snp.ref, "1":alt1, "2":alt2, "3":alt3, ".":"N"}

#                    print(alleles)
#                    print(snp.alt)
                    # in the - strand we will modify the seq in the revcompl
                    # when all SNPs are added rewill out put the revcompl of the revcompl again to have the proper ORF orientation
                    if (gene.strand=="-"):
                        seq = Seq(gene.nt, IUPAC.unambiguous_dna)
                        gene.set_revcompl(seq.reverse_complement())
                        offset=snp.pos-gene.spos

                        # go through every sample in this group
                        for i,samp in enumerate(snp.samples):
                            if (new_pop==True): # add all ORF seq into ORF_edits
                                ORF_edits[samp+"_A"]=gene.revcompl.tomutable()
                                ORF_edits[samp+"_B"]=gene.revcompl.tomutable()
                            
#                            print(len(ORF_edits[samp+"_A"]))
#                            print(offset)
                            # index for A and B allele
                            allele_i=list(set(snp.values[i].split(":")[0].split("/")))
                            if "." in allele_i:
                                allele_i.remove(".")
                                
                            if len(allele_i)==0: # genotype unknown > only "."
                                ORF_edits[samp+"_A"][offset]="N"
                                ORF_edits[samp+"_B"][offset]="N" 
                            elif len(allele_i)==1: # homozygous genotype
                                ORF_edits[samp+"_A"][offset]=alleles[allele_i[0]]
                                ORF_edits[samp+"_B"][offset]=alleles[allele_i[0]] 
                            elif len(allele_i)==2: # heterozygous genotype
                                ORF_edits[samp+"_A"][offset]=alleles[allele_i[0]]
                                ORF_edits[samp+"_B"][offset]=alleles[allele_i[1]]
                            else:
                                print("more than 2 alleles, case to be added")
#                            if snp.alt.count(",") == 1:
#                                print(ORF_edits[samp+"_A"][offset])
#                                print(ORF_edits[samp+"_B"][offset])
#                            i_A=snp.values[i].split(":")[0].split("/")[0] 
#                            i_B=snp.values[i].split(":")[0].split("/")[1]
#                            ORF_edits[samp+"_A"][offset]=alleles[i_A]
#                            ORF_edits[samp+"_B"][offset]=alleles[i_B]
                        new_pop=False
                    else:
                        offset=snp.pos-gene.spos
#                        print(offset)
                        # go through every sample in this group
                        for i,samp in enumerate(snp.samples):
                            if (new_pop==True): # add all ORF seq into ORF_edits
                                ORF_edits[samp+"_A"]=Seq(gene.nt).tomutable()
                                ORF_edits[samp+"_B"]=Seq(gene.nt).tomutable()

#                            print(len(ORF_edits[samp+"_A"]))
#                            print(offset)
                            
                            
                            # if position to edit is in the cds
                            # of course should alway be but N have been added to found ORFs when stiop codon incomplete or missing
                            if offset<= len(ORF_edits[samp+"_A"]) :
                                # index for A and B allele
#                                print(snp.values[i])
                                allele_i=list(set(snp.values[i].split(":")[0].split("/")))
                                if "." in allele_i:
                                    allele_i.remove(".")
                                
#                                print(allele_i)
#                                print(len(allele_i))
                                if len(allele_i)==0: # genotype unknown > only "."
                                    ORF_edits[samp+"_A"][offset]="N"
                                    ORF_edits[samp+"_B"][offset]="N" 
                                elif len(allele_i)==1: # homozygous genotype
#                                    print(allele_i[0])
#                                    print(alleles)
                                    ORF_edits[samp+"_A"][offset]=alleles[allele_i[0]]
                                    ORF_edits[samp+"_B"][offset]=alleles[allele_i[0]] 
                                elif len(allele_i)==2: # heterozygous genotype
                                    ORF_edits[samp+"_A"][offset]=alleles[allele_i[0]]
                                    ORF_edits[samp+"_B"][offset]=alleles[allele_i[1]]
                                else:
                                    print("more than 2 alleles, case to be added")
#                               i_A=snp.values[i].split(":")[0].split("/")[0] 
#                               i_B=snp.values[i].split(":")[0].split("/")[1]
#                               ORF_edits[samp+"_A"][offset]=alleles[i_A]
#                               ORF_edits[samp+"_B"][offset]=alleles[i_B]
                        new_pop=False



            # print out edited ORF fasta         
            for i_samp in ORF_edits.keys():
#                out=">"+pops_transl[pop]+"_"+str(count)+"; "+i_samp+"\n"
                out=">"+pops_transl[pop]+"; "+i_samp+"\n"
                count=count+1
                f_gene.write(out)
                seq_out=Seq(str(ORF_edits[i_samp]))
                if (strand=="-"):
                    out=str(seq_out.reverse_complement())+"\n"
                else:
                    out=str(seq_out)+"\n"
                f_gene.write(out)
            if (SNPpres==False): # no SNP is present in this population, use ref sequence
#                out=">"+pops_transl[pop]+"_"+str(count)+"; consensus\n"
                out=">"+pops_transl[pop]+"; consensus\n"
                f_gene.write(out)
                out=gene.nt+"\n"
                f_gene.write(out)
                    

        # end for pop in pops
        f_gene.close()
f.close()        


   


# go through all genes again and generate fasta files to be used as input for the MK test script 
# i.e. generate all "unphased haplotypes"

MK_folder=[]
for pair in range(0,len(ingroups)): # for each in-/out-group pairing
    in_a=ingroups[pair].split(",") # ingroup array
    out_a=outgroups[pair].split(",") # outgroup array
    MK_folder.append("MK_in_"+"_".join(in_a)+"_out_"+"_".join(out_a)+"/")
    os.system("mkdir "+MK_folder[pair])


f=open(args.ORF,"r")
for l in f: # gor each gene in ORF the input file
    if re.match('>',l):
        # headers are expected to contain the gene name, the chr name and the start pos
        gname,chr,spos,epos,strand=l.rstrip('\n').lstrip('>').split("_")
    else:
#        nuc = l.replace('\n','')

#        print(gname)
        for pair in range(0,len(ingroups)): # for each in-/out-group pairing
            in_a=ingroups[pair].split(",") # ingroup array
            out_a=outgroups[pair].split(",") # outgroup array

            f_gene=open(gname+"_tmp/"+gname+".fa",'r')
            f_gene_res=open(MK_folder[pair]+gname+"_in_"+"_".join(in_a)+"_out_"+"_".join(out_a)+".fa",'w')
            
            ## also only create to haplotypes harbouring all the variation for the ingroup
            nuc_x=[] # both ingroup "unphased" sequences
            nuc_y=[] # both ingroup "unphased" sequences
            
            nuc_a=[] # both outgroup "unphased" sequences
            nuc_b=[] # both outgroup "unphased" sequences
            for r in f_gene:
                if re.match('>',r):
                    # headers are expected this format: >Ldon1; BD12_A
                    group,samp_haplo=r.lstrip(">").rstrip("\n").split("; ")
                    
                else:
                    nuc = r.replace('\n','')
                    nuc_c=[c for c in nuc] # current nuc sequence
                    if group in in_a: # all haplotype from the ingroup
                        ## also only create to haplotypes harbouring all the variation for the ingroup
                        if len(nuc_x)==0: # first outgroup seq
                            nuc_x=[c for c in nuc]
                            nuc_y=[c for c in nuc]
                        else:
                            for n in range(len(nuc_x)):
                                if nuc_c[n] != nuc_x[n] and nuc_c[n] != "N": # if curent nuc unequal to a copy
                                    nuc_y[n]=nuc_c[n]   # put other variant in b copy                        
                        ## when two haplotypes are printed for each ingroup sample
#                        out=">"+"_".join(in_a)+"; "+"_".join([group,samp_haplo])+"\n"
#                        f_gene_res.write(out)
#                        out=nuc+"\n"
#                        f_gene_res.write(out)
                    else: # outgroup "haplo"
                        if len(nuc_a)==0: # first outgroup seq
                            nuc_a=[c for c in nuc]
                            nuc_b=[c for c in nuc]
                        else:
                            for n in range(len(nuc_a)):
                                if nuc_c[n] != nuc_a[n] and nuc_c[n] != "N": # if curent nuc unequal to a copy
                                    nuc_b[n]=nuc_c[n]   # put other variant in b copy
            ## ingroup
            ing=">"+"_".join(in_a)+"; consensus\n"
            ing=ing+"".join(nuc_x)+"\n"
            f_gene_res.write(ing)
            ing=">"+"_".join(in_a)+"; consensus\n"
            ing=ing+"".join(nuc_y)+"\n"
            f_gene_res.write(ing)
            # outgroup
            out=">"+"_".join(out_a)+"_A"+"; consensus\n"
            out=out+"".join(nuc_a)+"\n"
            f_gene_res.write(out)
            out=">"+"_".join(out_a)+"_B"+"; consensus\n"
            out=out+"".join(nuc_b)+"\n"
            f_gene_res.write(out)
                        
            f_gene.close()
            f_gene_res.close() # files to be opened and closed for each pairing
            
        # end for a gene
                    
f.close()



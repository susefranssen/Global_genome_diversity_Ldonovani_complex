#!/usr/bin/env bash

# alignment was converted to phylip by
# http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_phylip.php
# fasta was taken as alignment already as the have al the same positions so should be an alignment already
# using this file to convert
# /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/05variants/variants.ldon_complex.maxi.120.final.fa



queue=normal




# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 120 samples
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T1_SNP.mask
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T1_SNP.mask
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T1_SNP.mask" 
name1=raxml_T1   
bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.phylip -n T1_SNP.mask

# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 67 samples, removing samples with identical sequences
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T2_SNP.mask.red
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T2_SNP.mask.red
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T2_SNP.mask.red" 
name2=raxml_T2   
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.phylip.reduced -n T2_SNP.mask.red




# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 120 samples
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T3_SNP.nomask
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T3_SNP.nomask
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T3_SNP.nomask" 
name3=raxml_T3   
bsub -q $queue -o ${name3}.o -e ${name3}.e -J ${name3} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.phylip -n T3_SNP.nomask

# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 79 samples, removing samples with identical sequences
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T4_SNP.nomask.red
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T4_SNP.nomask.red
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T4_SNP.nomask.red" 
name4=raxml_T4   
bsub -q $queue -o ${name4}.o -e ${name4}.e -J ${name4} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.phylip.reduced -n T4_SNP.nomask.red



# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 120 samples
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T5_SNP.nomask.indel
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T5_SNP.nomask.indel
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T5_SNP.nomask.indel" 
name5=raxml_T5   
bsub -q $queue -o ${name5}.o -e ${name5}.e -J ${name5} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.indel.phylip -n T5_SNP.nomask.indel


# rapid bootstrapping (ML search + Bootstrapping) in one single step by typing
# 100 rapid Bootstrap searches, 20 ML searches and return the best ML tree with support values
# 84 samples, removing samples with identical sequences
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/
mkdir T6_SNP.nomask.indel.red
cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/T6_SNP.nomask.indel.red
# # draw bipartitions on the best ML tree as follows: 
echo "raxml rapid bootstrapping with best tree T6_SNP.nomask.indel.red" 
name6=raxml_T6   
bsub -q $queue -o ${name6}.o -e ${name6}.e -J ${name6} -R"select[mem>4000] rusage[mem=4000]" -M4000 raxmlHPC -f a -m GTRGAMMA -p 12345 -x 12345 -# 100 -s /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.indel.phylip.reduced -n T6_SNP.nomask.indel.red




exit 0
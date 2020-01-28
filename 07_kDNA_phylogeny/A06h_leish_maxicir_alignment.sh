#!/usr/bin/env bash


queue=normal

# 120 samples masking and high coverage region, SNPs no indels
#---------- Alignement 
# concatenate all fasta seqs of all 120 samples that were chosen because they have good maxicircle coverage
# ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.fa \
# | wc -l
# 120 --> correct
ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.txt
#
cat /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.fa
#
# alignment
echo "muscle alignment" 
name=muscle_A    
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} -R"select[mem>4000] rusage[mem=4000]" -M4000 muscle -in /lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.fa -out /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.fa   


# 120 samples no masking, high coverage region, SNPs no indels
#---------- Alignement 
# concatenate all fasta seqs of all 120 samples that were chosen because they have good maxicircle coverage
# ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.mask.fa \
# | wc -l
# 120 --> correct
ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.txt
#
cat /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.fa
#
# alignment
echo "muscle alignment" 
name1=muscle_B    
bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} -R"select[mem>4000] rusage[mem=4000]" -M4000 muscle -in /lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.fa -out /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.fa  


# 120 samples no masking, high coverage region, SNPs AND indels
#---------- Alignement 
# concatenate all fasta seqs of all 120 samples that were chosen because they have good maxicircle coverage
# ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.indel.fa \
# | wc -l
# 120 --> correct
ls /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.indel.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.indel.txt
#
cat /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi*final.nomask.indel.fa > \
/lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.indel.fa
#
# alignment
echo "muscle alignment" 
name2=muscle_C    
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} -R"select[mem>4000] rusage[mem=4000]" -M4000 muscle -in /lustre/scratch118/infgen/team133/sf18/05variants/variants.ldon_complex.maxi.120.final.nomask.indel.fa -out /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A06_maxicircle/alignment.ldon_complex.maxi.120.final.nomask.indel.fa  





exit 0
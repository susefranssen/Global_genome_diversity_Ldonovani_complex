




cd /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex

for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 
do
grep -e"^CHR" -e"LinJ."$chr  snps.filt.leish_global.linj.complex.vcflike.txt > snps.filt.leish_global.linj.complex.LinJ.${chr}.vcflike.txt
gzip snps.filt.leish_global.linj.complex.LinJ.${chr}.vcflike.txt &
done



queue=normal
memory=4000 # max seen 1991 MB

# formating data and removing NA and non polymorphic SNPs
name1=R_A01_StAMPP_format_1-3
bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/A01_leish_donovaniComplex_StAMPP.R format

memory=2000
# plotting chr trees
name2=R_A01_StAMPP_chrtree
# -w "done(${name1})"
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/A01_leish_donovaniComplex_StAMPP.R chrtree dphC

memory=4000
# plotting cons trees
name3=R_A01_StAMPP_constree
# -w "done(${name1})"
bsub -q $queue -o ${name3}.o -e ${name3}.e -J ${name3} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/A01_leish_donovaniComplex_StAMPP.R constree dphC


# bootstrap on constree with windows
memory=4000
# plotting cons trees
name=R_A01_StAMPP_constree_bs102
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/A01_leish_donovaniComplex_StAMPP_window_bootstrap.R 







# grep -L "Successfully completed." *.o

cd ~/leish_donovaniComplex/A18_variation_known_candidate_loci

# check variation in previously known candidate loci


####################
# AQP1 (Imamura et al 2016, elife) 
cd ~/leish_donovaniComplex/A18_variation_known_candidate_loci/AQP1

# get AQP1 coordinates
grep '^LinJ.31' /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep 'Aquaglyceroporin' > TriTrypDB-38_LinfantumJPCM5_APQ1.gff
# LinJ.31	EuPathDB	gene	10497	11441	.	-	.	ID=LinJ.31.0030;description=Aquaglyceroporin 1


# get AQP1 genotypes
# vcflike
awk '{if ($2 >= 10497 && $2 <=11441 || $1=="CHR") print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.31.vcflike.txt > AQP1_snps.filt.leish_global.linj.complex.LinJ.31.vcflike.txt
# up and downstream 200 bp
awk '{if ($2 >= (10497-200) && $2 <=(11441+200) || $1=="CHR") print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.31.vcflike.txt > AQP1.ud200_snps.filt.leish_global.linj.complex.LinJ.31.vcflike.txt
#
# vcf
bedtools intersect \
-header \
-a /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf \
-b  TriTrypDB-38_LinfantumJPCM5_APQ1.gff | sort -u > AQP1_snps.filt.leish_global.linj.complex.diploid.vcf
# up and downstream 200 bp
awk '{if ($1=="LinJ.31" && $2 >= (10497-200) && $2 <=(11441+200) || $1 ~ /^ *#/) print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf > AQP1.ud200_snps.filt.leish_global.linj.complex.diploid.vcf

# SNPeff annotation of SNPs
queue=normal
memory=4000
name2=A18_SNPeff
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
"/software/jdk1.8.0_74/bin/java -Xmx4g -jar \
/lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.jar ann \
-ud 0 \
-c /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.config \
-i vcf \
-csvStats SNPeff_AQP1_stats -v \
TT9_Linf_JPCM5_UTR AQP1_snps.filt.leish_global.linj.complex.diploid.vcf \
> SNPeff_APQP1.vcf"
# up and downstream 200 bp
name2=A18_SNPeff_ud
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
"/software/jdk1.8.0_74/bin/java -Xmx4g -jar \
/lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.jar ann \
-ud 200 \
-c /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.config \
-i vcf \
-csvStats SNPeff_AQP1.ud200_stats -v \
TT9_Linf_JPCM5_UTR AQP1.ud200_snps.filt.leish_global.linj.complex.diploid.vcf \
> SNPeff_APQP1.ud200.vcf"


#####
#APQP1
# extract gene fasta
bedtools getfasta -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa -bed TriTrypDB-38_LinfantumJPCM5_APQ1.gff -fo TriTrypDB-38_LinfantumJPCM5_APQ1_stranded.fasta -s
#
bedtools getfasta -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa -bed TriTrypDB-38_LinfantumJPCM5_APQ1.gff -fo TriTrypDB-38_LinfantumJPCM5_APQ1_ignoreStrand.fasta



#####
# Miltefosine transporter LdMT
cd ~/leish_donovaniComplex/A18_variation_known_candidate_loci/LdMT
#
# get LdMT coordinates
grep '^LinJ.13' /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep 'Miltefosine transporter' | grep 'gene' > TriTrypDB-38_LinfantumJPCM5_LdMT.gff
#
# extract gene fasta
bedtools getfasta -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa -bed TriTrypDB-38_LinfantumJPCM5_LdMT.gff -fo TriTrypDB-38_LinfantumJPCM5_LdMT_stranded.fasta -s
#
bedtools getfasta -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa -bed TriTrypDB-38_LinfantumJPCM5_LdMT.gff -fo TriTrypDB-38_LinfantumJPCM5_LdMT_ignoreStrand.fasta
#
# SNPs putatively involved in Miltefosine resistance
#
# LdMT gene length 3,294 bp
# genome pos 617319 - 620612
# gene on - strand
#
# E197D
#    gene pos nt (-strand)						589(G) 590(A) 591(G) (checked in jalview)
# to gene pos nt (+strand) (gene length - pos + 1) 2704(C) 2705(T) 2706(C)
# to genome pos (gene spos + gene pos nt -1)	620022(C) 620023(T) 620024(C)
# 
# E encoded by GAG or GAA
# D encoded by GAT and GAC
# so variants leading to E197D can have either of the genomic position SNP variants:
# genome pos	620022(A) 620023(T) 620024(C)
# genome pos	620022(G) 620023(T) 620024(C)
# --> look for SNPs at chr 13 pos 620022 (to A or G)
awk '{if ($1=="LinJ.13" && $2 >= (620022) && $2 <=(620024) || $1 ~ /^ *#/) print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf 
# --> there was no SNP called
#
# A691P
#    gene pos nt (-strand)						2071(G) 2072(C) 2073(T) (jalview)
# to gene pos nt (+strand) (gene length - pos + 1) 1222(A) 1223(G) 1224(C)
# to genome pos (gene spos + gene pos nt -1)	618540(A) 618541(G) 618542(C)
# A encoded by GCT, GCC, GCA, GCG
# P encoded by CCT, CCC, CCA, CCG
# so variants leading to A691P can have either of the genomic position SNP variants:
# genome pos	618542(G)
awk '{if ($1=="LinJ.13" && $2 >= (618540) && $2 <=(618542) || $1 ~ /^ *#/) print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf 
# --> there was no SNP called
#
awk '{if ($1=="LinJ.13" && $2 >= (617319) && $2 <=(620612) || $1 ~ /^ *#/) print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf > LdMT_snps.filt.leish_global.linj.complex.diploid.vcf



# run R script for gene CNV in known candidate genes
#
queue=normal
memory=17000

name=R_Hlocus
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
/software/R-3.5.0/bin/Rscript ~/00scripts/leish_donovaniComplex/a18_variation_known_candidate_loci.r 





# compare breakpoints of MSL with those from Carnielli paper
bwa mem /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta MSL_breakpoints_Juliana.fa > MSL_breakpoints_Juliana.sam

# offset for MSL, DR and IR location that Juliana sent for v5 is 262
# 1135445-1135183
# [1] 262
# 1122737-1122475
# [1] 262


bwa mem /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta MSL_large_region_v38_dotplot_1.fa > MSL_large_region_v38_dotplot_1.sam





# check if AQP1 gene is in a region with good mappability

grep -e "LinJ.31.0030" /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep gene
# LinJ.31	EuPathDB	gene	10497	11441	.	-	.	ID=LinJ.31.0030;description=Aquaglyceroporin 1

grep -e "LinJ.31.0010" /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep gene
# LinJ.31	EuPathDB	gene	6025	8337	.	-	.	ID=LinJ.31.0010;description=5-methyltetrahydropteroyltriglutamate-homocysteine S-methyltransferase%2C putative

grep -e "LinJ.31.0040" /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep gene
# LinJ.31	EuPathDB	gene	12415	15573	.	-	.	ID=LinJ.31.0040;description=RING-variant domain containing protein%2C putative

awk '{if ($1=="LinJ.31" && $3<20000) print $0}' ~/leish_donovaniComplex/A20_masked_ref_stats/mask.LinJ.old.match95.bed
LinJ.31	1	2510
LinJ.31	2524	2624
LinJ.31	2650	3447
LinJ.31	3469	3865
LinJ.31	4047	4530
LinJ.31	4532	4672
LinJ.31	4734	4952
LinJ.31	4953	5093
# 6025	8337 adjacent closer to telomere
LinJ.31	8447	8548
# 10497	11441 AQP1
# 12415	15573 further from telomere
LinJ.31	#13181	13281
LinJ.31	16595	16906
LinJ.31	17509	17609
LinJ.31	19143	19243
LinJ.31	19668	19768
LinJ.31	19776	19876
LinJ.31	19898	19998

# check intersection of the AQP1 gene and the two adjacent genes with masked regions
grep -e "LinJ.31.0030" -e "LinJ.31.0010" -e "LinJ.31.0040" /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | grep gene > LinJ.31.0030_AQP1_plus.gff 
#
bedtools intersect -wa -wb -a LinJ.31.0030_AQP1_plus.gff  -b ~/leish_donovaniComplex/A20_masked_ref_stats/mask.LinJ.old.match95.bed

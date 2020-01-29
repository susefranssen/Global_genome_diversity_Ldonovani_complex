
# grep -L "Successfully completed." *.o


cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats






###############################################
#
# A doing SNPeff for subsets of fixed SNPs (markers) of interest
#
# 1) subsetting input vcf
for marker in Linf_vs_Ldon CH.Linf.nohy_3 CUK.Linf_11 Ldon1star_44 Ldon2_7 Ldon3star_18 Ldon4_4 Ldon5star_7 Linf1star_43 Linf_vs_Ldon
do
awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print f1[$1$2]}' \
~/leish_donovaniComplex/A02_SNPeff/all_SNPs_vcf/allSNPs_NA.frac0.2_poly.vcf \
~/leish_donovaniComplex/A12_group_freqs/dat.mainX.fixed_${marker}.txt \
> allSNPs_NA.frac0.2_poly_${marker}.vcf
done
#
#
#
# 2)do SNPeff annot for marker sets of interest
#
# 2a) add additional annotation to for SNPeff to use
# add additional annotation
# 
# # add 
# # Leishmania infantum JPCM5_TT38
# TT38_Linf_JPCM5.genome : Leishmania
# to config file: /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.config 
cd /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/data/
mkdir TT38_Linf_JPCM5
# insert gff.gz file and rename to genes.gff.gz
cp /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff TT38_Linf_JPCM5/genes.gff
gzip TT38_Linf_JPCM5/genes.gff
cd ..
# add corresponding reference sequence into folder
cp /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/data/genomes/TT38_Linf_JPCM5.fa
#
# build new reference
queue=yesterday
memory=4000
# 
name1=A02_SNPeff_build
bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
"/software/jdk1.8.0_74/bin/java -Xmx4g -jar \
snpEff.jar build -gff3 -v TT38_Linf_JPCM5"
# 
# 2b) calculate SNPeff annotation for the different marker group marker sets
queue=normal
memory=4000
cd ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/#
for marker in  CH.Linf.nohy_3 CUK.Linf_11 Ldon1star_44 Ldon2_7 Ldon3star_18 Ldon4_4 Ldon5star_7 Linf1star_43 Linf_vs_Ldon
do
# marker=Linf_vs_Ldon
cd ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/
mkdir SNPeff_${marker}
name2=A13_SNPeff_${marker}
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
"/software/jdk1.8.0_74/bin/java -Xmx4g -jar \
/lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.jar ann \
-ud 0 \
-c /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.config \
-i vcf \
-csvStats SNPeff_${marker}/SNPeff_${marker}_stats -v \
TT38_Linf_JPCM5 /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats/allSNPs_NA.frac0.2_poly_${marker}.vcf \
> SNPeff_${marker}/SNPeff_${marker}.vcf"
done
# 
#
# 
# 3) doing GO enrichment of genes that have at least one moderate of high effect variant for the different marker SNP sets each 
cd ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/
#
# mapping info from SNPeff type to impact
grep -v '^#' ~/leish_donovaniComplex/A02_SNPeff/res_all_SNPs_vcf/allSNPs_NA.frac0.2_poly_annTT9_UTR.vcf | awk -F '[\t|]' '{print $9"\t"$10}' | sort -u > SNPeff_types_impacts.txt
#
# formatting on GO enrichment
for marker in  CH.Linf.nohy_3 CUK.Linf_11 Ldon1star_44 Ldon2_7 Ldon3star_18 Ldon4_4 Ldon5star_7 Linf1star_43 Linf_vs_Ldon
do
cat SNPeff_${marker}/SNPeff_${marker}_stats.genes.txt | \
awk '$5!="0" || $7!="0" {print $1}' | grep -v '#'  > SNPeff_${marker}/SNPeff_${marker}_stats.genes_high_mod.id
# awk '$5!="0" || $7!="0" {print $1}' | grep -v '#' | awk -F '[_-]' '{print $2}' > SNPeff_${marker}/SNPeff_${marker}_stats.genes_high_mod.id
done


# gene enrichemnt done from http://tritrypdb.org/tritrypdb/showQuestion.do?questionFullName=GeneQuestions.GeneByLocusTag
# including computational annotation


#-------------------
# test genes fixed for moderate to high effect SNPs between both species that are responsible for GO enrichment "pathogenesis"
#
grep -e 'LinJ.04.0710' -e 'LinJ.34.2950' -e 'LinJ.13.0930' -e 'LinJ.13.0940' -e 'LinJ.35.4250'  /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff |\
grep 'gene' > ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/Linf_Ldon_GOpathogenesis_TrypDB-38_LinfantumJPCM5.gff
#
# get gene sequences from the v38 version (will be compared to the cds of the new version v41)
/software/bedtools-2.22.0/bin/bedtools getfasta \
-fi /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta \
-bed ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/Linf_Ldon_GOpathogenesis_TrypDB-38_LinfantumJPCM5.gff \
-fo ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/Linf_Ldon_GOpathogenesis_TriTrypDB-38_LinfantumJPCM5_TriTrypDB-38_LinfantumJPCM5_Genome.fa
#
# # get gene sequences from the Leishmania_infantum_JPCM5_27July2011 version (will be compared to the cds of the new version v41)
/software/bedtools-2.22.0/bin/bedtools getfasta \
-fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa \
-bed ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/Linf_Ldon_GOpathogenesis_TrypDB-38_LinfantumJPCM5.gff \
-fo ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/Linf_Ldon_GOpathogenesis_TriTrypDB-38_LinfantumJPCM5_Leishmania_infantum_JPCM5_27July2011.fa








cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats


###############################################
#
# B doing SNPeff polymorphic loci among groups
# 
#
# 1) subsetting input vcf
awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print f1[$1$2]}' \
~/leish_donovaniComplex/A02_SNPeff/all_SNPs_vcf/allSNPs_NA.frac0.2_poly.vcf \
~/leish_donovaniComplex/A12_group_freqs/Segr_main8_maxshared.txt \
> allSNPs_NA.frac0.2_poly_Segr_main8_maxshared.vcf
#
#
#
# 2)do SNPeff annot for marker sets of interest
#
# 2a) add additional annotation to for SNPeff to use
# ... done above already
# 
# 2b) calculate SNPeff annotation for the different marker group marker sets
queue=normal
memory=4000
cd ~/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats_v38/
mkdir SNPeff_Segr_main8_maxshared
name2=A13_SNPeff_Segr_main8_maxshared
bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
"/software/jdk1.8.0_74/bin/java -Xmx4g -jar \
/lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.jar ann \
-ud 0 \
-c /lustre/scratch118/infgen/team133/sf18/software/snpEff_4.2/snpEff.config \
-i vcf \
-csvStats SNPeff_Segr_main8_maxshared/SNPeff_Segr_main8_maxshared_stats -v \
TT38_Linf_JPCM5 /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A13_SNPeff_fixed_segr_cats/allSNPs_NA.frac0.2_poly_Segr_main8_maxshared.vcf \
> SNPeff_Segr_main8_maxshared/SNPeff_Segr_main8_maxshared.vcf"
#
#
# 3) combine frequency and SNPeff annotation data
awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print f1[$1$2], $8}' \
~/leish_donovaniComplex/A12_group_freqs/Segr_main8_maxshared.txt \
SNPeff_Segr_main8_maxshared/SNPeff_Segr_main8_maxshared.vcf \
> SNPeff_Segr_main8_maxshared/SNPeff_Segr_main8_maxshared_group_freqs.vcf







# grep -L "Successfully completed." *.o


cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38

queue=normal


# get gff for genes only
awk '{if($3=="gene") print $0}' /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff | sort | grep -v 'LinJ.00' > TriTrypDB-38_LinfantumJPCM5.genes.gff
# change to bed format
awk -F '[\t=;]' '{print $1,$4-1,$5,$10}' TriTrypDB-38_LinfantumJPCM5.genes.gff > TriTrypDB-38_LinfantumJPCM5.genes.bed


sample=EP


for sample in EP Peking D_2 STRAIN_B STRAIN_A RACOON_DOG SKIN DOG_STRAIN IECDR1 BD09 BD12 BD14 BD15 BD17 BD21 BD22 BD24 BD25 BD27 Ldon282cl2 BPK029A1 BPK035A1 BPK067A1 BPK077A1 BPK156A1 BPK157A1 BPK164A1 BPK282I9 BPK294A1 BPK406A1 BPK413A1 BPK471A1 BPK512A1 BPK562A1 BPK612A1 BPK623A1 BPK648A1 BPK649A1 L60b OVN3 CL-SL LRC-L51p Chowd5 STL2-78 STL2-79 DD8 Nandi AG83 BHU41 BHU220A1 BHU1062-4 Don201 BHU816A1 BHU824A1 BHU931A1 BHU1064A1 BHU1065A1 BHU1137A1 BHU1139A1 LRC-L1311 LRC-L1313 Inf206 BUMM3 SUKKAR2 Don081 LdonLV9 GEBRE1 Don038 356WTV AM560WTI 363SKWTI 383WTI 364SPWTII AM563WTI LRC-L53 LRC-L57 MRC74 NLB-323 LRC-L445 LRC-L61 1S GILANI GE 38-UMK 45-UMK 452BM 597-2 597LN 762L LEM3472 1026-8 855-9 SUDAN1 Inf007 LRC-L699 LRC-L740 LRC-L1275 TH4 TH5 TH6 NT10 NT16 LRC-L1312 LRC-L1296 LRC_L1303 CH32 CH33 CH34 CH35 CH36 Inf152 CUK2 CUK3 CUK4 CUK5 CUK6 CUK7 CUK8 CUK9 CUK10 CUK11 CUK12 ISS174 ISS2420 ISS2426 ISS2429 ISS2508 Malta33 BUCK Inf001 LRC_L47 RM1 LEM1985 LPN114 LEM3278 Inf045 Inf004 Inf055 BCN83 BCN87 LinJPCM5 IMT260 IMT373cl1 ITMAP26 Cha001 MAM ARL WC WR285 HN167 HN336

do


# queue=yesterday
queue=normal

# mq20PP: -q 20 -f 0x0002 -F 0x0004 -F 0x0008 
# when this data was used regions with e.g. tandem repeated genes were removed and caused problems with gene coverage, therefore excluding the mq20 criterion
# 
# 1024: read is PCR or optical duplicate
# 0x0002: proper pair
# 0x0004: read unmapped
# 0x0008: mate unmapped
#
name0=a14_01a_rmdup_${sample} # map -> filter out unmapped reads
bsub  -q $queue -o ${name0}.o -e ${name0}.e -J ${name0} \
-R"select[mem>1000] rusage[mem=1000]" -M1000 \
"samtools-1.3 view -b -F 1024 -f 0x0002 -F 0x0004 -F 0x0008 /lustre/scratch118/infgen/team133/cd16/BAMfiles/leish_global/merged/linj.${sample}.sorted.markdup.realigned.bam > /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex_rmdup/linj.${sample}.sorted.markdup.realigned.PP.rmdup.bam"


name0a=index_${sample} 
bsub -w "done(${name0})" -q $queue -o ${name0a}.o -e ${name0a}.e -J ${name0a} \
-R"select[mem>1000] rusage[mem=1000]" -M1000 \
samtools-1.3 index /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex_rmdup/linj.${sample}.sorted.markdup.realigned.PP.rmdup.bam

# -d	Report the depth at each position in each B feature.
# -split	Treat "split" BAM or BED12 entries as distinct BED intervals.
name1=a14_02_bed_gene_cov_${sample}
# bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
bsub -w "done(${name0a})" -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
-R"select[mem>5000] rusage[mem=5000]" -M5000 \
"bedtools coverage -d -split -abam /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex_rmdup/linj.${sample}.sorted.markdup.realigned.PP.rmdup.bam -b TriTrypDB-38_LinfantumJPCM5.genes.gff > /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/linj.${sample}.sorted.markdup.realigned.PP.rmdup.gene.cov"

name2=a14_03_gene_cov_red_${sample}
# bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
bsub -w "done(${name1})" -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>1000] rusage[mem=1000]" -M1000 \
bash ~/00scripts/leish_donovaniComplex/a14_gene_cov_line.sh $sample


done






memory=7000

# run r script on server
name=R_A14
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a14_gene_cov.r



####
#
# genome-wide coverage
#

for sample in Peking D_2 STRAIN_B STRAIN_A RACOON_DOG SKIN DOG_STRAIN IECDR1 BD09 BD12 BD14 BD15 BD17 BD21 BD22 BD24 BD25 BD27 Ldon282cl2 BPK029A1 BPK035A1 BPK067A1 BPK077A1 BPK156A1 BPK157A1 BPK164A1 BPK282I9 BPK294A1 BPK406A1 BPK413A1 BPK471A1 BPK512A1 BPK562A1 BPK612A1 BPK623A1 BPK648A1 BPK649A1 L60b OVN3 CL-SL LRC-L51p Chowd5 STL2-78 STL2-79 DD8 Nandi AG83 BHU41 BHU220A1 BHU1062-4 Don201 BHU816A1 BHU824A1 BHU931A1 BHU1064A1 BHU1065A1 BHU1137A1 BHU1139A1 LRC-L1311 LRC-L1313 Inf206 BUMM3 SUKKAR2 Don081 LdonLV9 GEBRE1 Don038 356WTV AM560WTI 363SKWTI 383WTI 364SPWTII AM563WTI LRC-L53 LRC-L57 MRC74 NLB-323 LRC-L445 LRC-L61 1S GILANI GE 38-UMK 45-UMK 452BM 597-2 597LN 762L LEM3472 1026-8 855-9 SUDAN1 Inf007 LRC-L699 LRC-L740 LRC-L1275 TH4 TH5 TH6 NT10 NT16 LRC-L1312 LRC-L1296 LRC_L1303 CH32 CH33 CH34 CH35 CH36 EP Inf152 CUK2 CUK3 CUK4 CUK5 CUK6 CUK7 CUK8 CUK9 CUK10 CUK11 CUK12 ISS174 ISS2420 ISS2426 ISS2429 ISS2508 Malta33 BUCK Inf001 LRC_L47 RM1 LEM1985 LPN114 LEM3278 Inf045 Inf004 Inf055 BCN83 BCN87 LinJPCM5 IMT260 IMT373cl1 ITMAP26 Cha001 MAM ARL WC WR285 HN167 HN336

do

# sample=EP
# queue=yesterday
queue=normal

name1=a14_bed_genome_cov_${sample}
bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
-R"select[mem>5000] rusage[mem=5000]" -M5000 \
"bedtools genomecov -d -split -ibam /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex_rmdup/linj.${sample}.sorted.markdup.realigned.PP.rmdup.bam > /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/genomecov/linj.${sample}.sorted.markdup.realigned.PP.rmdup.genome.cov"



done

queue=yesterday
memory=7000

# run r script on server
name=R_A14_genome1
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a14_genome_cov.r




# inspect some indels of interest
cd ~/leish_donovaniComplex/A14_gene_cov/genomecov

bedtools intersect -a /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff -b Chr35_largestIndels_insertion_spos368.5_epos414.5.bed | grep 'gene' > Chr35_largestIndels_insertion_spos368.5_epos414.5_genes.gff


# run GO enrichment analysis for all indels shared between both species that have been marked with blue / red in dir ~/leish_donovaniComplex/A14_gene_cov/genomecov/indels_both_species
# these are indels in the middle of a chromosome and supported by several samples
for ID in 3 12 18 33 48 54 58 63 65 83 102 118 150 161 174 185 198 208 215 221 232 233
for ID in 168 169 170 171
do

echo $ID
grep 'id'$ID'_' indels_in_both_species.bed > 'id'$ID'.bed'

python3 ~/00scripts/general_scripts/python3/topGO_Fisher/topGO_Fisher.py \
--bed 'id'$ID'.bed' \
--gff /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-38_LinfantumJPCM5.gff \
--dir 'id'$ID \
--spath ~/00scripts/general_scripts/python3/topGO_Fisher/ \
--cpath ~/leish_donovaniComplex/A14_gene_cov/genomecov

done





bsub -R"select[mem>5000] rusage[mem=5000]" -M5000 -I java -jar /software/pathogen/external/apps/usr/local/IGV-2.3.8/igv.jar



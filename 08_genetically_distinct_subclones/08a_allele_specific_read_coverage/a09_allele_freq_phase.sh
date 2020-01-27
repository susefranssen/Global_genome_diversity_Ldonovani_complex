
# grep -L "Successfully completed." *.o


cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex

mkdir A09_allele_freq_phase

queue=normal


sample=EP

for sample in Peking D_2 STRAIN_B STRAIN_A RACOON_DOG SKIN DOG_STRAIN IECDR1 BD09 BD12 BD14 BD15 BD17 BD21 BD22 BD24 BD25 BD27 Ldon282cl2 BPK029A1 BPK035A1 BPK067A1 BPK077A1 BPK156A1 BPK157A1 BPK164A1 BPK282I9 BPK294A1 BPK406A1 BPK413A1 BPK471A1 BPK512A1 BPK562A1 BPK612A1 BPK623A1 BPK648A1 BPK649A1 L60b OVN3 CL-SL LRC-L51p Chowd5 STL2-78 STL2-79 DD8 Nandi AG83 BHU41 BHU220A1 BHU1062-4 Don201 BHU816A1 BHU824A1 BHU931A1 BHU1064A1 BHU1065A1 BHU1137A1 BHU1139A1 LRC-L1311 LRC-L1313 Inf206 BUMM3 SUKKAR2 Don081 LdonLV9 GEBRE1 Don038 356WTV AM560WTI 363SKWTI 383WTI 364SPWTII AM563WTI LRC-L53 LRC-L57 MRC74 NLB-323 LRC-L445 LRC-L61 1S GILANI GE 38-UMK 45-UMK 452BM 597-2 597LN 762L LEM3472 1026-8 855-9 SUDAN1 Inf007 LRC-L699 LRC-L740 LRC-L1275 TH4 TH5 TH6 NT10 NT16 LRC-L1312 LRC-L1296 LRC_L1303 CH32 CH33 CH34 CH35 CH36 EP Inf152 CUK2 CUK3 CUK4 CUK5 CUK6 CUK7 CUK8 CUK9 CUK10 CUK11 CUK12 ISS174 ISS2420 ISS2426 ISS2429 ISS2508 Malta33 BUCK Inf001 LRC_L47 RM1 LEM1985 LPN114 LEM3278 Inf045 Inf004 Inf055 BCN83 BCN87 LinJPCM5 IMT260 IMT373cl1 ITMAP26 Cha001 MAM ARL WC WR285 HN167 HN336
do

queue=yesterday

# nomq20PP: -b
# mq20PP: -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b
name0=a09_00_view_${sample}
bsub  -q $queue -o ${name0}.o -e ${name0}.e -J ${name0} \
-R"select[mem>1000] rusage[mem=1000]" -M1000 "samtools-1.3 view -q 20 -f 0x0002 -F 0x0004 -F 0x0008 -b /lustre/scratch118/infgen/team133/cd16/BAMfiles/leish_global/merged/linj.${sample}.sorted.markdup.realigned.bam > /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.bam"


name1=a09_01_index_${sample}
#bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
bsub -w "done(${name0})" -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
-R"select[mem>1000] rusage[mem=1000]" -M1000 "samtools-1.3 index /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.bam"

queue=normal

name2=a09_02_quality_${sample} # test for mapping quality
#bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
bsub -w "done(${name1})" -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} \
-R"select[mem>4000] rusage[mem=4000]" -M4000 "/software/bin/java -Xmx4g -jar ~/software/Readtools/ReadTools.jar QualityEncodingDetector -I /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.bam > /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.bam.qual"


name5=a09_05_pileup_${sample} # NOTE: marked duplicates are ignored by mpileup, so here
#bsub -q $queue -o ${name5}.o -e ${name5}.e -J ${name5} \
bsub -w "done(${name2})" -q $queue -o ${name5}.o -e ${name5}.e -J ${name5} \
-R"select[mem>4000] rusage[mem=4000]" -M4000 \
"samtools-1.3 mpileup -d 3500 -B -Q 10 -f /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011.fa -o /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.pileup /lustre/scratch118/infgen/team133/sf18/02X_BAM/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.bam "


name6=a09_06_sort_${sample}
#bsub -q $queue -o ${name6}.o -e ${name6}.e -J ${name6} \
bsub -w "done(${name5})" -q $queue -o ${name6}.o -e ${name6}.e -J ${name6} \
-R"select[mem>4000] rusage[mem=4000]" -M4000 \
"sort -k1,1 -k2n,2 /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.pileup > /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.sort.pileup"


name7=a09_07_pileup_SNP_${sample}
#bsub -q $queue -o ${name7}.o -e ${name7}.e -J ${name7} \
bsub -w "done(${name6})" -q $queue -o ${name7}.o -e ${name7}.e -J ${name7} \
-R"select[mem>4000] rusage[mem=4000]" -M4000 \
bash ~/00scripts/leish_donovaniComplex/a09_allele_freq_phase_07line.sh ${sample}




name8=a09_08_sync_${sample}
#bsub -q $queue -o ${name8}.o -e ${name8}.e -J ${name8} \
bsub -w "done(${name7})" -q $queue -o ${name8}.o -e ${name8}.e -J ${name8} \
-R"select[mem>5000] rusage[mem=5000]" -M5000 "/software/bin/java -Xmx5g -jar ~/software/popoolation2_1201/mpileup2sync.jar --input /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.sort.SNPs.pileup --output /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.sort.SNPs.pileup.sync --fastq-type sanger --min-qual 20 --threads 1"

done

# name9=a09_09_sync_${sample}
# #bsub -q $queue -o ${name9}.o -e ${name9}.e -J ${name9} \
# bsub -w "done(${name8})" -q $queue -o ${name9}.o -e ${name9}.e -J ${name9} \
# -R"select[mem>500] rusage[mem=500]" -M500 rm /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.*pileup


#######
# check some positions

# MAM
cat 06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike.txt | cut -f1,2,154 | grep 'HET' |head
# LinJ.01	753	0/1:TC:HET:2:63:99
# LinJ.01	1081	0/1:AG:HET:2:56:99
# LinJ.01	1405	0/1:TG:HET:2:68:89
# LinJ.01	2182	0/1:GA:HET:2:62:99
# LinJ.01	3406	0/1:TC:HET:2:49:99
# LinJ.01	5623	0/1:AG:HET:2:51:99
# LinJ.01	5659	0/1:TC:HET:2:43:77
# LinJ.01	6999	0/1:CG:HET:2:77:99
# LinJ.01	7003	0/1:CT:HET:2:77:99
# LinJ.01	8487	0/1:AC:HET:2:67:99
#
# dat[pos %in% c(753,1081,1405,2182,3406,5623,5659,6999,7003,8487)]
#         chr  pos ref  a  t  c  g sum ref.count ref.freq
#  1: LinJ.01  753   t  0 52 11  0  63        52     0.83
#  2: LinJ.01 1081   a 48  0  0  5  53        48     0.91
#  3: LinJ.01 1405   t  0 62  0  7  69        62     0.90
#  4: LinJ.01 2182   g  7  0  0 53  60        53     0.88
#  5: LinJ.01 3406   t  0 42  6  0  48        42     0.88
#  6: LinJ.01 5623   a 41  0  0  6  47        41     0.87
#  7: LinJ.01 5659   t  0 36  4  0  40        36     0.90
#  8: LinJ.01 6999   c  0  0 65 11  76        65     0.86
#  9: LinJ.01 7003   c  0 10 65  0  75        65     0.87
# 10: LinJ.01 8487   a 57  0  7  0  64        57     0.89
# 11: LinJ.05 5659   c  0  5 53  0  58        53     0.91
# 12: LinJ.10 3406   c  0  1 79  0  80        79     0.99
# 13: LinJ.20 6999   c  0  2 74  0  76        74     0.97



# was just run like this for speed but has to be adopted within the script
#
# queue=normal
# memory=4000
# 
# for sample in MAM EP GE CH34 GILANI 762L GEBRE1 CH32 45-UMK LRC-L61 452BM 1S 364SPWTII 363SKWTI LEM3472 383WTI 38-UMK 1026-8 855-9 SUDAN1 Inf152 Malta33 LRC-L740 597LN 597-2 SUKKAR2 LRC-L53 BUMM3 Inf055 CUK10 CUK9 356WTV CUK2 ISS2429 CUK4 ISS2426 CUK6 CUK7 CUK8 CUK11 LdonLV9 CUK5 CUK12 ISS174 BPK157A1 CUK3
# do
# 
# name=het_$sample
# bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
# -R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
# Rscript ~/00scripts/leish_donovaniComplex/a09_allele_freq_phase.r $sample
# 
# done


name=het
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a09_allele_freq_phase.r


name=tree
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/A01_leish_donovaniComplex_StAMPP_phylo_addhaplo.R

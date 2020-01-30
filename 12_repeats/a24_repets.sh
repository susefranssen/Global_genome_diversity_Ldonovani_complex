

#!/usr/bin/env bash

# grep -L "Successfully completed." *.o



module load bedtools/2.29.0--hc088bd4_3
module load samtools/1.9--h91753b0_8
module load mummer4/4.0.0beta2--pl526hf484d3e_4
# module load bwa/0.7.17=pl5.22.0_2

cd ~/leish_donovaniComplex/a24_repeats

###############################
#
# get bed file for repetive sequences from Ubeda et al. 2014
#
cat Ubeda_2014_ST3.txt | grep -v '^id' |\
awk 'BEGIN{OFS="\t"}{print $2,$3-1,$4-1,"RAG"$10"_seq1_SIDER"$11"_sposSeq1-"$3"_id"$1,$7,$3-1,$4-1}' > Ubeda_v3_seq1.bed
#
cat Ubeda_2014_ST3.txt | grep -v '^id' |\
awk 'BEGIN{OFS="\t"}{print $2,$5-1,$6-1,"RAG"$10"_seq2_SIDER"$11"_sposSeq1-"$3"_id"$1,$7,$5-1,$6-1}' > Ubeda_v3_seq2.bed
#
# combine and get each repeat (per position) only once
cat Ubeda_v3_seq1.bed Ubeda_v3_seq2.bed | sort -k1,1 -k2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"NA","NA",$6,$7}' | sort -k1,1 -k2n -u > Ubeda_v3_uniq.bed
#
# extract repeat sequences based on coordinates for repeats from Ubeda et al. 2014 and the genome version of TriTryp (no changes of chr sequences between version v1-38, change of additional contigs in v25, new assembly from version v39)
#
bedtools getfasta -s -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-9.0_LinfantumJPCM5_Genome.fasta \
-bed Ubeda_v3_uniq.bed -fo Ubeda_v3_uniq_TriTrypDB-9.0.fa 2> log_bedtools_TriTrypDB-9.txt
#
bedtools getfasta -s -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/Leishmania_infantum_JPCM5_27July2011_cap_sort_mainchr.fa \
-bed Ubeda_v3_uniq.bed -fo Ubeda_v3_uniq_27July2011.fa 2> log_bedtools_27July2011.txt
#
bedtools getfasta -s -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta \
-bed Ubeda_v3_uniq.bed -fo Ubeda_v3_uniq_TriTrypDB-38.fa 2> log_bedtools_TriTrypDB-38.txt
# NOTE: several sequences were not present in the chromosomes any more and are consequently not included



###############################
#
# use reference obtained from Frédéric Raymond <frederic.raymond@fsaa.ulaval.ca>, Ubeda et al. 2014 (as version v3 from GeneDB is not available any longer)
#
# these should now give of the correct sequences of the repeats as described in Ubeda et al. 2014
#
cat reference_genomes/20080508_RepeatAnalysisGenomes/LinJ*.fasta | sed 's/>LinJ/\n>LinJ/' | sed 's/LinJ/LinJ./' | sed 's/_20070420_V3.artemis//' | sed 's/_20070420_V3.0.artemis//' | fold > reference_genomes/20080508_RepeatAnalysisGenomes_format.fa
#
cat ~/leish_donovaniComplex/a24_repeats/reference_genomes/20080508_RepeatAnalysisGenomes_format.fa | sed 's/LinJ/Ubeda_LinJ/' | grep -e Lin -e a -e c -e g -e t > ~/leish_donovaniComplex/a24_repeats/reference_genomes/20080508_RepeatAnalysisGenomes_format_header.fa
#
samtools faidx reference_genomes/20080508_RepeatAnalysisGenomes_format.fa
# extract repeat sequences
bedtools getfasta -s -fi reference_genomes/20080508_RepeatAnalysisGenomes_format.fa -bed Ubeda_v3_uniq.bed -fo Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.fa 2> log_bedtools_20080508_RepeatAnalysisGenomes_format.txt
#
# map repeat based on Ubeda ref to my genome version
# bwa mem /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.fa -o Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.sam
# #
# samtools view -b Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.sam > Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.bam
# samtools sort Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.bam -o Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format_sort.bam
# samtools index Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format_sort.bam


# chr 27 of the reference used in the current study
samtools faidx /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta LinJ.27 > TriTrypDB-38_LinfantumJPCM5_Genome_LinJ.27.fasta





###############################
#
# nucmer comparison between genomes and to locate repeats based on older genome reference
# 
# current genome version used here TriTryp v38
#
# genome TT38 vs Ubeda -- whole genome comparison between current version and Ubeda
ref_qry=genomeTT38_genomeUbeda
ref_fasta=/lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta
qry_fasta=~/leish_donovaniComplex/a24_repeats/reference_genomes/20080508_RepeatAnalysisGenomes_format_header.fa
#
#
# LinJ.27 TT38 vs Ubeda -- chr 27 comparison between current version and Ubeda
ref_qry=LinJ27.TT38_LinJ27.Ubeda
ref_fasta=~/leish_donovaniComplex/a24_repeats/TriTrypDB-38_LinfantumJPCM5_Genome_LinJ.27.fasta
qry_fasta=~/leish_donovaniComplex/a24_repeats/reference_genomes/20080508_RepeatAnalysisGenomes/LinJ27_20070420_V3.artemis.fasta
#
# genome TT38 vs Ubeda_all repeats -- mapping Ubeda repeats to current genome version
ref_qry=genomeTT38_RepeatsUbeda
ref_fasta=/lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta
qry_fasta=~/leish_donovaniComplex/a24_repeats/Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.fa
#
#
# run for the above three different parameter settings above
mkdir nucmer_${ref_qry}
cd nucmer_${ref_qry}
#
nucmer -p ${ref_qry} ${ref_fasta} ${qry_fasta}
show-coords -r ${ref_qry}.delta > ${ref_qry}.coords # -lTH
show-coords -lTH -r ${ref_qry}.delta > ${ref_qry}_tab.coords # -lTH
#
delta-filter -q -r ${ref_qry}.delta > ${ref_qry}.filter
# following cmd run on farm3
mummerplot ${ref_qry}.filter -R ${ref_fasta} -Q ${qry_fasta} -t png
#
# inspect by dropping delta file here: http://assemblytics.com


#
# genome TT38 vs Ubeda_all repeats -- mapping Ubeda repeats to current genome version
# ... and getting coordinates of those sequences in the currently used genome version
ref_qry=genomeTT38_RepeatsUbeda
ref_fasta=/lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta
qry_fasta=~/leish_donovaniComplex/a24_repeats/Ubeda_v3_uniq_20080508_RepeatAnalysisGenomes_format.fa
cd nucmer_${ref_qry}
#
# get coordinates of repeats on LinJ.27 for TT38 ref
cat genomeTT38_RepeatsUbeda_tab.coords | awk 'BEGIN{OFS="\t"}{if ($10=="LinJ.27") print $10,$1-1,$2,$11"_"$7,"NA",$1-1,$2}' > genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ27.coords.bed
# get coordinates of repeats on LinJ.35 for TT38 ref
cat genomeTT38_RepeatsUbeda_tab.coords | awk 'BEGIN{OFS="\t"}{if ($10=="LinJ.35") print $10,$1-1,$2,$11"_"$7,"NA",$1-1,$2}' > genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ35.coords.bed





###############################
#
# finding repeats "de novo"
#
# LinJ.27_reg TT38 vs TT38 - LinJ.27:190000-300000
# this region contains the region that was found to be missing in the assembly used for repeat identification in Ubeda et al 2014 and is also of particular interest to us as it contain the common deletion on chr 27 (Fig. 7)
samtools faidx /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-38_LinfantumJPCM5_Genome.fasta LinJ.27:190000-300000 > TriTrypDB-38_LinfantumJPCM5_Genome_LinJ.27_reg.190-300k.fasta
#
ref_qry=LinJ27regTT38_LinJ27regTT38
ref_fasta=~/leish_donovaniComplex/a24_repeats/TriTrypDB-38_LinfantumJPCM5_Genome_LinJ.27_reg.190-300k.fasta
#
#
mkdir nucmer_${ref_qry}
cd nucmer_${ref_qry}
#
nucmer --maxmatch --nosimplify -p ${ref_qry} ${ref_fasta} ${ref_fasta}
show-coords -r ${ref_qry}.delta > ${ref_qry}.coords # -lTH
show-coords -lTH -r ${ref_qry}.delta > ${ref_qry}_tab.coords # -lTH
#
delta-filter -q -r ${ref_qry}.delta > ${ref_qry}.filter
# following cmd run on farm3
mummerplot ${ref_qry}.filter -R ${ref_fasta} -Q ${qry_fasta} -t png
#
#
nucmer --maxmatch --nosimplify --mincluster 30 --minmatch 7 -p ${ref_qry}_mincl30_minm7 ${ref_fasta} ${ref_fasta} 
show-coords -lTH -r ${ref_qry}_mincl30_minm7.delta > ${ref_qry}_mincl30_minm7_tab.coords # -lTH
#
#
# queue=normal
# memory=4000
# name1=log_nucmer_mincl20_minm7
# bsub -n1 -R"span[hosts=1]" -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
# -R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
# nucmer --maxmatch --nosimplify --mincluster 20 --minmatch 7 -p ${ref_qry}_mincl20_minm7 ${ref_fasta} ${ref_fasta} 
# show-coords -lTH -r ${ref_qry}_mincl20_minm7.delta > ${ref_qry}_mincl20_minm7_tab.coords # -lTH
# #
# #
# queue=normal
# memory=4000
# name1=log_nucmer_mincl20_minm5
# bsub -n1 -R"span[hosts=1]" -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
# -R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
# nucmer --maxmatch --nosimplify --mincluster 20 --minmatch 5 -p ${ref_qry}_mincl20_minm5 ${ref_fasta} ${ref_fasta} 
# show-coords -lTH -r ${ref_qry}_mincl20_minm5.delta > ${ref_qry}_mincl20_minm5_tab.coords # -lTH
# 
# repeat-match -n 50 ${ref_fasta} > ${ref_qry}.repeats
# # exact-tandems ${ref_fasta} 50 > ${ref_qry}.tandems



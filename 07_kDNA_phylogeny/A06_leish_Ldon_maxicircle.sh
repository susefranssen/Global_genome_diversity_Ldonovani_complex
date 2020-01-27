#!/usr/bin/env bash
# grep -L "Successfully completed." *.o

# a
# look up paths for all samples and write into text file
bash ~/00scripts/leish_donovaniComplex/A06a_leish_Ldon_maxicircle_fastq.sh > sampIDs.leish.global.all_paths.txt &

# b
# create simlinks to fastq files
bash ~/00scripts/leish_donovaniComplex/A06b_leish_Ldon_maxicircle_fastq.sh 

# c
# map samples --> maybe start in several jobs due to resulting load on the farm
bash ~/00scripts/leish_donovaniComplex/A06c_leish_Ldon_maxicircle_smalt_mapping.sh 1026-8,1S,356WTV,363SKWTI,364SPWTII,38-UMK,383WTI,45-UMK,452BM,597-2,597LN,762L,855-9,AG83,AM560WTI,AM563WTI,ARL,BCN83,BCN87,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,BHU1062-4,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,BHU220A1,BHU41,BHU816A1,BHU824A1,BHU931A1,BPK029A1_a,BPK029A1_b,BPK029A1_c,BPK029A1_d,BPK035A1_a,BPK035A1_b,BPK067A1_a,BPK067A1_b,BPK067A1_c,BPK077A1_a,BPK077A1_b,BPK077A1_c,BPK077A1_d,BPK156A1,BPK157A1,BPK164A1_a,BPK164A1_b,BPK164A1_c,BPK282I9_a,BPK282I9_b,BPK294A1,BPK406A1_a,BPK406A1_b,BPK413A1,BPK471A1,BPK512A1_a,BPK512A1_b,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,BUCK,BUMM3_a,BUMM3_b,CH32,CH33_a,CH33_b,CH34_a,CH34_b,CH35,CH36_a,CH36_b,CL-SL,CL1_aethiopica,CL4_aethiopica,CUK10,CUK11,CUK12,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,Cha001,Chowd5,DD8_a,DD8_b,DOG_STRAIN,D_2,Don038,Don081,Don201,EP,GE,GEBRE1,GILANI_a,GILANI_b,HN167_a,HN167_b,HN336,IECDR1,IMT260,IMT373cl1,ISS174_a,ISS174_b,ISS2420,ISS2426,ISS2429,ISS2508,ITMAP26,Inf001,Inf004,Inf007,Inf045,Inf055,Inf152,Inf206,L60b,LEM1985,LEM3278,LEM3472,LPN114,LRC-L1275_a,LRC-L1275_b,LRC-L1296,LRC-L1311,LRC-L1312,LRC-L1313,LRC-L445,LRC-L51p,LRC-L53,LRC-L57,LRC-L61,LRC-L699,LRC-L740,LRC_L1303,LRC_L47,Ldon282cl2,LdonLV9_a,LdonLV9_b,LinJPCM5,LmexU1103_v1,LmjFried,MAM,MRC74_a,MRC74_b,Malta33,NLB-323,NT10,NT16,Nandi,OVN3,P283,Peking_a,Peking_b,RACOON_DOG,RM1,SKIN,STL2-78,STL2-79,STRAIN_A,STRAIN_B,SUDAN1_a,SUDAN1_b,SUKKAR2,TH4,TH5,TH6,WC,WR285_a,WR285_b 0.8
# OVN3,CL-SL missing
bash ~/00scripts/leish_donovaniComplex/A06c_leish_Ldon_maxicircle_smalt_mapping.sh OVN3,CL-SL 0.8
# Lamazonensis samples
bash ~/00scripts/leish_donovaniComplex/A06c_leish_Ldon_maxicircle_smalt_mapping_Lamazonensis.sh Lamazonensis_A,Lamazonensis_B 0.8


# d 
# index reference
bash ~/00scripts/leish_donovaniComplex/A06d_leish_maxicir_gatk_index.sh


# e
# realignment and quality filtering
bash ~/00scripts/leish_donovaniComplex/A06e_leish_maxicir_gatk_realign.sh 1026-8,1S,356WTV,363SKWTI,364SPWTII,38-UMK,383WTI,45-UMK,452BM,597-2,597LN,762L,855-9,AG83,AM560WTI,AM563WTI,ARL,BCN83,BCN87,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,BHU1062-4,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,BHU220A1,BHU41,BHU816A1,BHU824A1,BHU931A1,BPK156A1,BPK157A1,BPK294A1,BPK413A1,BPK471A1,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,BUCK,CH32,CH35,CL1_aethiopica,CL4_aethiopica,CUK10,CUK11,CUK12,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,Cha001,Chowd5,DOG_STRAIN,D_2,Don038,Don081,Don201,EP,GE,GEBRE1,HN336,IECDR1,IMT260,IMT373cl1,ISS2420,ISS2426,ISS2429,ISS2508,ITMAP26,Inf001,Inf004,Inf007,Inf045,Inf055,Inf152,Inf206,L60b,LEM1985,LEM3278,LEM3472,LPN114,LRC-L1296,LRC-L1311,LRC-L1312,LRC-L1313,LRC-L445,LRC-L51p,LRC-L53,LRC-L57,LRC-L61,LRC-L699,LRC-L740,LRC_L1303,LRC_L47,Ldon282cl2,LinJPCM5,LmexU1103_v1,LmjFried,MAM,Malta33,NLB-323,NT10,NT16,Nandi,P283,RACOON_DOG,RM1,SKIN,STL2-78,STL2-79,STRAIN_A,STRAIN_B,SUKKAR2,TH4,TH5,TH6,WC 0.8
#### 
bash ~/00scripts/leish_donovaniComplex/A06e_leish_maxicir_gatk_realign.sh Lamazonensis_A,Lamazonensis_B 0.8
#
# samples with several lanes sequenced where bam files need to be merged first
bash ~/00scripts/leish_donovaniComplex/A06e1_leish_maxicir_gatk_merge.sh
# index merged bams
for sample in BPK029A1 BPK077A1 BPK067A1 BPK164A1 BPK035A1 BPK282I9 BPK406A1 BPK512A1 BUMM3 CH33 CH34 CH36 DD8 GILANI HN167 ISS174 LdonLV9 LRC-L1275 MRC74 Peking SUDAN1 WR285
do
	queue=yesterday
	name1=index_${sample}
	bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
	-R"select[mem>4000] rusage[mem=4000]" -M4000 \
	samtools-1.3 index /lustre/scratch118/infgen/team133/sf18/021smalt/ldon_complex_${sample}/leth0.8/out.sorted.markdup.bam
done
bash ~/00scripts/leish_donovaniComplex/A06e_leish_maxicir_gatk_realign.sh BPK029A1,BPK077A1,BPK067A1,BPK164A1,BPK035A1,BPK282I9,BPK406A1,BPK512A1,BUMM3,CH33,CH34,CH36,DD8,GILANI,HN167,ISS174,LdonLV9,LRC-L1275,MRC74,Peking,SUDAN1,WR285 0.8
# 


# f
# bam index, varCalls, 
bash ~/00scripts/leish_donovaniComplex/A06f_leish_maxicir_callvariants.sh 1026-8,1S,356WTV,363SKWTI,364SPWTII,38-UMK,383WTI,45-UMK,452BM,597-2,597LN,762L,855-9,AG83,AM560WTI,AM563WTI,ARL,BCN83,BCN87,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,BHU1062-4,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,BHU220A1,BHU41,BHU816A1,BHU824A1,BHU931A1,BPK156A1,BPK157A1,BPK294A1,BPK413A1,BPK471A1,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,BUCK,CH32,CH35,CL1_aethiopica,CL4_aethiopica,CUK10,CUK11,CUK12,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,Cha001,Chowd5,DOG_STRAIN,D_2,Don038,Don081,Don201,EP,GE,GEBRE1,HN336,IECDR1,IMT260,IMT373cl1,ISS2420,ISS2426,ISS2429,ISS2508,ITMAP26,Inf001,Inf004,Inf007,Inf045,Inf055,Inf152,Inf206,L60b,LEM1985,LEM3278,LEM3472,LPN114,LRC-L1296,LRC-L1311,LRC-L1312,LRC-L1313,LRC-L445,LRC-L51p,LRC-L53,LRC-L57,LRC-L61,LRC-L699,LRC-L740,LRC_L1303,LRC_L47,Ldon282cl2,LinJPCM5,LmexU1103_v1,LmjFried,MAM,Malta33,NLB-323,NT10,NT16,Nandi,P283,RACOON_DOG,RM1,SKIN,STL2-78,STL2-79,STRAIN_A,STRAIN_B,SUKKAR2,TH4,TH5,TH6,WC 1
#
bash ~/00scripts/leish_donovaniComplex/A06f_leish_maxicir_callvariants.sh BPK029A1,BPK077A1,BPK067A1,BPK164A1,BPK035A1,BPK282I9,BPK406A1,BPK512A1,BUMM3,CH33,CH34,CH36,DD8,GILANI,HN167,ISS174,LdonLV9,LRC-L1275,MRC74,Peking,SUDAN1,WR285 1
# 
bash ~/00scripts/leish_donovaniComplex/A06f_leish_maxicir_callvariants.sh Lamazonensis_A,Lamazonensis_B 1



# do var calling stats in between with R script
# --> A06f1_leish_maxicir_var_metric.r
# --> A06f1_leish_maxicir_covdepth.r   <--- !!!

# g
# SNP hard filter, to fasta
# only take the 120 samples that were chosen based on proper coverage:
# --> ~/leish_donovaniComplex/A06_var_SNP_caling/takesamples.minmed20.txt
bash ~/00scripts/leish_donovaniComplex/A06g_leish_maxicir_filtvar.sh 1026-8,1S,356WTV,363SKWTI,364SPWTII,38-UMK,383WTI,45-UMK,452BM,597-2,597LN,762L,AG83,AM560WTI,AM563WTI,ARL,BCN83,BCN87,BD22,BHU1062-4,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,BHU220A1,BHU816A1,BHU824A1,BHU931A1,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK156A1,BPK157A1,BPK164A1,BPK282I9,BPK294A1,BPK406A1,BPK413A1,BPK471A1,BPK512A1,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,BUMM3,CH32,CH33,CH34,CH35,CH36,CL1_aethiopica,CL4_aethiopica,CUK10,CUK11,CUK12,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,Chowd5,DD8,DOG_STRAIN,D_2,EP,GE,GEBRE1,GILANI,HN167,IMT260,IMT373cl1,ISS174,ISS2420,ISS2426,ISS2429,ISS2508,ITMAP26,L60b,LEM1985,LEM3278,LPN114,LRC-L1275,LRC-L1296,LRC-L1311,LRC-L1312,LRC-L1313,LRC-L445,LRC-L57,LRC-L61,LRC_L1303,LRC_L47,Ldon282cl2,LdonLV9,LinJPCM5,LmexU1103_v1,MAM,MRC74,Malta33,NT10,NT16,Nandi,P283,Peking,RACOON_DOG,SKIN,STL2-78,STL2-79,STRAIN_A,STRAIN_B,SUDAN1,SUKKAR2,TH4,TH6,WC,WR285
#
bash ~/00scripts/leish_donovaniComplex/A06g_leish_maxicir_filtvar.sh Lamazonensis_A,Lamazonensis_B
#running

# h
# alignment  
# to be run BEFORE Lamazonensis samples were added                                                                                                                                  
bash ~/00scripts/leish_donovaniComplex/A06h_leish_maxicir_alignment.sh                                                                                                                                                                            
# to be run AFTER Lamazonensis samples were added                                                                                                                                  
bash ~/00scripts/leish_donovaniComplex/A06h_leish_maxicir_alignment_includeLamazon.sh

# AAAAAAAA


A06i_leish_maxicir_raxml.sh
A06i_leish_maxicir_raxml_inclLamazon.sh

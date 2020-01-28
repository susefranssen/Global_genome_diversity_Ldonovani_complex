

# grep -L "Successfully completed." *.o

# transform Caroline's vcflike file to a vcf format accepted by vcftools for genotype
# correlation calculations (LD)


cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A04_vcf_LD




# names indicates the group ans how many samples are present for the group
#
# for the first 5 'groups' no LD was calculated, they were just run through this pipeline for Venn diagrams later
for name in "Ldon4_4" "EP_1" "other.Ldon_3" "CH_hy_2" "CH_nohy_3"
for name in "infantum_47" "donovani1_52" "donovani2a_19" "donovani2b_7" "donovani3_8" "TurkeyH_11" "donovani1_main42" "Ldon1_new_45" "Ldon5_7" "Ldon2_new_4" "Ldon1_noBPK157A1_44" "Linf1_Mon1_31" "Linf1_China_7" "Linf1_ChinaUzb_10" "Linf1_nonMon1_5" "Ldon3_noLRCL53_7" "Ldon4_noGILANI_18" "Linf1_Mon1_noISSm_31" "CH" "Linf_noISSextr_noInf152_43" "Ldon4_noGILANI_18_s7R1" "Ldon1_noBPK157A1_44_s7R1" "Linf_noISSextr_noInf152_43_s7R1" "TurkeyH_11_s7R1" "Linf1_Mon1_31_s5R1" "Linf_noISSextr_noInf152_43_s18R1" "Linf_noISSextr_noInf152_43_s18R2" "Linf_noISSextr_noInf152_43_s18R3" "Ldon1_noBPK157A1_44_s18R1" "Ldon1_noBPK157A1_44_s18R2" "Ldon1_noBPK157A1_44_s18R3"
for name in "Ldon1_noBPK157A1_44_s7R2" "Ldon1_noBPK157A1_44_s7R3" "Ldon4_noGILANI_18_s7R2" "Ldon4_noGILANI_18_s7R3" "Linf_noISSextr_noInf152_43_s7R2" "Linf_noISSextr_noInf152_43_s7R3" "TurkeyH_11_s7R2" "TurkeyH_11_s7R3"
for name in "Ldon1_noBPK157A1_44" "Ldon5_7" "Ldon4_noGILANI_18" "Ldon3_noLRCL53_7" "Linf_noISSextr_noInf152_43" "Ldon4_4" "CH" "TurkeyH_11"
do
	mkdir $name
	
	if [ "$name" = "infantum_47" ]; then
		group="Peking,D_2,STRAIN_B,STRAIN_A,RACOON_DOG,SKIN,DOG_STRAIN,LRC-L1311,LRC-L1313,Inf007,LRC-L699,LRC-L1275,TH4,TH5,TH6,NT10,NT16,LRC-L1312,LRC-L1296,LRC_L1303,Inf152,ISS174,ISS2420,ISS2426,ISS2429,ISS2508,Inf001,LRC_L47,RM1,LEM1985,LPN114,LEM3278,Inf045,Inf004,Inf055,BCN83,BCN87,LinJPCM5,IMT260,IMT373cl1,ITMAP26,Cha001,ARL,WC,WR285,HN167,HN336"
	elif [ "$name" = "Linf_noISSextr_noInf152_43" ]; then
		group="Peking,D_2,STRAIN_B,STRAIN_A,RACOON_DOG,SKIN,DOG_STRAIN,LRC-L1311,LRC-L1313,Inf007,LRC-L699,LRC-L1275,TH4,TH5,TH6,NT10,NT16,LRC-L1312,LRC-L1296,LRC_L1303,ISS2420,ISS2508,Inf001,LRC_L47,RM1,LEM1985,LPN114,LEM3278,Inf045,Inf004,Inf055,BCN83,BCN87,LinJPCM5,IMT260,IMT373cl1,ITMAP26,Cha001,ARL,WC,WR285,HN167,HN336"
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s18R1" ]; then
		group="ARL,BCN83,Cha001,HN336,IMT260,Inf001,Inf004,Inf045,ISS2508,ITMAP26,LEM3278,LPN114,LRC_L47,LRC-L699,Peking,RM1,SKIN,WR285"
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s18R2" ]; then
		group="BCN83,Cha001,DOG_STRAIN,Inf045,Inf055,LEM1985,LinJPCM5,LRC-L1275,LRC-L1296,LRC-L1312,LRC-L1313,NT10,RM1,SKIN,STRAIN_B,TH5,TH6,WR285"
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s18R3" ]; then
		group="D_2,DOG_STRAIN,HN167,IMT373cl1,Inf001,Inf004,Inf045,Inf055,ISS2420,ISS2508,LinJPCM5,LRC-L1296,LRC-L1312,NT10,RM1,STRAIN_A,TH4,TH6"
	elif [ "$name" = "Linf1_Mon1_31" ]; then
	group="Inf007,LRC-L699,LRC-L1275,TH4,TH5,TH6,NT10,NT16,LRC-L1296,LRC_L1303,ISS174,ISS2420,ISS2426,ISS2429,ISS2508,Inf001,RM1,LEM1985,LPN114,LEM3278,BCN83,BCN87,LinJPCM5,IMT260,IMT373cl1,ITMAP26,Cha001,ARL,WC,HN167,HN336"
	elif [ "$name" = "Linf1_Mon1_31_s5R1" ]; then
	group="ARL,BCN87,Inf001,ISS2508,TH5"
	elif [ "$name" = "Linf1_Mon1_noISSm_31" ]; then
	group="Inf007,LRC-L699,LRC-L1275,TH4,TH5,TH6,NT10,NT16,LRC-L1296,LRC_L1303,ISS2420,ISS2508,Inf001,RM1,LEM1985,LPN114,LEM3278,BCN83,BCN87,LinJPCM5,IMT260,IMT373cl1,ITMAP26,Cha001,ARL,WC,HN167,HN336"	
	elif [ "$name" = "Linf1_China_7" ]; then
	group="Peking,RACOON_DOG,DOG_STRAIN,D_2,STRAIN_B,STRAIN_A,SKIN"
	elif [ "$name" = "Linf1_ChinaUzb_10" ]; then
	group="Peking,RACOON_DOG,DOG_STRAIN,D_2,STRAIN_B,STRAIN_A,SKIN,LRC-L1311,LRC-L1313,LRC-L1312"	
	elif [ "$name" = "Linf1_nonMon1_5" ]; then
	group="Inf004,Inf055,Inf045,LRC_L47,WR285"	
	elif [ "$name" = "donovani1_52" ]; then
		group="IECDR1,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,Ldon282cl2,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK156A1,BPK157A1,BPK164A1,BPK282I9,BPK294A1,BPK406A1,BPK413A1,BPK471A1,BPK512A1,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,L60b,OVN3,CL-SL,Chowd5,STL2-78,STL2-79,DD8,Nandi,AG83,BHU41,BHU220A1,BHU1062-4,Don201,BHU816A1,BHU824A1,BHU931A1,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,Inf206,BUCK"
	elif [ "$name" = "donovani1_main42" ]; then
		group="IECDR1,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,Ldon282cl2,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK157A1,BPK164A1,BPK282I9,BPK294A1,BPK471A1,BPK562A1,BPK649A1,Chowd5,STL2-78,STL2-79,DD8,Nandi,AG83,BHU41,BHU220A1,BHU1062-4,Don201,BHU816A1,BHU824A1,BHU931A1,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,Inf206,BUCK"		
	elif [ "$name" = "donovani2a_19" ]; then
		group="LdonLV9,GEBRE1,356WTV,363SKWTI,383WTI,364SPWTII,LRC-L61,1S,GILANI,38-UMK,45-UMK,452BM,597-2,597LN,762L,1026-8,855-9,SUDAN1,Malta33"
	elif [ "$name" = "Ldon4_noGILANI_18" ]; then
		group="LdonLV9,GEBRE1,356WTV,363SKWTI,383WTI,364SPWTII,LRC-L61,1S,38-UMK,45-UMK,452BM,597-2,597LN,762L,1026-8,855-9,SUDAN1,Malta33"
	elif [ "$name" = "donovani2b_7" ]; then
		group="BUMM3,SUKKAR2,Don081,Don038,GE,LEM3472,LRC-L740"
	elif [ "$name" = "donovani3_8" ]; then
		group="LRC-L51p,AM560WTI,AM563WTI,LRC-L53,LRC-L57,MRC74,NLB-323,LRC-L445"
	elif [ "$name" = "Ldon3_noLRCL53_7" ]; then
		group="LRC-L51p,AM560WTI,AM563WTI,LRC-L57,MRC74,NLB-323,LRC-L445"
	elif [ "$name" = "TurkeyH_11" ]; then
		group="CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,CUK10,CUK11,CUK12"
	# downsampled
	elif [ "$name" = "infantum_7" ]; then
		group="Cha001,ISS2508,LEM1985,LRC_L1303,LRC-L1296,TH4,TH6"
	elif [ "$name" = "donovani1_7" ]; then
		group="BD27,BHU1062-4,BPK413A1,BPK623A1,BUCK,Inf206,OVN3" 
	elif [ "$name" = "donovani2a_7" ]; then
		group="1026-8,1S,356WTV,38-UMK,45-UMK,597LN,SUDAN1"
	elif [ "$name" = "donovani3_7" ]; then
		group="AM560WTI,AM563WTI,LRC-L445,LRC-L51p,LRC-L57,MRC74,NLB-323"		
	elif [ "$name" = "Ldon1_new_45" ]; then
		group="IECDR1,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,Ldon282cl2,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK157A1,BPK164A1,BPK282I9,BPK294A1,BPK471A1,BPK562A1,BPK649A1,L60b,OVN3,CL-SL,Chowd5,STL2-78,STL2-79,DD8,Nandi,AG83,BHU41,BHU220A1,BHU1062-4,Don201,BHU816A1,BHU824A1,BHU931A1,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,Inf206,BUCK"
	elif [ "$name" = "Ldon1_noBPK157A1_44" ]; then
	 group="IECDR1,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,Ldon282cl2,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK164A1,BPK282I9,BPK294A1,BPK471A1,BPK562A1,BPK649A1,L60b,OVN3,CL-SL,Chowd5,STL2-78,STL2-79,DD8,Nandi,AG83,BHU41,BHU220A1,BHU1062-4,Don201,BHU816A1,BHU824A1,BHU931A1,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,Inf206,BUCK"
	elif [ "$name" = "Ldon1_noBPK157A1_44_s18R1" ]; then
	 group="BD09,BD15,BD22,BD27,BHU1062-4,BHU1065A1,BHU1137A1,BHU220A1,BHU816A1,BPK077A1,BPK471A1,BPK562A1,BUCK,CL-SL,DD8,IECDR1,Ldon282cl2,STL2-79"
	elif [ "$name" = "Ldon1_noBPK157A1_44_s18R2" ]; then
	 group="AG83,BD12,BD17,BD24,BD25,BD27,BHU1137A1,BHU1139A1,BHU824A1,BPK035A1,BPK282I9,BPK562A1,BPK649A1,BUCK,DD8,L60b,OVN3,STL2-78"
	elif [ "$name" = "Ldon1_noBPK157A1_44_s18R3" ]; then
	 group="BD12,BD15,BD22,BD25,BD27,BHU1139A1,BHU220A1,BHU41,BHU931A1,BPK035A1,BPK067A1,BPK077A1,BPK164A1,BPK562A1,BPK649A1,BUCK,Nandi,STL2-79"		
	elif [ "$name" = "Ldon5_7" ]; then
		group="BPK512A1,BPK406A1,BPK623A1,BPK156A1,BPK413A1,BPK648A1,BPK612A1"
	elif [ "$name" = "Ldon2_new_4" ]; then
		group="SUKKAR2,BUMM3,Don038,Don081"
	elif [ "$name" = "CH" ]; then
		group="CH32,CH33,CH34,CH35,CH36"
	#######	
	elif [ "$name" = "Ldon1_noBPK157A1_44_s7R1" ]; then
	 group="BD15,BD27,BPK164A1,BPK649A1,BUCK,Chowd5,DD8"
	elif [ "$name" = "Ldon1_noBPK157A1_44_s7R2" ]; then
	 group="BD21,BPK029A1,BPK649A1,CL-SL,Don201,L60b,STL2-78"
	elif [ "$name" = "Ldon1_noBPK157A1_44_s7R3" ]; then
	 group="AG83,BD12,BD25,BPK282I9,BPK471A1,BPK649A1,Chowd5" 
	#
	elif [ "$name" = "Ldon4_noGILANI_18_s7R1" ]; then
		group="1S,364SPWTII,452BM,597LN,762L,GEBRE1,Malta33"
	elif [ "$name" = "Ldon4_noGILANI_18_s7R2" ]; then
		group="1S,45-UMK,597-2,597LN,855-9,LdonLV9,Malta33"
	elif [ "$name" = "Ldon4_noGILANI_18_s7R3" ]; then
		group="356WTV,383WTI,452BM,597LN,855-9,GEBRE1,LRC-L61"
	#
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s7R1" ]; then
		group="HN336,Inf045,ITMAP26,NT16,STRAIN_B,TH5,TH6"
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s7R2" ]; then
		group="ARL,HN167,Inf007,ITMAP26,LEM1985,LRC-L699,TH6"
	elif [ "$name" = "Linf_noISSextr_noInf152_43_s7R3" ]; then
		group="D_2,LRC_L1303,LRC_L47,LRC-L1311,Peking,RACOON_DOG,SKIN"
	#
	elif [ "$name" = "TurkeyH_11_s7R1" ]; then
		group="CUK2,CUK4,CUK5,CUK6,CUK7,CUK9,CUK11,CUK12"
	elif [ "$name" = "TurkeyH_11_s7R2" ]; then
		group="CUK10,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7"
	elif [ "$name" = "TurkeyH_11_s7R3" ]; then
		group="CUK10,CUK11,CUK12,CUK2,CUK3,CUK4,CUK9"
		
	elif [ "$name" = "Ldon4_4" ]; then
		group="SUKKAR2,BUMM3,Don038,Don081"
	elif [ "$name" = "EP_1" ]; then
		group="EP"
	elif [ "$name" = "other.Ldon_3" ]; then
		group="LEM3472,GE,LRC-L740"
	elif [ "$name" = "CH_hy_2" ]; then
		group="CH32,CH34"
	elif [ "$name" = "CH_nohy_3" ]; then
		group="CH33,CH35,CH36"
	fi;

	echo $name
	echo $group

	#for ((chr=1; chr <= 36 ; chr++))
	for chr in 00 
# 	for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
	do
		echo $chr
	
		# transform to vcf format
		python ~/00scripts/leish_donovaniComplex/a04_vcflike_to_vcf.py \
		--vcflike /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}.vcflike.txt \
		--sample_names ${group} >\
		/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}.vcf 
 
# 		# transform to diploids
# 		python ~/00scripts/general_scripts/vcf_transform_to_diploid_GT-DP-GQ.py \
# 		--vcf /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}.vcf >\
# 		/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf
# 	
# 		# calc genotype correlations
# 		maxdist=1000000
# 		queue=yesterday
# 		memory=2000
# 		name1=A04_vcf_genor2_${maxdist}_${name}_${chr}
# 		bsub -q $queue -o ${name1}.o -e ${name1}.e -J ${name1} \
# 		-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
# 		/software/vcftools-0.1.14/bin/vcftools \
# 		--vcf /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf \
# 		--geno-r2 --ld-window-bp ${maxdist} \
# 		--out ${name}/LinJ.${chr}_${name}_diploid_genor2_maxd${maxdist}.txt
# 		
		
		
	done

done

# all data in one file
# to vcf
group="Peking,D_2,STRAIN_B,STRAIN_A,RACOON_DOG,SKIN,DOG_STRAIN,IECDR1,BD09,BD12,BD14,BD15,BD17,BD21,BD22,BD24,BD25,BD27,Ldon282cl2,BPK029A1,BPK035A1,BPK067A1,BPK077A1,BPK156A1,BPK157A1,BPK164A1,BPK282I9,BPK294A1,BPK406A1,BPK413A1,BPK471A1,BPK512A1,BPK562A1,BPK612A1,BPK623A1,BPK648A1,BPK649A1,L60b,OVN3,CL-SL,LRC-L51p,Chowd5,STL2-78,STL2-79,DD8,Nandi,AG83,BHU41,BHU220A1,BHU1062-4,Don201,BHU816A1,BHU824A1,BHU931A1,BHU1064A1,BHU1065A1,BHU1137A1,BHU1139A1,LRC-L1311,LRC-L1313,Inf206,BUMM3,SUKKAR2,Don081,LdonLV9,GEBRE1,Don038,356WTV,AM560WTI,363SKWTI,383WTI,364SPWTII,AM563WTI,LRC-L53,LRC-L57,MRC74,NLB-323,LRC-L445,LRC-L61,1S,GILANI,GE,38-UMK,45-UMK,452BM,597-2,597LN,762L,LEM3472,1026-8,855-9,SUDAN1,Inf007,LRC-L699,LRC-L740,LRC-L1275,TH4,TH5,TH6,NT10,NT16,LRC-L1312,LRC-L1296,LRC_L1303,CH32,CH33,CH34,CH35,CH36,EP,Inf152,CUK2,CUK3,CUK4,CUK5,CUK6,CUK7,CUK8,CUK9,CUK10,CUK11,CUK12,ISS174,ISS2420,ISS2426,ISS2429,ISS2508,Malta33,BUCK,Inf001,LRC_L47,RM1,LEM1985,LPN114,LEM3278,Inf045,Inf004,Inf055,BCN83,BCN87,LinJPCM5,IMT260,IMT373cl1,ITMAP26,Cha001,MAM,ARL,WC,WR285,HN167,HN336"
python ~/00scripts/leish_donovaniComplex/a04_vcflike_to_vcf.py \
--vcflike /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike.txt \
--sample_names $group > /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike.vcf
# to diploid vcf
python ~/00scripts/general_scripts/vcf_transform_to_diploid_GT-DP-GQ.py \
--vcf /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike.vcf >\
/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike_diploid.vcf



ls ../../06snps/leish_donovaniComplex/






#-----------------------------------------------------
#
# R2 - dist --> boxplots
# summary across chr
#
queue=yesterday
memory=4000


#
#
#
memory=6500
name=R_A04_vcf_LD_box_1l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r donovani1_52 1000 40 200 100000
#
name=R_A04_vcf_LD_box_3l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r donovani2a_19 1000 7 200 100000
#
memory=7500
name=R_A04_vcf_LD_box_4l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r donovani2b_7 1000 7 200 100000
#
memory=6500
name=R_A04_vcf_LD_box_5l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r donovani3_8 1000 7 200 100000
#
name=R_A04_vcf_LD_box_2l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r infantum_47 1000 40 200 100000
#

#
name=R_A04_vcf_LD_box_6l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon1_new_45 1000 40 200 100000
#
name=R_A04_vcf_LD_box_7l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon2_new_4 1000 4 200 100000
#
name=R_A04_vcf_LD_box_8l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon5_7 1000 7 200 100000
#
name=R_A04_vcf_LD_box_9l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon1_noBPK157A1_44 1000 7 200 100000
#
name=R_A04_vcf_LD_box_10l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Linf1_Mon1_31 1000 28 200 100000
#
name=R_A04_vcf_LD_box_11l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Linf1_China_7 1000 6 200 100000
#
name=R_A04_vcf_LD_box_12l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Linf1_ChinaUzb_10 1000 8 200 100000
#
name=R_A04_vcf_LD_box_13l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Linf1_nonMon1_5 1000 5 200 100000
#
name=R_A04_vcf_LD_box_14l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon3_noLRCL53_7 1000 6 200 100000
#
name=R_A04_vcf_LD_box_15l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Ldon4_noGILANI_18 1000 16 200 100000
#
name=R_A04_vcf_LD_box_16l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r Linf1_Mon1_noISSm_31 1000 28 200 100000
#
name=R_A04_vcf_LD_box_17l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r CH 1000 5 200 100000
#
name=R_A04_vcf_LD_box_18l
# -w "done(${name1})"
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
Rscript ~/00scripts/leish_donovaniComplex/a04_LD_dist_02.r TurkeyH_11 1000 9 200 100000





#---------------------
# R2 between chr

cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A04_vcf_LD


repl=1

# for 100 randomly chosen SNPs for each chr
#
for name in "infantum_47" "donovani1_52" "donovani2a_19" "donovani2b_7" "donovani3_8" "TurkeyH_11" "donovani1_main42" "Ldon1_new_45" "Ldon5_7" "Ldon2_new_4" "Ldon1_noBPK157A1_44" "Linf1_Mon1_31" "Linf1_China_7" "Linf1_ChinaUzb_10" "Linf1_nonMon1_5" "Ldon3_noLRCL53_7" "Ldon4_noGILANI_18" "Linf1_Mon1_noISSm_31" "CH" "Linf_noISSextr_noInf152_43"
for name in "Ldon1_noBPK157A1_44" "Ldon5_7" "Ldon4_noGILANI_18" "Ldon3_noLRCL53_7" "Linf_noISSextr_noInf152_43" "TurkeyH_11"
for name in "Ldon1_noBPK157A1_44_s18R1" "Ldon1_noBPK157A1_44_s18R2" "Ldon1_noBPK157A1_44_s18R3" "Linf_noISSextr_noInf152_43_s18R1" "Linf_noISSextr_noInf152_43_s18R2" "Linf_noISSextr_noInf152_43_s18R3"
for name in "Ldon1_noBPK157A1_44_s7R2" "Ldon1_noBPK157A1_44_s7R3" "Ldon4_noGILANI_18_s7R2" "Ldon4_noGILANI_18_s7R3" "Linf_noISSextr_noInf152_43_s7R2" "Linf_noISSextr_noInf152_43_s7R3" "TurkeyH_11_s7R2" "TurkeyH_11_s7R3"
for name in "Ldon1_noBPK157A1_44_s7R1" "Ldon4_noGILANI_18_s7R1" "Linf_noISSextr_noInf152_43_s7R1" "TurkeyH_11_s7R1"
do
	grep '^\#' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.01_${name}_diploid.vcf > header.txt

	for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
	do
		grep -e '1/0' -e '0/1' -e '1/1' -e '2/0' -e '0/2' -e '2/2' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf | grep '0/0' | sort -R | head -n 100 | sort -k1,1 -k2,2n > /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf.tmp
		
		# echo $name $chr
# 		grep -e '1/0' -e '0/1' -e '1/1' -e '2/0' -e '0/2' -e '2/2' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf | grep '0/0' | wc -l > SNPnum_${name}_${chr}.txt
# 		grep -e '1/0' -e '0/1' -e '1/1' -e '2/0' -e '0/2' -e '2/2' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.${chr}_${name}_diploid.vcf | grep '0/0' | head -n 100 | wc -l >> SNPnum_${name}_${chr}.txt
	done

	cat header.txt /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.*_${name}_diploid.vcf.tmp > /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.all_${name}_diploid_random100_${repl}.vcf
	
	rm /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.*_${name}_diploid.vcf.tmp
	

	# calc genotype correlations between chr
	queue=yesterday
	name2=A04_vcf_genor2_interchrom_${name}_${repl}
	bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} -R "select[mem>4000] rusage[mem=4000]" -M 4000 /software/vcftools-0.1.14/bin/vcftools --vcf /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.LinJ.all_${name}_diploid_random100_${repl}.vcf --interchrom-geno-r2 --out ${name}/LinJ.all_${name}_diploid_genor2_r100_${repl}
		
done








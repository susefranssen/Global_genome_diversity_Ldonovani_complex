

cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A12_group_freqs

path=/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/


# "Ldon4_4" "EP_1" "other.Ldon_3" "CH" "CH_hy_2" "CH_nohy_3" are not needed for group freqs but for SNP venn diagrams
for group in "Ldon4_4" "EP_1" "other.Ldon_3" "CH" "CH_hy_2" "CH_nohy_3"
for group in "Ldon4_noGILANI_18_s7R1"  "Ldon1_noBPK157A1_44_s7R1" "Linf_noISSextr_noInf152_43_s7R1" "TurkeyH_11_s7R1" "Ldon4_noGILANI_18_s7R2" "Ldon1_noBPK157A1_44_s7R2" "Linf_noISSextr_noInf152_43_s7R2"   "TurkeyH_11_s7R2" "Ldon4_noGILANI_18_s7R3" "Ldon1_noBPK157A1_44_s7R3"   "Linf_noISSextr_noInf152_43_s7R3" "TurkeyH_11_s7R3"
for group in "Ldon1_noBPK157A1_44_s18R1" "Ldon1_noBPK157A1_44_s18R2" "Ldon1_noBPK157A1_44_s18R3" "Linf_noISSextr_noInf152_43_s18R1" "Linf_noISSextr_noInf152_43_s18R2" "Linf_noISSextr_noInf152_43_s18R3" "Ldon1_noBPK157A1_44" "Linf_noISSextr_noInf152_43" "Ldon4_noGILANI_18" "TurkeyH_11" "Ldon5_7" "donovani3_8" "Ldon3_noLRCL53_7" 
do
mkdir ${group}
for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
echo "${group}_${chr}"
python3 ~/00scripts/general_scripts/python3/vcf_get_SFS/vcf_get_SFS_aneuploid.py \
--vcf ${path}snps.filt.leish_global.linj.complex.LinJ.${chr}_${group}.vcf \
> ${group}/group_freqs_${group}_${chr}_aneuploid.txt
done
cat ${group}/group_freqs_${group}_*_aneuploid.txt > group_freqs_${group}_all_aneuploid.tmp
# remove multiple header within file
awk ' /^chr\tpos\tref\talt\treffreq\talt1freq\talt2frq\tmajallelefreq\tcov/ && FNR > 1 {next} {print $0} ' group_freqs_${group}_all_aneuploid.tmp \
> group_freqs_${group}.all_aneuploid.txt
rm -r ${group}
rm group_freqs_${group}_all_aneuploid.tmp
done

# this error occurs when no genotypes are called for this position, na s are written in this case and the warning can be ignored
#
#   freqs=[str(round(ac_row_sum[i]/sum(ac_row_sum),4)) for i in [0,1,2]] # freqs for ref and both alternate alleles



# temp copy to read into R
cp /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A12_group_freqs/group_freqs_*.all_aneuploid.txt ~/tmp/.

cp /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-9.0_LinfantumJPCM5_Genome.fasta.fai ~/leish_donovaniComplex/A12_group_freqs/.


# subset of interest
#    72909 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44.all.txt
#      684 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R1.all.txt
#    23304 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R2.all.txt
#      686 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R3.all.txt
#
#     9682 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon3_noLRCL53_7.all.txt
#
#    47008 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18.all.txt
#    38768 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R1.all.txt
#    35822 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R2.all.txt
#    38623 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R3.all.txt
#
#    45588 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon5_7.all.txt
#
#   116269 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43.all.txt
#    14525 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R1.all.txt
#     2882 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R2.all.txt
#    15510 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R3.all.txt
#
#    10523 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11.all.txt
#    10347 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R1.all.txt
#    10362 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R2.all.txt
#    10295 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R3.all.txt


# all
#    99658 /nfs/users/nfs_s/sf18/tmp/group_freqs_donovani1_52.all.txt
#    47683 /nfs/users/nfs_s/sf18/tmp/group_freqs_donovani2a_19.all.txt
#    16974 /nfs/users/nfs_s/sf18/tmp/group_freqs_donovani3_8.all.txt
#    52962 /nfs/users/nfs_s/sf18/tmp/group_freqs_infantum_47.all.txt
#    27208 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_new_45.all.txt
#    72909 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44.all.txt
#    23050 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s18R1.all.txt
#    23848 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s18R2.all.txt
#     1183 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s18R3.all.txt
#      684 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R1.all.txt
#    23304 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R2.all.txt
#      686 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon1_noBPK157A1_44_s7R3.all.txt
#     9682 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon3_noLRCL53_7.all.txt
#    47008 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18.all.txt
#    38768 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R1.all.txt
#    35822 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R2.all.txt
#    38623 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon4_noGILANI_18_s7R3.all.txt
#    45588 /nfs/users/nfs_s/sf18/tmp/group_freqs_Ldon5_7.all.txt
#    24590 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf1_Mon1_31.all.txt
#     3054 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf1_Mon1_31_s5R1.all.txt
#     5673 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf1_Mon1_noISSm_31.all.txt
#    20957 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf1_nonMon1_5.all.txt
#   116269 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43.all.txt
#    27332 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s18R1.all.txt
#    27941 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s18R2.all.txt
#    27234 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s18R3.all.txt
#    14525 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R1.all.txt
#     2882 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R2.all.txt
#    15510 /nfs/users/nfs_s/sf18/tmp/group_freqs_Linf_noISSextr_noInf152_43_s7R3.all.txt
#    10523 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11.all.txt
#    10347 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R1.all.txt
#    10362 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R2.all.txt
#    10295 /nfs/users/nfs_s/sf18/tmp/group_freqs_TurkeyH_11_s7R3.all.txt



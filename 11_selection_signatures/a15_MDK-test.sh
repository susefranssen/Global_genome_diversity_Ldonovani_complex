
# grep -L "Successfully completed." *.o


cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A15_MDK-test


# This approach is ignoring introns as we know from the newest annotation version v39, August 2018 that there are non- multi exon genes any more

###################
# overview
# 1. extract mRNA seq (ignores intron annotation)
# 2a find longest ORF (of multiexon genes only) 
# 3 compare ORF of previous multiexon genes (ignoring introns) with ORF of new annotation
# 2b find longest ORF od all genes
# 4 run MK test (MK.pl)
###################



###################
# 1
# extract mRNA sequences
###################
# extract nt seqs ignoring introns (as seem to be obsolete in the newest annotation version)
awk '($3=="mRNA") {OFS="\t"; print $0}' /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-9.0_LinfantumJPCM5.gff | awk -F '[\t;=_]' 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$11"_"$1"_"$4"_"$5"_"$7,1,$7}' > TriTrypDB-9.0_LinfantumJPCM5_mRNA.bed

# extract mRNA from gff file
bedtools getfasta -name -s -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-9.0_LinfantumJPCM5_Genome.fasta -bed TriTrypDB-9.0_LinfantumJPCM5_mRNA.bed -fo TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa

# for extracting ID for multi exon genes only
#
# get bed file containing CDS seqs
awk '($3=="CDS") {OFS="\t"; print $0}' /lustre/scratch118/infgen/team133/sf18/refgenomes/annotation/TriTrypDB-9.0_LinfantumJPCM5.gff | awk -F '[\t;=_]' 'BEGIN{OFS="\t"}{print $1,$4-1,$5,$11"_"$1"_"$4"_"$5"_"$7,1,$7}' > TriTrypDB-9.0_LinfantumJPCM5_CDS.bed
#
# extract CDS from gff file
bedtools getfasta -name -s -fi /lustre/scratch118/infgen/team133/sf18/refgenomes/TriTrypDB-9.0_LinfantumJPCM5_Genome.fasta -bed TriTrypDB-9.0_LinfantumJPCM5_CDS.bed -fo TriTrypDB-9.0_LinfantumJPCM5_CDS.fa
#
grep '0-2' TriTrypDB-9.0_LinfantumJPCM5_CDS.fa | cut -d'-' -f1 | cut -d'>' -f2 | sort -u > multi_exon_genes.txt
#
# get longest ORF from within the mRNA seq
python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/get_subset_fasta.py \
--fasta TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa \
--out TriTrypDB-9.0_LinfantumJPCM5_mRNA_multiexon.fa \
--tag multi_exon_genes.txt



###################
# 2a
# find longest ORF (of multiexon genes only)

python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_a_find_ORF.py --CDS TriTrypDB-9.0_LinfantumJPCM5_mRNA_multiexon.fa --log TriTrypDB-9.0_LinfantumJPCM5_mRNA_multiexon.log > TriTrypDB-9.0_LinfantumJPCM5_mRNA_multiexon_ORF.fa



###################
# 3
# compare ORF of previous multiexon genes (ignoring introns) with ORF of new annotation 

# get mapping from new to old gene annotation
grep '^LINF' ~/leish_donovaniComplex/A15_MDK-test/geneid_mapping_new_old_annot.txt | wc -l
# 8748		--> number of new annotated genes on chr scaffolds
grep '^LINF' ~/leish_donovaniComplex/A15_MDK-test/geneid_mapping_new_old_annot.txt | grep 'LinJ'  | wc -l
#    8614   
#
grep '^LINF' ~/leish_donovaniComplex/A15_MDK-test/geneid_mapping_new_old_annot.txt | grep 'LinJ' | awk 'BEGIN{OFS="\t"}{print $1,$5}' > geneid_mapping_new_old_annot.tmp
#
sed 's/,//g' geneid_mapping_new_old_annot.tmp > geneid_mapping_new_old_annot.txt
rm geneid_mapping_new_old_annot.tmp
#
cut -f1 geneid_mapping_new_old_annot.txt | sort -u |wc -l
# 8614		--> number of new annotated genes that have orthologs in the old annotated gens on the chr scaffolds
cut -f2 geneid_mapping_new_old_annot.txt | sort -u |wc -l
# 8594		--> number of old annotated genes on chr scaffolds

python3 ~/00scripts/general_scripts/python3/fasta_to_1line.py --fasta ../../refgenomes/TriTrypDB-39_LinfantumJPCM5_AnnotatedCDSs.fasta > ../../refgenomes/TriTrypDB-39_LinfantumJPCM5_AnnotatedCDSs_single-line.fasta



# get a matching new version gene ids for all previous intron genes
echo "" > multi_exon_genes_match_new_geneids.txt
while read gene
do
# new name of respective gene
LINF=$(grep $gene geneid_mapping_new_old_annot.txt | cut -f1 )
echo -e $gene '\t' $LINF >> multi_exon_genes_new_geneids.txt
done < multi_exon_genes.txt
cat multi_exon_genes_new_geneids.txt | awk 'BEGIN{OFS="\t"}{if($2!="") print $1,$2}' > multi_exon_genes_new_geneids_matched.txt
wc -l multi_exon_genes_new_geneids*
#   35 multi_exon_genes_new_geneids_matched.txt
#   43 multi_exon_genes_new_geneids.txt
# --> 8 of the previous intron genes do not exist in the new genome and annotation version any more
# only old gene name of all previous multi exon genes that have a match to the new ones
cat multi_exon_genes_new_geneids_matched.txt | cut -f1 > multi_exon_genes_matched.txt



echo "" > old_intron_gene_comp.txt
while read gene
do
# gene=LinJ.01.0850
echo $gene
# combine CDS from new annotation version with extracted longest ORFs from old version annotation mRNA
# new name of respective gene
LINF=$(grep $gene multi_exon_genes_new_geneids_matched.txt | cut -f2 )
echo -e $gene '\t' $LINF
#
# for every old "intron" gene get the CDS of the same gene in the new annotation version
grep '>'$LINF -A1 ../../refgenomes/TriTrypDB-39_LinfantumJPCM5_AnnotatedCDSs_single-line.fasta > 	tmp_${LINF}_CDS.fa
#
# get the respective gene sequence from the old version where the longest ORF was searched in the annotated mRNA region (ignoring introns)
grep '>'$gene -A1 TriTrypDB-9.0_LinfantumJPCM5_mRNA_multiexon_ORF.fa > tmp_${gene}_mRNA_ORF.fa
#
# combine both CDS seqs in one file (from the new annotation and the "Suse-curated" one using the old version)
cat tmp_${gene}_mRNA_ORF.fa tmp_${LINF}_CDS.fa > ${gene}_${LINF}_ORF.fa
#
# generate alignment of both
muscle -in ${gene}_${LINF}_ORF.fa -out tmp_${gene}_${LINF}_ORF.aln
#
# get stats how different both seqs are reading in the alignment of both
python3 ~/00scripts/general_scripts/python3/compare_fasta_aln.py --fasta tmp_${gene}_${LINF}_ORF.aln >> old_intron_gene_comp.txt
done < multi_exon_genes_matched.txt

rm tmp_L*
cat LinJ.*_ORF.fa > all_matched_intron_genes.fa
bwa mem ../../refgenomes/TriTrypDB-9.0_LinfantumJPCM5_Genome.fasta all_matched_intron_genes.fa > all_matched_intron_genes.sam
samtools view -bS all_matched_intron_genes.sam > all_matched_intron_genes.bam
samtools sort all_matched_intron_genes.bam all_matched_intron_genes_sort
samtools index all_matched_intron_genes_sort.bam

# !!! muscle alignments to not do the same kind of alignment as expected from mapping against the genome and results are therefore not suitable

# --> conclusion, probably do support more
# I will use the strategy of getting mRNA seqs of each gene, extracting the longest ORF (ignoring introns) and useing this one for each gene

mkdir tmp_multiexon_genes
mv LinJ.* tmp_multiexon_genes/.
mv multi_exon_genes* tmp_multiexon_genes/.
mv all_matched_intron_genes* tmp_multiexon_genes/.
mv old_intron_gene_comp.txt tmp_multiexon_genes/.
mv geneid_mapping_new_old_annot.txt tmp_multiexon_genes/.
mv *_mRNA_multiexon* tmp_multiexon_genes/
cp tmp_multiexon_genes/all_matched_intron_genes*bam* ~/tmp/.


###################
# 2b
# find longest ORF of all genes 

python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_a_find_ORF.py --CDS TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa  --log TriTrypDB-9.0_LinfantumJPCM5_mRNA.log > TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF.fa


# test gene of interest
# grep '^>LinJ.33.3220' -A1 TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa > aa.fa
# python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_a_find_ORF.py --CDS aa.fa --log aa.log > aa_ORF.fa


###################
# 3
# edit SNPs in ORF sequences

cd /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/
# Ldon1_noBPK157A1_44			Ldon1* 44/40
# Ldon5_7						Ldon2 7/6
# Ldon4_noGILANI_18				Ldon3* 18/16
# Ldon4_4						Ldon4
# Ldon3_noLRCL53_7				Ldon5* 7/6
# Linf_noISSextr_noInf152_43	Linf1* 43/40
# TurkeyH_11					CUK.Linf
# CH							CH.Linf
for pop in Ldon1_noBPK157A1_44 Ldon5_7 Ldon4_noGILANI_18 Ldon3_noLRCL53_7 Linf_noISSextr_noInf152_43 Ldon4_4 CH TurkeyH_11
do
for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
bgzip snps.filt.leish_global.linj.complex.LinJ.${chr}_${pop}.vcf
tabix -p vcf snps.filt.leish_global.linj.complex.LinJ.${chr}_${pop}.vcf.gz
done
done

cd /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A15_MDK-test

# python3 ~/00scripts/general_scripts/python3/MDK-test_Ldon_A15/MDK_b_add_vcf_SNPs.py \
# --ORF aplus.fa \
# --ingroup CHLinf,CUKLinf,Linf1_Ldon1,Ldon2,Ldon3,Ldon4,Ldon5 \
# --outgroup Ldon1,Ldon2,Ldon3,Ldon4,Ldon5_CHLinf,CUKLinf,Linf1
# 
# python3 ~/00scripts/general_scripts/python3/MDK-test_Ldon_A15/MDK_b_add_vcf_SNPs.py \
# --ORF aminus.fa \
# --ingroup CHLinf,CUKLinf,Linf1_Ldon1,Ldon2,Ldon3,Ldon4,Ldon5 \
# --outgroup Ldon1,Ldon2,Ldon3,Ldon4,Ldon5_CHLinf,CUKLinf,Linf1

queue=normal
memory=4000

for chr in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
do
grep 'LinJ.'${chr} -A1 TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF.fa > TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF_LinJ.${chr}.fa
#
name=chr_${chr}
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_b_add_vcf_SNPs.py \
--ORF TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF_LinJ.${chr}.fa  \
--ingroup CHLinf,CUKLinf,Linf1_Ldon1,Ldon2,Ldon3,Ldon4,Ldon5 \
--outgroup Ldon1,Ldon2,Ldon3,Ldon4,Ldon5_CHLinf,CUKLinf,Linf1
done





# number of genes
grep '>' TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF.fa | wc -l
# 8234
grep '>' TriTrypDB-9.0_LinfantumJPCM5_mRNA_ORF.fa | sort -u | wc -l
# 8234 
# --> number of genes with ORF
#
grep '>' TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa | wc -l
# 8239
grep '>' TriTrypDB-9.0_LinfantumJPCM5_mRNA.fa | sort -u | wc -l
# 8239
# number of genes with mRNA (with and without ORF)

tail -2 TriTrypDB-9.0_LinfantumJPCM5_mRNA.log 
# Number of genes with no ORF 5
# Number of genes with ORF 8234

mv LinJ*_tmp tmp/.
ls tmp/LinJ*/LinJ* | wc -l
# 8234
ls MK_in_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_out_CHLinf_CUKLinf_Linf1/* | wc -l
# 8234
ls MK_in_CHLinf_CUKLinf_Linf1_out_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5/* | wc -l
# 8234

mkdir tmp
mv LinJ.*_tmp tmp/.



###################
# 4
# run MK test (MK.pl from Holloway et al. 2007)
###################
# focal one species and polarised by all alleles in the other species
queue=normal
memory=4000
# 
name=ingroup_inf
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
perl ~/software/McDonald_Kreitman_script/MK.pl -outfile MKres_ingroup_mRNA_ORF_inf.txt -pol pol \
-ingroup CHLinf_CUKLinf_Linf1 \
-outgroup Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_A \
-outgroup Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_B \
-dir ./MK_in_CHLinf_CUKLinf_Linf1_out_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5/

name=ingroup_don
bsub -q $queue -o ${name}.o -e ${name}.e -J ${name} \
-R"select[mem>${memory}] rusage[mem=${memory}]" -M${memory} \
perl ~/software/McDonald_Kreitman_script/MK.pl -outfile MKres_ingroup_mRNA_ORF_don.txt -pol pol \
-ingroup Ldon1_Ldon2_Ldon3_Ldon4_Ldon5 \
-outgroup CHLinf_CUKLinf_Linf1_A \
-outgroup CHLinf_CUKLinf_Linf1_B \
-dir ./MK_in_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_out_CHLinf_CUKLinf_Linf1/


mv MK_in_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_out_CHLinf_CUKLinf_Linf1/ MK_mRNA_ORF_in_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5_out_CHLinf_CUKLinf_Linf1/
mv MK_in_CHLinf_CUKLinf_Linf1_out_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5/ MK_mRNA_ORF_in_CHLinf_CUKLinf_Linf1_out_Ldon1_Ldon2_Ldon3_Ldon4_Ldon5/




###################
# 6
# add results parameters to MK.pl results
###################
python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_c_complete_MKout.py --MKres MKres_ingroup_mRNA_ORF_inf.txt 

python3 ~/00scripts/general_scripts/python3/leish_donovaniComplex_A15/MDK-test_Ldon_A15/MDK_c_complete_MKout.py --MKres MKres_ingroup_mRNA_ORF_don.txt 

cp MKres_ingroup_mRNA_ORF_* ~/leish_donovaniComplex/A15_MDK-test/.


###################
# 7
# analyse / plot MK results in R
###################



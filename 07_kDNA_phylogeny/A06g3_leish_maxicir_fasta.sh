#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide the name of the sample"';
    exit 1;
fi


RDsub="ldon_complex_maxi"
sample=$1

echo "5 SNPs only to fasta"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.fa\n"
/software/bin/java -Xmx4g -jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar \
-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
-T FastaAlternateReferenceMaker \
-o /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.fa \
--variant /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.vcf \
-IUPAC 1
 
echo "6 SNPs only to fasta, masking fasta for non covered positions across samples"
echo "lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.fa\n"
bedtools maskfasta \
-fi /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.fa \
-bed ~/leish_donovaniComplex/A06_var_SNP_caling/maxicircle_mask_pos.bed \
-fo /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.fa


echo "7 SNPs only to fasta, masked, cutting out 'essentially' the coding region"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.fa\n"
bedtools getfasta \
-fi /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.fa \
-bed ~/leish_donovaniComplex/A06_var_SNP_caling/maxicircle_take_pos.bed \
-fo /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.fa 

	
echo "8 add header to fasta file"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.name.fa\n"
awk -v s=$sample 'BEGIN{OFS=""} {if($1 ~ /^>/) {print ">",s} else {print $0}}' \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.fa \
| fold -w60 > /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.name.fa
#
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.fa\n"
cp -s /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.mask.reg.name.fa \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.fa
	
	
echo "9 cut out 'essentially' the coding region taking UNmasked sampes"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.fa\n"
bedtools getfasta \
-fi /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.fa \
-bed ~/leish_donovaniComplex/A06_var_SNP_caling/maxicircle_take_pos.bed \
-fo /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.fa 

	
echo "10 add header to fasta file taking UNmasked samples"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.name.fa\n"
awk -v s=$sample 'BEGIN{OFS=""} {if($1 ~ /^>/) {print ">",s} else {print $0}}' \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.fa \
| fold -w60 > /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.name.fa
#
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.nomask.fa\n"
cp -s /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.reg.name.fa \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.nomask.fa
	
	
# include indels
echo "11 call indels"
echo "/lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.raw.vcf\n"
/software/bin/java -Xmx4g \
-jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T SelectVariants \
-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
-V /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.vcf -selectType INDEL \
-o /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.raw.vcf 

echo "12 filter indels"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.filtered.vcf\n"
/software/bin/java -Xmx4g \
-jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T VariantFiltration \
-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
-V /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.raw.vcf \
--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "indel_filter" \
-o /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.filtered.vcf
	
	
# add contig name to SNP modified fasta
echo "12a change on SNP .fa reference"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa\n"
echo '>Contig1' |\
cat - /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.fa |\
grep -v '^>1' > /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa 

head /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa
echo "\n"
	
# create dictionary
echo "12b dict on SNP .fa reference"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.dict\n"
samtools-1.3 dict \
-o /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.dict \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa
	

# add contig name to SNP modified fasta
echo "12c fai on SNP .fa reference"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa.fai\n"
samtools-1.3 faidx /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa 
	
	
echo "13 indels to SNP fasta"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.fa\n"
/software/bin/java -Xmx4g -jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar \
-R /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.noindel.var.n.fa \
-T FastaAlternateReferenceMaker \
-o /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.fa \
--variant /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.indels.filtered.vcf \
-IUPAC 1
 
echo "14 cut out 'essentially' the coding region from SNP & indel fasta"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.fa\n"
bedtools getfasta \
-fi /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.fa \
-bed ~/leish_donovaniComplex/A06_var_SNP_caling/maxicircle_take_pos.bed \
-fo /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.fa 

echo "15 change header from from SNP & indel & reg fasta"
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.name.fa\n"
echo '>'${sample} |\
cat /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.fa |\
grep -v '^>1:' | fold -w60 \
> /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.name.fa
#
echo "/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.nomask.indel.fa\n"
cp -s /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.filter.var.indels.filtered.reg.name.fa \
/lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.final.nomask.indel.fa	

rm /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.*var.fa*
rm /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.*mask.fa*
rm /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.*reg.fa*
rm /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.*.n.fa*
rm /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.*red.fa*


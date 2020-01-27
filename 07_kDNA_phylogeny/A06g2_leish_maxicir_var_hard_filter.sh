#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide the name of the sample"';
    exit 1;
fi


RDsub="ldon_complex_maxi"
sample=$1
somy=1

medcov=$(cat ~/leish_donovaniComplex/A06_var_SNP_caling/DP.med.allSamples.txt | awk -v s="$sample" '$1==s {print $2}')

echo 'here'
echo $medcov

echo "applying hard filters ${sample}, median coverage ${medcov}, somy ${somy}"





echo "gatk hard filter variants...."
DPmax=$(echo "1.9 * $medcov" | bc)
DPmin=$(echo "$somy * 4" | bc)
echo "median cov $medcov, DPmax $DPmax, DPmin $DPmin"
/software/bin/java -Xmx4g -jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
--variant /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.vcf \
--filterExpression "QD < 2.0" \
--filterName "QDs2" \
--filterExpression "MQ < 40.0" \
--filterName "MQs40" \
--filterExpression "FS > 13.0" \
--filterName "FSl13" \
--filterExpression "SOR > 4" \
--filterName "SORl4" \
--filterExpression "BaseQRankSum > 3.1 || BaseQRankSum <-3.1" \
--filterName "BQRSls3.1" \
--filterExpression "ClippingRankSum > 3.1 || ClippingRankSum <-3.1" \
--filterName "CRSls3.1" \
--filterExpression "MQRankSum > 3.1 || MQRankSum <-3.1" \
--filterName "MQRSls3.1" \
--filterExpression "ReadPosRankSum > 3.1 || ReadPosRankSum <-3.1" \
--filterName "ReadPosRSls3.1" \
--filterExpression "DP > ${DPmax}" \
--filterName "DPmax1.9medcov" \
--filterExpression "DP < ${DPmin}" \
--filterName "DPmin4somy" \
-o /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.vcf
e1=$?
echo "exit $e1"



echo "exclude indels, non-variant loci and filtered loci (trim remaining alleles by default)"
/software/bin/java -Xmx4g -jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.6/GenomeAnalysisTK.jar \
-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
-T SelectVariants \
--variant /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.vcf \
-o /lustre/scratch118/infgen/team133/sf18//05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.filter.noindel.var.vcf \
--selectTypeToExclude INDEL \
-env \
-ef
e2=$?
echo "exit $e2"


eall=$(($e1 + $e2))
echo "eall $eall"

exit $eall



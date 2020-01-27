#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide arrary of sample names, e.g. "sample1,sample2,sample3"';
    exit 1;
fi
if [ -z $2 ]; then
    echo 'Please provide the somy, e.g. 2"';
    exit 1;
fi

echo "sample $1"
echo "y parameter $2"
# create sample array
IFS=',' read -r -a array <<< "$1"
somy=$2

#queue=yesterday
queue=normal
memory=2000

RDsub="ldon_complex_maxi"




# do gatk realignment
for index in "${!array[@]}"
do
	sample=${array[index]}
	echo $sample
	
	echo "1 index bam"
	samtools-1.3 index /lustre/scratch118/infgen/team133/sf18/022gatk_realign/ldon_complex_${sample}/lmjf.0.8.sorted.markdup.realigned.mq20.PP.bam

	echo "2 gatk var calling creating variants, sample ${sample}"
	name2=2varCalling_${sample}
	bsub -q $queue -o ${name2}.o -e ${name2}.e -J ${name2} -R"select[mem>4000] rusage[mem=4000]" -M4000 \
	/software/bin/java -Xmx4g -jar /lustre/scratch118/infgen/team133/sf18/software/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar -T HaplotypeCaller \
	-R /lustre/scratch118/infgen/team133/sf18/refgenomes/maxicircles/maxicircle_Ldonovani.fa \
	-I /lustre/scratch118/infgen/team133/sf18/022gatk_realign/ldon_complex_${sample}/lmjf.0.8.sorted.markdup.realigned.mq20.PP.bam \
	-o /lustre/scratch118/infgen/team133/sf18/05variants/${RDsub}/variants.ldon_complex.maxi.${sample}.vcf --sample_ploidy $somy -dt NONE --annotateNDA


	echo "3 creation of one col file per sample"
	name3=3vcf_inspection_${sample}
	# 	bsub -w "done(${name2})" -q $queue -o ${name3}.o -e ${name3}.e -J ${name3} 
	bsub -w "done(${name2})" -q $queue -o ${name3}.o -e ${name3}.e -J ${name3} -R"select[mem>4000] rusage[mem=4000]" -M4000 \
	bash ~/00scripts/leish_donovaniComplex/A06f1_aa.sh ${sample}
	

	echo "4 get depth" 
	name4=4depth_${sample}    
	#bsub -w "done(${name3})" -q $queue -o ${name4}.o -e ${name4}.e -J ${name4} \
	bsub -w "done(${name3})" -q $queue -o ${name4}.o -e ${name4}.e -J ${name4} \
	-R"select[mem>4000] rusage[mem=4000]" -M4000 \
	'samtools-1.3 depth /lustre/scratch118/infgen/team133/sf18/022gatk_realign/ldon_complex_'${sample}'/lmjf.0.8.sorted.markdup.realigned.mq20.PP.bam > /lustre/scratch118/infgen/team133/sf18/022gatk_realign/ldon_complex_'${sample}'/ldon.maxi.0.8.sorted.markdup.realigned.mq20.PP.cov'


	

done
	
exit 0


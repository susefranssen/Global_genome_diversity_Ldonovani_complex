#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide arrary of sample names, e.g. "sample1,sample2,sample3"';
    exit 1;
fi


echo "sample $1"

# create sample array
IFS=',' read -r -a array <<< "$1"


#queue=yesterday
queue=normal
memory=2000

RDsub="ldon_complex_maxi"




# do
for index in "${!array[@]}"
do
	sample=${array[index]}
	echo $sample
	
	# echo "4 hard filter SNPs" 
	name4=4hF_${sample}    
	bsub -q $queue -o ${name4}.o -e ${name4}.e -J ${name4} \
	-R"select[mem>4000] rusage[mem=4000]" -M4000 \
	bash ~/00scripts/leish_donovaniComplex/A06g2_leish_maxicir_var_hard_filter.sh $sample

	echo "5 var2fasta"
	name5=5_${sample}
	# -w "done(${name4})"
	bsub -w "done(${name4})" -q $queue -o ${name5}.o -e ${name5}.e -J ${name5} \
	-R"select[mem>4000] rusage[mem=4000]" -M4000 \
	bash ~/00scripts/leish_donovaniComplex/A06g3_leish_maxicir_fasta.sh $sample
	

done
	
exit 0


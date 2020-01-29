

#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide the name of the sample"';
    exit 1;
fi

sample=$1



awk -F '[\t;=]' '{print $1,$4,$10,$14}' /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/linj.${sample}.sorted.markdup.realigned.PP.rmdup.gene.cov > /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/linj.${sample}.sorted.markdup.realigned.PP.rmdup.gene.cov.red








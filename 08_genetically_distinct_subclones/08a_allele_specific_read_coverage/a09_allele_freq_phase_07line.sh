
#!/usr/bin/env bash
if [ -z $1 ]; then
    echo 'Please provide the name of the sample"';
    exit 1;
fi


sample=$1


awk 'BEGIN{OFS="\t"} NR==FNR {f1[$1$2] = $0; next} ($1$2 in f1) {print $0}' /lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/snps.filt.leish_global.linj.complex.vcflike.txt /lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.sort.pileup > \
/lustre/scratch118/infgen/team133/sf18/02X_pileup/leish_donovaniComplex/linj.${sample}.sorted.markdup.realigned.sort.SNPs.pileup
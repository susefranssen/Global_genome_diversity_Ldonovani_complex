#!/usr/bin/env bash

if [ -z $1 ]; then
    echo 'sample"';
    exit 1;
fi

sample=$1

bcftools-1.2 query -f '%CHROM\t%POS\t%AC\t%AF\t%DP\t%QD\t%MQ\t%FS\t%SOR\t%BaseQRankSum\t%ClippingRankSum\t%MQRankSum\t%ReadPosRankSum\n' /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.vcf > /lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.${sample}.col
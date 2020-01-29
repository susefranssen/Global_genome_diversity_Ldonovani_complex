

# Large CNVs and gene CNVs

".sh" scripts on the console were run first but are dependent on the file structure of the server they were run:
**a14_gene_cov.sh** 
which requires **a14_gene_cov_line.sh**
and the R script:
**a10_raw_somies.r**

Those file are required as input with the following R script, however, intermediate results files can also be loaded from within the r script:
**a14_genome_cov.r**  identification and plotting of large indel
**a14_gene_cov.r**    identification and summary stats of gene copy number variation, including for candidate genes


# Genetic variation associated with previously described drug resistance loci

".sh" scripts on the console were run first but are dependent on the file structure of the server they were run:
**a18_variation_known_candidate_loci.sh**

Those file are required as input with the following R script, however, intermediate results files can also be loaded from within the r script:
**a18_variation_known_candidate_loci.r**  plotting of coverage around drug resistance loci in isolates harbouring variation


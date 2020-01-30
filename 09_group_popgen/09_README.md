
# LD by group

The ".sh" script was run first but is dependent on the file structure of the server it was run:

**A04_vcf_LD.sh** using the script **a04_vcflike_to_vcf.py**.

It creates files that are required as input in the following R script, however, intermediate result files are also provided and can also be loaded from within the R script:

**a04_LD_dist_02_groups.r**  within chromosomes, data by group

**a04_LD_dist_03_groups.r**  within chromosomes, merge group data

**a04_LD_dist_04_groups.r**  between chromsomes, plotting including between and within chr LD



# Site frequency spectrum by group and marker SNPs (Venn diagrams)

The ".sh" script was run first but is dependent on the file structure of the server it was run:

**a12_group_freqs_aneuploid_vcf.sh** 

Examples files of this type are provided in:data/SFS_file_examples


Those file are required as input with the following R script:

**a12_group_freqs.r**  








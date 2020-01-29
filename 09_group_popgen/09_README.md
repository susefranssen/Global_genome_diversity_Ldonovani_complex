
# LD by group

".sh" scripts on the console were run first but are dependent on the file structure of the server they were run:

**A04_vcf_LD.sh** using the script **a04_vcflike_to_vcf.py**.

This creates out put files in the following form

"/LinJ.all_",name,"_diploid_genor2_r100_",repl,".interchrom.geno.ld"

Examples files of this type are provided in:data/LD_file_examples


Those file are required as input with the following R script:

**a04_LD_dist_02_groups.r**  within chromosomes, data by group

**a04_LD_dist_03_groups.r**  within chromosomes, merge group data

**a04_LD_dist_04_groups.r**  between chromsomes, plotting including between and within chr LD


# Site frequency spectrum by group and marker SNPs (Venn diagrams)

".sh" scripts on the console were run first but are dependent on the file structure of the server they were run:

**a12_group_freqs_aneuploid_vcf.sh** 

Examples files of this type are provided in:data/SFS_file_examples

Those file are required as input with the following R script:

**a12_group_freqs.r**  








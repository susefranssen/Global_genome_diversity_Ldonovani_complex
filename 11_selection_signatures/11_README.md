

# McDonald Kreitman test

The ".sh" scripts were run first but are dependent on the file structure of the server they were run. These are:

**a15_MDK-test.sh**     which requires several python scripts and a perl script (MK.pl, Holloway et al. 2007)

**a15_MDK_python_scripts/MDK_a_find_ORF.py**

**a15_MDK_python_scripts/MDK_b_add_vcf_SNPs.py**

**a15_MDK_python_scripts/MDK_c_complete_MKout.py**

**a15_MDK_python_scripts/compare_fasta_aln.py**

**a15_MDK_python_scripts/fasta_to_1line.py**

**a15_MDK_python_scripts/get_subset_fasta.py**

**a15_MDK_python_scripts/MK.pl**


They create files that are required as input in the following R script, however, intermediate results files can also be loaded from within the R script:

**a15_MDK-test.r**  



# Annotation and GO enrichment of group-specific marker SNPs

The ".sh" script was run first but is dependent on the file structure of the server it was run:

**a13_SNPeff_fixed_segr_cats_v38.sh**

It creates files that are required as input in the following R script, however, intermediate result files are also provided and can also be loaded from within the R script:

**a13_SNPeff_fixed_segr_cats_v38_GOenrichment.r**  GO enrichment test

**A13_all_REVIGO_treemap_p0.05_topGOweight_v38.r** visualisation of GO enrichment results using REVIGO


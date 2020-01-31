

#step 1 - generate Q and P matrix output files with admixture
#analyis used admixture Version 1.23, from http://software.genetics.ucla.edu/admixture 
admixture --cv=10 input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.cv10.ped 8
admixture --cv=10 input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.cv10.ped 11
admixture --cv=10 input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.cv10.ped 13


#step 2 - generate SVG file from Q matrices and tree
perl tree_to_SVGcircletree.pl --qmatrix=input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.8.Q --qmatrix=input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.11.Q --qmatrix=input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.13.Q --ldonglobal_fix_spp_names --annularsectors --annularsectorgap_p=0 --space_for_bars=40 --pedfile=input.ADMIXTURE.filt.leish.global.complex.linj.diploid.biallelic.gt3.thinned.ped LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_newick.2.mpt --species_name_colors=LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_newick.leaf_names_to_group_and_colours.3.txt --species_name_colors_column=2 --learncolors --also_color_dots --root_branch_length=0.02 --no_negative_branches LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_newick.2.mpt.K8andK11andK13_annular.LEARNCOLS.sortcols.svg


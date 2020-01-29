


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/11_selection_signatures/")
dir.create("REVIGO_GO_res_visualisation", showWarnings = F)
setwd("REVIGO_GO_res_visualisation")

library(treemap)# treemap package by Martijn Tennekes
library(data.table)
library(foreach)

#################
#
# Visualisation of GO results using REVIGO:
# 
# revigo was used with the following parameters:
#   
#   http://revigo.irb.hr/
#   
#   all GO terms with a value for weightGO <0.05 (p-vals given to revive)
# 
# for revigo standard parameters were used
# allowed similarity: medium (0.7)
# show: abs log10 pval
#
#
# input files that were used for REVIGO are:
# 
# imput used for REVIGO were the GO terms with pval<0.05 (weightFish)
# "SNPeff_",comp,"/SNPeff_",comp,"_stats.genes_high_mod_topGO_classic.weight.0.05.txt"
# e.g. "SNPeff_Linf_vs_Ldon/SNPeff_Linf_vs_Ldon_stats.genes_high_mod_topGO_classic.weight.0.05.txt"
#
#################
#

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

# arranging data for revigo first

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0007018","microtubule-based movement",0.287,4.2676,0.617,0.000,"microtubule-based movement"),
                     c("GO:0007059","chromosome segregation",0.476,1.7605,0.717,0.167,"microtubule-based movement"),
                     c("GO:0015698","inorganic anion transport",0.872,1.8380,0.698,0.000,"inorganic anion transport"),
                     c("GO:0055085","transmembrane transport",8.916,2.3261,0.700,0.400,"inorganic anion transport"),
                     c("GO:0051656","establishment of organelle localization",0.180,1.5177,0.709,0.257,"inorganic anion transport"),
                     c("GO:0006974","cellular response to DNA damage stimulus",2.360,2.0650,0.693,0.034,"cellular response to DNA damage stimulus"),
                     c("GO:1902531","regulation of intracellular signal transduction",0.547,1.6770,0.674,0.474,"cellular response to DNA damage stimulus"),
                     c("GO:0006310","DNA recombination",1.641,1.3788,0.703,0.041,"DNA recombination"),
                     c("GO:0006464","cellular protein modification process",7.726,3.1938,0.664,0.510,"DNA recombination"),
                     c("GO:0001522","pseudouridine synthesis",0.350,1.3080,0.656,0.248,"DNA recombination"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="CH.Linf.nohy_3"]
all<-revigo.data

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006468","protein phosphorylation",4.137,4.0132,0.246,0.000,"protein phosphorylation"),
                     c("GO:0052106","quorum sensing involved in interaction with host",0.000,1.6198,0.021,0.021,"quorum sensing involved in interaction with host"),
                     c("GO:0048874","homeostasis of number of cells in a free-living population",0.029,1.4437,0.020,0.692,"quorum sensing involved in interaction with host"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="CUK.Linf_11"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006468","protein phosphorylation",4.137,6.4685,0.388,0.000,"protein phosphorylation"),
                     c("GO:0071897","DNA biosynthetic process",0.676,1.7399,0.374,0.172,"protein phosphorylation"),
                     c("GO:0043087","regulation of GTPase activity",0.510,1.6289,0.447,0.000,"regulation of GTPase activity"),
                     c("GO:0052564","response to immune response of other organism involved in symbiotic interaction",0.010,2.4559,0.444,0.000,"response to immune response of other organism involved in symbiotic interaction"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Ldon1star_44"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006261","DNA-dependent DNA replication",0.576,2.6198,0.624,0.000,"DNA-dependent DNA replication"),
                     c("GO:0006486","protein glycosylation",0.317,1.5058,0.563,0.243,"DNA-dependent DNA replication"),
                     c("GO:0009190","cyclic nucleotide biosynthetic process",0.182,2.4559,0.441,0.170,"DNA-dependent DNA replication"),
                     c("GO:0016310","phosphorylation",7.764,2.1427,0.583,0.398,"DNA-dependent DNA replication"),
                     c("GO:0006555","methionine metabolic process",0.358,2.1367,0.443,0.268,"DNA-dependent DNA replication"),
                     c("GO:0046854","phosphatidylinositol phosphorylation",0.173,1.7423,0.504,0.413,"DNA-dependent DNA replication"),
                     c("GO:0007018","microtubule-based movement",0.287,2.0555,0.606,0.030,"microtubule-based movement"),
                     c("GO:0007205","protein kinase C-activating G-protein coupled receptor signaling pathway",0.018,1.7423,0.645,0.057,"protein kinase C-activating G-protein coupled receptor signaling pathway"),
                     c("GO:0009966","regulation of signal transduction",0.857,1.4377,0.621,0.380,"protein kinase C-activating G-protein coupled receptor signaling pathway"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Ldon2_7"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006468","protein phosphorylation",4.137,6.6198,0.515,0.000,"protein phosphorylation"),
                     c("GO:0006259","DNA metabolic process",5.607,2.9208,0.570,0.301,"protein phosphorylation"),
                     c("GO:0006261","DNA-dependent DNA replication",0.576,1.8356,0.584,0.169,"protein phosphorylation"),
                     c("GO:0016311","dephosphorylation",1.250,1.6440,0.566,0.467,"protein phosphorylation"),
                     c("GO:1902531","regulation of intracellular signal transduction",0.547,2.0809,0.569,0.039,"regulation of intracellular signal transduction"),
                     c("GO:0070887","cellular response to chemical stimulus",1.007,1.6946,0.586,0.433,"regulation of intracellular signal transduction"),
                     c("GO:0032940","secretion by cell",0.763,1.6946,0.638,0.081,"secretion by cell"),
                     c("GO:0007006","mitochondrial membrane organization",0.065,1.6946,0.649,0.153,"secretion by cell"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Ldon3star_18"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006468","protein phosphorylation",4.137,6.7696,0.354,0.000,"protein phosphorylation"),
                     c("GO:0006298","mismatch repair",0.165,1.9101,0.431,0.147,"protein phosphorylation"),
                     c("GO:0001522","pseudouridine synthesis",0.350,1.4237,0.330,0.474,"protein phosphorylation"),
                     c("GO:0015698","inorganic anion transport",0.872,1.4486,0.566,0.000,"inorganic anion transport"),
                     c("GO:0006325","chromatin organization",0.668,2.0044,0.528,0.040,"chromatin organization"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Ldon4_4"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0007018","microtubule-based movement",0.287,4.5528,0.734,0.000,"microtubule-based movement"),
                     c("GO:0006643","membrane lipid metabolic process",0.382,2.0706,0.779,0.171,"microtubule-based movement"),
                     c("GO:0007059","chromosome segregation",0.476,2.3468,0.796,0.167,"microtubule-based movement"),
                     c("GO:0032259","methylation",3.103,1.7773,0.876,0.000,"methylation"),
                     c("GO:0048193","Golgi vesicle transport",0.297,2.0410,0.819,0.000,"Golgi vesicle transport"),
                     c("GO:0055085","transmembrane transport",8.916,2.9586,0.815,0.396,"Golgi vesicle transport"),
                     c("GO:0046903","secretion",0.810,1.5143,0.783,0.269,"Golgi vesicle transport"),
                     c("GO:0048874","homeostasis of number of cells in a free-living population",0.029,2.8539,0.682,0.023,"homeostasis of number of cells in a free-living population"),
                     c("GO:0000723","telomere maintenance",0.133,1.9136,0.554,0.545,"homeostasis of number of cells in a free-living population"),
                     c("GO:0052106","quorum sensing involved in interaction with host",0.000,2.3872,0.695,0.692,"homeostasis of number of cells in a free-living population"),
                     c("GO:0044145","modulation of development of symbiont involved in interaction with host",0.000,2.7447,0.765,0.477,"homeostasis of number of cells in a free-living population"),
                     c("GO:0051276","chromosome organization",1.477,1.4976,0.729,0.592,"homeostasis of number of cells in a free-living population"),
                     c("GO:0006302","double-strand break repair",0.211,1.9318,0.747,0.027,"double-strand break repair"),
                     c("GO:0034470","ncRNA processing",2.222,1.3851,0.743,0.326,"double-strand break repair"),
                     c("GO:0006464","cellular protein modification process",7.726,7.4318,0.741,0.510,"double-strand break repair"),
                     c("GO:0001522","pseudouridine synthesis",0.350,1.7100,0.725,0.205,"double-strand break repair"),
                     c("GO:0016311","dephosphorylation",1.250,1.8041,0.824,0.056,"dephosphorylation"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Ldon5star_7"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006928","movement of cell or subcellular component",0.973,3.5850,0.492,0.000,"movement of cell or subcellular component"),
                     c("GO:0009405","pathogenesis",0.095,2.3990,0.640,0.000,"pathogenesis"),
                     c("GO:0006468","protein phosphorylation",4.137,1.8105,0.493,0.042,"protein phosphorylation"),
                     c("GO:0006650","glycerophospholipid metabolic process",0.543,1.5596,0.284,0.420,"protein phosphorylation"),
                     c("GO:0006643","membrane lipid metabolic process",0.382,1.4372,0.365,0.651,"protein phosphorylation"),
                     c("GO:0046903","secretion",0.810,1.8207,0.566,0.086,"secretion"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Linf_vs_Ldon"]
all<-rbind(all, revigo.data)

revigo.names <- c("term_ID","description","freqInDbPercent","abslog10pvalue","uniqueness","dispensability","representative");
revigo.data <- rbind(c("GO:0006298","mismatch repair",0.165,2.9208,0.373,0.000,"mismatch repair"),
                     c("GO:0006979","response to oxidative stress",0.575,2.5229,0.390,0.509,"mismatch repair"),
                     c("GO:0006928","movement of cell or subcellular component",0.973,1.5482,0.507,0.030,"movement of cell or subcellular component"),
                     c("GO:0055085","transmembrane transport",8.916,2.6990,0.512,0.228,"movement of cell or subcellular component"),
                     c("GO:0060285","cilium-dependent cell motility",0.006,1.3575,0.429,0.131,"movement of cell or subcellular component"))
revigo.data <-data.table(revigo.data); setnames(revigo.data,colnames(revigo.data),revigo.names)
revigo.data[,group:="Linf1star_43"]
all<-rbind(all, revigo.data)
all<-all[order(representative)]
write.table(all, file="A13_all_REVIGO_treemap_p0.05_topGOweight_v38.txt", row.names = F, quote = T)

#--------------------

# this table is equivalent to the table saved just above but two additional columns were added manually 
# that specify the colour that should be used 
dat<-data.table(read.table("A13_all_REVIGO_treemap_p0.05_topGOweight_v38_col_mod.txt", header = T))
# dat[,abslog10pvalue:=abs(log10pvalue)]
dat[,color:=paste0("#",color)]
dat[,representative:=as.character(representative)]
dat[,match:=(term_ID==term_ID.1)]

for (mm in unique(dat$group))
{
  # mm<-"Linf_vs_Ldon"
  print(mm)
  # print(unique(dat[group==mm, color]))
  # print(unique(dat[group==mm, representative]))
  # 
  pdf( file=paste0("revigo_treemap_",mm,".pdf"), width=9, height=5 ) # width and height are in inches
  treemap(
    dat[group==mm],
    index = c("representative","description"),
    vSize = "abslog10pvalue",
    type = "categorical",
    vColor = "representative",
    title = mm,
    inflate.labels = F,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
    # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
    position.legend = "none",
    palette=unique(dat[group==mm, color]),
    fontsize.labels = 18,
    fontsize.title = 26
  )
  dev.off()
}



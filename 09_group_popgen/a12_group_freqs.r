

# cp /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A12_group_freqs/group_freqs_*.all.txt ~/tmp/.

library(data.table)
library(foreach)
library(ggplot2)
library(UpSetR)

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/09_group_popgen/")
dir.create("SFS")
setwd("SFS/")

load("../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")



# NOTE the following if requires input files generated as described in a12_group_freqs_aneuploid_vcf.sh
# an example input is provided for group<-"Ldon1_noBPK157A1_44" in data/SFS_file_examples
if (F)
{
  groups<-c("Ldon1_noBPK157A1_44_s18R1","Ldon1_noBPK157A1_44_s18R2","Ldon1_noBPK157A1_44_s18R3",
            "Linf_noISSextr_noInf152_43_s18R1", "Linf_noISSextr_noInf152_43_s18R2","Linf_noISSextr_noInf152_43_s18R3",
            "Ldon1_noBPK157A1_44",
            "Linf_noISSextr_noInf152_43",
            "Ldon4_noGILANI_18",
            "TurkeyH_11",
            "Ldon5_7",
            "donovani3_8","Ldon3_noLRCL53_7",
            "Ldon4_noGILANI_18_s7R1", "Ldon1_noBPK157A1_44_s7R1", "Linf_noISSextr_noInf152_43_s7R1", "TurkeyH_11_s7R1",
            "Ldon4_noGILANI_18_s7R2",  "Ldon1_noBPK157A1_44_s7R2", "Linf_noISSextr_noInf152_43_s7R2",  "TurkeyH_11_s7R2",
            "Ldon4_noGILANI_18_s7R3",  "Ldon1_noBPK157A1_44_s7R3",  "Linf_noISSextr_noInf152_43_s7R3",  "TurkeyH_11_s7R3",
            "Ldon4_4", "EP_1", "other.Ldon_3", "CH", "CH_hy_2", "CH_nohy_3")
  # groups<-"Ldon1_noBPK157A1_44"
  
  dat.all<-foreach (group = groups, .combine=rbind) %do%
  {
    # group<-"Ldon1_noBPK157A1_44"
    print(group)
    dat<-data.table(read.table(paste0("../data/SFS_file_examples/group_freqs_",group,".all_aneuploid.txt"), header = T))
    dat[, refcount:=round(reffreq*cov)]
    
    if (group=="donovani1_52" | group=="donovani1_7" | group=="donovani1_main42" | group=="Ldon1_new_45" 
        | group=="Ldon1_noBPK157A1_44" | group=="Ldon1_noBPK157A1_44_s7R1" | group=="Ldon1_noBPK157A1_44_s7R2" | group=="Ldon1_noBPK157A1_44_s7R3" 
        | group=="Ldon1_noBPK157A1_44_s18R1"| group=="Ldon1_noBPK157A1_44_s18R2"| group=="Ldon1_noBPK157A1_44_s18R3") {dat[,maingroup:="Ldon1"]}
    if (group=="donovani2a_19" | group=="donovani2a_7" | group=="Ldon4_noGILANI_18" 
        | group=="Ldon4_noGILANI_18_s7R1" | group=="Ldon4_noGILANI_18_s7R2"| group=="Ldon4_noGILANI_18_s7R3") {dat[,maingroup:="Ldon3"]}
    if (group=="donovani2b_7" | group=="Ldon2_new_4") {dat[,maingroup:="Ldon4"]}
    if (group=="donovani3_8" | group=="donovani3_7" | group=="Ldon3_noLRCL53_7") {dat[,maingroup:="Ldon5"]}
    if (group=="infantum_47" | group=="infantum_7" | group=="Linf1_Mon1_31" | group=="Linf1_China_7" 
        | group=="Linf1_ChinaUzb_10" | group=="Linf1_nonMon1_5" | group=="Linf1_Mon1_noISSm_31"
        | group=="Linf_noISSextr_noInf152_43" | group=="Linf1_Mon1_31_s5R1" | group=="Linf_noISSextr_noInf152_43_s7R1" | group=="Linf_noISSextr_noInf152_43_s7R2"
        | group=="Linf_noISSextr_noInf152_43_s7R3"
        | group=="Linf_noISSextr_noInf152_43_s18R1"| group=="Linf_noISSextr_noInf152_43_s18R2"| group=="Linf_noISSextr_noInf152_43_s18R3") {dat[,maingroup:="Linf1"]}
    if (group=="Ldon5_7") {dat[,maingroup:="Ldon2"]}
    if (group=="TurkeyH_11" | group=="TurkeyH_11_s7R1" | group=="TurkeyH_11_s7R2"| group=="TurkeyH_11_s7R3") {dat[,maingroup:="CUK_Linf"]}
    if (group=="Ldon4_4") {dat[,maingroup:="Ldon4"]}
    if (group=="EP_1") {dat[,maingroup:="EP"]}
    if (group=="other.Ldon_3") {dat[,maingroup:="other_Ldon"]}
    if (group=="CH" | group=="CH_hy_2" | group=="CH_nohy_3") {dat[,maingroup:="CH_Linf"]}
    
    dat[,group.id:=group]
    
  }
  
  # renaming groups
  dat.all[group.id=="Ldon1_new_45", group:="Ldon1_45"]
  dat.all[group.id=="Ldon1_noBPK157A1_44", group:="Ldon1*_44"]
  dat.all[group.id=="Ldon1_noBPK157A1_44_s18R1", group:="Ldon1*_18R1"]
  dat.all[group.id=="Ldon1_noBPK157A1_44_s18R2", group:="Ldon1*_18R2"]
  dat.all[group.id=="Ldon1_noBPK157A1_44_s18R3", group:="Ldon1*_18R3"]
  #
  dat.all[group.id=="infantum_47", group:="Linf1_47"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43", group:="Linf1*_43"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s18R1", group:="Linf1*_18R1"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s18R2", group:="Linf1*_18R2"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s18R3", group:="Linf1*_18R3"]
  #
  dat.all[group.id=="donovani2a_19", group:="Ldon3_19"]
  dat.all[group.id=="Ldon4_noGILANI_18", group:="Ldon3*_18"]
  #
  dat.all[group.id=="TurkeyH_11", group:="CUK.Linf_11"]
  #
  dat.all[group.id=="Ldon5_7", group:="Ldon2_7"]
  #
  dat.all[group.id=="donovani3_8", group:="Ldon5_8"]
  dat.all[group.id=="Ldon3_noLRCL53_7", group:="Ldon5*_7"]
  #
  dat.all[group.id=="Ldon1_noBPK157A1_44_s7R1", group:="Ldon1*_7R1"]
  dat.all[group.id=="Ldon1_noBPK157A1_44_s7R2", group:="Ldon1*_7R2"]
  dat.all[group.id=="Ldon1_noBPK157A1_44_s7R3", group:="Ldon1*_7R3"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s7R1", group:="Linf1*_7R1"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s7R2", group:="Linf1*_7R2"]
  dat.all[group.id=="Linf_noISSextr_noInf152_43_s7R3", group:="Linf1*_7R3"]
  dat.all[group.id=="Ldon4_noGILANI_18_s7R1", group:="Ldon3*_7R1"]
  dat.all[group.id=="Ldon4_noGILANI_18_s7R2", group:="Ldon3*_7R2"]
  dat.all[group.id=="Ldon4_noGILANI_18_s7R3", group:="Ldon3*_7R3"]
  dat.all[group.id=="TurkeyH_11_s7R1", group:="CUK.Linf_7R1"]
  dat.all[group.id=="TurkeyH_11_s7R2", group:="CUK.Linf_7R2"]
  dat.all[group.id=="TurkeyH_11_s7R3", group:="CUK.Linf_7R3"]
  #
  dat.all[group.id=="Ldon4_4", group:="Ldon4_4"]
  dat.all[group.id=="EP_1", group:="other.Linf_1"]
  dat.all[group.id=="other.Ldon_3", group:="other.Ldon_3"]
  dat.all[group.id=="CH", group:="CH.Linf_5"]
  dat.all[group.id=="CH_hy_2", group:="CH.Linf.hy_2"]
  dat.all[group.id=="CH_nohy_3", group:="CH.Linf.nohy_3"]
  dat.all[majallelefreq<0.5, majallelefreq:=1-majallelefreq]
  
  #####
  
  save(dat.all, file="SFS.dat.all.RData")
}

load(file="SFS.dat.all.RData") # dat.all
dat.all[,SNP.id:=paste0(chr,"_",pos)]

dat.all.poly<-dat.all[majallelefreq!=1] # only polymorphic sites in the respective population

dir.create("../SFS_res", showWarnings = F)

pdf(paste0("../SFS_res/Minfreq_6main_groups.pdf"), width=7, height=5)
ggplot(dat.all.poly[group.id %in% c("Ldon1_noBPK157A1_44","Linf_noISSextr_noInf152_43","Ldon4_noGILANI_18","TurkeyH_11","Ldon5_7","Ldon3_noLRCL53_7")],
       aes(1-majallelefreq)) + 
  geom_histogram(aes(fill=group), bins = 25) +
  facet_wrap(~group, scales = "free") +
  xlim(0,0.52) +
  scale_fill_manual(values=unique(sample.info[,.(groups, group.col)])[order(groups),group.col][c(2,3,4,5,7,8)]) +
  guides(fill=FALSE) + labs(x="Minor allele frequency", y="Count")
dev.off()

pdf(paste0("../SFS_res/Minfreq_18sreps_maingroups.pdf"))
ggplot(dat.all.poly[group.id %in% c("Ldon1_noBPK157A1_44_s18R1","Ldon1_noBPK157A1_44_s18R2","Ldon1_noBPK157A1_44_s18R3",
                               "Linf_noISSextr_noInf152_43_s18R1","Linf_noISSextr_noInf152_43_s18R2","Linf_noISSextr_noInf152_43_s18R3")],
       aes(1-majallelefreq)) + 
  geom_histogram(aes(fill=maingroup), bins = 25) +
  facet_wrap(~group, scales = "free") +
  xlim(0,0.52) +
  scale_fill_manual(values=unique(sample.info[,.(groups, group.col)])[order(groups),group.col][c(3,8)]) +
  guides(fill=FALSE) + labs(x="Minor allele frequency", y="Count")
dev.off()

pdf(paste0("../SFS_res/Minfreq_7sreps_maingroups.pdf"), width=7, height=6)
aa<-dat.all.poly[group.id %in% c("Ldon4_noGILANI_18_s7R1", "Ldon1_noBPK157A1_44_s7R1", "Linf_noISSextr_noInf152_43_s7R1", "TurkeyH_11_s7R1", 
                            "Ldon4_noGILANI_18_s7R2",  "Ldon1_noBPK157A1_44_s7R2", "Linf_noISSextr_noInf152_43_s7R2",  "TurkeyH_11_s7R2",
                            "Ldon4_noGILANI_18_s7R3",  "Ldon1_noBPK157A1_44_s7R3",  "Linf_noISSextr_noInf152_43_s7R3",  "TurkeyH_11_s7R3")]
ggplot(aa, aes(1-majallelefreq)) + 
  geom_histogram(aes(fill=maingroup), bins = 25) +
  facet_wrap(~group, scales = "free", ncol = 3) +
  xlim(0,0.52) +
  scale_fill_manual(values=unique(sample.info[,.(groups, group.col)])[order(groups),group.col][c(2,3,5,8)]) +
  guides(fill=FALSE) + labs(x="Minor allele frequency", y="Count")
dev.off()


pdf(paste0("../SFS_res/Minfreq_all_groups.pdf"), width=7, height=5)
ggplot(dat.all.poly[group.id %in% c("Ldon1_noBPK157A1_44","Linf_noISSextr_noInf152_43","Ldon4_noGILANI_18","TurkeyH_11","Ldon5_7","Ldon3_noLRCL53_7",
                               "Ldon4_4", "EP_1", "other.Ldon_3", "CH","CH_hy_2","CH_nohy_3")],
       aes(1-majallelefreq)) + 
  geom_histogram(aes(fill=group), bins = 25) +
  facet_wrap(~group, scales = "free") +
  xlim(0,0.52) +
  scale_fill_manual(values=unique(sample.info[,.(groups, group.col)])[order(groups),group.col][c(1,1,1,2,3,4,5,6,7,8,9,10)]) +
  guides(fill=FALSE) + labs(x="Minor allele frequency", y="Count")
dev.off()


# getting minor allele counts corresponding to the frequency plots
dat.all.poly[,mincount:=ifelse(reffreq<0.5,refcount,cov-refcount)]
pdf("../SFS_res/dat.all.poly_Ldon1*_7RX_mincounts.pdf")
par(mfrow=c(3,1))
hist(dat.all.poly[group=="Ldon1*_7R1" ,mincount], nclass=100)
hist(dat.all.poly[group=="Ldon1*_7R2" ,mincount], nclass=100)
hist(dat.all.poly[group=="Ldon1*_7R3" ,mincount], nclass=100)
dev.off()



#########################
# SNPs venn diagram

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/09_group_popgen/")
dir.create("SNP_Venns", showWarnings = F)
setwd("SNP_Venns/")
dir.create("../SNP_Venns_res", showWarnings = F)


dat.mainX<-dat.all[group.id %in% c("Ldon1_noBPK157A1_44","Linf_noISSextr_noInf152_43","Ldon4_noGILANI_18","TurkeyH_11","Ldon5_7",
                                   "Ldon3_noLRCL53_7", "Ldon4_4", "CH_nohy_3")][,.(chr,pos,group,reffreq)]
dat.mainX[,SNP.id:=paste0(chr,"_",pos)]
dat.mainX
# get markers that are fixed for different alleles between groups
dat.mainX.wide<-data.table(dcast(dat.mainX, formula = chr + pos + SNP.id ~ group, value.var = "reffreq"))
# NOTE: SNPS with NA in any group are excluded
dat.mainX.fixed<-dat.mainX.wide[CH.Linf.nohy_3 %in% 0:1 & CUK.Linf_11 %in% 0:1 & get("Ldon1*_44") %in% 0:1 & Ldon2_7 %in% 0:1 & get("Ldon3*_18")
               %in% 0:1 & Ldon4_4 %in% 0:1 & get("Ldon5*_7") %in% 0:1 & get("Linf1*_43") %in% 0:1]
dat.mainX.fixed[,group.sum:=CH.Linf.nohy_3 +CUK.Linf_11 +get("Ldon1*_44") +Ldon2_7 +get("Ldon3*_18") +Ldon4_4 +get("Ldon5*_7") +get("Linf1*_43")]
#
# number of SNPs fixed in each group but polymorphic between at least two groups
nrow(dat.mainX.fixed[group.sum!=0 & group.sum!=8,])
# [1] 154,422
#
# SNPs fixed in each group, where at least two groups share the same allele
dat.mainX.fixed.aa<-dat.mainX.fixed[group.sum!=0 & group.sum!=1 & group.sum!=7 & group.sum!=8] # removing invariant SNPs and those specific to one group
png("../SNP_Venns_res/Fixed_main8_groupSmarker.png", width=1100, height=800)
upset(dat.mainX.fixed.aa[,3:11], nsets = 8, order.by = "freq", nintersects = 10, text.scale=2) # pdf saved manually
dev.off()
# [1] 57,232
#
# save SNP ids of markers for Linf vs Ldon
setnames(dat.mainX.fixed.aa,c("chr","pos"),c("#CHROM","POS"))
write.table(dat.mainX.fixed.aa[group.sum==3 & CH.Linf.nohy_3==1 & CUK.Linf_11==1 & `Linf1*_43`==1], file="../SNP_Venns_res/dat.mainX.fixed_Linf_vs_Ldon.txt",
            quote = F, row.names = F)
#
# marker SNPs for exactly one group (SNPs fixd in each group only)
dat.mainX.fixed.bb<-dat.mainX.fixed[group.sum==7 | group.sum==1] 
png("../SNP_Venns_res/Fixed_main8_groupmarker.png", width=1100, height=800)
upset(dat.mainX.fixed.bb[,2:11], nsets = 9, order.by = "freq", nintersects = 16, text.scale=2) # pdf saved manually
dev.off()
nrow(dat.mainX.fixed.bb)
# [1] 97,190
# 
# save SNP ids of fixed markers for individual groups 
setnames(dat.mainX.fixed.bb,c("chr","pos"),c("#CHROM","POS"))
write.table(dat.mainX.fixed.bb[`Ldon5*_7`==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_Ldon5star_7.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[CH.Linf.nohy_3==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_CH.Linf.nohy_3.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[`Ldon3*_18`==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_Ldon3star_18.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[`Ldon1*_44`==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_Ldon1star_44.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[Ldon2_7==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_Ldon2_7.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[`Linf1*_43`==1 & group.sum==1], file="../SNP_Venns_res/dat.mainX.fixed_Linf1star_43.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[Ldon4_4==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_Ldon4_4.txt", quote = F, row.names = F)
write.table(dat.mainX.fixed.bb[CUK.Linf_11==0 & group.sum==7], file="../SNP_Venns_res/dat.mainX.fixed_CUK.Linf_11.txt", quote = F, row.names = F)





# look which SNPs are segregating in which populations
dat.all.mainX.poly<-dat.all[group.id %in% c("Ldon1_noBPK157A1_44","Linf_noISSextr_noInf152_43","Ldon4_noGILANI_18","TurkeyH_11","Ldon5_7",
                                            "Ldon3_noLRCL53_7", "Ldon4_4", "CH_nohy_3")][majallelefreq!=1] # only polymorphic sites in the respective population
dat.all.mainX.poly[,pres:=1] # presence indicator 
dat.mainX.poly.wide<-data.table(dcast(dat.all.mainX.poly, formula = chr + pos + SNP.id ~ group, value.var = "pres")) # not presense combinations will be indicated as NA
dat.mainX.poly.wide[is.na(dat.mainX.poly.wide)]<-0 # NAs replaced by 0
# sum in how many groups an allele is segregating
dat.mainX.poly.wide[,group.sum:=CH.Linf.nohy_3 +CUK.Linf_11 +get("Ldon1*_44") +Ldon2_7 +get("Ldon3*_18") +Ldon4_4 +get("Ldon5*_7") +get("Linf1*_43")]
# in how many groups is an allele seregating?
table(dat.mainX.poly.wide$group.sum)
# 1      2      3      4      5      6 
# 184625  17882   1551    131     21      4 
# number of SNPs segregating in at last one group
nrow(dat.mainX.poly.wide)
# [1] 204214
# SNPs segregating in the maximum number of groups (5-6)
png("../SNP_Venns_res/Segr_main8_maxshared.png", width=1100, height=800)
upset(dat.mainX.poly.wide[,3:11], nsets = 8, order.by = c("degree"), nintersects = 17, text.scale=2) # pdf saved manually
dev.off()

segr.ids<-dat.mainX.poly.wide[group.sum>=5,SNP.id]
write.table(dat.mainX.wide[SNP.id %in% segr.ids], file="../SNP_Venns_res/Segr_main8_maxshared.txt", quote = F, row.names = F)

#############################
# SNP counts only

table(dat.all.poly$group)
# aneuploid genotypes
# CH.Linf_5   CH.Linf.hy_2 CH.Linf.nohy_3    CUK.Linf_11   CUK.Linf_7R1   CUK.Linf_7R2   CUK.Linf_7R3    Ldon1*_18R1    Ldon1*_18R2 
# 57063          29341           1675          10524          10348          10363          10296          23052          23851 
# Ldon1*_18R3      Ldon1*_44     Ldon1*_7R1     Ldon1*_7R2     Ldon1*_7R3        Ldon2_7      Ldon3*_18     Ldon3*_7R1     Ldon3*_7R2 
# 1182          24837            683          23305            685          45595          47027          38780          35832 
# Ldon3*_7R3        Ldon4_4        Ldon5_8       Ldon5*_7    Linf1*_18R1    Linf1*_18R2    Linf1*_18R3      Linf1*_43     Linf1*_7R1 
# 38638          52582          16978           9681          27336          27949          27237          33774          14525 
# Linf1*_7R2     Linf1*_7R3   other.Ldon_3   other.Linf_1 
# 2881          15515          87845          30345
#
# # diploid
# # CH.Linf_5   CH.Linf.hy_2 CH.Linf.nohy_3    CUK.Linf_11   CUK.Linf_7R1   CUK.Linf_7R2   CUK.Linf_7R3    Ldon1*_18R1    Ldon1*_18R2 
# # 57058          29341           1674          10522          10346          10361          10294          23049          23847 
# # Ldon1*_18R3      Ldon1*_44     Ldon1*_7R1     Ldon1*_7R2     Ldon1*_7R3        Ldon2_7      Ldon3*_18     Ldon3*_7R1     Ldon3*_7R2 
# # 1182          24830            683          23303            685          45587          47007          38767          35821 
# # Ldon3*_7R3        Ldon4_4        Ldon5_8       Ldon5*_7    Linf1*_18R1    Linf1*_18R2    Linf1*_18R3      Linf1*_43     Linf1*_7R1 
# # 38622          52561          16973           9681          27331          27940          27233          33764          14524 
# # Linf1*_7R2     Linf1*_7R3   other.Ldon_3   other.Linf_1 
# # 2881          15509          87785          30345


snp.summary<-data.table(melt(table(dat.all.poly$group)))
setnames(snp.summary,colnames(snp.summary),c("group","n.SNP"))
write.table(snp.summary, file="../SNP_Venns_res/SNP.summary.txt", quote=FALSE, row.names=F)



snp.summary[,Group:=unlist(lapply(strsplit(as.character(snp.summary$group), split = "_"), function(x) x[1]))]
a<-unlist(lapply(strsplit(as.character(snp.summary$group), split = "_"), function(x) x[2]))
snp.summary[,size:=unlist(lapply(strsplit(as.character(a), split = "R"), function(x) x[1]))]
snp.summary[,rep:=unlist(lapply(strsplit(as.character(a), split = "R"), function(x) x[2]))]

snp.sub<-snp.summary[size %in% c(1,3,4,5,7,11,18,43,44)]# & !(size==18 & rep %in% 1:3)]
snp.sub[,n.SNP.mean:=round(mean(n.SNP)), by=.(Group,size)]
snp.sub[,n.SNP.sd:=round(sd(n.SNP)), by=.(Group,size)]
# get genome size
genome<-data.table(read.table("../TriTrypDB-38_LinfantumJPCM5_Genome.fasta.fai"))
sum(genome[V1!="LinJ.00", V2])
# [1] 31,924,954 <- size of the used genome
snp.sub[,d.SNP.mean_1Mb:=n.SNP.mean/31.924954]
snp.sub[,d.SNP.mean_10kb:=n.SNP.mean/3192.4954]
snp.sub[,d.SNP.mean_1kb:=n.SNP.mean/31924.954]
snp.sub[,size:=as.numeric(size)]


snp.sub.1<-unique(snp.sub[,.(Group,size,n.SNP.mean,n.SNP.sd,d.SNP.mean_1Mb,d.SNP.mean_10kb,d.SNP.mean_1kb)])
write.table(snp.sub.1, file="../SNP_Venns_res/SNP.summary.sub.1.txt", quote=FALSE, row.names=F)

pdf(paste0("../SNP_Venns_res/SNP.stats_perGroup_all.pdf"), width=6.5, height=5)
ggplot(snp.sub.1, aes(x=size, y=n.SNP.mean, color=Group)) +
  geom_errorbar(aes(ymin=n.SNP.mean-n.SNP.sd, ymax=n.SNP.mean+n.SNP.sd), width=12, position=position_dodge(width=1), cex=0.5) +
  geom_line(linetype="dashed") +
  geom_point(size=3, position=position_dodge(width=1)) +
  # scale_y_continuous(breaks=seq(0,120000,10000)) +
  scale_colour_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(1,1,2,3,4,5,6,7,8,9,10), group.col]) +
  labs(x="Sample size", y="SNP number [mean across three replicates]") +
  scale_y_continuous(sec.axis = sec_axis(~./3192.4954, name = "Mean SNP number per 10 kb]")) +
  theme_bw()
dev.off()  



pdf(paste0("../SNP_Venns_res/SNP.stats_perGroup_main8.pdf"), width=6.5, height=5)
ggplot(snp.sub.1[size %in% c(4,7,11,18,43,44) | Group=="CH.Linf.nohy"], aes(x=size, y=n.SNP.mean, color=Group)) +
  geom_errorbar(aes(ymin=n.SNP.mean-n.SNP.sd, ymax=n.SNP.mean+n.SNP.sd), width=12, position=position_dodge(width=1), cex=0.5) +
  geom_line(linetype="dashed") +
  geom_point(size=3, position=position_dodge(width=1)) +
  # scale_y_continuous(breaks=seq(0,120000,10000)) +
  scale_colour_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(1,2,3,4,5,6,7,8), group.col]) +
  labs(x="Sample size", y="SNP number [mean across three replicates]") +
  scale_y_continuous(sec.axis = sec_axis(~./3192.4954, name = "Mean SNP number per 10 kb]")) +
  theme_bw()
dev.off()  


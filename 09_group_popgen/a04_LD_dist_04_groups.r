

# LD between chromosomes



setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/09_group_popgen/")
dir.create("LD_between_chr")
setwd("LD_between_chr/")

library(data.table)
library(ggplot2)
library(foreach)
library(igraph)


# calc window means
# NOTE: this condition can only be run for samples, where the raw data is also provided: c("Ldon1_noBPK157A1_44_s7R1","Ldon1_noBPK157A1_44_s7R2","Ldon1_noBPK157A1_44_s7R3")
# However, the resulting full data produced by this step can be loaded after the if condition
if(F) 
{
  samples1=c("Ldon1_noBPK157A1_44","Ldon5_7","Ldon4_noGILANI_18","Ldon3_noLRCL53_7","Linf_noISSextr_noInf152_43","TurkeyH_11")
  samples2=c("Ldon1_noBPK157A1_44_s18R1","Ldon1_noBPK157A1_44_s18R2","Ldon1_noBPK157A1_44_s18R3",
            "Linf_noISSextr_noInf152_43_s18R1","Linf_noISSextr_noInf152_43_s18R2","Linf_noISSextr_noInf152_43_s18R3")
  samples3=c("Ldon4_noGILANI_18_s7R1", "Ldon1_noBPK157A1_44_s7R1", "Linf_noISSextr_noInf152_43_s7R1", "TurkeyH_11_s7R1",
             "Ldon4_noGILANI_18_s7R2",  "Ldon1_noBPK157A1_44_s7R2", "Linf_noISSextr_noInf152_43_s7R2",  "TurkeyH_11_s7R2",
             "Ldon4_noGILANI_18_s7R3",  "Ldon1_noBPK157A1_44_s7R3",  "Linf_noISSextr_noInf152_43_s7R3",  "TurkeyH_11_s7R3")
  # samples<-c(samples1,samples2, samples3)
  samples <- c("Ldon1_noBPK157A1_44_s7R1","Ldon1_noBPK157A1_44_s7R2","Ldon1_noBPK157A1_44_s7R3")
  
  LD.inter.chr <-foreach (repl = 1:2, .combine=rbind) %do%
  {
    print(repl)
    LD.inter.chr <-foreach (name = samples, .combine=rbind) %do%
    {
      print(name)
      dat<-data.table(read.table(paste0("../../data/LD_file_examples/LinJ.all_",name,"_diploid_genor2_r100_",repl,".interchrom.geno.ld"), head=T))
      
      dat[,meanR2:=mean(R.2, na.rm = T), by=c("CHR1","CHR2")]
      dat[,sd:=sd(R.2, na.rm = T), by=c("CHR1","CHR2")]
      dat[,meanR2.all:=mean(R.2, na.rm = T)]
      dat[,sd.all:=sd(R.2, na.rm = T)]
      dat.mean<-unique(dat[,.(CHR1,CHR2,meanR2,sd,meanR2.all,sd.all)])
      dat.mean[,group:=name]
      dat.mean[,rep:=repl]
    }
  }
  
  # update the group name to the latest version
  LD.inter.chr[group=="Ldon1_noBPK157A1_44", Group:="Ldon1*"] # * star indicating when putatively mixed samples are removed
  LD.inter.chr[group=="Ldon5_7", Group:="Ldon2"]
  LD.inter.chr[group=="Ldon4_noGILANI_18", Group:="Ldon3*"] # * star indicating when putatively mixed samples are removed
  LD.inter.chr[group=="Ldon3_noLRCL53_7", Group:="Ldon5*"] # * star indicating when putatively mixed samples are removed
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43", Group:="Linf1*"] # * star indicating when putatively mixed samples are removed
  LD.inter.chr[group=="TurkeyH_11", Group:="CUK.Linf"]
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s18R1", Group:="Ldon1*_18R1"]
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s18R2", Group:="Ldon1*_18R2"]
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s18R3", Group:="Ldon1*_18R3"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s18R1", Group:="Linf1*_18R1"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s18R2", Group:="Linf1*_18R2"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s18R3", Group:="Linf1*_18R3"]
  #
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s7R1", Group:="Ldon1*_7R1"]
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s7R2", Group:="Ldon1*_7R2"]
  LD.inter.chr[group=="Ldon1_noBPK157A1_44_s7R3", Group:="Ldon1*_7R3"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s7R1", Group:="Linf1*_7R1"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s7R2", Group:="Linf1*_7R2"]
  LD.inter.chr[group=="Linf_noISSextr_noInf152_43_s7R3", Group:="Linf1*_7R3"]
  LD.inter.chr[group=="Ldon4_noGILANI_18_s7R1", Group:="Ldon3*_7R1"]
  LD.inter.chr[group=="Ldon4_noGILANI_18_s7R2", Group:="Ldon3*_7R2"]
  LD.inter.chr[group=="Ldon4_noGILANI_18_s7R3", Group:="Ldon3*_7R3"]
  LD.inter.chr[group=="TurkeyH_11_s7R1", Group:="CUK.Linf_7R1"]
  LD.inter.chr[group=="TurkeyH_11_s7R2", Group:="CUK.Linf_7R2"]
  LD.inter.chr[group=="TurkeyH_11_s7R3", Group:="CUK.Linf_7R3"]
  LD.inter.chr[,Group.rep:=paste0(Group,"_",rep)]
  LD.inter.chr<-LD.inter.chr[order(Group,rep, CHR1, CHR2)]
  
  save(LD.inter.chr, file="LD.inter.chr.RData")
}


load(file="LD.inter.chr.RData") # LD.inter.chr
#
load("../../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData") # sample.info


means<-unique(LD.inter.chr[,.(meanR2.all, sd.all, group, Group, rep, Group.rep)])
write.table( means, file="LD_between_chr_mean.txt",
             col.names = T, row.names = F, quote = F)


pdf("LD_between_chr_mean.pdf",width=7, height=6)
ggplot(means,
       aes(x=Group, y=meanR2.all, color=Group, shape=as.factor(rep))) + 
  geom_point(size=3, position = position_dodge(width=0.6)) +
  geom_errorbar(aes(ymin=meanR2.all-sd.all, ymax=meanR2.all+sd.all), width=.2, position = position_dodge(width=0.6)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(2,2,2,2,3,3,3,3,3,3,3,4,5,5,5,5,7,8,8,8,8,8,8,8), group.col]) +
  geom_hline(yintercept=c(0,1))
dev.off()



##################
# LD within a chromosome


dir.create("../LD_res", showWarnings = F)
rep=1

means.lr <-data.table(read.table("../LD_within_chr/longrangeLD.all.txt", header=T))
means.lr[group=="Ldon1_noBPK157A1_44", Group:="Ldon1*"] # * star indicating when putatively mixed samples are removed
means.lr[group=="Ldon5_7", Group:="Ldon2"]
means.lr[group=="Ldon4_noGILANI_18", Group:="Ldon3*"] # * star indicating when putatively mixed samples are removed
means.lr[group=="Ldon3_noLRCL53_7", Group:="Ldon5*"] # * star indicating when putatively mixed samples are removed
means.lr[group=="Linf_noISSextr_noInf152_43", Group:="Linf1*"] # * star indicating when putatively mixed samples are removed
means.lr[group=="TurkeyH_11", Group:="CUK.Linf"]
means.lr[group=="Ldon1_noBPK157A1_44_s18R1", Group:="Ldon1*_18R1"]
means.lr[group=="Ldon1_noBPK157A1_44_s18R2", Group:="Ldon1*_18R2"]
means.lr[group=="Ldon1_noBPK157A1_44_s18R3", Group:="Ldon1*_18R3"]
means.lr[group=="Linf_noISSextr_noInf152_43_s18R1", Group:="Linf1*_18R1"]
means.lr[group=="Linf_noISSextr_noInf152_43_s18R2", Group:="Linf1*_18R2"]
means.lr[group=="Linf_noISSextr_noInf152_43_s18R3", Group:="Linf1*_18R3"]
#
means.lr[group=="Ldon1_noBPK157A1_44_s7R1", Group:="Ldon1*_7R1"]
means.lr[group=="Ldon1_noBPK157A1_44_s7R2", Group:="Ldon1*_7R2"]
means.lr[group=="Ldon1_noBPK157A1_44_s7R3", Group:="Ldon1*_7R3"]
means.lr[group=="Linf_noISSextr_noInf152_43_s7R1", Group:="Linf1*_7R1"]
means.lr[group=="Linf_noISSextr_noInf152_43_s7R2", Group:="Linf1*_7R2"]
means.lr[group=="Linf_noISSextr_noInf152_43_s7R3", Group:="Linf1*_7R3"]
means.lr[group=="Ldon4_noGILANI_18_s7R1", Group:="Ldon3*_7R1"]
means.lr[group=="Ldon4_noGILANI_18_s7R2", Group:="Ldon3*_7R2"]
means.lr[group=="Ldon4_noGILANI_18_s7R3", Group:="Ldon3*_7R3"]
means.lr[group=="TurkeyH_11_s7R1", Group:="CUK.Linf_7R1"]
means.lr[group=="TurkeyH_11_s7R2", Group:="CUK.Linf_7R2"]
means.lr[group=="TurkeyH_11_s7R3", Group:="CUK.Linf_7R3"]
means.lr[,Group.rep:=paste0(Group,"_",rep)]


#
means[,distbinkb:=125]
setnames(means, c("meanR2.all","sd.all"), c("mean","sd"))

means.all<-rbind(means[rep==1,.(Group,mean,sd,distbinkb)], means.lr[,.(Group,mean,sd,distbinkb)])
means.all[,maingroup:=Group]
for (gg in c("Ldon1","Ldon2","Ldon3","Ldon5","Linf1","CUK.Linf"))
{
  means.all[grep(gg,Group),maingroup:=gg]
}
means.all[,n.samp:=0] # number of samples
for (ss in c(18,7))
{
  means.all[grep(paste0("_",ss),Group),n.samp:=ss]
}
# adding samples size for the remaining ones
means.all[Group=="CUK.Linf",n.samp:=11]
means.all[Group=="Ldon1*",n.samp:=44]
means.all[Group=="Ldon2",n.samp:=7]
means.all[Group=="Ldon3*",n.samp:=18]
means.all[Group=="Ldon5*",n.samp:=7]
means.all[Group=="Linf1*",n.samp:=43] 
means.all[,rep:=unlist(lapply(strsplit(means.all$Group,"_"), function(x) ifelse(length(x)==2, substr(x[[2]], nchar(x[[2]])-1, nchar(x[[2]])), "all")))] # this one contains 18 samples
# means.all[rep=="nf", rep:="all"]
means.all[,rep_size:=paste0(rep,"_",n.samp)]
means.all[,Group_size:=paste0(Group,"_",n.samp)]


pdf("../LD_res/LD_between_chr_longrange_mean_orig_noerrorbars.pdf",width=8, height=4)
ggplot(means.all[Group %in% c("Ldon1*","Ldon2","Ldon3*","Ldon5*","Linf1*","CUK.Linf")],
       aes(x=distbinkb, y=mean, shape=rep_size, color=Group)) + 
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=12, position=position_dodge(width=2), cex=0.05) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  theme_bw() +
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(2,3,4,5,7,8), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean for 1 kb windows]") +
  scale_y_continuous(breaks = 0:5 *0.2) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(values=c(15:18,0:2)) 
dev.off()

pdf("../LD_res/LD_between_chr_longrange_mean_s18_noerrorbars.pdf",width=8, height=4)
ggplot(means.all[n.samp>=18],
       aes(x=distbinkb, y=mean, shape=rep_size, color=maingroup)) + 
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=12, position=position_dodge(width=2), cex=0.05) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  theme_bw() +
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(3,5,8), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean for 1 kb windows]") +
  scale_y_continuous(breaks = 0:5 *0.2) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(values=c(15:17,0:2)) 
dev.off()


pdf("../LD_res/LD_between_chr_longrange_mean_s7_noerrorbars.pdf",width=8, height=4)
ggplot(means.all[n.samp>18 | n.samp==11| n.samp==7 | Group=="Ldon3*"],
       aes(x=distbinkb, y=mean, shape=rep_size, color=maingroup)) + 
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=12, position=position_dodge(width=2), cex=0.05) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  theme_bw() + 
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(2,3,4,5,7,8), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean for 1 kb windows]") +
  scale_y_continuous(breaks = 0:5 *0.2) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(values=c(15:18,25,0:2)) 
dev.off()




# calculate mean and sd of the means between the 3 independent replicates
means.all[,acReps.count:=length(rep), by=.(maingroup,n.samp,distbinkb)] # number of reps in the specific group / across reps (as a test)
means.all[,acReps.mean:=mean(mean), by=.(maingroup,n.samp,distbinkb)] # window mean of means across reps
means.all[,acReps.sd:=sd(mean), by=.(maingroup,n.samp,distbinkb)] # window sd of means across reps
means.all[,group:=maingroup]
means.all[maingroup %in% c("Ldon1","Ldon3","Ldon5","Linf1"),group:=paste0(maingroup,"*")]

# means across replicates data table
means.all.acReps<-unique(means.all[,.(distbinkb,maingroup,group,n.samp,acReps.count,acReps.mean,acReps.sd)])
means.all.acReps[is.na(acReps.sd), acReps.sd:=0]

write.table(means.all.acReps, file="../LD_res/means.all.acReps.txt", col.names = T, row.names = F, quote = F)

pdf("../LD_res/LD_between_chr_longrange_mean_s7_noerrorbars_Ldon1*.pdf",width=8, height=4)
ggplot(means.all[n.samp==7 & group=="Ldon1*"],
       aes(x=distbinkb, y=mean, shape=rep_size, color=group, linetype=rep_size)) + 
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=12, position=position_dodge(width=2), cex=0.05) +
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  theme_bw() + 
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(3), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean for 1 kb windows]") +
  scale_y_continuous(breaks = 0:5 *0.2) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(values=c(15:18,25,0:2)) 
dev.off()

pdf("../LD_res/LD_between_chr_longrange_mean_s7_acReps.pdf",width=8, height=4)
ggplot(means.all.acReps[n.samp==7],
       aes(x=distbinkb, y=acReps.mean, color=group, shape=as.factor(acReps.count))) + 
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  geom_errorbar(aes(ymin=acReps.mean-acReps.sd, ymax=acReps.mean+acReps.sd), width=12, position=position_dodge(width=2), cex=0.1) +
  theme_bw() + 
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(2,3,4,5,7,8), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean +/-sd across 3 replicates]") +
  scale_y_continuous(breaks = 0:5 *0.2, limits = c(-.2,1.15)) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(name  ="# subsampled \nreps of size 7",values=c(15:18,25,0:2)) 
dev.off()



pdf("../LD_res/LD_between_chr_longrange_mean_s18_acReps.pdf",width=8, height=4)
ggplot(means.all.acReps[n.samp==18],
       aes(x=distbinkb, y=acReps.mean, color=group, shape=as.factor(acReps.count))) + 
  geom_line(position=position_dodge(width=2)) +
  geom_point(size=4,position=position_dodge(width=2)) +
  geom_errorbar(aes(ymin=acReps.mean-acReps.sd, ymax=acReps.mean+acReps.sd), width=12, position=position_dodge(width=2), cex=0.1) +
  theme_bw() + 
  scale_color_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][c(3,5,8), group.col]) +
  xlab("Pair-wise SNP distance [kb]") + ylab("LD estimate [mean +/-sd across 3 replicates]") +
  scale_y_continuous(breaks = 0:5 *0.2, limits = c(-.2,1.15)) + scale_x_continuous(breaks = c(0:10 *10, 125)) +
  geom_hline(yintercept = c(0,1), size=0.5) +
  scale_shape_manual(name  ="# subsampled \nreps of size 18",values=c(15:18,25,0:2)) 
dev.off()

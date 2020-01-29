

library(data.table)
library(ggplot2)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/11_selection_signatures/")
dir.create("MDK_test", showWarnings = F)
setwd("MDK_test")



for (comp in c("inf","don"))
{
  dat<-data.table(read.table(paste0("MKres_ingroup_mRNA_ORF_",comp,"_table.txt"), header = T))
  dat[,neglogFETpval:=-log(FETpval)]
  
  dat[,P:=Pn+Ps] # number of polymorphisms
  dat[,D:=Dn+Ds] # number of fixed differences
  # from Molecular Population Genetics , Mathew D. Hahn
  # you need at least 4 P and 4 D to be able to get pval  <0.05
  dat[,test:=P >=2 & D>=2]
  dat[test==T, FDR:=p.adjust(FETpval)]
  
  # pdf(paste0("MKres_ingroup_ORF_",comp,"_NI_FET.pdf"),width=5,height=4)
  # gg<-ggplot(dat, aes(neglogNI, neglogFETpval, col=test)) + 
  #   geom_point(shape=21) +
  #   geom_hline(yintercept = -log(c(0.05))) +
  #   theme_bw() +
  #   scale_color_manual(values=c("gray", "black")) 
  # plot(gg)
  # dev.off()
  
  pdf(paste0("MKres_ingroup_ORF_",comp,"_NI_FET_all.pdf"),width=4,height=4)
  gg<-ggplot(dat, aes(neglogNI, neglogFETpval)) + 
    geom_point(shape=21) +
    geom_hline(yintercept = -log(c(0.05))) +
    theme_bw() +
    scale_color_manual(values=c("gray", "black")) 
  plot(gg)
  dev.off()
  
  # pdf(paste0("MKres_ingroup_ORF_",comp,"_NI_FET_filt.pdf"),width=4,height=4)
  # gg<-ggplot(dat[test==T], aes(neglogNI, neglogFETpval)) + 
  #   geom_point(shape=21) +
  #   geom_hline(yintercept = -log(c(0.05))) +
  #   theme_bw() +
  #   scale_color_manual(values=c("gray", "black")) #+
  #   # geom_text(aes(label=ifelse(dat[test==T]$FETpval<= 0.05,as.character(dat[test==T]$gene),""), hjust=0.2, vjust=-1)) 
  # plot(gg)
  # dev.off()
  
  info<-paste("Number of gene with at least 2 fixed and 2 polymorphic sites:",nrow(dat[test==T]), sep=" ")
  write.table(info, file=paste0("MKres_ingroup_ORF_",comp,"_cand_info.txt"), quote = F, row.names = F)
  
  dat.sig<-dat[FETpval<0.05, .(gene,Pn,Dn,Ps,Ds,P,D,MKcodons,NI,neglogNI,alpha,FETpval,neglogFETpval,test,FDR)]
  dat.sig<-dat.sig[order(-alpha)]
  
  write.table(dat.sig, file=paste0("MKres_ingroup_ORF_",comp,"_cand.txt"), quote = F, row.names = F)
  
  comp
  nrow(dat[NI==1])
  nrow(dat[NI>1])
  nrow(dat[NI<1])
  nrow(dat)
  nrow(dat[NI==1])/nrow(dat)
  nrow(dat[NI>1])/nrow(dat)
  nrow(dat[NI<1])/nrow(dat)

  # > comp
  # [1] "don"
  # >   nrow(dat[NI==1])
  # [1] 1513
  # >   nrow(dat[NI>1])
  # [1] 3542
  # >   nrow(dat[NI<1])
  # [1] 3179
  # >   nrow(dat)
  # [1] 8234
  # >   nrow(dat[NI==1])/nrow(dat)
  # [1] 0.1837503
  # >   nrow(dat[NI>1])/nrow(dat)
  # [1] 0.4301676
  # >   nrow(dat[NI<1])/nrow(dat)
  # [1] 0.3860821
  
  # > comp
  # [1] "inf"
  # >   nrow(dat[NI==1])
  # [1] 2447
  # >   nrow(dat[NI>1])
  # [1] 3550
  # >   nrow(dat[NI<1])
  # [1] 2237
  # >   nrow(dat)
  # [1] 8234
  # >   nrow(dat[NI==1])/nrow(dat)
  # [1] 0.2971824
  # >   nrow(dat[NI>1])/nrow(dat)
  # [1] 0.4311392
  # >   nrow(dat[NI<1])/nrow(dat)
  # [1] 0.2716784
  # 4/8234
  # [1] 0.0004857906
  # 12/8234
  # [1] 0.001457372
}


# all genes with mod and high effect variants differentially fixed between species
# rom: TabS4_SNPeff_high_mod_GO_bp_all.xlsx
diff_species<-c("LinJ.05.0760","LinJ.09.0150","LinJ.13.1390","LinJ.14.1130","LinJ.14.1190","LinJ.16.1680","LinJ.17.0890","LinJ.20.0700","LinJ.21.1280","LinJ.22.0780","LinJ.22.0930","LinJ.23.1570","LinJ.24.1470","LinJ.27.2460","LinJ.28.3110","LinJ.31.0310","LinJ.33.2690","LinJ.33.3170","LinJ.34.4090","LinJ.04.0710","LinJ.13.0930","LinJ.13.0940","LinJ.34.2950","LinJ.35.4250","LinJ.05.0760","LinJ.09.0150","LinJ.13.1390","LinJ.14.1130","LinJ.14.1190","LinJ.16.1680","LinJ.17.0890","LinJ.20.0700","LinJ.21.1280","LinJ.22.0780","LinJ.22.0930","LinJ.23.1570","LinJ.24.1470","LinJ.27.1530","LinJ.27.2460","LinJ.28.3110","LinJ.31.0310","LinJ.33.2690","LinJ.33.3170","LinJ.34.4090","LinJ.04.0710","LinJ.29.1320","LinJ.34.2950")

xdon<-data.table(read.table("MKres_ingroup_ORF_don_cand.txt",header = T))
don<-unlist(lapply(strsplit(as.character(xdon[,gene]), split = "-"), function(x) x[1]))

xinf<-data.table(read.table("MKres_ingroup_ORF_inf_cand.txt",header = T))
inf<-unlist(lapply(strsplit(as.character(xinf[,gene]), split = "-"), function(x) x[1]))

intersect(don,inf)
# [1] "LinJ.36.4770" "LinJ.23.1090"
intersect(don,diff_species)
# character(0)
intersect(inf,diff_species)
# character(0)




args <- commandArgs(TRUE)

# comma separated list of functions to perform
name <- args[1]
#"infantum_47" "donovani1_52" "donovani2a_19" "donovani2b_7" "donovani3_8" "TurkeyH_11"
#name="donovani1_52"
binsize <- as.numeric(args[2])
# 100 1000 5000
#binsize=1000
minind <- as.numeric(args[3])
#minind=40
mincountsbin=0
mincountsbin <- as.numeric(args[4])
maxplotdist <- as.numeric(args[5])


# name="Ldon1_noBPK157A1_44"; minind=40
binsize=1000; mincountsbin=200; maxplotdist=100000


maxdist="1000000"


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/09_group_popgen/LD_within_chr")

library(data.table)
library(ggplot2)
library(foreach)



# note intermediate data files are only available for sample "Ldon1_noBPK157A1_44", 
# to run this part for all files LD values for all groups firs have to be generated as detailed in A04_vcf_LD.sh
# otherwise the results file for this section is already provided 09_group_popgen/LD_within_chr/longrangeLD.all.txt
if (F)  
{
  samples1=c("Ldon1_noBPK157A1_44", 
             "Ldon5_7", 
             "Ldon4_noGILANI_18", 
             "Ldon3_noLRCL53_7", 
             "Linf_noISSextr_noInf152_43",
             "TurkeyH_11")
  samples2=c("Ldon1_noBPK157A1_44_s18R1","Ldon1_noBPK157A1_44_s18R2","Ldon1_noBPK157A1_44_s18R3",
             "Linf_noISSextr_noInf152_43_s18R1","Linf_noISSextr_noInf152_43_s18R2","Linf_noISSextr_noInf152_43_s18R3")
  samples3=c("Ldon4_noGILANI_18_s7R1", "Ldon1_noBPK157A1_44_s7R1", "Linf_noISSextr_noInf152_43_s7R1", "TurkeyH_11_s7R1",
             "Ldon4_noGILANI_18_s7R2",  "Ldon1_noBPK157A1_44_s7R2", "Linf_noISSextr_noInf152_43_s7R2",  "TurkeyH_11_s7R2",
             "Ldon4_noGILANI_18_s7R3",  "Ldon1_noBPK157A1_44_s7R3",  "Linf_noISSextr_noInf152_43_s7R3",  "TurkeyH_11_s7R3")
  samples<-c(samples1,samples2, samples3)
  # samples<-"Ldon1_noBPK157A1_44"
  
  longrangeLD.all <-foreach (name =samples, .combine=rbind) %do%
  {
    print(name)
    
    if (name %in% c("Ldon1_noBPK157A1_44","Linf_noISSextr_noInf152_43")) {
      minind=40
    } else if (name %in% c("Ldon5_7", "Ldon3_noLRCL53_7")) {
      minind=6
    } else if (name %in% c("Ldon4_noGILANI_18")) {
      minind=16
    } else if (name %in% c("TurkeyH_11")) {
      minind=9
    } else if (length(grep("s7R", name)) ==1) {
      minind=6
    } else {
      minind=16
    }
    
    
    load(file=paste0("Boxplot_",name,"_diploid_genor2_meanR2_bin",binsize,"_minind",minind,"_maxdist",maxplotdist,"_mincount",mincountsbin,".RData"))
    dat.all<-data.table(dat.all)
    
    # do manually when all dat.all data has been saved
    longrangeLD <- foreach (ddistbinkb =c(1, 1:10*10), .combine=rbind) %do%
    {
      print(ddistbinkb)
      
      distXkb<-dat.all[distbinkb==ddistbinkb]
      distXkb[,mean:=mean(R.2)]
      distXkb[,sd:=sd(R.2)]
      distXkb<-unique(distXkb[,.(distbinkb,mean,sd)])
      distXkb[,group:=name]
      distXkb[,width:=1000]
      
      distXkb
    }
    longrangeLD
  }
  write.table(longrangeLD.all, file="longrangeLD.all.txt", quote=FALSE, col.names=T, row.names=F)
} else {
  longrangeLD.all<-data.table(read.table("longrangeLD.all.txt", header=T))
}




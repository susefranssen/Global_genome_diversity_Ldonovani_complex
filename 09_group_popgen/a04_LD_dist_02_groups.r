
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



# name="Ldon1_noBPK157A1_44"; minind=40;   binsize=1000; mincountsbin=200; maxplotdist=100000

maxdist="1000000"

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/09_group_popgen/")
dir.create("LD_within_chr", showWarnings = F)
setwd("LD_within_chr")

library(data.table)
library(ggplot2)
library(foreach)


snpnum_group <- c()

#dat.all <- foreach(chr=c(paste0("0",1:3)), .combine=rbind) %do%
dat.all <- foreach(chr=c(paste0("0",1:9),10:36), .combine=rbind) %do%
{
  # chr="01"
  print(chr)
  dat <- NULL
  dat <- data.table(read.table(paste0("../../data/LD_file_examples/LinJ.",chr,"_",name,"_diploid_genor2_maxd",maxdist,".txt.geno.ld"), header=T))
  print(nrow(dat))

  snpnum_group <- rbind( snpnum_group, c( name, chr, length(sort(unique(union(dat$POS1,dat$POS2)))) ))

  dat[,dist:=POS2-POS1]
  dat[,distbin:=dist %/% binsize]
  dat[,distbinkb:=round(distbin*binsize/1000,2)]
  dat <- dat[N_INDV>=minind]
  dat <- dat[!is.na(R.2)]

  dat[,binsampnum:=length(R.2),by=distbinkb]

  dat[distbinkb<=maxplotdist/1000,]
}

f.name=paste0("Boxplot_",name,"_diploid_genor2_meanR2_bin",binsize,"_minind",minind,"_maxdist",maxplotdist,"_mincount",mincountsbin,".RData")
save(dat.all, file=f.name)
# load(file=f.name)



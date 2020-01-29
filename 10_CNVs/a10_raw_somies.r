
library(data.table)
library(foreach)
library(ggplot2)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/10_CNVs/")

#---------------------
# functions


# get deviance (error) based on frequency data
# peakx: frequencies where the peaks are located
getdeviance <- function(peakx, somy)
{
  freqs<-1/somy * 1:(somy-1)
  
  deviance <-foreach (peakx.f=peakx, .combine="+") %do%
  {
    # peakx.f=peakx[1]
    peakx.f.i <-which.min(abs(peakx.f - freqs))
    sqrt(abs(freqs[peakx.f.i] - peakx.f))
  }

  deviance / length(peakx)
}


#-------------------------------------
load("../05_heterozygosities/hets_sample_heterozygosities.RData")
sample.het<-names(sort(hets[hets>0.004], decreasing = T))
sample.het<-sub("X","", sample.het)
sample.het<-sub(".","-", sample.het, fixed=T)

# 
# # load all somy data estimated by Caroline
# all<-foreach (s = names(hets), .combine=rbind) %do% # go through all samples
# {
#   s.mod<-sub("X","", s)
#   s.mod<-sub(".","-", s.mod, fixed=T)
#   print(paste(s.mod,s))
#   dat<-data.table(read.table(paste0("/lustre/scratch118/infgen/team133/cd16/ploidy/leish_global/Somy.summarystats.linj.",s.mod,".maxsomy10.txt"),
#                              head=T))
#   # dat<-data.table(read.table(paste0(path,"/Somy.summarystats.linj.",s.mod,".maxsomy10.txt"),head=T))
#   if(ncol(dat)==6) {dat[,Somy.mod:=0]}
#   dat[,sample:=s]
#   dat[,sample_noX:=s.mod]
#   
#   dat
# }
# 
# all[Somy==3 & sample %in% c("EP","CH32","CH34","GE"),] # AAAAAAA
# 
# 
# all[,Deviance.sum:=sum(Deviance), by=.(sample)]
# all[,Deviance.max:=max(Deviance), by=.(sample)]
# all[,Somy.mod.T:=sum(Somy.mod)>1, by=.(sample)]
# all[, het:=F]
# 
# # sample.het<-names(sort(hets[hets>0.004], decreasing = T))
# all[sample_noX %in% sample.het, het:=T]
# 
# save(all, file="somy_all.RData")
load(file="somy_all.RData")

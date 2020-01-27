


setwd("~/leish_donovaniComplex/A06_var_SNP_caling/")

library(data.table)
library(ggplot2)
library(foreach)
library(plyr)

lay_out = function(...) {    
  x <- list(...)
  n <- max(sapply(x, function(x) max(x[[2]])))
  p <- max(sapply(x, function(x) max(x[[3]])))
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(n, p)))    
  
  for (i in seq_len(length(x))) {
    print(x[[i]][[1]], vp = grid::viewport(layout.pos.row = x[[i]][[2]], 
                                           layout.pos.col = x[[i]][[3]]))
  }
} 




# samples.all <-as.character(unlist(read.table("/lustre/scratch118/infgen/team133/sf18/01fastq/leish_donovaniComplex/sampIDs.leish.global.all.txt")))
# samples.all <- setdiff(samples.all,c("OVN3","CL-SL"))
samples.all <-as.character(unlist(read.table("~/leish_donovaniComplex/A06_var_SNP_caling/takesamples.minmed20.txt")))
samples.all<-c(samples.all,"Lamazonensis_A","Lamazonensis_B")


DP.max=1.9 # multiplied by median
QD.min=10 #3
MQ.min=40
FS.max=13
SOR.max=4
RS.t=3.1 # RankSum threshold applied to all 4 RankSum tests

DP.med.all <- data.table(chr=paste0("Contig",1))

dat.all <- foreach (sample=samples.all, .combine=rbind) %do%
{
  print(sample)
  dat<-NULL
  dat=data.table(read.table(paste0("/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.",sample,".col")))
  setnames(dat, old=names(dat), new=c("chr","pos","AC","AF","DP","QD","MQ","FS","SOR","BaseQRS","ClippingRS","MQRS","ReadPosRS"))
  print(head(dat))
  xx <- ddply(dat, "chr", summarize, med = median(DP))
  DP.med.all[, samp:=xx$med]; setnames(DP.med.all, old="samp", new=sample)
  
  dat[,sample.id:=sample]
}

write.table(t(DP.med.all), file=paste0("DP.med.allSamples.txt"), col.names = T, row.names = T, quote = F)

#save(dat.all, file = "dat.all.RData")
#load("dat.all.RData")

#### density


# MQ
g1 <- ggplot(dat.all, aes(MQ, colour=sample.id)) + 
  geom_density(adjust=20) +
  #geom_histogram(binwidth=2) +
  ggtitle(paste0("MQ (",MQ.min,", red)")) +
  geom_vline(xintercept=MQ.min, col="red") + guides(colour=FALSE)
#print(g1)

# FS
g2 <- ggplot(dat.all, aes(FS, colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.5)+
  ggtitle(paste0("FS (",FS.max,", red)")) +
  geom_vline(xintercept=FS.max, col="red") + guides(colour=FALSE)
#print(g2)

# SOR
g3 <- ggplot(dat.all, aes(SOR, colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.5)+
  ggtitle(paste0("SOR (",SOR.max,", red)")) +
  geom_vline(xintercept=SOR.max, col="red")
#print(g3)

pdf(paste0("MQ_FS_SOR_allSamples.pdf"),width=24,height=5)
lay_out(list(g1, 1, 1:3), list(g2, 1, 4:6), list(g3, 1, 7:10))
dev.off()



# BaseQRS
g4 <- ggplot(dat.all, aes(as.numeric(as.character(dat.all$BaseQRS)), colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.25) +
  geom_vline(xintercept=c(-RS.t,RS.t), col="red") + guides(colour=FALSE)
#print(g4)

# ClippingRS
g5 <- ggplot(dat.all, aes(as.numeric(as.character(dat.all$ClippingRS)), colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.25) +
  geom_vline(xintercept=c(-RS.t,RS.t), col="red") + guides(colour=FALSE)
#print(g5)

# MQRS
g6 <- ggplot(dat.all, aes(as.numeric(as.character(dat.all$MQRS)), colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.25) +
  geom_vline(xintercept=c(-RS.t,RS.t), col="red") + guides(colour=FALSE)
#print(g6)

# ReADPOSRS
g7 <- ggplot(dat.all, aes(as.numeric(as.character(dat.all$ReadPosRS)), colour=sample.id)) + 
  geom_density() +
  #geom_histogram(binwidth=0.25) +
  geom_vline(xintercept=c(-RS.t,RS.t), col="red")
#print(g7)

pdf(paste0("BaseQ_Clipping_MQ_ReadPos.RS_allsamples.pdf"),width=34,height=5)
lay_out(list(g4, 1, 1:3), list(g5, 1, 4:6), list(g6, 1, 7:9), list(g7, 1, 10:13))
dev.off()


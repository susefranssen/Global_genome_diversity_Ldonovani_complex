


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/07_kDNA_phylogeny/")

library(data.table)
library(ggplot2)
library(foreach)
library(plyr)


samples.all <-as.character(unlist(read.table("/lustre/scratch118/infgen/team133/sf18/01fastq/leish_donovaniComplex/sampIDs.leish.global.all.txt")))
samples.all <- setdiff(samples.all,c("OVN3","CL-SL")) # exclude samples not yet processed
samples.all <- setdiff(samples.all,c("TH5")) # exclude samples with no mappiing

tag=""
tag="_addLamazon"
if (tag=="_addLamazon") {samples.all<-c(samples.all,"Lamazonensis_A","Lamazonensis_B")}

# dat.all <- foreach (sample=samples.all, .combine=rbind) %do%
# {
#   print(sample)
#   dat=data.table(read.table(paste0("/lustre/scratch118/infgen/team133/sf18/022gatk_realign/ldon_complex_",sample,"/ldon.maxi.0.8.sorted.markdup.realigned.mq20.PP.cov")))
#   setnames(dat,colnames(dat),c("contig","pos","cov"))
#   dat[,samp:=sample]
# #   print(head(dat))
#   
#   dat
# }
# 
# save(dat.all, file=paste0("ldon.maxi.0.8.sorted.markdup.realigned.mq20.PP_cov",tag,".RData"))
load(file=paste0("/Volumes/sf18/leish_donovaniComplex/A06_var_SNP_caling/ldon.maxi.0.8.sorted.markdup.realigned.mq20.PP_cov",tag,".RData"))

dat.all[,type:="ingroup"]
dat.all[samp %in% c("LmexU1103_v1","Lamazonensis_A","Lamazonensis_B","P283","CL1_aethiopica","CL4_aethiopica"),type:="outgroup"]
dat.all<-dat.all[order(type,samp)]

pdf(paste0("maxicicle.cov.boxplot",tag,".pdf"),width=40,height=6)
ggplot(dat.all, aes(x=samp, y=cov)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

pdf(paste0("maxicicle.cov.boxplot",tag,"_1row.pdf"),width=20,height=10)
ggplot(dat.all[samp %in% unique(dat.all[,samp])[1:77]], aes(x=samp, y=cov)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,6000)
dev.off()
pdf(paste0("maxicicle.cov.boxplot",tag,"_2row.pdf"),width=20,height=10)
ggplot(dat.all[samp %in% unique(dat.all[!samp %in% c("CL1_aethiopica","CL4_aethiopica"),samp])[77:153]], aes(x=samp, y=cov)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0,6000)
dev.off()



dat.all[,cov:=as.double(cov)]
dat.all[,samp:=as.character(samp)]
dat.all[,med.cov:=median(cov),by=samp]
med.cov<-unique(dat.all[,.(samp,med.cov)])[order(samp)]




png(paste0("maxicicle.cov.alongchr.1-60",tag,".png"),width=2000,height=600)
p <- ggplot(dat.all[samp %in% ss[1:60]], aes(pos, cov, col=samp)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()

png(paste0("maxicicle.cov.alongchr.60-120",tag,".png"),width=2000,height=600)
p <- ggplot(dat.all[samp %in% ss[61:length(ss)]], aes(pos, cov, col=samp)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p)
dev.off()


# take only samples where the medium coverage is at least 20
ss<-sort(med.cov[med.cov>=20]$samp)
write.table(ss,file="Samples_medCov_atleast20.txt", quote = F, row.names = F)

length(unique(dat.all[samp %in% ss & type =="ingroup", samp]))
# [1] 116

# take only samples where the medium coverage is at least 20
dat.ss<-dat.all[samp %in% ss]
dat.ss[,min.pos.cov:=min(cov),by=pos]

min.cov.ss<-unique(dat.ss[,.(pos,min.pos.cov)])


pdf(paste0("maxicicle.cov.alongchr.takesamples.mincov",tag,".pdf"),width=16,height=5)
p <- ggplot(min.cov.ss, aes(pos, min.pos.cov)) +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_vline(xintercept = c(min(min.cov.ss[min.pos.cov>=10]$pos), max(min.cov.ss[min.pos.cov>=10]$pos)), col="red") # 985, 17204
print(p)
dev.off()

min.cov.ss
write.table(min.cov.ss[min.pos.cov>=10]$pos, file=paste0("positions.takesamples.mincov10",tag,".txt"), quote=FALSE, col.names=F, row.names=F)
write.table(ss, file=paste0("takesamples.minmed20",tag,".txt"), quote=FALSE, col.names=F, row.names=F)

max(min.cov.ss[min.pos.cov>=20]$pos)
# [1] 17204 --> last position to use

mask <-min.cov.ss[min.pos.cov<10]
mask[,contig:="1"]
mask[,start:=pos-1]
mask[,end:=pos]

region<- data.table(contig="1", start=min(min.cov.ss[min.pos.cov>=10]$pos)-1, end=max(min.cov.ss[min.pos.cov>=10]$pos)-1)

write.table(mask[,.(contig,start,end)], file=paste0("maxicircle_mask_pos",tag,".bed"), quote=FALSE, col.names=F, row.names=F, sep = "\t")
write.table(region, file=paste0("maxicircle_take_pos",tag,".bed"), quote=FALSE, col.names=F, row.names=F, sep = "\t")




library(foreach)
library(data.table)
library(reshape)
library(ggplot2)
library(stringr)

setwd("~/leish_donovaniComplex/A06_var_SNP_caling/")
load(file="~/leish_donovaniComplex/A01_StAMPP/addoutgroup/sample.info.RData")

sample.info$sample<-sub("X","", sample.info$sample)
sample.info$sample<-sub(".","-", sample.info$sample, fixed=T)

samp.use<-c()
for (s in sample.info$sample) # of the other samples the filtered vcf do not exist ... need to check why
{
  if(file.exists(paste0("/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.",
                        s,".filter.noindel.var.vcf")) )
  {
    x <- try(read.table(paste0("/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.",
                               s,".filter.noindel.var.vcf")) )
    if(inherits(x, "try-error"))
    {
      print(paste("no SNPs for sample",s))
    } else {
      samp.use<-c(samp.use,s)
    } 
  }
}


all <-foreach (samp = samp.use, .combine=rbind) %do%
{
  #   samp="MAM"
  print(samp)
  dat <-data.table(read.table(paste0("/lustre/scratch118/infgen/team133/sf18/05variants/ldon_complex_maxi/variants.ldon_complex.maxi.",samp,".filter.noindel.var.vcf")))
  setnames(dat, colnames(dat),c("chr","pos","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","data"))
  dat[,data:=as.character(data)]
  
  dat[,sample:=samp]
  dat[,geno.str:=unlist(lapply(strsplit(dat[,data], split = ":"), function(x) x[1]))]
  somy=str_count(dat$geno.str[1],"/")+1
  dat[,ref.g.c:=str_count(geno.str,"0")] # reference allel count in genotype 
  dat[,alt.g.c:=str_count(geno.str,"1")] # alternate allel count in genotype 
  dat[ref.g.c==somy,geno.code:=0]
  dat[alt.g.c==somy,geno.code:=2]
  dat[ref.g.c!=somy & alt.g.c!=somy,geno.code:=1]
  
  dat[,AD:=unlist(lapply(strsplit(dat[,data], split = ":"), function(x) x[2]))] # allele depth
  dat[,ref.count:=unlist(lapply(strsplit(dat[,AD], split = ","), function(x) x[1]))] # ref allele count
  dat[,alt.count:=unlist(lapply(strsplit(dat[,AD], split = ","), function(x) x[2]))] # alternate allele count
  dat[,cov:=unlist(lapply(strsplit(dat[,data], split = ":"), function(x) x[3]))] # coverage
  dat[,ref.count:=as.numeric(ref.count)]
  dat[,alt.count:=as.numeric(alt.count)]
  dat[,cov:=as.numeric(cov)]
  dat[,GQ:=unlist(lapply(strsplit(dat[,data], split = ":"), function(x) x[4]))] # genotype quality
  dat[,ID:=NULL]
  dat[,REF:=NULL]
  dat[,ALT:=NULL]
  dat[,QUAL:=NULL]
  dat[,FILTER:=NULL]
  dat[,INFO:=NULL]
  dat[,FORMAT:=NULL]
  dat[,AD:=NULL]

  dat
}
all[,ref.freq:=ref.count/cov]



pdf("SNPfreq2_all.pdf", width=11,height=8)
ggplot(all, aes(x=ref.count/cov)) + geom_histogram() +
  facet_wrap( ~ sample, scales = "free_y") + xlim(-0.05,1.05) +
  geom_vline(xintercept = 0.5, colour="blue", size=0.4) +
  geom_vline(xintercept = 0.25, colour="orange", size=0.4) +
  geom_vline(xintercept = 0.75, colour="orange", size=0.4) +
  geom_vline(xintercept = c(1/3), colour="green", size=0.4) +
  geom_vline(xintercept = c(2/3), colour="green", size=0.4) +
  xlab("Reference allele frequency")
dev.off()

aa<-all[ref.freq>0.1]

pdf("SNPfreq2_all_polyloci.pdf", width=11,height=8)
ggplot(aa, aes(x=ref.count/cov)) + geom_histogram() +
  facet_wrap( ~ sample, scales = "free_y") + xlim(0,1) +
  geom_vline(xintercept = 0.5, colour="blue", size=0.4) +
  geom_vline(xintercept = 0.25, colour="orange", size=0.4) +
  geom_vline(xintercept = 0.75, colour="orange", size=0.4) +
  geom_vline(xintercept = c(1/3), colour="green", size=0.4) +
  geom_vline(xintercept = c(2/3), colour="green", size=0.4) +
  xlab("Reference allele frequency")
dev.off()


write.table(all, file="SNP_sample_stats.txt", quote = F, row.names = F)



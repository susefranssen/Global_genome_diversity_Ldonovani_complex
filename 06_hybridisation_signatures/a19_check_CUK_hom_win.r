

library(data.table)
library(foreach)
library(ggplot2)
library(reshape2)
library(patchwork)
library(StAMPP)
# library(stringr)


#----------------------------------



setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/06_hybridisation_signatures/")
dir.create("CUK_hom_win", showWarnings = F)
setwd("CUK_hom_win/")


# data table specifiying homozygous and heterozygous fractions of the CUK samples (polarised on reference JPCM5)
# information / file obtained from the original paper on the CUK samples, i.e. Imamura et al. 2016
dat<-data.table(read.table("ALL.5000.windows..skip.txt.NEW.2"))
setnames(dat, colnames(dat), c("sample","chr","win","homFrac","hetFrac"))
dat[,spos:=as.numeric(unlist(lapply(strsplit(as.character(dat[,win]), split=":"), function(x) x[1])))]
dat[,epos:=as.numeric(unlist(lapply(strsplit(as.character(dat[,win]), split=":"), function(x) x[2])))]
dat<-dat[!sample %in% c("BP206","BPK282","LV9")]
mdat<-melt(dat, id.vars = c("sample","chr","win","spos","epos"))
mdat[,variable:=as.character(variable)]

# plotting of hom / het regions by chrom and sample
dir.create("Hom_het_regions_bysample", showWarnings = F)
# for (chrom in c(paste0("0",1:9),10:36))
# {
#   pdf(paste0("Hom_het_regions_bysample/Hom_het_regions_bysample_chr",chrom,".pdf"), width = 8, height=5)
#   gg<-ggplot(data=mdat[!sample %in% c("BP206","BPK282","LV9i") & chr==paste0("LinJ.",chrom)], aes(x=spos, y=value, col=sample)) +
#     geom_line() +
#     geom_point() +
#     facet_wrap(.~variable, ncol=1)
#   print(gg)
#   dev.off()
# }

# get the mean fractions across all CUK samples
dat[, homFrac_mean:=mean(homFrac), by=.(chr,spos)]
dat[, hetFrac_mean:=mean(hetFrac), by=.(chr,spos)]  
# identifies regions that are hom in all samples
dat[,p_JPCM5:=0]
dat[,p_other:=0]
fJPCM5=0.0004
dat[homFrac_mean<fJPCM5 & hetFrac_mean < 0.0002, p_JPCM5:=1] # indicates when heterozygosity is low and regions/ windows very CLOSE to the reference
dat[homFrac_mean>0.001 & hetFrac_mean < 0.0002, p_other:=1] # indicates when heterozygosity is low and regions/ windows very DISTANT from the reference
# 0.0002*5000
# [1] 1 --> 1 SNP per 5000 kb window
# 
# 0.001*5000
# [1] 5 --> 5 SNPs per 5000 kb window
dat<-dat[order(chr,spos)]
dat[,chr:=as.character(chr)]
mdat_mean<-melt(dat[,.(sample,chr,win,spos,epos,homFrac_mean,hetFrac_mean,p_JPCM5,p_other)], id.vars = c("sample","chr","win","spos","epos","p_JPCM5","p_other"))
mdat_mean<-mdat_mean[,.(chr,win,spos,epos,variable,value,p_JPCM5,p_other)]

# plotting of hom / het regions by chrom and mean across samples
# dir.create("Hom_het_regions_mean", showWarnings = F)
# for (chrom in c(paste0("0",1:9),10:36))
# {
#   pdf(paste0("Hom_het_regions_mean/Hom_het_regions_mean_chr",chrom,".pdf"), width = 8, height=5)
#   gg<-ggplot(data=mdat_mean[chr==paste0("LinJ.",chrom)], aes(x=spos, y=value, col=variable)) +
#     geom_line() +
#     geom_point()
#   print(gg)
#   dev.off()
# }

# decide on the cutoff when a hom region might be from the JPCM5-LIKE parent or the other clearly more distant parent
pdf("Hist_mean_fracs_hom_het.pdf")
ggplot(mdat_mean, aes(x=value, color=variable)) + 
  # geom_histogram(fill="white", position="dodge") 
  geom_histogram( alpha=0.5, position="identity", bins=100) +
  geom_vline(xintercept=fJPCM5) +
  geom_vline(xintercept=0.004) +
  geom_vline(xintercept=0.001)
dev.off()


dat_red<-unique(dat[,.(chr,spos,epos,homFrac_mean,hetFrac_mean,p_JPCM5,p_other)])



######################
#
# identification of larger region that are either as the JPCM5-LIKE parent or the other clearly more distant parent
#
# regions resembling (identical) to JPCM5
p1=0 # indicating parent 1 region or not of previous row
p1chr="LinJ" # indicating chr of previous row
p1epos=0 # indicating spos of previous row
p1_res=data.table(matrix(c(rep(NA,4)), ncol = 4))
for (i in 1:nrow(dat_red))
{
  # start of p1 region
  if(dat_red[i,p_JPCM5] > p1)
  {
    # print(i)
    p1_reg<-c(dat_red[i,chr], dat_red[i,spos]) 
  }
  
  # end of p1 region within a chr
  if(dat_red[i,p_JPCM5] < p1 & p1chr== dat_red[i,chr])
  {
    p1_res=rbind(p1_res, data.table(matrix(c(p1_reg, p1chr, p1epos), ncol = 4)))
  } 
  
  # end of chromosome ending with p1 region
  if(p1chr!= dat_red[i,chr] & p1==1)
  {
    p1_res=rbind(p1_res, data.table(matrix(c(p1_reg, p1chr, p1epos), ncol = 4)))
    if(dat_red[i,p_JPCM5]==1)
    {
      p1_reg<-c(dat_red[i,chr], dat_red[i,spos]) 
    }
  } 
  
  p1=dat_red[i,p_JPCM5]
  p1chr=dat_red[i,chr]
  p1epos=dat_red[i,epos]
}
#
setnames(p1_res,colnames(p1_res),c("chrA","spos","chrB","epos"))
p1_res<-p1_res[!is.na(chrA)]
p1_res[,spos:=as.numeric(spos)]
p1_res[,epos:=as.numeric(epos)]
p1_res[,len:=epos-spos]
p1_res[,parent:="p1"]


######################
#
# identification of larger region that are either as the other clearly more distant parent
#
# regions describing the other parent
p2=0 # indicating parent 1 region or not of previous row
p2chr="LinJ" # indicating chr of previous row
p2epos=0 # indicating spos of previous row
p2_res=data.table(matrix(c(rep(NA,4)), ncol = 4))
for (i in 1:nrow(dat_red))
{
  # start of p2 region
  if(dat_red[i,p_other] > p2)
  {
    # print(i)
    p2_reg<-c(dat_red[i,chr], dat_red[i,spos]) 
  }
  
  # end of p2 region within a chr
  if(dat_red[i,p_other] < p2 & p2chr== dat_red[i,chr])
  {
    p2_res=rbind(p2_res, data.table(matrix(c(p2_reg, p2chr, p2epos), ncol = 4)))
  } 
  
  # end of chromosome ending with p2 region
  if(p2chr!= dat_red[i,chr] & p2==1)
  {
    p2_res=rbind(p2_res, data.table(matrix(c(p2_reg, p2chr, p2epos), ncol = 4)))
    if(dat_red[i,p_other]==1)
    {
      p2_reg<-c(dat_red[i,chr], dat_red[i,spos]) 
    }
  } 
  
  p2=dat_red[i,p_other]
  p2chr=dat_red[i,chr]
  p2epos=dat_red[i,epos]
}
#
setnames(p2_res,colnames(p2_res),c("chrA","spos","chrB","epos"))
p2_res<-p2_res[!is.na(chrA)]
p2_res[,spos:=as.numeric(spos)]
p2_res[,epos:=as.numeric(epos)]
p2_res[,len:=epos-spos]
p2_res[,parent:="p2"]
#
#
pRes<-rbind(p1_res,p2_res)
pRes[,chr:=as.numeric(unlist(lapply(strsplit(pRes$chrA, split="J."), function(x) x[2])))] # add chr column as number


pdf("Hist_parent_blocks.pdf")
ggplot(pRes, aes(x=len, color=parent)) + 
  geom_histogram(fill="white", position="dodge", bins=100) 
# geom_histogram( alpha=0.5, position="identity", bins=100) 
dev.off()

pdf("parent_blocks_JPCM5.pdf", width = 10, height=5)
ggplot(dat_red, aes(x=spos, y=as.factor(p_JPCM5), col=as.factor(p_JPCM5))) +
  geom_point(size=0.1) +
  facet_grid(chr~.)
dev.off()
pdf("parent_blocks_other.pdf", width = 10, height=5)
ggplot(dat_red, aes(x=spos, y=as.factor(p_other), col=as.factor(p_other))) +
  geom_point(size=0.1) +
  facet_grid(chr~.)
dev.off()


#--------------------------
# make phylogenetic tress from obtained parent-like regions

#
# data without outgroup
# load(file = "../A01_StAMPP/l.gtStampp_cov0.2_poly.RData") # l.gtStampp.freq, list of simple tables for each chromosome with SNP data in StaMPP format
#
# data with outgroup
# load(file = "../A01_StAMPP/addoutgroup/l.gtStampp_cov0.2_poly1_20.RData") # l.gtStampp.freq, list of simple tables for each chromosome with SNP data in StaMPP format
# l1_20<-l.gtStampp.freq
# load(file = "../A01_StAMPP/addoutgroup/l.gtStampp_cov0.2_poly21_31.RData") # l.gtStampp.freq, list of simple tables for each chromosome with SNP data in StaMPP format
# l21_31<-l.gtStampp.freq
# load(file = "../A01_StAMPP/addoutgroup/l.gtStampp_cov0.2_poly32_35.RData") # l.gtStampp.freq, list of simple tables for each chromosome with SNP data in StaMPP format
# l32_35<-l.gtStampp.freq
# load(file = "../A01_StAMPP/addoutgroup/l.gtStampp_cov0.2_poly36.RData") # l.gtStampp.freq, list of simple tables for each chromosome with SNP data in StaMPP format
# l36<-l.gtStampp.freq
# l.gtStampp.freq<-c(l1_20,l21_31,l32_35,l36)
# save(l.gtStampp.freq, file="../A01_StAMPP/addoutgroup/l.gtStampp_cov0.2_poly.RData")
load(file="../../01_phylogenetic_reconstruction/addoutgroup/l.gtStampp_cov0.2_poly.RData") # l.gtStampp.freq
allSamp<-l.gtStampp.freq[[1]]$Sample

pRes[,maxlenChr:=max(len), by=.(chr, parent)] # to identify longest block for each parent by chromosome
pResLong<-pRes[len==maxlenChr,]
pResLong<-pResLong[order(-len,chr,parent)] # longest blocks per parent and chr

######################
#
# plot hom het regions for the largest one parent blocks found
for (chrom in c(paste0("0",c(5,7)),c(13,21,24,32,34,35)))
{
  pdf(paste0("Hom_het_regions_bysample/Hom_het_regions_bysample_chr",chrom,"_blocks.pdf"), width = 8, height=5)
  a<-pRes[chr==as.numeric(chrom)]
  a<-cbind(a, y1=rep(0,nrow(a)), y2=rep(  max(mdat[!sample %in% c("BP206","BPK282","LV9i") & chr==paste0("LinJ.",chrom), value]) ,nrow(a)))
  
  gg<-ggplot(data=mdat[!sample %in% c("BP206","BPK282","LV9i") & chr==paste0("LinJ.",chrom)], aes(x=spos, y=value, col=sample)) +
    annotate("rect", xmin=a[parent=="p1"]$spos, xmax=a[parent=="p1"]$epos, ymin=a[parent=="p1"]$y1 , ymax=a[parent=="p1"]$y2, 
             alpha=0.2,  fill="blue") +
    annotate("rect", xmin=a[parent=="p2"]$spos, xmax=a[parent=="p2"]$epos, ymin=a[parent=="p2"]$y1 , ymax=a[parent=="p2"]$y2, 
             alpha=0.2,  fill="orange") +
    annotate("rect", xmin=a[parent=="p1" & len==maxlenChr]$spos, xmax=a[parent=="p1" & len==maxlenChr]$epos, ymin=a[parent=="p1" & len==maxlenChr]$y1 , ymax=a[parent=="p1" & len==maxlenChr]$y2, 
             alpha=0.2,  color="blue") +
    annotate("rect", xmin=a[parent=="p2" & len==maxlenChr]$spos, xmax=a[parent=="p2" & len==maxlenChr]$epos, ymin=a[parent=="p2" & len==maxlenChr]$y1 , ymax=a[parent=="p2" & len==maxlenChr]$y2, 
             alpha=0.2,  color="orange") +
    geom_line() +
    geom_point() +
    facet_wrap(.~variable, ncol=1) +
    theme_bw()
  print(gg)
  dev.off()
}

######################
#
# figs with 4 largest blocks for each parent
gg_Hom_het_regions_bysample_chr_block <- function(chrom, parent)
{
  a<-pRes[chr==as.numeric(chrom)]
  a<-cbind(a, y1=rep(0,nrow(a)), y2=rep(  max(mdat[!sample %in% c("BP206","BPK282","LV9i") & chr==paste0("LinJ.",chrom), value]) ,nrow(a)))
  
  gg<-ggplot(data=mdat[!sample %in% c("BP206","BPK282","LV9i") & chr==paste0("LinJ.",chrom)], aes(x=spos, y=value, col=sample)) +
    annotate("rect", xmin=a[parent=="p1"]$spos, xmax=a[parent=="p1"]$epos, ymin=a[parent=="p1"]$y1 , ymax=a[parent=="p1"]$y2, 
             alpha=0.2,  fill="blue") +
    annotate("rect", xmin=a[parent=="p2"]$spos, xmax=a[parent=="p2"]$epos, ymin=a[parent=="p2"]$y1 , ymax=a[parent=="p2"]$y2, 
             alpha=0.2,  fill="orange") +
    geom_line() +
    geom_point() +
    facet_grid(variable~chr) +
    theme_bw() +
    guides(col=FALSE) + labs(x="Chromosomal position", y="Variant fraction")
  
  if (parent=="p1"){gg<-gg+annotate("rect", xmin=a[parent=="p1" & len==maxlenChr]$spos, xmax=a[parent=="p1" & len==maxlenChr]$epos, ymin=a[parent=="p1" & len==maxlenChr]$y1 , ymax=a[parent=="p1" & len==maxlenChr]$y2, 
                                    alpha=0.2,  color="blue")}
  if (parent=="p2"){gg<-gg+annotate("rect", xmin=a[parent=="p2" & len==maxlenChr]$spos, xmax=a[parent=="p2" & len==maxlenChr]$epos, ymin=a[parent=="p2" & len==maxlenChr]$y1 , ymax=a[parent=="p2" & len==maxlenChr]$y2, 
                                    alpha=0.2,  color="orange")}
  gg
  
}
#
g1<-gg_Hom_het_regions_bysample_chr_block(chrom="05",parent="p1")
g2<-gg_Hom_het_regions_bysample_chr_block(chrom="07",parent="p1")
g3<-gg_Hom_het_regions_bysample_chr_block(chrom="13",parent="p1")
g4<-gg_Hom_het_regions_bysample_chr_block(chrom="35",parent="p1")
pdf(paste0("Hom_het_regions_bysample/Hom_het_regions_bysample_p1_blocks.pdf"), width = 10, height=8)
g1+g2+g3+g4
dev.off()
#
g1<-gg_Hom_het_regions_bysample_chr_block(chrom="21",parent="p2")
g2<-gg_Hom_het_regions_bysample_chr_block(chrom="24",parent="p2")
g3<-gg_Hom_het_regions_bysample_chr_block(chrom="32",parent="p2")
g4<-gg_Hom_het_regions_bysample_chr_block(chrom="34",parent="p2")
pdf(paste0("Hom_het_regions_bysample/Hom_het_regions_bysample_p2_blocks.pdf"), width = 10, height=8)
g1+g2+g3+g4
dev.off()


# chrom=24;spos=670000;epos=835000

dir.create("block_trees", showWarnings = F)

get_block_phylo <-function(chrom,spos,epos)
{
  library(StAMPP)
  
  # SNP position without 5 trailing header "columns"
  pos<-as.numeric(unlist(lapply(strsplit(names(l.gtStampp.freq[[chrom]]), split="_"), function(x) x[2])))
  pos<-pos[6:length(pos)]
  posBoo<-c(rep(TRUE,5),pos>=spos & pos<=epos) # boolean vector indicating which positions to take
  
  freq_red<-l.gtStampp.freq[[chrom]][,posBoo, with=F]
  freq_red<-freq_red[!freq_red$Sample %in% c("CL1_ethiopica","CL4_ethiopica"),]
  # freq_red<-freq_red[!freq_red$Sample %in% c("CH32","CH33","CH34","CH35","CH36"),]
  # summary( as.numeric(unlist(lapply(strsplit(names(freq_red), split="_"), function(x) x[2])))[6:length(pos)] ) # just to double check correct positions have been chosen
  
  NeisD.ind<-stamppNeisD(data.frame(freq_red), pop=F)
  tree <- nj(NeisD.ind)
  
  load("../../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData") # sample.info
  sample.info <-sample.info[order(sample)]
  if (nrow(freq_red)>151) # add sample colors for outgroup samples
  {
    outg<-allSamp[!allSamp %in% sample.info$sample]
    sample.info<-rbind(sample.info, list(outg[1], "NA", "black"))
    sample.info<-rbind(sample.info, list(outg[2], "NA", "black"))
    sample.info<-rbind(sample.info, list(outg[3], "NA", "black"))
    # sample.info<-rbind(sample.info, list(outg[4], "NA", "black"))
    # sample.info<-rbind(sample.info, list(outg[5], "NA", "black"))
  }
  sample.info<-sample.info[order(match(sample,tree$tip.label))]
  tip.label.col<-sample.info$group.col

  sspos<-spos
  pdf(paste0("block_trees/",pRes[chr==chrom & spos==sspos, parent],"_chr",chrom,"_spos",spos,"_len",pRes[chr==chrom & spos==sspos, len]/1000,"k.pdf"),
      width=25,height=25)
  plot(root(tree, which(tree$tip.label == "LmexU1103_v1")), tip.color=tip.label.col)
  # plot(tree, tip.color=tip.label.col)
  dev.off()
  write.tree(tree, file = paste0("block_trees/",pRes[chr==chrom & spos==sspos, parent],"_chr",chrom,"_spos",spos,"_len",pRes[chr==chrom & spos==sspos, len]/1000,"k.txt"), append = FALSE, digits = 10, tree.names = FALSE)
}

# table containing the information of the longest blocks specific to either parent
pResLong[1:10,]
# chrA    spos    chrB    epos    len parent chr maxlenChr
# 1: LinJ.35       0 LinJ.35  215000 215000     p1  35    215000
# 2: LinJ.32  525000 LinJ.32  710000 185000     p2  32    185000
# 3: LinJ.05       0 LinJ.05  180000 180000     p1   5    180000
# 4: LinJ.24  670000 LinJ.24  835000 165000     p2  24    165000
# 5: LinJ.07  435000 LinJ.07  590000 155000     p1   7    155000
# 6: LinJ.13  275000 LinJ.13  425000 150000     p1  13    150000
# 7: LinJ.21  605000 LinJ.21  755000 150000     p2  21    150000
# 8: LinJ.34 1330000 LinJ.34 1480000 150000     p2  34    150000
# 9: LinJ.35  335000 LinJ.35  485000 150000     p2  35    150000
# 10: LinJ.30  985000 LinJ.30 1125000 140000     p1  30    140000

get_block_phylo(chrom=35,spos=0,epos=215000)
get_block_phylo(chrom=32,spos=525000,epos=710000)
get_block_phylo(chrom=5,spos=0,epos=180000)
get_block_phylo(chrom=24,spos=670000,epos=835000)
get_block_phylo(chrom=7,spos=435000,epos=590000)
get_block_phylo(chrom=13,spos=275000,epos=425000)
get_block_phylo(chrom=21,spos=605000,epos=755000)
get_block_phylo(chrom=34,spos=1330000,epos=1480000)




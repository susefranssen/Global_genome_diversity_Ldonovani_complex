
library(foreach)
library(data.table)
library(ggplot2)
library(GGally)
library(gridExtra)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/04_aneuploidy")

# # data tables tht will be created 
# ploidy        # sample ploidies in different columns
# ploidy.a      # sample ploidies in one column along with meta data
# sample.info   # sample info including groups and colours
# sd.g          # sd in somy by chromosome and group, in different columns for different groups
# sd.g.format   # sd in somy by chromosome and group along with metadata, i.e. group


# generate and safe file with sample metadata information
#
# load(file=paste0(path,"/leish_donovaniComplex/A05_heterozygosities/sample.info.RData"))
# sample.info[groups=="infantum1", groups:="Linf1"]
# sample.info[groups=="donovani1", groups:="Ldon1"]
# sample.info[groups=="donovani3", groups:="Ldon5"]
# sample.info[groups=="donovani2", groups:="Ldon4"]
# sample.info[groups=="donovani4", groups:="Ldon3"]
# sample.info[grep("CUK", sample), groups:="CUK_Linf"]
# sample.info[grep("CUK", sample), group.col:="green"]
# sample.info[grep("CH", sample), groups:="CH_Linf"]
# sample.info[grep("CH", sample), group.col:="gray34"]
# sample.info[sample %in% c("EP","MAM"), groups:="other_Linf"]
# sample.info[sample %in% c("EP","MAM"), group.col:="black"]
# sample.info[sample %in% c("GE","LEM3472","LRC.L740"), groups:="other_Ldon"]
# sample.info[sample %in% c("GE","LEM3472","LRC.L740"), group.col:="#ff0080"]
# sample.info[sample %in% c("BPK156A1","BPK406A1","BPK413A1","BPK512A1","BPK612A1","BPK623A1","BPK648A1"), groups:="Ldon2"]
# sample.info[sample %in% c("BPK156A1","BPK406A1","BPK413A1","BPK512A1","BPK612A1","BPK623A1","BPK648A1"), group.col:="#7300da"]
# 
# table(sample.info$groups)
# # CH        CUK      Ldon1      Ldon2      Ldon3      Ldon4      Ldon5      Linf1 other_Ldon other_Linf
# # 5         11         45          7         19          4          8         47          3          2
# 
# save(sample.info, file="../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")
load(file="../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData") # sample.info
ploidy<-data.table(read.table("somies_updated.txt", header = T))

######################
# functions


# reorder correlation matrix
reorder.cor.mat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  cormat
}

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[lower.tri(cormat)] <- NA
  return(cormat)
}

get_upper_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}

# The F statistic on the last line is telling you whether the regression as a whole
# is performing 'better than random' - any set of random predictors will have some
# relationship with the response, so it's seeing whether your model fits better than you'd
# expect if all your predictors had no relationship with the response (beyond what would
# be explained by that randomness). This is used for a test of whether the model outperforms
# 'noise' as a predictor. The p-value in the last row is the p-value for that test, 
# essentially comparing the full model you fitted with an intercept-only model.
#
# method for getting the overall pval of the regression model
lmp <- function (modelobject) {
  if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
  f <- summary(modelobject)$fstatistic
  p <- pf(f[1],f[2],f[3],lower.tail=F)
  attributes(p) <- NULL
  return(p)
}

#######################


# 
# 
# # samples we want to phase for chromosomes with uneven somies
# samp<-c("MAM","EP",sample.info[groups=="Ldon4",sample], colnames(ploidy)[grep("CUK",colnames(ploidy))])
# samp<-c("MAM","EP",sample.info[groups=="Ldon4",sample], colnames(ploidy)[grep("CUK",colnames(ploidy))], "LRC.L53", sample.info[groups=="Ldon3",sample])
# 
# 
# # get chr with uneven somies
# ploidy.sub<-ploidy[,colnames(ploidy) %in% samp, with=F]
# 
# chroms<-c() # chrom with uneven somies in our samples
# for (chr in 1:36){ chroms<- c(chroms, ifelse(sum((ploidy.sub[chr,]%% 2 ==1))>0,chr,0)) }
# chroms<-chroms [! chroms %in% 0] # remove zeros
# # number of samples with uneven somy at respective chromosome
# n.chrodd<-foreach (chr=chroms, .combine=c) %do% { sum((ploidy.sub[chr,]%% 2 ==1)) }
# 
# # with most unven somies across the sample set of interest
# cc<-chroms[n.chrodd==max(n.chrodd)]
# # > chroms[n.chrodd==12]
# # [1] 26
# # > chroms[n.chrodd==11]
# # [1] 9
# # > chroms[n.chrodd==10]
# # [1] 23
# 
# # > chroms[n.chrodd==19]
# # [1] 23
# # > chroms[n.chrodd==18]
# # [1] 26
# # > chroms[n.chrodd==17]
# # [1] 9
# 
# 


#############################
# look at ploidy across populations

ploidy.a<-cbind(ploidy)
ploidy.a<-melt(ploidy.a, measure.vars = colnames(ploidy))
ploidy.a[,chr:=rep(1:36, ncol(ploidy))]
setnames(ploidy.a, colnames(ploidy.a)[1:2], c("sample","somy"))
for (i in unique(sample.info$groups)) {
  ploidy.a[sample %in% sample.info[groups==i, sample],group:=i]
}
ploidy.a[,size:=length(sample)/36, by=(group)] # size of the group
ploidy.a[,c:=1]
ploidy.a[,c_somy.chr.g:=sum(c), by=.(group, somy, chr)] # counts for each somy per chr and per group
ploidy.a[,sd_somy.chr.g:=sd(somy), by=.(group, chr)]
ploidy.a[,sd_somy.chr:=sd(somy), by=.(chr)]
ploidy.a[,med_somy.chr:=median(somy), by=.(chr)]
ploidy.a[,mean_somy.chr:=mean(somy), by=.(chr)]
ploidy.a[,mean_somy.chr.g:=mean(somy), by=.(group,chr)]
ploidy.a[size>=9,med_somy.chr.g:=median(somy), by=.(group,chr)]

color.palette  <- c(colorRampPalette(c("yellow","orange","green4","blue","red"))(6) ,c("purple","pink"))

pdf("Somy_range_median_all.pdf", width=9, height=4)
ggplot(data=ploidy.a, aes(x=as.factor(chr), y=somy, color=as.factor(med_somy.chr))) + 
  geom_boxplot(fatten = 5) + 
  labs(x = "Chromosome", y = "Somy", color = "Median\nsomy") +
  scale_color_manual(values=color.palette[sort(unique(ploidy.a$med_somy.chr))])
dev.off()
#
pdf("Somy_range_median_groups.pdf", width=9, height=4)
ggplot(data=ploidy.a[size>=9], aes(x=as.factor(chr), y=somy, color=as.factor(med_somy.chr.g))) + 
  geom_boxplot(fatten = 5) + 
  labs(x = "Chromosome", y = "Somy", color = "Median\nsomy") +
  facet_wrap(~group, ncol = 2) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values=color.palette[sort(unique(ploidy.a$med_somy.chr.g))])
dev.off()
#
table(ploidy.a[,.(chr,somy)])
ploidy.summary<-table(ploidy.a[,.(chr,somy)])
ploidy.summary.t<-data.table(ploidy.summary)
# for all chromosomes at least 2 samples were trisomic
summary(ploidy.summary.t[somy==3, .(N)])
# N       
# Min.   : 2.0  
# 1st Qu.: 8.5  
# Median :28.0  
# Mean   :29.0  
# 3rd Qu.:46.0  
# Max.   :76.0
ploidy.summary.t[chr==31]
# chr somy   N
# 1:  31    1   0
# 2:  31    2   0
# 3:  31    3   9
# 4:  31    4 121
# 5:  31    5  19
# 6:  31    6   1
# 7:  31    8   1

ploidy.a[,.(chr, group,mean_somy.chr,mean_somy.chr.g)][chr==1 & group=="Linf1"]

# somy variability by chromosome
#
xx<-unique(ploidy.a[,.(chr,group,size,med_somy.chr,mean_somy.chr.g,sd_somy.chr.g)])[size>=9]
xx[,med_sd_somy.chr.g:=median(sd_somy.chr.g), by=.(chr)]
xx[,med_mean_somy.chr.g:=median(mean_somy.chr.g), by=.(chr)]
#
xx[,chr:=as.factor(chr)]
xx$chr <- factor(xx$chr, levels = unique(xx[,.(chr,med_sd_somy.chr.g)])[order(med_sd_somy.chr.g)]$chr) # sorting chr levels accoring to the median of the group sd in somy
#
pdf("Sdsomy_bychr_groups.pdf", width=8, height=3.5)
ggplot(data=xx[size>=9], aes(as.factor(chr), y=sd_somy.chr.g, fill=as.factor(med_somy.chr))) +
  geom_boxplot() +
  labs(x="Chromosomes [ordered by median sd of group somies]", 
       y="Sd of somies of the four largest groups") + 
  guides(fill=guide_legend(title="Median\nsomy across\n151 samples")) +
  theme(legend.position="none") +
  scale_fill_manual(values=color.palette[unique(ploidy.a$med_somy.chr)])
dev.off()
xx.red<-xx[size>=9][, sd_somy.chr.g.med:=median(sd_somy.chr.g), by=chr]
xx.red<-unique(xx.red[,.(chr, sd_somy.chr.g.med)])
write.table(xx.red, file="Sdsomy_bychr_groups_MedianChrVariablityByChr_groupSizeL8.txt", quote=F, row.names = F)

xx$chr <- factor(xx$chr, levels = unique(xx[,.(chr,med_mean_somy.chr.g)])[order(med_mean_somy.chr.g)]$chr)
#
# pdf("Meansomy_bychr_groups.pdf", width=8, height=3.5)
# ggplot(data=xx[size>=9], aes(as.factor(chr), y=mean_somy.chr.g, fill=as.factor(med_somy.chr))) +
#   geom_boxplot() +
#   labs(x="Chromosomes [ordered by median of mean group somies]", 
#        y="Mean somies of the four largest groups") + 
#   guides(fill=guide_legend(title="Median\nsomy across\n151 samples")) +
#   theme(legend.position="none") +
#   scale_fill_manual(values=color.palette[unique(ploidy.a$med_somy.chr)])
# dev.off()


# correlation between sd values of different groups
#
xx$chr <- factor(xx$chr, levels = 1:36) # sorting chr levels by chr number again

# correlation between sd in population somy
sd.g.format<-unique(xx[size>=9][,.(chr, group, sd_somy.chr.g, med_sd_somy.chr.g)])
sd.g<-dcast(sd.g.format, chr ~ group, value.var = "sd_somy.chr.g")
plot(sd.g[,2:5, with=F])
cor(sd.g[,2:5, with=F], method = "spearman")
#          CUK_Linf     Ldon1     Ldon3     Linf1
# CUK_Linf 1.0000000 0.4485932 0.4276383 0.4778165
# Ldon1    0.4485932 1.0000000 0.8065094 0.7775469
# Ldon3    0.4276383 0.8065094 1.0000000 0.8184087
# Linf1    0.4778165 0.7775469 0.8184087 1.0000000
pvals<-cor.test(unlist(sd.g[,2, with=F]), unlist(sd.g[,3, with=F]))$p.value
pvals<-c(pvals, cor.test(unlist(sd.g[,2, with=F]), unlist(sd.g[,4, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(sd.g[,2, with=F]), unlist(sd.g[,5, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(sd.g[,3, with=F]), unlist(sd.g[,4, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(sd.g[,3, with=F]), unlist(sd.g[,5, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(sd.g[,4, with=F]), unlist(sd.g[,5, with=F]))$p.value)
p.adjust(pvals)
# [1] 1.203509e-02 1.203509e-02 5.252020e-03 8.415838e-09 2.139609e-06 2.139609e-06
#
pdf("Sdsomy_bychr_groups_cor.pdf", width=4, height=4)
ggcorr(sd.g[,2:5, with=F], palette = "RdBu", label = TRUE, method=c("pairwise", "spearman"),
       digits = 4, name="Spearman\ncorrelation\nfor sd in somy")
dev.off()

# correlation between mean population somy
mean.g.format<-unique(xx[size>=9][,.(chr, group, mean_somy.chr.g)])
mean.g<-dcast(mean.g.format, chr ~ group, value.var = "mean_somy.chr.g")
plot(mean.g[,2:5, with=F])
cor(mean.g[,2:5, with=F], method = "spearman")
#          CUK_Linf     Ldon1     Ldon3     Linf1
# CUK_Linf 1.0000000 0.5472025 0.4972307 0.4909716
# Ldon1    0.5472025 1.0000000 0.8277376 0.8047534
# Ldon3    0.4972307 0.8277376 1.0000000 0.7928307
# Linf1    0.4909716 0.8047534 0.7928307 1.0000000
pvals<-cor.test(unlist(mean.g[,2, with=F]), unlist(mean.g[,3, with=F]))$p.value
pvals<-c(pvals, cor.test(unlist(mean.g[,2, with=F]), unlist(mean.g[,4, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(mean.g[,2, with=F]), unlist(mean.g[,5, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(mean.g[,3, with=F]), unlist(mean.g[,4, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(mean.g[,3, with=F]), unlist(mean.g[,5, with=F]))$p.value)
pvals<-c(pvals, cor.test(unlist(mean.g[,4, with=F]), unlist(mean.g[,5, with=F]))$p.value)
p.adjust(pvals)
# [1] 1.242714e-10 2.399375e-10 1.242714e-10 3.619844e-14 4.621638e-14 1.899877e-13
#
pdf("Meansomy_bychr_groups_cor.pdf", width=4, height=4)
ggcorr(mean.g[,2:5, with=F], digits=4,  palette = "RdBu", label = TRUE, method=c("pairwise", "spearman"),
       name="Spearman\ncorrelation\nfor mean group somy")
dev.off()

# correlation between mean and sd of group ploidy
pvals<-c(cor.test(mean.g$Linf1,sd.g$Linf1, method="spearman")$p.val)
pvals<-c(pvals, cor.test(mean.g$Ldon3,sd.g$Ldon3, method="spearman")$p.val)
pvals<-c(pvals, cor.test(mean.g$Ldon1,sd.g$Ldon1, method="spearman")$p.val)
pvals<-c(pvals, cor.test(mean.g$CUK_Linf,sd.g$CUK_Linf, method="spearman")$p.val)
p.adjust(pvals)
# [1] 1.810863e-06 1.923409e-05 4.827512e-04 4.827512e-04
rsq<-c(cor.test(mean.g$Linf1,sd.g$Linf1, method="spearman")$estimate)
rsq<-c(rsq, cor.test(mean.g$Ldon3,sd.g$Ldon3, method="spearman")$estimate)
rsq<-c(rsq, cor.test(mean.g$Ldon1,sd.g$Ldon1, method="spearman")$estimate)
rsq<-c(rsq, cor.test(mean.g$CUK_Linf,sd.g$CUK_Linf, method="spearman")$estimate)
rsq
# rho       rho       rho       rho 
# 0.9566017 0.9594311 0.9160066 0.9634732 


write.table(ploidy.a, file="group.ploidy.sd.txt", quote = F, row.names = F)



# 
# somy counts per group for each different somy observed
ploidy.a.g<-unique(ploidy.a[,.(somy,chr,group,size,c_somy.chr.g)])
ploidy.a.g[,c_somy.chr.g.f:=round(c_somy.chr.g/size,4)]
#
# pdf("Chr.pop.frac.polysomy.pdf", width=9, height=6)
# ggplot(data=ploidy.a.g, aes(group, y=c_somy.chr.g.f, fill=as.factor(somy))) +
#   geom_bar(stat="identity") +
#   facet_wrap(~ chr) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 9, hjust = 1)) +
#   ylab("Group somies [fraction]") + 
#   scale_fill_manual(values=color.palette[sort(unique(ploidy.a.g$somy))])
# dev.off()
#
pdf("Chr.pop.frac.polysomy_min5.pdf", width=11, height=7)
ggplot(data=ploidy.a.g[size>=5], aes(group, y=c_somy.chr.g.f, fill=as.factor(somy))) +
  geom_bar(stat="identity") +
  facet_wrap(~ chr) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 13, hjust = 1),
        axis.text = element_text(size=13),
        axis.title = element_text(size=13),
        strip.text.x = element_text(size = 15)) +
  ylab("Group somies [fraction]") + 
  scale_fill_manual(values=color.palette[sort(unique(ploidy.a.g$somy))])
dev.off()



#
#
# test if on the chromosome level low heterozygosity is associated with high chromosome copy umber variability
#
# data used:
sd.g.format
load("../05_heterozygosities/chr.hets.RData") # chr.hets
setnames(chr.hets,colnames(chr.hets),c("samp","chr","het"))
chr.hets[,chr:=gsub("chr", "", chr)]

for (ss in sample.info$sample)
{
  chr.hets[samp==ss, group:=sample.info[sample==ss, groups]]
}

for (chrom in 1:36) # sd on chr variability per group
{
  print(chrom)
  for (gg in unique(chr.hets$group)) 
  {
    chr.hets[chr==chrom & group==gg, sd_somy.chr.g:= unique(ploidy.a[chr==chrom & group==gg, sd_somy.chr.g])]
  }
}

for (chrom in 1:36) # median across 4 major groups of # sd on chr variability per group
{
  chr.hets[chr==chrom, med_sd_somy.chr.g:=unique(sd.g.format[chr==chrom, med_sd_somy.chr.g])] # median across 4 major groups of # sd on chr variability per group
}

aa<-chr.hets[chr!=31][group!="other_Linf"]
bb<-aa[group %in% c("Ldon1","Ldon2","Ldon3","Ldon5","Linf1","CH_Linf","CUK_Linf")]


# with chr variability estimate from variability estimate per group
ggplot(data=bb, aes(x=sd_somy.chr.g, y=het)) +
  geom_point() + facet_wrap(~group) +
  geom_smooth(method = "lm", se = FALSE)

ggplot(data=bb, aes(x=sd_somy.chr.g, y=het)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

cor.test(bb$sd_somy.chr.g, bb$het, method="spearman")
# S = 2.1637e+10, p-value = 5.018e-05
#   rho -0.05748314 

pdf("Cor_Sdsomy_bygroups_het.pdf")
hist_top <- ggplot(bb, aes(x=sd_somy.chr.g, fill=group))+geom_histogram() + 
  theme(legend.position='none') + xlab("") +
  scale_fill_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff75ff","#0000ff")) 
  # scale_fill_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff7f00","#ff75ff","#0000ff","#ff0080","black"))
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(), 
        panel.background=element_blank(), 
        axis.text.x=element_blank(), axis.text.y=element_blank(),           
        axis.title.x=element_blank(), axis.title.y=element_blank())
scatter <- ggplot(bb, aes(x=sd_somy.chr.g, y=het, color=group))+geom_point(size=0.7) + theme(legend.position='none') + 
  geom_smooth(method = "lm", se = F, aes(color=group, fill=group), fullrange=T, size=0.6) +
  scale_color_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff75ff","#0000ff")) +
  scale_fill_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff75ff","#0000ff")) 
  # scale_color_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff7f00","#ff75ff","#0000ff","#ff0080","black")) 
hist_right <- ggplot(bb, aes(x=het, fill=group))+geom_histogram()+coord_flip() + theme(legend.position='none') + ylab("") +
  scale_fill_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff75ff","#0000ff")) 
  # scale_fill_manual(values=c("#737373","#00ff00","#ff0000","#7300da","#8b1c62","#ff7f00","#ff75ff","#0000ff","#ff0080","black")) 
  
grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
dev.off()







res.cor<-foreach (gg = unique(bb$group), .combine=rbind) %do%
{
  # gg="Linf1"
  tmp<-cor.test(bb[group==gg]$het , bb[group==gg]$sd_somy.chr.g, method="spearman")
  xx<-lm(data=bb[group==gg], formula= het ~ sd_somy.chr.g)
  c(gg,round(tmp$estimate,2), round(tmp$p.value,6), 
    ifelse(tmp$p.value<0.001, "***", ifelse(tmp$p.value<0.01, "**", ifelse(tmp$p.value<0.05, "*", "-"))),
    round(xx$coefficients[1],5),round(xx$coefficients[2],5),
    round(summary(xx)$coefficients[1,4],5),round(summary(xx)$coefficients[2,4],5),
    lmp(xx))
}
res.cor<-data.table(res.cor)
res.cor[order(V3, rho)]
setnames(res.cor,c("V1","V3","V4","(Intercept)","sd_somy.chr.g","V7","V8","V9"),
         c("group","cor.pval","cor.sig","lm.inter","lm.slope","lm.inter.pval","lm.slope.pval","lm.overall.pval"))
res.cor[,lm.slope.pval.sig:=ifelse(lm.slope.pval<0.001,"***", ifelse(lm.slope.pval<0.01,"**", ifelse(lm.slope.pval<0.05,"*", "-")))]
res.cor[,lm.slope.padj:=p.adjust(lm.slope.pval)]
res.cor[,lm.slope.padj.sig:=ifelse(lm.slope.padj<0.001,"***", ifelse(lm.slope.padj<0.01,"**", ifelse(lm.slope.padj<0.05,"*", "-")))]
res.cor[, lm.overall.pval:=as.numeric(lm.overall.pval)]
res.cor[,lm.overall.pval.sig:=ifelse(lm.overall.pval<0.001,"***", ifelse(lm.overall.pval<0.01,"**", ifelse(lm.overall.pval<0.05,"*", "-")))]
res.cor[,lm.overall.padj:=p.adjust(lm.overall.pval)]
res.cor[,lm.overall.padj.sig:=ifelse(lm.overall.padj<0.001,"***", ifelse(lm.overall.padj<0.01,"**", ifelse(lm.overall.padj<0.05,"*", "-")))]
res.cor[,cor.padj:=p.adjust(cor.pval)]
res.cor[,cor.padj.sig:=ifelse(cor.padj<0.001,"***", ifelse(cor.padj<0.01,"**", ifelse(cor.padj<0.05,"*", "-")))]

res.cor
# group   rho cor.pval cor.sig lm.inter lm.slope lm.inter.pval lm.slope.pval lm.overall.pval lm.slope.pval.sig lm.slope.padj
# 1:    Linf1 -0.05    0.031       *    9e-04  0.00025             0       0.53468    5.346823e-01                 -       1.00000
# 2:    Ldon1 -0.07   0.0068      **  0.00019    1e-05             0       0.90075    9.007469e-01                 -       1.00000
# 3:    Ldon2  0.04   0.4941       -  0.00028  0.00015             0       0.13563    1.356325e-01                 -       0.54252
# 4:    Ldon5  0.06   0.3428       -  0.00144  0.00077         1e-05       0.16703    1.670334e-01                 -       0.54252
# 5:    Ldon3 -0.17        0     ***  0.01327 -0.00538             0             0    4.110694e-06               ***       0.00000
# 6:  CH_Linf -0.07   0.3477       -   0.0126 -0.00726             0       0.03959    3.959247e-02                 *       0.19795
# 7: CUK_Linf  0.02   0.7158       -  0.00929 -0.00655             0        0.0099    9.896757e-03                **       0.05940
# lm.slope.padj.sig lm.overall.pval.sig lm.overall.padj lm.overall.padj.sig cor.padj cor.padj.sig
# 1:                 -                   -    1.000000e+00                   -   0.1550            -
#   2:                 -                   -    1.000000e+00                   -   0.0408            *
#   3:                 -                   -    5.425301e-01                   -   1.0000            -
#   4:                 -                   -    5.425301e-01                   -   1.0000            -
#   5:               ***                 ***    2.877486e-05                 ***   0.0000          ***
#   6:                 -                   *    1.979624e-01                   -   1.0000            -
#   7:                 -                  **    5.938054e-02                   -   1.0000            -
  
 


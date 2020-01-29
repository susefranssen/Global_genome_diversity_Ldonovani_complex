

library(data.table)
library(foreach)
library(ggplot2)
library(reshape2)
library(ape)
library(vegan)
library(topGO)

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/10_CNVs/")
dir.create("CNV_gene", showWarnings = F)
setwd("CNV_gene")



#----------------------------------
# functions

cormat_heatmap <-function(cormat, reorder=T, cor.stat="Pearson\nCorrelation")
{
  library(ggplot2)
  library(reshape2)
  
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  # Reorder the correlation matrix
  if (reorder) {cormat <- reorder_cormat(cormat)}
  upper_tri <- get_upper_tri(cormat)
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  melted_cormat$value<-round(melted_cormat$value,2)
  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                         midpoint = 0.5, limit = c(0,1), space = "Lab",
                         name=cor.stat) +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                     size = 12, hjust = 1))+
    coord_fixed() +
    # geom_text(aes(Var2, Var1, label = value), color = "black", size = 4) + ######
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          legend.justification = c(1, 0),
          legend.position = c(0.6, 0.7),
          legend.direction = "horizontal") +
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  # Print the heatmap
  print(ggheatmap)
  
}

#----

# GOid: to search
# GO2geneID: GO to gene mapping
get_genes_for_GO <- function(GOid, GO2geneID)
{
  if (GOid %in% names(GO2geneID))
  {
    ii<- GOid == names(GO2geneID)
    GO2geneID[ii][[1]]
  } else {
    NA
  }
}



#----------------------------------


# load somy (median chr coverage data not to be used from this file)
load(file="../somy_all.RData")
setnames(all, c("chr","sample"), c("chrom","samp"))
# sample.info of groups
load(file="../../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData") 


#####################
#
# get coverages for samples from the various study origins
#
#####################
#
all[,RDmedianHaploid:=RDmedian/Somy]
med.med.cov<-all[,.(samp.median.haploid=median(RDmedianHaploid)), by=samp]
med.med.cov[, source:="thisStudy"]
med.med.cov[samp %in% c("BHU1065A1","BHU1137A1","BHU1139A1","BHU220A1","BHU816A1","BHU824A1","BHU931A1","BHU1064A1",
                        "BPK164A1","BPK649A1","BPK029A1","BPK035A1","BPK067A1","BPK077A1","BPK157A1",
                        "BPK294A1","BPK471A1","BPK562A1","BPK612A1","BPK623A1","BPK648A1","BPK156A1",
                        "BPK413A1","BPK406A1","BPK512A1",
                        "BD09","BD12","BD14","BD15","BD17","BD21","BD24","BD25"), source:="Imamura"]
# "Ldon282cl2","BPK282I9","BD27"
med.med.cov[samp %in% c("LinJPCM5"), source:="Peacock"]
med.med.cov[samp %in% c("CUK2","CUK3","CUK4","CUK5","CUK6","CUK7","CUK8","CUK10",
                        "CUK11","CUK12","CUK9"), source:="Rogers"]
med.med.cov[samp %in% c("X356WTV","X363SKWTI","X364SPWTII","X383WTI","AM560WTI","AM563WTI"), source:="Zackay"]
med.med.cov[samp %in% c("OVN3","CL.SL"), source:="Zhang"]

table(med.med.cov[,source])
# Imamura   Peacock    Rogers thisStudy    Zackay     Zhang 
# 33         1        11        98         6         2 
summary(med.med.cov[source=="thisStudy", samp.median.haploid])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 9.75   20.12   26.73   28.38   33.75   87.58 


#####################
#
# calculate gene copy number
#
#####################
#
# median chromosome coverage based on median chr coverage across medians of 5kb windows
# the code how this was calaculated is provided in script a14_genome_cov.r
winsize=5000
load(file=paste0("../all.genome.cov.median_",winsize/1000,"kb.RData"))
genome.cov[,chr.medcov_across_win:=median(cov.median.win), by=.(sample,chr)]
# Note: this can only be run when the intermediate files as described in a14_gene_cov.sh have been generated before
# otherwise the results table that was generated in the following if statement is provided and can be loaded  with the load cmd after the if statement.
if (F)
{
  all.gene.cp<-foreach(sample.c=sample.info$sample, .combine=rbind) %do%
  # all.gene.cp<-foreach(sample.c=c("Peking","D_2"), .combine=rbind) %do%
  {
    # sample.c="AG83"
    
    somy<-all[samp==sample.c, .(chrom, RDmedian, Somy)]
    
    # files in /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38
    ss<-sample.c
    ss<-gsub("[.]", "-", ss)
    ss<-gsub("X", "", ss)
    
    if (server)
    {
      file=paste0("/lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/linj.",ss,".sorted.markdup.realigned.PP.rmdup.gene.cov.red")
      
    } else {
      file=paste0("/Volumes/sf18/tmp/linj.",ss,".sorted.markdup.realigned.PP.rmdup.gene.cov.red")
    }
    #
    if (!file.exists(file))
    {
      print(paste("file", file, " does not exist!"))
    } else {
      print(paste("sample", file))
    }
    
    dat<-data.table(read.table(file))
    setnames(dat,colnames(dat),c("chr","spos","gene","cov"))
    
    dat[,cov.mean:=mean(cov), by=gene]
    dat[,cov.sd:=sd(cov), by=gene]
    dat[,cov.med:=as.double(median(as.double(cov))), by=gene]
    dat.red<-unique(dat[,.(chr,spos,gene,cov.mean,cov.sd,cov.med)])
    # dat.red<-dat.red[chr !="LinJ.00"]
    dat.red[,chr.ori:=chr] # original chr naming including LinJ.00
    dat.red[,gene:=as.character(gene)]
    dat.red[,chr:=substr(dat.red[,gene],1,7)] # also assign LinJ.00 genes to chr
    dat.red[,spos.ori:=spos]
    dat.red[chr.ori=="LinJ.00",spos:=-spos]
    
    # add somy and chr median depth information
    for (cc in as.character(somy[,chrom]))
    {
      # cc="LinJ.01"
      dat.red[chr==cc, chr.medcov_across_win:=unique(genome.cov[sample==sample.c & chr==as.numeric(substr(cc,6,7)), chr.medcov_across_win])]
      # dat.red[chr==cc, chr.RDmed:=somy[chrom==cc, RDmedian]] # RDmedian calcuated by caroline
      dat.red[chr==cc, chr.somy:=somy[chrom==cc, Somy]]
    }
    #  gen cov normalized by haplid chr depth
    dat.red[,cov.mean.norm:=cov.mean/(chr.medcov_across_win/chr.somy)]
    dat.red[,cov.sd.norm:=cov.sd/(chr.medcov_across_win/chr.somy)]
    dat.red[,absdiff:=abs(cov.mean.norm-chr.somy)]
    dat.red[,diff:=round(cov.mean.norm-chr.somy, 2)]
    # dat.red[, copych:=round(absdiff)] #round up if larger than 0.5
    # dat.red[diff<0, copych:=copych*(-1)] #round up if larger than 0.5
    dat.red[,sample:=sample.c]
    dat.red[,group:=sample.info[sample==sample.c,groups]]
    
    # dat.red[,.(sample, group, chr, spos, gene, chr.somy, cov.mean,cov.sd,cov.med,chr.medcov_across_win, cov.mean.norm, cov.sd.norm, absdiff, copych)]
    dat.red[,.(sample, group, chr, spos, gene, chr.somy, cov.mean,cov.sd,cov.med,chr.medcov_across_win, cov.mean.norm, cov.sd.norm, absdiff, diff)]
  }
  save(all.gene.cp,file="all.gene.cp.RData")
}
load(file="all.gene.cp.RData")
#
#
aa<-all.gene.cp
aa[, CN.abs:=chr.somy+diff]
aa[,round.diff:=round(diff)]
aa[,round.diff_gene.median:=median(round.diff), by=.(gene)]
aa[,round.diff_gene.mean:=round(mean(round.diff),2), by=.(gene)]
aa[,round.diff_gene.sd:=round(sd(round.diff),2), by=.(gene)]
#
aa[,round.diff_gene.zeroCount:=sum(round.diff==0), by=.(gene)] # by gene: number of samples with diff.round == 0
aa[,round.diff_gene.num_g0:=sum(round.diff>0), by=.(gene)] # by gene: number of samples with diff.round greater 0
aa[,round.diff_gene.num_s0:=sum(round.diff<0), by=.(gene)] # by gene: number of samples with diff.round smaller 0
aa[,round.diff_gene.num_gs0:=sum(abs(round.diff)>0), by=.(gene)] # by gene: number of samples with diff.round greater and/or smaller 0
#
aa[,cov.mean:=round(cov.mean, 2)]
aa[,cov.sd:=round(cov.sd, 2)]
aa[,cov.mean.norm:=round(cov.mean.norm, 2)]
aa[,cov.sd.norm:=round(cov.sd.norm, 2)]
aa[,absdiff:=round(absdiff, 2)]
aa<-aa[order(group, sample, chr, spos)]
aa[,CN.abs.round:=chr.somy+round.diff]

write.table(aa[group %in% c("CH_Linf","CUK_Linf","Linf1","other_Linf"), 
               .(sample, group, chr, spos, gene, chr.somy, cov.med, cov.sd, chr.medcov_across_win, cov.mean.norm, diff, round.diff, round.diff_gene.median, round.diff_gene.sd, round.diff_gene.num_gs0)], 
            file="Sample_gene_copyn_change_Linf.txt", quote=F, row.names=F)
write.table(aa[group %in% c("Ldon1","Ldon2","Ldon3","Ldon4","Ldon5","other_Ldon"), 
               .(sample, group, chr, spos, gene, chr.somy, cov.med, cov.sd, chr.medcov_across_win, cov.mean.norm, diff, round.diff, round.diff_gene.median, round.diff_gene.sd, round.diff_gene.num_gs0)], 
            file="Sample_gene_copyn_change_Ldon.txt", quote=F, row.names=F)


##########


#####################
#
# copy number changes in genes associated with drug resistance
#
#####################
#
if(F)
{
  aa[,round.diff.frac:=round.diff/chr.somy]
  
  # Miltefosine transporter
  #
  table(aa[gene=="LinJ.13.1590",.(round.diff)])
  # 0   1 
  # 136  15 
  table(aa[gene=="LinJ.13.1590",.(group,round.diff)])
  # round.diff
  # group         0  1
  # CH_Linf     5  0
  # CUK_Linf    3  8
  # Ldon1      43  2
  # Ldon2       5  2
  # Ldon3      17  2
  # Ldon4       4  0
  # Ldon5       7  1
  # Linf1      47  0
  # other_Ldon  3  0
  # other_Linf  2  0
  #
  # gene deleted with Miltefosine transporter in one experiment
  table(aa[gene=="LinJ.13.1600",.(round.diff)])
  # 0   1 
  # 148   3 
  table(aa[gene=="LinJ.13.1600",.(group,round.diff)])
  # round.diff
  # group         0  1
  # CH_Linf     5  0
  # CUK_Linf   11  0
  # Ldon1      44  1
  # Ldon2       7  0
  # Ldon3      19  0
  # Ldon4       4  0
  # Ldon5       8  0
  # Linf1      45  2
  # other_Ldon  3  0
  # other_Linf  2  0
  #
  # associated with mMiltefosine reistance
  table(aa[gene=="LinJ.32.1040",.(round.diff)])
  # 0   1 
  # 150   1 
  table(aa[gene=="LinJ.32.1040",.(group,round.diff)])
  # round.diff
  # group         0  1
  # CH_Linf     5  0
  # CUK_Linf   11  0
  # Ldon1      44  1
  # Ldon2       7  0
  # Ldon3      19  0
  # Ldon4       4  0
  # Ldon5       8  0
  # Linf1      47  0
  # other_Ldon  3  0
  # other_Linf  2  0
  #
  #
  #
  # Miltefosine sensitivity locus (4 genes)
  #
  xx<-aa[ gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400") ]
  table(xx[,round.diff])
  # -8  -4  -3  -2  -1   0   1 
  # 4  11   6  26  96 439  22
  22/604
  # [1] 0.03642384
  (4+1+6+226+96)/604
  # [1] 0.5513245
  #
  # samples with loss of locus
  aa[CN.abs.round==0 & gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400"),.(gene,sample, group, chr.somy, round.diff)]
  length(unique(aa[round.diff.frac==-1 & gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400"), sample]))
  # [1] 4
  # number of genes that have a loss in these samples
  table(table(aa[round.diff.frac==-1 & gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400"), sample]))
  # 3 4 
  # 1 3
  #
  # samples with reduction of locus
  xx<-aa[round.diff.frac<0 & round.diff.frac!=-1 & gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400"), ]
  length(unique(xx[, sample]))
  # [1] 105
  # number of genes that have a reduction in these samples
  table(table(xx[, sample]))
  # 1  2  4 
  # 86 17  2
  # samples that have a reduction in 4 genes
  table(xx[,sample])[table(xx[,sample])==4]
  # CH35 IMT373cl1 
  # 4         4 
  # samples that have a reduction in 2 genes
  table(xx[,sample])[table(xx[,sample])==2]
  # AM563WTI BHU1139A1  BPK156A1  BPK648A1      CH32      CH33      CH34      CH36     CUK12      CUK7      CUK8 LRC_L1303  LRC.L445   LRC.L53 
  # 2         2         2         2         2         2         2         2         2         2         2         2         2         2 
  # LRC.L740     MRC74   X45.UMK 
  # 2         2         2  
  #
  #
  #
  # AQP1 -> LinJ.31.0030
  table(aa[gene=="LinJ.31.0030", .(round.diff)])
  # -2 -1  0  1  2  3 
  # 1  8 87 48  6  1 
  table(aa[gene=="LinJ.31.0030", .(group,round.diff)])
  # round.diff
  # group        -2 -1  0  1  2  3
  # CH_Linf     0  0  2  3  0  0
  # CUK_Linf    0  0  2  8  1  0
  # Ldon1       0  0 24 19  2  0
  # Ldon2       1  0  3  2  1  0
  # Ldon3       0  0 17  2  0  0
  # Ldon4       0  1  2  1  0  0
  # Ldon5       0  3  3  1  0  1
  # Linf1       0  3 31 11  2  0
  # other_Ldon  0  1  2  0  0  0
  # other_Linf  0  0  1  1  0  0
  9/151
  # [1] 0.05960265
  55/151
  # [1] 0.3642384
  #
  #
  # H-locus: YIP1, MRPA, PTR1
  #
  xx<-aa[ gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") ]
  table(xx[,round.diff])
  # -1   0   1   2   3   4   5   6   7   9  37  38 
  # 51 281  15  27  27  28  12   7   1   1   2   1 
  51/453
  # [1] 0.1125828
  (15 + 27 + 27 + 28 + 12 +  7  + 1  + 1  + 2  + 1)/453
  # 0.2671082
  #
  xx<-aa[ gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0300","LinJ.23.0310") ]
  table(xx[,round.diff])
  # -1   0   1   2   3   4   5   6   7   9  11  37  38  44 
  # 55 367  23  37  39  43  19  14   1   1   1   2   1   1
  55/604
  # [1] 0.0910596
  (23  +37 + 39 + 43 + 19 + 14 +  1  + 1  + 1  + 2 +  1 +  1)/604
  # 0.3013245
  #
  # 
  # CN deletions across genes in a sample
  table(aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff<0, sample])
  table(table(aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff<0, sample]))
  # 1  2  3 
  # 36  6  1 --> 1 sample that has all 3 genes reduced
  xx<-aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff<0,]
  table(xx[,.(sample)])[table(xx[,.(sample)])==3]
  # X1026.8 -> sample with all 3 genes reduced
  # 3 
  table(xx[,.(sample)])[table(xx[,.(sample)])==2]
  # BPK406A1     CH36     CUK4     CUK9       EP  LRC.L57 -> sample2 with 2 genes reduced
  # 2        2        2        2        2        2
  #
  # all samples with at least one gene with reduced gene copy number
  length(unique((xx[,sample])))
  # [1] 43
  unique((xx[,sample]))
  # [1] "CH33"       "CH34"       "CH35"       "CH36"       "CUK10"      "CUK2"       "CUK3"       "CUK4"       "CUK8"      
  # [10] "CUK9"       "BD09"       "BD12"       "BD14"       "BD21"       "BD22"       "BD24"       "BD25"       "BPK164A1"  
  # [19] "L60b"       "OVN3"       "BPK406A1"   "BPK512A1"   "BPK612A1"   "BPK623A1"   "BPK648A1"   "X1026.8"    "X597LN"    
  # [28] "BUMM3"      "AM560WTI"   "LRC.L445"   "LRC.L53"    "LRC.L57"    "MRC74"      "NLB.323"    "Cha001"     "DOG_STRAIN"
  # [37] "ISS2426"    "LEM1985"    "LRC.L1296"  "LRC_L1303"  "SKIN"       "GE"         "EP" 
  #
  # gene copy numbers are only reduced by 1
  table(xx[,round.diff])
  # -1 
  # 51
  # gene tht have the deletion
  table(xx[,gene])
  # LinJ.23.0280 LinJ.23.0290 LinJ.23.0310 
  # 28            5           18 
  # 
  # group distribution of samples with at least one reduced gene copy number
  table(unique(xx[,.(sample,group)])[,group])
  # CH_Linf   CUK_Linf      Ldon1      Ldon2      Ldon3      Ldon4      Ldon5      Linf1 other_Ldon other_Linf 
  # 4          6         10          5          2          1          6          7          1          1
  #
  # CN insertions across genes in a sample
  table(aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff>0, sample])
  table(table(aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff>0, sample]))
  # 1  2  3 
  # 5 55  2 --> 2 samples that has all 3 genes increased
  57/151
  # [1] 0.3774834
  xx<-aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0310") & round.diff>0,]
  table(xx[,.(sample)])[table(xx[,.(sample)])==3]
  # BPK157A1 LRC.L51p -> samples with all 3 genes increased
  # 3        3
  table(xx[,.(sample)])[table(xx[,.(sample)])==2]
  names(table(xx[,.(sample)])[table(xx[,.(sample)])==2])
  # [1] "AG83"       "BD09"       "BD12"       "BD14"       "BD15"       "BD17"       "BD21"       "BD22"       "BD24"       "BD25"      
  # [11] "BD27"       "BHU1062.4"  "BHU1064A1"  "BHU1065A1"  "BHU1137A1"  "BHU1139A1"  "BHU220A1"   "BHU41"      "BHU816A1"   "BHU824A1"  
  # [21] "BHU931A1"   "BPK029A1"   "BPK035A1"   "BPK067A1"   "BPK077A1"   "BPK164A1"   "BPK282I9"   "BPK294A1"   "BPK471A1"   "BPK562A1"  
  # [31] "BPK649A1"   "BUCK"       "Chowd5"     "CL.SL"      "DD8"        "Don201"     "GEBRE1"     "IECDR1"     "Inf206"     "Ldon282cl2"
  # [41] "LdonLV9"    "LRC.L61"    "Malta33"    "Nandi"      "STL2.78"    "STL2.79"    "SUDAN1"     "X1S"        "X356WTV"    "X363SKWTI" 
  # [51] "X364SPWTII" "X38.UMK"    "X383WTI"    "X45.UMK"    "X762L"
  table(xx[,.(sample)])[table(xx[,.(sample)])==1]
  names(table(xx[,.(sample)])[table(xx[,.(sample)])==1])
  # [1] "Cha001"   "Inf001"   "Inf045"   "LRC.L740" "Peking" 
  #
  #
  table(table(aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0300","LinJ.23.0310") & round.diff>0, sample]))
  # 1  2  3  4 
  # 4  4 54  2 --> 2 samples that has all 4 genes increased
  56/151
  # [1] 0.3708609
  xx<-aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0300","LinJ.23.0310") & round.diff>0,]
  table(xx[,.(sample)])[table(xx[,.(sample)])==4]
  # BPK157A1 LRC.L51p -> samples with all 4 genes increased
  # 4        4
  table(xx[,.(sample)])[table(xx[,.(sample)])==3]
  names(table(xx[,.(sample)])[table(xx[,.(sample)])==3])
  # [1] "AG83"       "BD09"       "BD12"       "BD14"       "BD15"       "BD17"       "BD21"       "BD22"       "BD24"       "BD25"      
  # [11] "BD27"       "BHU1062.4"  "BHU1064A1"  "BHU1065A1"  "BHU1137A1"  "BHU1139A1"  "BHU220A1"   "BHU41"      "BHU816A1"   "BHU824A1"  
  # [21] "BHU931A1"   "BPK029A1"   "BPK035A1"   "BPK067A1"   "BPK077A1"   "BPK164A1"   "BPK282I9"   "BPK294A1"   "BPK471A1"   "BPK562A1"  
  # [31] "BPK649A1"   "BUCK"       "Chowd5"     "CL.SL"      "DD8"        "Don201"     "GEBRE1"     "IECDR1"     "Inf206"     "Ldon282cl2"
  # [41] "LdonLV9"    "LRC.L61"    "Malta33"    "Nandi"      "STL2.78"    "STL2.79"    "SUDAN1"     "X1S"        "X356WTV"    "X363SKWTI" 
  # [51] "X364SPWTII" "X38.UMK"    "X383WTI"    "X45.UMK"    "X762L"
  table(xx[,.(sample)])[table(xx[,.(sample)])==2]
  names(table(xx[,.(sample)])[table(xx[,.(sample)])==2])
  # [1] "Cha001"   "CL.SL"    "Inf001"   "LRC.L740"
  names(table(xx[,.(sample)])[table(xx[,.(sample)])==1])
  # [1] "Don081" "Inf007" "Inf045" "Peking"
  #
  table(xx[sample %in% names(table(xx[,.(sample)])[table(xx[,.(sample)])==3]), gene])
  # LinJ.23.0280 LinJ.23.0290 LinJ.23.0300 
  # 54           54           54
  #
  #
  #
  #
  # MAPK 1, LinJ.36.6760
  #
  table(aa[gene=="LinJ.36.6760", round.diff])
  # -1  0  1  2  3  4  5 12 13 14 18 19 21 22 23 24 28 29 31 33 34 37 41 
  # 1 82  3 11  3  4  1 10  1  1  1  1  3  4 11  2  1  1  1  2  4  1  2 
  table(aa[gene=="LinJ.36.6760", .(group,round.diff)])
  # round.diff
  # group        -1  0  1  2  3  4  5 12 13 14 18 19 21 22 23 24 28 29 31 33 34 37 41
  # CH_Linf     0  5  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 0
  # CUK_Linf    0 11  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 0
  # Ldon1       0  0  0  0  0  0  0 10  1  0  1  1  3  4 11  2  1  1  1  2  4  1  2 --> all increased (+12 to +41)
  # Ldon2       0  7  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 0
  # Ldon3       0  0  2  9  3  4  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> all increased (+1 to +5)
  # Ldon4       0  3  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 1/4 increased
  # Ldon5       0  7  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 1/8 increased
  # Linf1       1 46  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 1/47 decreased
  # other_Ldon  0  1  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 2/3 increased
  # other_Linf  0  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 --> 0
  table(aa[gene=="LinJ.36.6760" & round.diff>0, group])
  # Ldon1      Ldon3      Ldon4      Ldon5 other_Ldon 
  # 45         19          1          1          2    -->68
  table(aa[gene=="LinJ.36.6760" & round.diff==0, group])
  # CH_Linf   CUK_Linf      Ldon2      Ldon4      Ldon5      Linf1 other_Ldon other_Linf 
  # 5         11          7          3          7         46          1          2      --> 82
  table(aa[gene=="LinJ.36.6760" & round.diff<0, group])
  # Linf1 
  # 1 

  # ICS5 highly antimony resistant, Imamura 2016
  # BPK164A1
  # BPK649A1
  
  aa[gene=="LinJ.36.6760" & round.diff>0, .(group,sample,round.diff)]
  # group     sample round.diff
  # 1:      Ldon1       AG83         23
  # 2:      Ldon1       BD09         12
  # 3:      Ldon1       BD12         12
  # 4:      Ldon1       BD14         12
  # 5:      Ldon1       BD15         12
  # 6:      Ldon1       BD17         12
  # 7:      Ldon1       BD21         12
  # 8:      Ldon1       BD22         12
  # 9:      Ldon1       BD24         12
  # 10:      Ldon1       BD25         13
  # 11:      Ldon1       BD27         12
  # 12:      Ldon1  BHU1062.4         23
  # 13:      Ldon1  BHU1064A1         22
  # 14:      Ldon1  BHU1065A1         22
  # 15:      Ldon1  BHU1137A1         24
  # 16:      Ldon1  BHU1139A1         21
  # 17:      Ldon1   BHU220A1         22
  # 18:      Ldon1      BHU41         24
  # 19:      Ldon1   BHU816A1         23
  # 20:      Ldon1   BHU824A1         23
  # 21:      Ldon1   BHU931A1         21
  # 22:      Ldon1   BPK029A1         23
  # 23:      Ldon1   BPK035A1         23
  # 24:      Ldon1   BPK067A1         19
  # 25:      Ldon1   BPK077A1         23
  # 26:      Ldon1   BPK157A1         18
  # 27:      Ldon1   BPK164A1         34 !!!
  # 28:      Ldon1   BPK282I9         29
  # 29:      Ldon1   BPK294A1         23
  # 30:      Ldon1   BPK471A1         23
  # 31:      Ldon1   BPK562A1         31
  # 32:      Ldon1   BPK649A1         33 !!!
  # 33:      Ldon1       BUCK         37
  # 34:      Ldon1      CL.SL         34
  # 35:      Ldon1     Chowd5         34
  # 36:      Ldon1        DD8         23
  # 37:      Ldon1     Don201         22
  # 38:      Ldon1     IECDR1         12
  # 39:      Ldon1     Inf206         21
  # 40:      Ldon1       L60b         23
  # 41:      Ldon1 Ldon282cl2         28
  # 42:      Ldon1      Nandi         34
  # 43:      Ldon1       OVN3         33
  # 44:      Ldon1    STL2.78         41
  # 45:      Ldon1    STL2.79         41
  # 46:      Ldon3     GEBRE1          2
  # 47:      Ldon3     GILANI          2
  # 48:      Ldon3    LRC.L61          5
  # 49:      Ldon3    LdonLV9          2
  # 50:      Ldon3    Malta33          2
  # 51:      Ldon3     SUDAN1          4
  # 52:      Ldon3    X1026.8          2
  # 53:      Ldon3        X1S          2
  # 54:      Ldon3    X356WTV          2
  # 55:      Ldon3  X363SKWTI          3
  # 56:      Ldon3 X364SPWTII          2
  # 57:      Ldon3    X38.UMK          4
  # 58:      Ldon3    X383WTI          3
  # 59:      Ldon3    X45.UMK          4
  # 60:      Ldon3     X452BM          3
  # 61:      Ldon3     X597.2          1
  # 62:      Ldon3     X597LN          1
  # 63:      Ldon3      X762L          4
  # 64:      Ldon3     X855.9          2
  # 65:      Ldon4     Don081         14
  # 66:      Ldon5   LRC.L51p          1
  # 67: other_Ldon         GE          2
  # 68: other_Ldon    LEM3472          2
  # group     sample round.diff
  
  # MAPK 13, LinJ.35.5330
  #
  table(aa[gene=="LinJ.35.5330", round.diff])
  # 0   1   2   4   5 
  # 141   4   4   1   1 
  table(aa[gene=="LinJ.35.5330", .(group,round.diff)])
  # round.diff
  # group         0  1  2  4  5
  # CH_Linf     5  0  0  0  0
  # CUK_Linf   11  0  0  0  0
  # Ldon1      43  0  2  0  0 -> 2/(2+43)
  # Ldon2       7  0  0  0  0
  # Ldon3      14  2  1  1  1 -> 5/(5+14)
  # Ldon4       3  0  1  0  0 -> 1/(4)
  # Ldon5       8  0  0  0  0
  # Linf1      46  1  0  0  0 -> 1/(47)
  # other_Ldon  2  1  0  0  0 -> 1/(3)
  # other_Linf  2  0  0  0  0
  aa[gene=="LinJ.35.5330" & round.diff>0, .(group,sample,round.diff)]
  # group    sample round.diff
  # 1:      Ldon1     BHU41          2
  # 2:      Ldon1  BPK649A1          2
  # 3:      Ldon3   LRC.L61          1
  # 4:      Ldon3   Malta33          5
  # 5:      Ldon3       X1S          1
  # 6:      Ldon3 X363SKWTI          2
  # 7:      Ldon3   X383WTI          4
  # 8:      Ldon4   SUKKAR2          2
  # 9:      Linf1   LRC_L47          1
  # 10: other_Ldon   LEM3472          1
  
  # H-locus
  # LinJ.23.0280
  # LinJ.23.0290
  # LinJ.23.0300
  # LinJ.23.0310
  # MAPK
  # LinJ.35.5330
  # AQP1
  # LinJ.31.0030
  # MIltefosine transporter and associated genes
  # LinJ.13.1590
  # LinJ.13.1600
  # LinJ.32.1040
  # Miltefosine sensitivity locus
  # LinJ.31.2370
  # LinJ.31.2380
  # LinJ.31.2390
  # LinJ.31.2400
  
  zz<-aa[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0300","LinJ.23.0310","LinJ.35.5330","LinJ.31.0030",
                 "LinJ.13.1590","LinJ.13.1600","LinJ.32.1040","LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400","LinJ.36.6760")]
  zz[gene %in% c("LinJ.23.0280","LinJ.23.0290","LinJ.23.0300","LinJ.23.0310"), locus:="H-locus"]
  zz[gene %in% "LinJ.35.5330", locus:="MAPK13"]
  zz[gene %in% "LinJ.36.6760", locus:="MAPK1"]
  zz[gene %in% "LinJ.31.0030", locus:="AQP1"]
  zz[gene %in% c("LinJ.13.1590","LinJ.13.1600","LinJ.32.1040"), locus:="Miltefosine_transporter_and_associated"]
  zz[gene %in% c("LinJ.31.2370","LinJ.31.2380","LinJ.31.2390","LinJ.31.2400"), locus:="MSL"]
  
  dir.create("CN_drug_resis_genes", showWarnings = F)
  
  for (l in unique(zz$locus))
  {
    pdf(paste0("CN_drug_resis_genes/locus_CNV_",l,".pdf"), 
        width=4, height=length(unique(zz[locus ==l, gene]))*1.1+2 )
    gg <-ggplot(zz[locus ==l], aes(x=group, y=round.diff, col=group)) + 
      # geom_violin() + geom_jitter(height = 0, width = 0.1) +
      geom_boxplot() + geom_jitter(shape=16, height=0.1, width=0.2) +#position=position_jitter(0.1)) +
      facet_grid(gene~locus, scale="free_y") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_manual(values=unique(sample.info[order(groups),group.col])) 
    if (l %in% unique(zz$locus)[c(3,5)]) {gg <- gg+ ylim(-2,+6)}
    if (l %in% unique(zz$locus)[c(1)]) {gg <- gg+ ylim(-2,+3)}
    # if (l %in% unique(zz$locus)[c(2)]) {gg <- gg+ scale_y_continuous(trans='log2')}
    if (l %in% unique(zz$locus)[c(2)]) {gg <- gg+ ylim(-2,+7)}
    plot(gg)
    dev.off()
  }
}



#####################
#
# summary stats on gene copy number for all genes
#
#####################
#
length(unique(aa$gene))
# [1] 8330
# for each sample and gene the copynumber change has been determined --> round.diff 
# this is now summarised for each gene by the mean, median and number of samples with CN >0 across samples
aa_round.diff_gene.stat<-unique(aa[,.(gene, round.diff_gene.median, round.diff_gene.mean, round.diff_gene.sd, round.diff_gene.zeroCount,
                                      round.diff_gene.num_g0, round.diff_gene.num_s0, round.diff_gene.num_gs0)])[order(gene)]


# histogram for frequency of genes with no CNV to the reference
aa_round.diff_gene.stat[round.diff_gene.num_gs0==0, num_gs0:="no samples"]
aa_round.diff_gene.stat[round.diff_gene.num_gs0>0 & round.diff_gene.num_gs0 <76, num_gs0:="< half of samples"]
aa_round.diff_gene.stat[round.diff_gene.num_gs0<151 & round.diff_gene.num_gs0 >=76, num_gs0:=">= half of samples"]
aa_round.diff_gene.stat[round.diff_gene.num_gs0==151, num_gs0:="all samples"]
pdf("Gene_CN_hist_min1_change.pdf", width=6.5, height=4)
ggplot(data=aa_round.diff_gene.stat, aes(round.diff_gene.num_gs0), col=zeroCount) + 
  # scale_colour_manual("Density", values = c("red", "black")) +
  annotate("rect", xmin=-0.7, xmax=0.6, ymin=0, ymax=710, alpha=0.5, fill="red") +
  annotate("rect", xmin=0.6, xmax=75.6, ymin=0, ymax=710, alpha=0.25, fill="orange") +
  annotate("rect", xmin=75.6, xmax=150.4, ymin=0, ymax=710, alpha=0.25, fill="green") +
  annotate("rect", xmin=150.4, xmax=151.6, ymin=0, ymax=710, alpha=0.5, fill="blue") +
  geom_histogram( binwidth = 1, col="black", fill="white", alpha = 0.1) + theme_bw() +
  labs(x="Number of isolates with gene CNV by gene", y="Gene count") 
dev.off()
# number of genes in the 4 different categories
table(aa_round.diff_gene.stat[,num_gs0])
# < half of samples >= at least half of samples                 all samples                  no samples 
# 7232                         332                          61                         705 
table(aa_round.diff_gene.stat[,num_gs0])/8330
# < half of samples >= at least half of samples                 all samples                  no samples 
# 0.868187275                 0.039855942                 0.007322929                 0.084633854
#
# look at genes that have copy number changes in all isolates
table(aa_round.diff_gene.stat[round.diff_gene.num_gs0==151, .(round.diff_gene.num_g0, round.diff_gene.num_s0)])
#                     round.diff_gene.num_s0
# round.diff_gene.num_g0   0  1  3 86 151
#                     0    0  0  0  0   5
#                     65   0  0  0  1   0
#                     148  0  0  2  0   0
#                     150  0  2  0  0   0
#                     151 51  0  0  0   0
# genes that are of diffeent CN from the reference in all isolates
nrow(aa_round.diff_gene.stat[round.diff_gene.num_gs0==151])
# [1] 61
nrow(aa_round.diff_gene.stat[round.diff_gene.num_gs0==151 & round.diff_gene.num_g0!=151 & round.diff_gene.num_s0!=151])
# [1] 5  -> number of genes , where the CN of all 151 samples differs but some >=+1 and some <=-1
nrow(aa_round.diff_gene.stat[round.diff_gene.num_g0==151])
# [1] 51 -> number of genes , where the CN of all 151 samples differs (>=+1) from the reference (assembly error?)
nrow(aa_round.diff_gene.stat[round.diff_gene.num_s0==151])
# [1] 5 -> number of genes , where the CN of all 151 samples differs (<=-1) from the reference (assembly error?)
#
#
#
# nrow(aa_round.diff_gene.stat[round.diff_gene.median!=0])
# # [1] 304   --> number of genes where the median across samples shows a copy number change
# nrow(aa_round.diff_gene.stat[round.diff_gene.median<0])
# # [1] 103
# nrow(aa_round.diff_gene.stat[round.diff_gene.median>0])
# # [1] 201
# nrow(aa_round.diff_gene.stat[round.diff_gene.median>3])
# # [1] 52
#
# looking at median CN change in different  frequence classes
aa_round.diff_gene.stat$num_gs0<-factor(aa_round.diff_gene.stat$num_gs0, 
                                        levels = c("no samples", "< half of samples", ">= half of samples", "all samples"))
aa_round.diff_gene.stat[round.diff_gene.median<0, med.size:="negmed"]
aa_round.diff_gene.stat[round.diff_gene.median==0, med.size:="zeromed"]
aa_round.diff_gene.stat[round.diff_gene.median>0, med.size:="posmed"]
pdf("Gene_CN_cat_medCHchange.pdf", width=6.5, height=4)
ggplot(aa_round.diff_gene.stat, aes(x=factor(num_gs0), y=round.diff_gene.median, colour=med.size)) +
  geom_boxplot(alpha=0.6) + theme_bw() +
  geom_point(position = position_jitterdodge(0.2)) +
  scale_color_manual(values=c("red", "darkgreen", "blue"), 
                    name="Median gene CN\nchange across\nisolates",
                    labels=c("negative", "positive", "zero")) +
  labs(x="Categories for different amounts of gene CNV across isolates", y="Median CN change by gene") 
dev.off()
#
table(aa_round.diff_gene.stat[med.size!="zeromed", .(med.size, num_gs0)])
#           num_gs0
# med.size no samples < half of samples >= half of samples all samples
# negmed          0                 0                 97           6
# posmed          0                 0                146          55


#####################
#
# median genes copy numbers within and among groups
#
#####################
#
# copy number change
aa[,copych.median.g:=median(diff), by=.(group, gene)]
aa[,copych.sd.g:=sd(diff), by=.(group, gene)]
aa[,copych.median:=median(diff), by=.(gene)]
aa[,copych.sd:=sd(diff), by=.(gene)]
aa[,copych.q10:=quantile(diff, probs=0.10), by=.(gene)]
aa[,copych.q25:=quantile(diff, probs=0.25), by=.(gene)]
aa[,copych.q75:=quantile(diff, probs=0.75), by=.(gene)]
aa[,copych.q90:=quantile(diff, probs=0.90), by=.(gene)]
# absolute CN
aa[,CN.abs.median.g:=median(CN.abs), by=.(group, gene)]
aa[,CN.abs.sd.g:=sd(CN.abs), by=.(group, gene)]
aa[,CN.abs.median:=median(CN.abs), by=.(gene)]
aa[,CN.abs.sd:=sd(CN.abs), by=.(gene)]
aa.medians<-unique(aa[,.(group, chr, spos, gene, copych.median.g, copych.sd.g, copych.median, copych.q10, copych.q25, copych.q75, copych.q90, copych.sd, CN.abs.median.g, CN.abs.sd.g, CN.abs.median, CN.abs.sd)])
#
pdf("Group_median_gene_copyn_change.pdf", width=7,height=5)
# png("AAGroup_median_gene_copyn_change.png", width=800,height=600)
ggplot(aa.medians[!group %in% c("other_Ldon", "other_Linf")], aes(copych.median, copych.median.g, col=group)) +
  geom_point() +
  geom_errorbar(aes(ymin=copych.median.g-copych.sd.g, ymax=copych.median.g+copych.sd.g), width=.15) +
  facet_wrap(~group, ncol=4) +
  geom_abline(intercept = 0, slope = 1, size=0.4) +
  geom_abline(intercept = 2, slope = 1, size=0.2, linetype="dashed") +
  geom_abline(intercept = -2, slope = 1, size=0.2, linetype="dashed") +
  labs(x="Median gene copy number change across all samples",
       y="Median gene copy number change across group samples") +
  scale_colour_manual(values=unique(sample.info[,.(groups,group.col)])[order(groups)][1:8, group.col]) +
  theme_bw()
dev.off()




#####################
#
# nj tree based on gene copy number differences of all samples
#
#####################
#
form.copych.samp<-dcast(aa, formula =gene ~ sample, value.var = "diff")
dist.copych.samp<-dist(form.copych.samp[,2:ncol(form.copych.samp)])
dist.copych.samp<-dist(t(form.copych.samp[,2:152]))
nj.copych.samp<-nj(dist.copych.samp)

tipcol<-foreach(i=nj.copych.samp$tip.label, .combine=c) %do% {sample.info[sample==i,group.col]}
tipgr<-foreach(i=nj.copych.samp$tip.label, .combine=c) %do% {sample.info[sample==i,groups]}
pdf("nj.copych.samp.pdf", width=6, height=12)
plot(nj.copych.samp,label.offset = 5, font = 1, cex=0.4)
tiplabels(pch = 19, col = tipcol, adj = 2.5, cex = 0.8)
dev.off()
#
nj.copych.samp$tip.label<-paste0(nj.copych.samp$tip.label,"_",tipgr)
write.nexus(nj.copych.samp, file="nj.copych.samp.newick")




#####################
#
# GSE - gene set enrichment analysis for genes that show a particular stregth of copy number change#
#####################
#

if(F)
{
  # setup
  #
  geneList<-unique(aa.medians[,.(gene,copych.median)])
  geneList[,Rank.median:=rank(-geneList$copych.median)]
  
  
  # perform enrichment analysis on the ranked list of number of gene copies
  geneID2GO <-readMappings(file = "TriTrypDB-38_LinfantumJPCM5_geneID_to_GOterm.txt")
  GO2geneID <-readMappings(file = "TriTrypDB-38_LinfantumJPCM5_GOterm_to_geneID.txt")
  # length(names(geneID2GO))
  # # [1] 4704 --> # geneIDs
  # length(names(GO2geneID))
  # # [1] 3445 --> # GOterms
  # GOn <-foreach (gID=geneID2GO, .combine=c) %do% {gID}
  # length(unique(GOn))
  # # [1] 3445
  # genen <-foreach (GO=GO2geneID, .combine=c) %do% {GO}
  # length(unique(genen))
  # # [1] 4704
  
  
  # testing for enrichment of genes where the majority of samples has an increased copy number change >=1
  # testing for enrichment of genes where the majority of samples has a decreased copy number change <=-1
  #
  CN_le_1 = function (score) { return(round(score) >0) }
  CN_le_4 = function (score) { return(round(score) >3) }
  CN_se_neg1 = function (score) { return(round(score) <0) }
  
  
  geneL<-geneList$copych.median
  names(geneL)<-geneList$gene
  
  GO_analysis <-function(geneSelFunc, geneL)
  {
    GOdata <- new("topGOdata", ontology = "BP", allGenes = geneL, 
                  annot = annFUN.gene2GO, gene2GO = geneID2GO,
                  geneSelectionFun=geneSelFunc, nodeSize=5)
    
    allnodes = length(GOdata@graph@nodes) #get number of tested nodes
    classic.fish = runTest(GOdata, algorithm = "classic", statistic = "fisher") #do classic fisher test
    weight.fish = runTest(GOdata, algorithm = "weight", statistic = "fisher") #do weight fisher test
    weight01.fish = runTest(GOdata, algorithm = "weight01", statistic = "fisher") #do weight01 fisher test
    elim.fish = runTest(GOdata, algorithm = "elim", statistic = "fisher") #do elim fisher test
    classic.ks = runTest(GOdata, algorithm="classic", statistic="KS") #do classic GSE test
    
    all.res = GenTable(GOdata, classicFish=classic.fish, weightFish=weight.fish,  weight01Fish=weight01.fish, classicKS=classic.ks, elimFish=elim.fish,
                       orderBy = "classicFish", topNodes = allnodes) #get all tested terms
    all.res<-data.table(all.res)
    all.res[,classicKS:=as.numeric(gsub(pattern = "< ", replacement = "", classicKS))]
    all.res[,classicFish:=as.numeric(gsub(pattern = "< ", replacement = "", classicFish))]
    all.res[,weightFish:=as.numeric(gsub(pattern = "< ", replacement = "", weightFish))]
    all.res[,weight01Fish:=as.numeric(gsub(pattern = "< ", replacement = "", weight01Fish))]
    all.res[,elimFish:=as.numeric(gsub(pattern = "< ", replacement = "", elimFish))]
    
    list(GOdata,all.res,classic.fish,weight.fish,weight01.fish,elim.fish)
  }
  
  
  get_res<-function(cand_function, geneList, nametag)
  {
    GOres<-GO_analysis(cand_function, geneL)
    GOres[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)]
    write.table(GOres[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)], file=paste0("GOres_",nametag,"_algocomp_classicAndweight0.05.txt"), 
                quote=T, row.names=F)
    #
    # showGroupDensity(GOres[[1]], GOres[[2]][weightFish<0.05][1,GO.ID], ranks = F, rm.one =F)
    genes<-foreach(id=GOres[[2]][weightFish<0.05, GO.ID], .combine=c) %do%
    {
      get_genes_for_GO(GOid = id, GO2geneID = GO2geneID)
    }
    xx<-geneList[gene %in% genes & cand_function(copych.median),.(gene,copych.median)]
    xx<-xx[,copych.median:=round(copych.median)]
    write.table(xx[order(copych.median)], file=paste0("GOres_",nametag,"_algocomp_enrichmentGenesForweight0.05.txt"), 
                quote=F, row.names=F)
    #
    pdf(paste0("GOres_",nametag,"_algocomp.pdf"),width=24,height=20)
    par(mfrow=c(2,2))
    showSigOfNodes(GOres[[1]], score(GOres[[3]]), firstSigNodes = nrow(GOres[[2]][classicFish<0.05]), useInfo = "all")
    showSigOfNodes(GOres[[1]], score(GOres[[4]]), firstSigNodes = nrow(GOres[[2]][weightFish<0.05]), useInfo = "all")
    showSigOfNodes(GOres[[1]], score(GOres[[5]]), firstSigNodes = nrow(GOres[[2]][weight01Fish<0.05]), useInfo = "all")
    showSigOfNodes(GOres[[1]], score(GOres[[6]]), firstSigNodes = nrow(GOres[[2]][elimFish<0.05]), useInfo = "all")
    dev.off()
  }
  
  # modify geneL >= half of samples and all sample categories
  aa_round.diff_gene.stat[num_gs0==">= half of samples", gene]
  aa_round.diff_gene.stat[num_gs0=="all samples", gene]
  
  
  
  #
  # testing for enrichment of genes where the majority of samples has an increased copy number change >=4
  #
  GOres_CN_le_4<-GO_analysis(CN_le_4, geneL)
  GOres_CN_le_4[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)]
  write.table(GOres_CN_le_4[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)], file="GOres_CN_le_4_algocomp_classicAndweight0.05.txt", 
              quote=T, row.names=F)
  #
  # showGroupDensity(GOres_CN_le_4[[1]], GOres_CN_le_4[[2]][weightFish<0.05][1,GO.ID], ranks = F, rm.one =F)
  genes<-foreach(id=GOres_CN_le_4[[2]][weightFish<0.05, GO.ID], .combine=c) %do%
  {
    get_genes_for_GO(GOid = id, GO2geneID = GO2geneID)
  }
  xx<-geneList[gene %in% genes & CN_le_4(copych.median),.(gene,copych.median)]
  xx<-xx[,copych.median:=round(copych.median)]
  write.table(xx[order(copych.median)], file="GOres_CN_le_4_algocomp_enrichmentGenesForweight0.05.txt", 
              quote=F, row.names=F)
  #
  pdf("GOres_CN_le_4_algocomp.pdf",width=24,height=20)
  par(mfrow=c(2,2))
  showSigOfNodes(GOres_CN_le_4[[1]], score(GOres_CN_le_4[[3]]), firstSigNodes = nrow(GOres_CN_le_4[[2]][classicFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_4[[1]], score(GOres_CN_le_4[[4]]), firstSigNodes = nrow(GOres_CN_le_4[[2]][weightFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_4[[1]], score(GOres_CN_le_4[[5]]), firstSigNodes = nrow(GOres_CN_le_4[[2]][weight01Fish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_4[[1]], score(GOres_CN_le_4[[6]]), firstSigNodes = nrow(GOres_CN_le_4[[2]][elimFish<0.05]), useInfo = "all")
  dev.off()
  
  
  
  #
  # testing for enrichment of genes where the majority of samples has an increased copy number change >=1
  #
  GOres_CN_le_1<-GO_analysis(CN_le_1, geneL)
  GOres_CN_le_1[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)]
  write.table(GOres_CN_le_1[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)], file="GOres_CN_le_1_algocomp_classicAndweight0.05.txt", 
              quote=T, row.names=F)
  #
  # showGroupDensity(GOres_CN_le_1[[1]], GOres_CN_le_1[[2]][weightFish<0.05][1,GO.ID], ranks = F, rm.one =F)
  genes<-foreach(id=GOres_CN_le_1[[2]][weightFish<0.05, GO.ID], .combine=c) %do%
  {
    get_genes_for_GO(GOid = id, GO2geneID = GO2geneID)
  }
  xx<-geneList[gene %in% genes & CN_le_1(copych.median),.(gene,copych.median)]
  xx<-xx[,copych.median:=round(copych.median)]
  write.table(xx[order(copych.median)], file="GOres_CN_le_1_algocomp_enrichmentGenesForweight0.05.txt", 
              quote=F, row.names=F)
  #
  pdf("GOres_CN_le_1_algocomp.pdf",width=24,height=20)
  par(mfrow=c(2,2))
  showSigOfNodes(GOres_CN_le_1[[1]], score(GOres_CN_le_1[[3]]), firstSigNodes = nrow(GOres_CN_le_1[[2]][classicFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_1[[1]], score(GOres_CN_le_1[[4]]), firstSigNodes = nrow(GOres_CN_le_1[[2]][weightFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_1[[1]], score(GOres_CN_le_1[[5]]), firstSigNodes = nrow(GOres_CN_le_1[[2]][weight01Fish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_le_1[[1]], score(GOres_CN_le_1[[6]]), firstSigNodes = nrow(GOres_CN_le_1[[2]][elimFish<0.05]), useInfo = "all")
  dev.off()
  
  
  #
  # testing for enrichment of genes where the majority of samples has a decreased copy number change <=-1
  GOres_CN_se_neg1<-GO_analysis(CN_se_neg1, geneL)
  GOres_CN_se_neg1[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)]
  write.table(GOres_CN_se_neg1[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)], file="GOres_CN_se_neg1_algocomp_classicAndweight0.05.txt", 
              quote=T, row.names=F)
  #
  showGroupDensity(GOres_CN_se_neg1[[1]], GOres_CN_se_neg1[[2]][weightFish<0.05][1,GO.ID], ranks = F, rm.one =F)
  genes<-foreach(id=GOres_CN_se_neg1[[2]][classicFish<0.05, GO.ID], .combine=c) %do%
  {
    get_genes_for_GO(GOid = id, GO2geneID = GO2geneID)
  }
  xx<-geneList[gene %in% genes & CN_se_neg1(copych.median),.(gene, copych.median)]
  xx<-xx[,copych.median:=round(copych.median)]
  write.table(xx[order(copych.median)], file="GOres_CN_se_neg1_algocomp_enrichmentGenesForweight0.05.txt", 
              quote=F, row.names=F)
  #
  pdf("GOres_CN_se_neg1_algocomp.pdf",width=24,height=20)
  par(mfrow=c(2,2))
  showSigOfNodes(GOres_CN_se_neg1[[1]], score(GOres_CN_se_neg1[[3]]), firstSigNodes = nrow(GOres_CN_se_neg1[[2]][classicFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_se_neg1[[1]], score(GOres_CN_se_neg1[[4]]), firstSigNodes = nrow(GOres_CN_se_neg1[[2]][weightFish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_se_neg1[[1]], score(GOres_CN_se_neg1[[5]]), firstSigNodes = nrow(GOres_CN_se_neg1[[2]][weight01Fish<0.05]), useInfo = "all")
  showSigOfNodes(GOres_CN_se_neg1[[1]], score(GOres_CN_se_neg1[[6]]), firstSigNodes = nrow(GOres_CN_se_neg1[[2]][elimFish<0.05]), useInfo = "all")
  dev.off()
  
  
  xx<-geneList[round(copych.median)!=0,.(gene,copych.median)]
  xx<-xx[,copych.median:=round(copych.median)][order(copych.median)]
  write.table(xx, file="GO_all_genes.txt", quote=F, row.names=F)
  
  
  # test gene CN souurounding of AQP1
  #
  Linf<-data.table(read.table(file="Sample_gene_copyn_change_Linf.txt", header = T))
  Ldon<-data.table(read.table(file="Sample_gene_copyn_change_Ldon.txt", header = T))
  
  sub1<-rbind(Linf[gene %in% c("LinJ.31.0010","LinJ.31.0030")], Linf[gene %in% c("LinJ.31.0010","LinJ.31.0030")])[,.(sample,gene,round.diff)]
  sub1[gene=="LinJ.31.0030",AQP1:=round.diff]
  sub1[,AQP1:=max(AQP1, na.rm = T), by=sample]
  # 
  # look at the CN change of LinJ.31.0010 when AQP1 is duplicated 
  nrow(sub1[AQP1>0 & round.diff<=0 & gene=="LinJ.31.0010"])
  # [1] 40
  nrow(sub1[AQP1>0 & round.diff>0 & gene=="LinJ.31.0010"])
  # [1] 12
  
  nrow(sub1[AQP1<0 & round.diff<0 & gene=="LinJ.31.0010"])
  # [1] 2
  nrow(sub1[AQP1<0 & round.diff>=0 & gene=="LinJ.31.0010"])
  # [1] 4
  
  
  
  sub2<-rbind(Linf[gene %in% c("LinJ.31.0040","LinJ.31.0030")], Linf[gene %in% c("LinJ.31.0040","LinJ.31.0030")])[,.(sample,gene,round.diff)]
  sub2[gene=="LinJ.31.0030",AQP1:=round.diff]
  sub2[,AQP1:=max(AQP1, na.rm = T), by=sample]
  
  
}

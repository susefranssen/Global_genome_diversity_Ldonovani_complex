

# generate somy heatmap

library(data.table)
library(ggplot2)
library(foreach)
library(gplots)



path="~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/github_scripts/Global_genome_diversity_Ldonovani_complex"
setwd(paste0(path,"/04_aneuploidy"))



# new groups and names
# 45, #ff0000
Ldon1=c("AG83","BD09","BD12","BD14","BD15","BD17","BD21","BD22","BD24","BD25","BD27","BHU1062.4","BHU1064A1","BHU1065A1","BHU1137A1",
        "BHU1139A1","BHU220A1","BHU41","BHU816A1","BHU824A1","BHU931A1","BPK029A1","BPK035A1","BPK067A1","BPK077A1","BPK157A1","BPK164A1",
        "BPK282I9","BPK294A1","BPK471A1","BPK562A1","BPK649A1","BUCK","CL.SL","Chowd5","DD8","Don201","IECDR1","Inf206","L60b","Ldon282cl2",
        "Nandi","OVN3","STL2.78","STL2.79")
# 7, #7300da
Ldon2=c("BPK156A1","BPK406A1","BPK413A1","BPK512A1","BPK612A1","BPK623A1","BPK648A1")
# 19, #8b1c62
Ldon3=c("LdonLV9","GEBRE1","X356WTV","X363SKWTI","X383WTI","X364SPWTII","LRC.L61","X1S",
        "GILANI","X38.UMK","X45.UMK","X452BM","X597.2","X597LN","X762L","X1026.8",
        "X855.9","SUDAN1","Malta33")
# 4, #ff7f00
Ldon4=c("BUMM3","Don038","Don081","SUKKAR2")
# 8, #ff75ff
Ldon5=c("AM560WTI","AM563WTI","LRC.L445","LRC.L51p","LRC.L53","LRC.L57","MRC74","NLB.323")
#
# 3, #ff0080
Ldon_other=c("GE","LEM3472","LRC.L740")
#
#
# 47, #0000ff
Linf1=c("Peking","D_2","STRAIN_B","STRAIN_A","RACOON_DOG","SKIN","DOG_STRAIN","LRC.L1311","LRC.L1313",
       "Inf007","LRC.L699","LRC.L1275","TH4","TH5","TH6","NT10","NT16","LRC.L1312","LRC.L1296","LRC_L1303",
       "Inf152","ISS174","ISS2420","ISS2426","ISS2429","ISS2508","Inf001","LRC_L47",
       "RM1","LEM1985","LPN114","LEM3278","Inf045","Inf004","Inf055","BCN83",
       "BCN87","LinJPCM5","IMT260","IMT373cl1","ITMAP26","Cha001","ARL","WC","WR285","HN167","HN336")
# 11, #00ff00
CUK=c("CUK10","CUK11","CUK12","CUK2","CUK3","CUK4","CUK5","CUK6","CUK7","CUK8","CUK9")
# 5, #737373
CH=c("CH32","CH33","CH34","CH35","CH36")
# 2, black
Linf_other=c("MAM","EP")

######



ploidy<-data.table(read.table(paste0(path,"/04_aneuploidy/somies_updated.txt"), header = T))
ploidy[,chr:=1:36]
ploidy.t <- melt(ploidy, id.vars = "chr")
setnames(ploidy.t, "variable", "sample")
ploidy.t[,sample:=as.character(sample)]
ploidy.t[sample %in% Ldon1, group:="Ldon1"]
ploidy.t[sample %in% Ldon2, group:="Ldon2"]
ploidy.t[sample %in% Ldon3, group:="Ldon3"]
ploidy.t[sample %in% Ldon4, group:="Ldon4"]
ploidy.t[sample %in% Ldon5, group:="Ldon5"]
ploidy.t[sample %in% Ldon_other, group:="Ldon_other"]
ploidy.t[sample %in% Linf1, group:="Linf1"]
ploidy.t[sample %in% CUK, group:="CUK"]
ploidy.t[sample %in% CH, group:="CH"]
ploidy.t[sample %in% Linf_other, group:="Linf_other"]
#
ploidy.t <- ploidy.t[order(group)]
# ploidy.t <- ploidy.t[!is.na(group)]
#
ploidy.t[sample %in% Ldon1, group.col:="#ff0000"]
ploidy.t[sample %in% Ldon2, group.col:="#7300da"]
ploidy.t[sample %in% Ldon3, group.col:="#8b1c62"]
ploidy.t[sample %in% Ldon4, group.col:="#ff7f00"]
ploidy.t[sample %in% Ldon5, group.col:="#ff75ff"]
ploidy.t[sample %in% Ldon_other, group.col:="#ff0080"]
ploidy.t[sample %in% Linf1, group.col:="#0000ff"]
ploidy.t[sample %in% CUK, group.col:="#00ff00"]
ploidy.t[sample %in% CH, group.col:="#737373"]
ploidy.t[sample %in% Linf_other, group.col:="black"]
#
# check correct number of sample is in each group
table(ploidy.t$group)/36
# CH        CUK Ldon_other      Ldon1      Ldon2      Ldon3      Ldon4      Ldon5 Linf_other      Linf1
# 5         11          3         45          7         19          4          8          2         47



######
# plot somies
hclust.ave <- function(x) hclust(x, method="average")



ploidy<-data.table(read.table(paste0(path,"/04_aneuploidy/somies_updated.txt"), header = T))

rowcols <- rep(NA,ncol(ploidy)-1)
rowcols[colnames(ploidy) %in% Ldon1] <- unique(ploidy.t[group=="Ldon1",group.col])
rowcols[colnames(ploidy) %in% Ldon2] <- unique(ploidy.t[group=="Ldon2",group.col])
rowcols[colnames(ploidy) %in% Ldon3] <- unique(ploidy.t[group=="Ldon3",group.col])
rowcols[colnames(ploidy) %in% Ldon4] <- unique(ploidy.t[group=="Ldon4",group.col])
rowcols[colnames(ploidy) %in% Ldon5] <- unique(ploidy.t[group=="Ldon5",group.col])
rowcols[colnames(ploidy) %in% Ldon_other] <- unique(ploidy.t[group=="Ldon_other",group.col])
rowcols[colnames(ploidy) %in% Linf1] <- unique(ploidy.t[group=="Linf1",group.col])
rowcols[colnames(ploidy) %in% Linf_other] <- unique(ploidy.t[group=="Linf_other",group.col])
rowcols[colnames(ploidy) %in% CH] <- unique(ploidy.t[group=="CH",group.col])
rowcols[colnames(ploidy) %in% CUK] <- unique(ploidy.t[group=="CUK",group.col])



color.palette  <- c(colorRampPalette(c("yellow","orange","green4","blue","red"))(6) ,c("purple","pink"))


pdf(paste0(path,"/04_aneuploidy/Somy2_dendro_groups.pdf"), width=11, height=12)
aa<-t(as.matrix(ploidy[,1:(ncol(ploidy))]))
heatmap.2(aa, dendrogram="row", col= color.palette,
          RowSideColors=rowcols,
          trace='none',
          hclustfun=hclust.ave,
          Colv = FALSE,
          key.xlab = "Ploidy", keysize = 0.6, key.title ="",
          xlab = "Chromosomes")
dev.off()












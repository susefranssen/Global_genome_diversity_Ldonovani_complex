
args <- commandArgs(TRUE)

# comma separated list of functions to perform
mode <- args[1] # c("format","chrtree","constree")
# grouping to choose
dph.tag <- args[2] # dphA, dphB






# mode="format"
# dph.tag="dphC"




print("chosen parameters:")
print(paste("mode", mode))
print(paste("dph.tag", dph.tag))


# SNP file origin:
setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/01_phylogenetic_reconstruction/")

library(data.table)
library(foreach)
library(gplots)
library(StAMPP)
library(ggplot2)
#library(help = StAMPP)



#-----------functions--------------

source("A01_leish_donovaniComplex_StAMPP_functions.R")


#-----------code--------------

dir.create("plots", showWarnings = F)
dir.create("plots/chr", showWarnings = F)
dir.create("NeisD", showWarnings = F)
dir.create("tmp", showWarnings = F)


# getting grouping vector and sample names
#
# location
chr <- "01"
# files are located in github folder data/snp_files
dat=data.table(head(read.table(paste0("../data/snp_files/snps.filt.leish_global.linj.complex.LinJ.",chr,".vcflike.txt"), header=T)))
samp.ori.order <- names(dat)[9:ncol(dat)]

write.table(samp.ori.order,"samp.ori.order.txt", quote=FALSE, row.names=F)
# to get location phenotype
samp.loc <- data.table(read.table("sample_location.txt"))
yy <- data.table(cbind(samp.ori.order,rep("vcf",length(samp.ori.order)))); setnames(yy, "samp.ori.order", "V1")
write.table(yy, file="yy.txt", quote=FALSE, row.names=F)
xx <- merge(yy, samp.loc, by="V1", all=T)
write.table(xx, file="samplebothnames_location.txt", quote=FALSE, row.names=F)
xx[is.na(V2.x)] # samples in the full sheet but shouldn't be included in the analysis

#
# bring to original order
ori.order <- xx[!is.na(V2.x)]$V1 # original order in sheet to be adjusted to order in vcflike
loc1 <- xx[!is.na(V2.x)]$V2.y # matching location order to be adjusted in the same way
loc2 <- xx[!is.na(V2.x)]$V3 # matching location order to be adjusted in the same way
loc3 <- xx[!is.na(V2.x)]$V4 # matching phylo group A order to be adjusted in the same way
loc4 <- xx[!is.na(V2.x)]$V5 # matching phylo group B order to be adjusted in the same way
loc5 <- xx[!is.na(V2.x)]$V6 # matching phylo group C order to be adjusted in the same way
aimed.order <- samp.ori.order # order to adjust to
samp.order.to.use <- as.character(ori.order[order(match(ori.order,aimed.order))]) # this vector should now be identical to samp.ori.order
dph1 <- as.character(loc1[order(match(ori.order,aimed.order))])
dph2 <- as.character(loc2[order(match(ori.order,aimed.order))])
dphA <- as.character(loc3[order(match(ori.order,aimed.order))])
dphB <- as.character(loc4[order(match(ori.order,aimed.order))])
dphC <- as.character(loc5[order(match(ori.order,aimed.order))])
write.table(cbind(samp.order.to.use ,dphA,dphB,dphC), file="samp.order.to.use.txt", quote=FALSE, row.names=F)
#
# order has to be adjusted once more as StAMPP changd it in its internal function
sample.order.StAMPP <- read.table(file="sample.order.StAMPP.txt")
sample.order.StAMPP <- as.character(unlist(sample.order.StAMPP))[2:dim(sample.order.StAMPP)[1]]
ori.order <- samp.order.to.use
aimed.order <- sample.order.StAMPP # order to adjust to
samp.reordered.test <- as.character(ori.order[order(match(ori.order,aimed.order))])
dphAr <- as.character(dphA[order(match(ori.order,aimed.order))])
dphBr <- as.character(dphB[order(match(ori.order,aimed.order))])
dphCr <- as.character(dphC[order(match(ori.order,aimed.order))])
print("end of prep")


# mode="format"
# dph.tag="dphC"
###############
#
# formatting of the data (has to be done prior to any other analysis)
#
if ("format" %in% mode)
{
  print("")
  print("formating and saving data")
  
  # format SNP data into gtStampp format and save a txt file
  for (chr in  c(paste0("0",1:9),10:36))
  # for (chr in  c(paste0("0",1:2)))
  {
    print(chr)
    print("Read dat")
    dat=data.table(read.table(paste0("../data/snp_files/snps.filt.leish_global.linj.complex.LinJ.",chr,".vcflike.txt"), header=T))
      
    print("format_StAMPP")
    gtStampp=format_StAMPP(dat)#, iloci = 1223) chr="02" 
    gtStampp[,Pop:=dph2]
    write.table(gtStampp,file=paste0("tmp/LmjF.",chr,".txt"), quote=FALSE, row.names=F)
  }
  
  # reduce to covered and polymorphic sites, ave stats as txt, save dat as .RData
  NA.frac=0.2 # maximal fraction of missing data alowed at one SNP pos
  #
  snp.stats=data.frame(chr=integer(), n_SNPs=integer(),n_SNPs_poly=integer(),n_SNPs_cov=integer(),SNPs_cov_poly=integer()); i=1
  l.gtStampp.freq <- foreach (chr=c(paste0("0",1:9),10:36)) %do% #c(paste0("0",1:9),10:36)) %do%
  # l.gtStampp.freq <- foreach (chr=c(paste0("0",1:2))) %do% #c(paste0("0",1:9),10:36)) %do%
  {
    print(chr)
    print("Read formatted")
    gtStampp=read.table(file=paste0("tmp/LmjF.",chr,".txt"), header=T) # data.frame object
    if (file.exists(paste0("tmp/LmjF.",chr,".txt")))
    {
      file.remove(paste0("tmp/LmjF.",chr,".txt"))
    }
    print("stampp Convert")
    gtStampp.freq <- stamppConvert(gtStampp, "r")
    n.SNPs=ncol(gtStampp.freq)-5
    n.samp=nrow(gtStampp.freq)
    
    print("reduce to polymorphic sites")
    n.SNPs.poly<-ncol(reduce.poly(gtStampp.freq))-5
    
    print("reduce to sufficiently called SNPs")
    gtStampp.freq<-reduce.cov(gtStampp.freq,NA.frac)
    n.SNPs.cov=ncol(gtStampp.freq)-5
    
    print("reduce to polymorphic sites")
    gtStampp.freq<-reduce.poly(gtStampp.freq)
    n.SNPs.cov.poly=ncol(gtStampp.freq)-5
    
    snp.stats[i,]=c(as.integer(chr), n.SNPs, n.SNPs.poly, n.SNPs.cov, n.SNPs.cov.poly); i=i+1
    
    gtStampp.freq
  }
  
  #####
  ### to order to aim to as previous was messed up by StAMPP coversion
  write.table(l.gtStampp.freq[[1]]$Sample,file="sample.order.StAMPP.txt", quote=FALSE, row.names=F)
  
  save(l.gtStampp.freq, file = paste0("l.gtStampp_cov",NA.frac,"_poly.RData"))
  write.table(snp.stats, file=paste0("plots/snp.stats_cov",NA.frac,"_poly.txt"), quote=FALSE, row.names=F)
  write.table(colSums(snp.stats), file=paste0("plots/snp.stats_cov",NA.frac,"_poly_sums.txt"), quote=FALSE, row.names=T)
}



# mode="chrtree"
# dph.tag="dphC"
###############
#
# calculate trees per chromosome
#
if ("chrtree" %in% mode)
{
  print("")
  print("chrtree")
  
  NA.frac=0.2 # maximal fraction of missing data allowed at one SNP pos
  
  #dph.tag <- "dphB"
  if (dph.tag=="dphC") 
  {
    print("dphC")
    dph.s <-dphC # phenotype / grouping to be used here
  } else
  {
    print("stop")
    stop("Please choose a grouping either dphA or dphB")
  }
 
  load(paste0("l.gtStampp_cov",NA.frac,"_poly.RData"))
  
  chr.trees = foreach (chr=c(1:36)) %do%
  # chr.trees = foreach (chr=c(1:2)) %do%
  {
      gtStampp.freq <- l.gtStampp.freq[[chr]]
      #gtStampp.freq$Pop <- dph.s
      NeisD.ind<-stamppNeisD(gtStampp.freq, pop=F)
      print(paste("chr",chr,"NeisD"))
      write.table(NeisD.ind, file=paste0("NeisD/LmjF.",chr,"_NeisD.ind_NA.frac",NA.frac,"_",dph.tag,".txt"), quote=FALSE, col.names=F, row.names=T)
      tree <- nj(NeisD.ind)
      
      print("Plot - plot neighbor joining tree based on Nei's distance")
      # plot bs neighbor joining tree based on Nei's distance
      # add phenotype group information, reordering dph vector
      samp.ori.order # original order
      aimed.order=tree$tip.label # order to adjust to
      dph.n1=dph1[order(match(samp.ori.order,aimed.order))]
      dph <- dph.s[order(match(samp.ori.order,aimed.order))] # chosen diseade phenotype for coloring
      
      tip.label.col <- data.table(cbind(tree$tip.label, rep("black",length(samp.ori.order)), as.character(dph)))
      tip.label.col[grep("Asia", tip.label.col$V3), V2:="red"]
      tip.label.col[grep("WAsia", tip.label.col$V3), V2:="orange"]
      tip.label.col[grep("Europe", tip.label.col$V3), V2:="blue"]
      tip.label.col[grep("SAmerica", tip.label.col$V3), V2:="yellow"]
      tip.label.col[grep("Africa", tip.label.col$V3), V2:="green"]
      #
      tip.label.col[grep("infantum", tip.label.col$V3), V2:="#0000CD"]
      tip.label.col[grep("donovani", tip.label.col$V3), V2:="#FF0000"]
      tip.label.col[grep("donovani1", tip.label.col$V3), V2:="#FF0000"]
      tip.label.col[grep("donovani2", tip.label.col$V3), V2:="#8B1C62"]
      tip.label.col[grep("donovani2b", tip.label.col$V3), V2:="#FF7F00"]
      tip.label.col[grep("donovani3", tip.label.col$V3), V2:="#FFBBFF"]
      
      tree$tip.label <- paste(tree$tip.label, dph.n1) # countries are added to the sample names
      
      # plot neighbor joining tree based on Nei's distance
      pdf(paste0("plots/chr/LmjF.",chr,"_NA.frac",NA.frac,"_poly_NeisD.pop_Nj_",dph.tag,".pdf"),width=12, height=26)
      plot(root(tree, which(tree$tip.label == "MAM Brazil")), tip.color=tip.label.col$V2)
      xxdim=par("usr") 
      add.scale.bar((xxdim[2]-xxdim[1])*2/3, (xxdim[4]-xxdim[3])/30*2)
      dev.off()
      
      # plot heatmap based on Nei's distance
      pdf(paste0("plots/chr/LmjF.",chr,"_NA.frac",NA.frac,"_poly_NeisD.pop_heat_",dph.tag,".pdf"), width=8, height=16)
      heatmap.2(NeisD.ind,trace="none", RowSideColors=tip.label.col$V2, ColSideColors=tip.label.col$V2)
      dev.off()
      
      write.tree(tree, file = paste0("plots/chr/LmjF.",chr,"_NA.frac",NA.frac,"_poly_NeisD.pop_Nj_newick_",dph.tag,".txt"), append = FALSE,
                 digits = 10, tree.names = FALSE)
  }
}



# mode="constree"
# dph.tag="dphC"
###############
#
# calculate consensus tree across chromosomes incl bootstraps
#
if ("constree" %in% mode)
{
  print("")
  print("constree")
  
  NA.frac=0.2 # maximal fraction of missing data allowed at one SNP pos
  
  if (dph.tag=="dphC") 
  {
    print("dphC")
    dph.s <-dphC # phenotype / grouping to be used here
  } else
  {
    print("stop")
    stop("Please choose a grouping either dphA or dphB")
  }
  #
  load(paste0("l.gtStampp_cov",NA.frac,"_poly.RData"))
  
  NeisD.allchr = foreach (chr=c(1:36)) %do% 
  # NeisD.allchr = foreach (chr=c(1:2)) %do%
  {
    print(paste("NeisD, chr",chr))
    gtStampp.freq <- l.gtStampp.freq[[chr]]
    NeisD.ind=stamppNeisD(gtStampp.freq, pop=FALSE) # Nei's genetic difference between individuals
  }
  print("NeisD.allchr finished")
  
  snp.stats=read.table("plots/snp.stats_cov0.2_poly.txt", header = T)
  snp.weights <- snp.stats$SNPs_cov_poly
  NeisD.chrweighted <- Reduce(`+`,Map(`*`, NeisD.allchr, snp.weights))/sum(snp.weights)
  write.table(NeisD.chrweighted, file=paste0("NeisD/LmjF.allchr.weighted_NeisD.ind_NA.frac",NA.frac,"_",dph.tag,".txt"), quote=FALSE, col.names=F, row.names=T)
  NeisD.chrweighted<-as.matrix(read.table(file=paste0("NeisD/LmjF.allchr.weighted_NeisD.ind_NA.frac",NA.frac,"_",dph.tag,".txt"), row.names = 1))
  tree.cons <- nj(NeisD.chrweighted)
  #
  # plot bs neighbor joining tree based on Nei's distance
  # add phenotype group information, reordering dph vector
  samp.ori.order # original order
  aimed.order=tree.cons$tip.label # order to adjust to
  dph.n1=dph1[order(match(samp.ori.order,aimed.order))]
  dph <- dph.s[order(match(samp.ori.order,aimed.order))] # chosen diseade phenotype for coloring
  
  tip.label.col <- data.table(cbind(tree.cons$tip.label, rep("black",length(samp.ori.order)), as.character(dph)))
  tip.label.col[grep("Asia", tip.label.col$V3), V2:="red"]
  tip.label.col[grep("WAsia", tip.label.col$V3), V2:="orange"]
  tip.label.col[grep("Europe", tip.label.col$V3), V2:="blue"]
  tip.label.col[grep("SAmerica", tip.label.col$V3), V2:="yellow"]
  tip.label.col[grep("Africa", tip.label.col$V3), V2:="green"]
  #
  tip.label.col[grep("infantum", tip.label.col$V3), V2:="#0000CD"]
  tip.label.col[grep("donovani", tip.label.col$V3), V2:="#FF0000"]
  tip.label.col[grep("donovani1", tip.label.col$V3), V2:="#FF0000"]
  tip.label.col[grep("donovani2", tip.label.col$V3), V2:="#8B1C62"]
  tip.label.col[grep("donovani2b", tip.label.col$V3), V2:="#FF7F00"]
  tip.label.col[grep("donovani3", tip.label.col$V3), V2:="#FFBBFF"]
  
  tree.cons$tip.label <- paste(tree.cons$tip.label, dph.n1) # countries are added to the sample names
  
  # plot neighbor joining tree based on Nei's distance
  pdf(paste0("plots/LmjF.allchr.weighted_NeisD.ind_NA.frac",NA.frac,"_poly_NeisD.pop_Nj.boot_",dph.tag,".pdf"),width=12, height=26)
  plot(root(tree.cons, which(tree.cons$tip.label == "MAM Brazil")), tip.color=tip.label.col$V2)
  # legend("bottomright", legend = c(paste0("BS Clades % (",bb," bs repl)")) #c("Bipartitions", "Clades")
  #        , pch = 22, pt.bg = c("lightblue"), pt.cex = 2.5)
  xxdim=par("usr") 
  add.scale.bar((xxdim[2]-xxdim[1])*2/3, (xxdim[4]-xxdim[3])/30*3)
  dev.off()
  
  # plot heatmap based on Nei's distance
  pdf(paste0("plots/LmjF.allchr.weighted_NeisD.ind_NA.frac",NA.frac,"_poly_NeisD.pop_Nj.boot_",dph.tag,"_heat.pdf"),width=12, height=26)
  heatmap.2(NeisD.chrweighted,trace="none", RowSideColors=tip.label.col$V2, ColSideColors=tip.label.col$V2)
  dev.off()
  
  write.tree(tree.cons, file = paste0("plots/LmjF.allchr.weighted_NeisD.ind_NA.frac",NA.frac,"_poly_NeisD.pop_Nj.boot_",dph.tag,"_newick.txt"),
             append = FALSE, digits = 10, tree.names = FALSE)
  
}
















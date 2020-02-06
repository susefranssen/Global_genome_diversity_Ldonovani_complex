

# get Neis distance matrix for windows along the chromosome

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/01_phylogenetic_reconstruction/")
dir.create("window", showWarnings=F)
setwd("window/")

library(data.table)
library(foreach)
library(StAMPP)
library(ggplot2)
library(circlize)
#library(help = StAMPP)


#-----------variables--------------
load("../sample.info_NEWsubgroups.RData")
NA.frac=0.2 # maximal fraction of missing data allowed at one SNP pos

#-----------functions--------------

source("../A01_leish_donovaniComplex_StAMPP_functions.R")


#-----------code--------------


###############
#
# formatting of the data (has to be done prior to any other analysis)
#

print("")
print("formating and saving data")

# format SNP data into gtStampp format for each window and save as a list for each chromosome
for (chr in  c(paste0("0",1:9),10:36))
#for (chr in  c(paste0("0",1:2)))
{
  # chr="01"
  print(chr)
  print("Read dat")
  dat=data.table(read.table(paste0("../../data/snp_files/snps.filt.leish_global.linj.complex.LinJ.",chr,".vcflike.txt"), header=T))
  dat[,w10kb:=POS.CHR%/%10000+1, by=CHR]

  samples<-names(dat)[9:ncol(dat)]
  if(sum(sample.info$sample!=samples[1:151])>0) {warning("sample orders do not match")}

  print("format_StAMPP for each window")
  chr.gtStampp<-foreach(win=1:max(dat$w10kb)) %do%  #sort(unique(dat$w10kb))) %do%
  {
#     print(win)
#     print(dim(dat[w10kb==win,1:(ncol(dat)-1),with=F]))
    if (nrow(dat[w10kb==win,1:(ncol(dat)-1),with=F]) > 0) # window is not empty
    {
      gtStampp=format_StAMPP(dat[w10kb==win,1:(ncol(dat)-1),with=F])
      gtStampp[,Pop:=sample.info$groups]
    } else {
      gtStampp<-NULL
    }
    gtStampp
  }
  save(chr.gtStampp, file=paste0(chr,"_chr.gtStmpp.RData"))
}


print("")
print("formating, reducing and saving data")

# reduce to covered and polymorphic sites, ave stats as txt, save dat as .RData
NA.frac=0.2 # maximal fraction of missing data alowed at one SNP pos
#
for (chr in c(paste0("0",1:9),10:36))
#for (chr in c(paste0("0",1:2)))
{
  # chr="01"
  print(chr)
  print("Read formatted")
  # load chr.gtStampp, list with Stampp formmatted 10kb windows
  load(file=paste0(chr,"_chr.gtStmpp.RData"))
#   if (file.exists(paste0("/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/tmp/",chr,"_chr.gtStmpp.RData")))
#   { file.remove(paste0("/lustre/scratch118/infgen/team133/sf18/06snps/leish_donovaniComplex/tmp/",chr,"_chr.gtStmpp.RData")) }
  print("stampp Convert")

  chr.gtStampp.freq <-foreach(w=1:length(chr.gtStampp)) %do%
  {
    win.gtStampp<-chr.gtStampp[[w]]
    if(!is.null(win.gtStampp))
    {
      print(paste("win",w))
      gtStampp.freq <- stamppConvert(win.gtStampp, "r")

      #     print("reduce to sufficiently called SNPs")
      gtStampp.freq<-reduce.cov(gtStampp.freq,NA.frac)

      #     print("reduce to polymorphic sites")
      gtStampp.freq<-reduce.poly(gtStampp.freq)
    } else {
      gtStampp.freq <-NULL
    }
    gtStampp.freq
  }
  save(chr.gtStampp.freq, file=paste0(chr,"_chr.gtStmpp.freq.RData"))
}





###############
#
# calculate trees per chromosome and windows and save
#
NA.frac=0.2 # maximal fraction of missing data allowed at one SNP pos
dir.create("NeisD", showWarnings = F)

win.count<-list()
# load(file="win.count.RData")
for (chr in c(paste0("0",1:9),10:36))
# for (chr in c(paste0("0",1:2)))
{
  # chr="01"
  dir.create(paste0("NeisD/chr_",chr), showWarnings = F)
  load( file=paste0(chr,"_chr.gtStmpp.freq.RData"))
  #
  win.count[as.integer(chr)]<-length(chr.gtStampp.freq)
  for (win in 1:length(chr.gtStampp.freq))
  {
    print(win)
    win.gtStampp.freq <-chr.gtStampp.freq[[win]]
    if(!is.null(win.gtStampp.freq)) # otherwise no SNPs in this window
    {
      if(ncol(win.gtStampp.freq) > 6) # otherwise only 1 SNPs in this window and stampNeisD does not work then
      {
        NeisD.ind<-stamppNeisD(data.frame(win.gtStampp.freq), pop=F)
        if(sum(rownames(NeisD.ind)!=win.gtStampp.freq$Sample) > 0){warning("cc")} # test sample order # needs to be loaded as was reordered by stamppConvert
        samples<-win.gtStampp.freq$Sample
        Pops<-win.gtStampp.freq$Pop
        NeisD<-data.table(NeisD.ind)
        setnames(NeisD, colnames(NeisD), samples)
        NeisD[,samples:=samples]
        NeisD[,pop:=Pops]
        
        # order sample.info color information by sample order in NeisD
        # sample.info[order(match(sample, NeisD$samples))]$sample == NeisD$samples
        NeisD[,group.col:=sample.info[order(match(sample, NeisD$samples))]$group.col]
        # NeisD[,group.col:=as.integer(unlist(phylocolc[NeisD$pop]))]

        #     tree <- nj(NeisD.ind)
        #
        #     write.tree(tree, file = paste0("NeisD/chr_",chr,"/LmjF.",chr,"_NA.frac",NA.frac,"_poly_NeisD.pop_Nj_newick.txt"), append = FALSE,
        #                digits = 10, tree.names = FALSE)
      } else {
        # as this window has no SNPs NeisD information with all 0s will be loaded
        load("NeisD_empty.RData")
      }
    } else {
      # as this window has no SNPs NeisD information with all 0s will be loaded
      load("NeisD_empty.RData")
    }
#     print(head(NeisD))
    print(paste("chr",chr,"NeisD, win",win))
    write.table(NeisD, file=paste0("NeisD/chr_",chr,"/LmjF.",chr,"_win",win,"_NeisD.ind_NA.frac",NA.frac,".txt"),
                quote=FALSE, col.names=F, row.names=F)
  }
}

save(win.count, file="win.count.RData")
load(file="win.count.RData")



#################
# stop here for phylogentic analysis and continue in script "A01_leish_donovaniComplex_StAMPP_window_bootstrap.R"
# 
# continue for calculating "hybrid circos plots"
#
# ###############
# #
# # get window based NeisD for sample of interest

# for (c.sample in sample.info$sample)
for (c.sample in c("EP","MAM","LdonLV9","CH32","CH34","GE","LEM3472","BUMM3","LRC.L740","SUKKAR2","BPK157A1","LRC.L53","BPK512A1","BPK612A1","ISS2429","ISS2426","ISS174","Inf055","Inf152","L60b","CL.SL","OVN3","LRC.L1311","LRC.L1312"))
{
  # dir.create(sample.info[sample==c.sample,]$groups, showWarnings = F)
  dir.create("dat", showWarnings = F)
  print(c.sample)
  load(file="win.count.RData")
  # all.chr.dists <-foreach (chr=c(paste0("0",1:2)), .combine=rbind) %do%
  all.chr.dists <-foreach (chr=c(paste0("0",1:9),10:36), .combine=rbind) %do%
  {
    print(chr)
    chr.dists <-foreach (win=1:win.count[[as.integer(chr)]], .combine=rbind) %do%
    {
      NeisD <-data.table(read.table(paste0("NeisD/chr_",chr,"/LmjF.",chr,"_win",win,"_NeisD.ind_NA.frac",NA.frac,".txt"), comment.char = ""))
      setnames(NeisD, colnames(NeisD), c(as.character(NeisD$V152),"sample", "group", "group.col"))
      aa<-NeisD[,.(get(c.sample),sample,group,group.col)]
      aa<-aa[order(V1)]
      aa<-aa[2:nrow(aa),]
      aa[,win:=win]
    }
    setnames(chr.dists,"V1","NeisD")
    chr.dists[,ordered.sample:=rep(1:150,max(chr.dists$win))]
    chr.dists[,chr:=as.integer(chr)]

    chr.dists
  }
  save(all.chr.dists, file=paste0("dat/all.chr.dists_",c.sample,".RData"))
}



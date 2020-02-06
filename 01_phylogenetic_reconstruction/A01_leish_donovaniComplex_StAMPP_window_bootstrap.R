

# get Neis distance matrix for windows along the chromosome

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/01_phylogenetic_reconstruction/window/")

library(data.table)
library(foreach)
library(StAMPP)
library(ggplot2)
# library(circlize)
#library(help = StAMPP)


#-----------variables--------------
load("../sample.info_NEWsubgroups.RData")
load(file="win.count.RData")
phylocol<-list("infantum_47"="0000CD", "donovani1_52"="FF0000", "donovani2a_19"="8B1C62", "donovani2b_7"="FF7F00", "donovani3_8"="FFBBFF", "other"="737373")
phylocolc<-list("infantum_47"=1, "donovani1_52"=2, "donovani2a_19"=3, "donovani2b_7"=4, "donovani3_8"=5, "other"=6)
NA.frac=0.2 # maximal fraction of missing data allowed at one SNP pos

#-----------functions--------------

source("../A01_leish_donovaniComplex_StAMPP_functions.R")


#-----------code--------------




###############
#
# A01_leish_donovaniComplex_StAMPP_window.R
# has to be run beforhand and save the neisD data for windows
#



print("start")
rep=1000
print(rep)

# ###############
# #
# # get window based NeisD for sample of interest

win.count # list with number of windows for each chromosome
win.sum <-sum(unlist(win.count)) # sum of all windows

# table for all windows
wins <-foreach (cc=1:36, .combine=rbind) %do%
# wins <-foreach (cc=1:2, .combine=rbind) %do%
{
  foreach (w=1:win.count[[as.numeric(cc)]], .combine=rbind) %do% 
  {
    c(cc,w)
  }
}
wins <-data.table(wins)
setnames(wins,colnames(wins),c("chr","win"))
wins[,id:=1:nrow(wins)]


# create a bootstrap replicate
#
bs.trees <-foreach (bs=1:rep) %do%
{
  print(paste("bs",bs))
  # pick windows for this bs replicate
  ids<-data.table(table(sort(sample(1:nrow(wins),nrow(wins),replace = T))))
  setnames(ids,colnames(ids),c("w.id","w.c")) # window id and window count
  count.wins<-matrix(0,ncol=151,nrow=151) # counts windows for which a secific comparison could not be done because of NA (will be substracted when mean over windows is calculated)
  
  print("calulate NeisD for all windows")
  NeisD.sum <-foreach (bb=1:nrow(ids), .combine="+") %do%
  #NeisD.sum <-foreach (bb=1:500, .combine="+") %do%
  {
    print(bb)
    w.id<-as.numeric(ids[bb,w.id]) # current window id
    w.c<-as.numeric(ids[bb,w.c]) # current window count (how often this window should be used for the current bs replicate)
    #   print(ids[bb,])
    xx <-wins[id==w.id] # get chr and win number for this current window id
    #   print(xx)
    
    if(xx$chr<10) {xx$chr<-paste0("0", xx$chr)}
    #   print(paste0("NeisD/chr_",xx$chr,"/LmjF.",xx$chr,"_win",xx$win,"_NeisD.ind_NA.frac",NA.frac,".txt"))
    NeisD <-data.table(read.table(paste0("NeisD/chr_",xx$chr,"/LmjF.",xx$chr,"_win",xx$win,"_NeisD.ind_NA.frac",NA.frac,".txt"), comment.char = ""))
    setnames(NeisD, colnames(NeisD), c(as.character(NeisD$V152),"sample", "group", "group.col"))
    
    foreach (cc=1:w.c, .combine="+") %do% # add up this matrix as often it is chosen
    {
      #     print(cc)
      #as.matrix(NeisD[1:3,1:5, with=F])
      #     print(as.matrix(NeisD[144:145,144:145, with=F]))
      aa <-as.matrix(NeisD[,1:151, with=F])
      
      # add count for respetive window comparison if window could not be added (due to NA)
      # add count for respetive window comparison if window could not be added (due to Inf)
      zz<-matrix(0,ncol=151,nrow=151)
      #   
      zz[is.na(aa)]<-1 # matrix indication all current NA positions marked with 1 else 0
      zz[is.infinite(aa)]<-1 # matrix indication all current Inf positions marked with 1 else 0
      count.wins <- count.wins+zz
      #
      aa[is.na(aa)] <- 0
      aa[is.infinite(aa)] <- 0
      
      aa
    }
  
  }
  
  div<-3198-count.wins
  NeisD.mean<-NeisD.sum/div
  
  tree.bs<-nj(NeisD.mean)
  
  tree.bs
}


dir.create("bs_trees", showWarnings = F)
save(bs.trees, file=paste0("bs_trees/bs.trees.",rep,".RData"))
load(file=paste0("bs_trees/bs.trees.",rep,".RData"))


# tree <-read.tree("../plots/LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_noCountry_newick.tx")
tree <-read.tree("../plots/LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_newick.txt")
# remove country names from tip labels
tips<-strsplit(tree$tip.label,split = "_")
tree$tip.label<-unlist(lapply(tips, function(x) ifelse(length(x)==2, x[1], paste0(x[1],"_", x[2])) ))

# tree$tip.label <-sapply(strsplit(tree$tip.label, split='_'),function(x) ifelse(length(x)==2,x[1],paste(x[1],x[2], sep = "_")))
clad <-round(prop.clades(tree, bs.trees, rooted = F)/rep*100)
clad[is.na(clad)]<-""
tree$node.label<- clad

pdf(paste0("../plots/LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_noCountry_newick_BS",rep,".pdf"),
    height=50, width=15)
plot(tree, show.node.label = T)
# nodelabels(clad)
dev.off()
write.tree(tree,file=paste0("../plots/LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_poly_NeisD.pop_Nj.boot_0_dphC_noCountry_newick_BS",rep,".txt"),
           digits = 4)






# get Neis distance matrix for windows along the chromosome

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/01_phylogenetic_reconstruction/")
dir.create("window", showWarnings=F)
setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/01_phylogenetic_reconstruction/window/")

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
# continue for "hybrid circos plots"



# ###############
# #
# # get window based NeisD for sample of interest

# for (c.sample in sample.info$sample)
# {
#   # dir.create(sample.info[sample==c.sample,]$groups, showWarnings = F)
#   dir.create("dat", showWarnings = F)
#   print(c.sample)
#   load(file="win.count.RData")
#   # all.chr.dists <-foreach (chr=c(paste0("0",1:9),10:36), .combine=rbind) %do%
#   all.chr.dists <-foreach (chr=c(paste0("0",1:2)), .combine=rbind) %do%
#   {
#     print(chr)
#     chr.dists <-foreach (win=1:win.count[[as.integer(chr)]], .combine=rbind) %do%
#     {
#       NeisD <-data.table(read.table(paste0("NeisD/chr_",chr,"/LmjF.",chr,"_win",win,"_NeisD.ind_NA.frac",NA.frac,".txt"), comment.char = ""))
#       setnames(NeisD, colnames(NeisD), c(as.character(NeisD$V152),"sample", "group", "group.col"))
#       aa<-NeisD[,.(get(c.sample),sample,group,group.col)]
#       aa<-aa[order(V1)]
#       aa<-aa[2:nrow(aa),]
#       aa[,win:=win]
#     }
#     setnames(chr.dists,"V1","NeisD")
#     chr.dists[,ordered.sample:=rep(1:150,max(chr.dists$win))]
#     chr.dists[,chr:=as.integer(chr)]
# 
#     chr.dists
#   }
#   save(all.chr.dists, file=paste0("dat/all.chr.dists_",c.sample,".RData"))
# }








####
#
# make circos plots for sample of interests for 10kb windows with
# 1. het count
# 2. heatmap with ordered NeisD coloured by groups
# 3. min and max Neisd present for the sample

load("../../A05_heterozygosities/het.wins.RData")
all.chr.dists
# all.chr.dists[,group.col:=as.character(group.col)]
# # update groups and colours
# for (gg in unique(sample.info[,groups]))
# {
#   all.chr.dists[sample %in% sample.info[groups==gg, sample], group:=gg]
#   all.chr.dists[sample %in% sample.info[groups==gg, sample], group.col:=unique(sample.info[groups==gg, group.col])]
# }
# save(all.chr.dists, file=paste0("dat/all.chr.dists_",c.sample,".RData"))

plot_circos <-function(c.sample)
{
  ploidy <-data.table(read.table(paste0(path,"/leish_donovaniComplex/A04_vcf_LD/somies_updated.txt"), head=T))
  # load("../../A05_heterozygosities/het.wins.RData")
  load(file=paste0("dat/all.chr.dists_",c.sample,".RData"))
  all.chr.dists[,chr:=as.character(chr)]
  all.chr.dists[chr %in% c(1,2,3,4,5,6,7,8,9), chr:=paste0("0",chr)]
  #
  
  dcir <- het.wins[,.(chr,w10kb,get(c.sample))]
  dcir[,chr:=chr]
  setnames(dcir, colnames(dcir)[3], c("y"))
  #
  png(paste0("Circos",c.sample,".png"),width=1200,height=1200)
  circos.clear()
  par(mar = c(1,1,1,1))#3,2,3,2)) 
  circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.1, start.degree=90)
  circos.initialize( factors=dcir$chr, x=dcir$w10kb )
  circos.trackPlotRegion(ylim = c(0, max(dcir$y)), factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
  {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.axis(
      h="top",                   # x axis on the inner or outer part of the track?
      labels=FALSE,              # show the labels of the axis?
      major.tick=TRUE,           # show ticks?
      labels.cex=0.5,            # labels size (higher=bigger)
      labels.font=2,             # labels font (1, 2, 3 , 4)
      direction="outside",       # ticks point to the outside or inside of the circle ?
      minor.ticks=4,             # Number of minor (=small) ticks
      major.tick.percentage=0.1, # The size of the ticks in percentage of the track height
      lwd=1.5                      # thickness of ticks and x axis.   
    )
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    aa = c(3.4, 0.5)
    if(theta < 90 || theta > 270)  aa = c(-2.4, 0.5)
    circos.text(x=mean(xlim), y=4, labels=name, facing = dd, cex=2,  adj = aa)
  })
  # add lines
  circos.trackLines(as.character(dcir$chr), dcir$w10kb, dcir$y, col="blue", lwd=2.5)
  #
  # add next tract
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.trackPlotRegion(ylim = c(0, 150), factors = dcir$chr, track.height=0.6,
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                         })
  # add samples
  for (os in 1:150) # os : ordered sample
  {
    aa<-all.chr.dists[ordered.sample==os]# & chr<10]
    circos.trackPoints(as.character(aa$chr), x=aa$win, y=rep(151-os,nrow(aa)), 
                       col=aa$group.col , pch = 20, cex = 0.1-(os*0.0005))
#                       col=paste0("#",unlist(phylocol[aa$group.col])) , pch = 20, cex = 0.1-(os*0.0005))
  }
  #
  # add track
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.trackPlotRegion(ylim = c(0, 1), factors = dcir$chr, track.height=0.1,
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                         })
  # add lines
  aa<-all.chr.dists[ordered.sample==1] # closest sample
  circos.trackLines(as.character(aa$chr), aa$win, aa$NeisD, col="darkgreen", lwd=2)
  aa<-all.chr.dists[ordered.sample==150] # furthest sample
  circos.trackLines(as.character(aa$chr), aa$win, aa$NeisD, col="orange", lwd=1.2)
  #
  # add somy tract
  circos.trackPlotRegion(ylim = c(0, 1), factors = as.character(aa$chr), track.height=0.1,
                         #panel.fun for each sector
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           
                           #plot main sector
                           somycol <- colorRampPalette(c("orange","yellow","green4","blue","red","pink"))(8)
                           circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                       col = somycol[ploidy[,get(c.sample)][i]], 
                                       border=somycol[ploidy[,get(c.sample)][i]])
                         })
  
  dev.off()
}

# n.plot.samp: indicates the number of nearest sample to be plotted inthe circos
plot_circos_part <-function(c.sample, n.plot.samp=50)
{
  ploidy <-data.table(read.table(paste0(path,"/leish_donovaniComplex/A04_vcf_LD/somies_updated.txt"), head=T))
  # load("../../A05_heterozygosities/het.wins.RData")
  load(file=paste0("dat/all.chr.dists_",c.sample,".RData"))
  all.chr.dists[,group.col:=as.character(group.col)]
  # update groups and colours
  for (gg in unique(sample.info[,groups]))
  {
    all.chr.dists[sample %in% sample.info[groups==gg, sample], group:=gg]
    all.chr.dists[sample %in% sample.info[groups==gg, sample], group.col:=unique(sample.info[groups==gg, group.col])]
  }
  save(all.chr.dists, file=paste0("dat/all.chr.dists_",c.sample,".RData"))
  all.chr.dists[,chr:=as.character(chr)]
  all.chr.dists[chr %in% c(1,2,3,4,5,6,7,8,9), chr:=paste0("0",chr)]
  #
  
  dcir <- het.wins[,.(chr,w10kb,get(c.sample))]
  dcir[,chr:=chr]
  setnames(dcir, colnames(dcir)[3], c("y"))
  summary(dcir[,y])
  #
  png(paste0("Circos_firstX/Circos",c.sample,"_first",n.plot.samp,"_maxhet",max(dcir$y),".png"),width=1200,height=1200)
  circos.clear()
  # par(mar = c(3,2,3,2)) #c(1,1,1,1))#
  par(mar = c(1,1,1,1))#
  circos.par(cell.padding = c(0, 0, 0, 0), "track.height" = 0.2, start.degree=90)
  circos.initialize( factors=dcir$chr, x=dcir$w10kb )
  circos.trackPlotRegion(ylim = c(0, max(dcir$y)), factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
  {
    name = get.cell.meta.data("sector.index")
    i = get.cell.meta.data("sector.numeric.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    circos.axis(
      h="top",                   # x axis on the inner or outer part of the track?
      labels=FALSE,              # show the labels of the axis?
      major.tick=TRUE,           # show ticks?
      labels.cex=0.5,            # labels size (higher=bigger)
      labels.font=2,             # labels font (1, 2, 3 , 4)
      direction="outside",       # ticks point to the outside or inside of the circle ?
      minor.ticks=4,             # Number of minor (=small) ticks
      major.tick.percentage=0.1, # The size of the ticks in percentage of the track height
      lwd=1.5                      # thickness of ticks and x axis.   
    )
    
    #text direction (dd) and adjusmtents (aa)
    theta = circlize(mean(xlim), 1.3)[1, 1] %% 360
    dd <- ifelse(theta < 90 || theta > 270, "clockwise", "reverse.clockwise")
    a1 <- ifelse(theta < 90 || theta > 270, c(-4.4), c(5.4))
    a2 <- 0.5; aa<-c(a1,a2) # 0.5 mean text willl be printed inthe middle of the cell
    circos.text(x=mean(xlim), y=0.2, labels=name, facing = dd, cex=2,  adj = aa)
    # circos.text(x=mean(xlim), y=4, labels=name, facing = dd, cex=2,  adj = aa)
  })
  # add lines
  circos.trackLines(as.character(dcir$chr), dcir$w10kb, dcir$y, col="blue", lwd=2.5)
  #
  # add next tract
  circos.par(cell.padding = c(0, 0, 0, 0))
#   n.plot.samp=50
  circos.trackPlotRegion(ylim = c(0, n.plot.samp), factors = dcir$chr, track.height=0.3,
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                         })
  # add samples
  for (os in 1:n.plot.samp) # os : ordered sample
  {
    aa<-all.chr.dists[ordered.sample==os]# & chr<10]
    circos.trackPoints(as.character(aa$chr), x=aa$win, y=rep((n.plot.samp+1)-os,nrow(aa)), 
                       col=aa$group.col , pch = 20, cex = 0.5-(os*0.0005))
    #                       col=paste0("#",unlist(phylocol[aa$group.col])) , pch = 20, cex = 0.1-(os*0.0005))
  }
  #
  # add track
  circos.par(cell.padding = c(0, 0, 0, 0))
  circos.trackPlotRegion(ylim = c(0, 1), factors = dcir$chr, track.height=0.15,
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                         })
  # add lines
  aa<-all.chr.dists[ordered.sample==1] # closest sample
  circos.trackLines(as.character(aa$chr), aa$win, aa$NeisD, col="darkgreen", lwd=2)
  aa<-all.chr.dists[ordered.sample==n.plot.samp]; aa[NeisD==Inf,NeisD:=NA] # furthest sample
  circos.trackLines(as.character(aa$chr), aa$win, aa$NeisD, col="orange", lwd=1.2)
  #
  # add somy tract
  circos.trackPlotRegion(ylim = c(0, 1), factors = as.character(aa$chr), track.height=0.15,
                         #panel.fun for each sector
                         panel.fun = function(x, y) {
                           #select details of current sector
                           name = get.cell.meta.data("sector.index")
                           i = get.cell.meta.data("sector.numeric.index")
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           
                           #plot main sector
                           somycol <- colorRampPalette(c("orange","yellow","green4","blue","red","pink"))(8)
                           circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2], ytop=ylim[2], 
                                       col = somycol[ploidy[,get(c.sample)][i]], 
                                       border=somycol[ploidy[,get(c.sample)][i]])
                         })
  
  dev.off()
}

dir.create("Circos_firstX", showWarnings = F)

for (i in c("EP")) {plot_circos(c.sample=i)}
for (i in c("MAM")) {plot_circos(c.sample=i)}
for (i in sample.info[groups=="donovani3_8"]$sample){plot_circos(c.sample=i)}
for (i in c("LdonLV9")) {plot_circos(c.sample=i)}
for (i in c("CH32","CH34")) {plot_circos(c.sample=i)}
for (i in c("ISS174","ISS2426","ISS2429")) {plot_circos(c.sample=i)}

for (i in c("MAM")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("MAM","EP")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("CH32","CH34")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("GE","LEM3472","BUMM3","LRC.L740","SUKKAR2")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("BPK157A1")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("LRC.L53")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in sample.info[groups=="Ldon3"]$sample){plot_circos_part(c.sample=i, n.plot.samp = 60)}

for (i in c("BPK512A1","BPK612A1")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("ISS2429","ISS2426","ISS174","Inf055","Inf152")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}

for (i in c("BPK512A1")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("L60b","CL.SL","OVN3")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}
for (i in c("LRC.L1311","LRC.L1312")) {plot_circos_part(c.sample=i, n.plot.samp = 60)}

##
# 
#
# for all samples in group donovani2a_19 (=donovani4)
foreach (samp = sample.info[groups %in% c("donovani2a_19")]$sample, .combine=cbind) %do%
{
  #samp="SUDAN1"
  load(file=paste0("dat/all.chr.dists_",samp,".RData"))
  
  aa <-all.chr.dists[chr==7 & group %in% c("infantum_47","donovani1_52")]
  
  ggplot(aa[order(win,ordered.sample)], aes(ordered.sample, NeisD, colour = factor(group))) + 
    geom_point() + facet_wrap(~win, ncol=5)


}


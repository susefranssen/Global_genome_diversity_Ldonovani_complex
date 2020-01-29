


library(data.table)
library(foreach)
library(ggplot2)
library(reshape2)
library(stringr)
library(gplots)
library(igraph)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/10_CNVs/")
dir.create("CNV_large", showWarnings = F)
setwd("CNV_large")

load("../../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")

#----------------------------------
# functions
#----------------------------------

# boVec: boolean vector indicating where an insert (or deletion, not both at the same time) is
# valVec: int vector indicating which copy number differenece the respective window has
# will output spos and epos of consecutive Trues
get_spos_epos<-function(boVec,valVec)
{
  res<-c(NA,NA,NA,NA,NA,NA)
  spos=NA
  epos=NA
  indelVals=c() # copy number vals of the windows in the current indel
  for (w in 1:length(boVec)) #go through all windows of this chromosome and sample
  {
    # middle pos in insert
    if (boVec[w]==T & !is.na(spos)) {indelVals=c(indelVals,valVec[w])}
    # first pos in insert
    if (boVec[w]==T & is.na(spos)) {spos=w; indelVals=valVec[w]}
    # last pos in insert
    if (boVec[w]==F & !is.na(spos)) {
      epos=w-1
      res<-rbind(res,c(spos, epos, epos-spos+1, mean(indelVals),  median(indelVals),  sd(indelVals)))
      spos=NA; epos=NA
    }
  }
  if (!is.na(spos)) {epos=w;  res<-rbind(res,c(spos, epos, epos-spos+1, mean(indelVals),  median(indelVals),  sd(indelVals)))}
  res<-data.table(res)
  setnames(res, colnames(res), c("spos","epos", "length","meanSomyDiff","medianSomyDiff","sdSomyDiff"))
  res[!is.na(epos)]
}

hclust.ave <- function(x) hclust(x, method="average")

# determine the list of unique indels for a given chromosome and indel type
# chrom chromosome number
# indel indel type either insertion or deletion
# win_dist number of windows the spos and epos can be different and still be counted as identical location
get_unique_indels<-function(chrom=27, indel="deletion", indels=indels, win_dist)
{
  if (indel=="deletion") {indel.sub<-indels[chr==chrom & deletion.pres==T]
  } else if (indel=="insertion") {indel.sub<-indels[chr==chrom & insert.pres==T]
  } else {stop()}
  
  
  if (nrow(indel.sub)==0) # no indels are present on this chr
  {
    zz<-data.table(rbind(c(indel, chrom, spos=NA, epos=NA, mean_dist=NA, samples.num=NA, groups.num=NA, samples.id=NA, groups.id=NA),
                         c(indel, chrom, spos=NA, epos=NA, mean_dist=NA, samples.num=NA, groups.num=NA, samples.id=NA, groups.id=NA)))
    setnames(zz,paste0("V",1:2), c("indel","chr"))
    zz
  } else 
  {
    indel.sub[, id:=paste0(name,":",spos)] # unique ID for each line
    edges<-c()
    if (nrow(indel.sub)>1) # only one indel in total
    {
      for (i in 1:(nrow(indel.sub)-1)) 
      {
        for (j in (i+1):nrow(indel.sub)) 
        {
          if ((abs(indel.sub[id==indel.sub[i,id], spos] - indel.sub[id==indel.sub[j,id], spos]) <=win_dist) & (abs(indel.sub[id==indel.sub[i,id], epos] - indel.sub[id==indel.sub[j,id], epos]) <=win_dist))
          {
            edges<-c(edges, unlist(c(indel.sub[i,id], indel.sub[j,id])))
          }
        }
      }
    } 
    
    
    if (is.null(edges)) # none of the indels are "identical"
    {
      zz<-foreach (l = 1:nrow(indel.sub), .combine=rbind) %do%
      {
        c(indel, chrom, spos=indel.sub[l,spos], epos=indel.sub[l,epos], mean_dist=1, samples.num=1, groups.num=1, samples.id=indel.sub[l,id], groups.id=as.character(indel.sub[l,group]))
      }
      if (nrow(indel.sub)==1) {zz<-rbind(zz,c(indel, chrom, spos=NA, epos=NA, mean_dist=NA, samples.num=NA, groups.num=NA, samples.id=NA, groups.id=NA))}
      zz<-data.table(zz)
      setnames(zz,paste0("V",1:2), c("indel","chr"))
      zz
    } else # identical indels exist
    {
      # make graph to find "identical" indels
      # all nodes that are unconnected, need to be added separately
      addnodes<-indel.sub[,id][!indel.sub[,id] %in% unique(edges)]
      edges<-c(edges,sort(rep(addnodes,times=2)))
      a<-make_graph(edges, directed = FALSE) # only includes nodes that are connected
      
      b<-decompose(a)
      
      clusters<-list()
      clusters.mean_dist<-c()
      for ( i in 1:length(b)) # go through all disconnected subgraphs
      {
        clusters[[i]]<-V(b[[i]])$name
        clusters.mean_dist<-c(clusters.mean_dist, round(mean_distance(b[[i]]),2))
        # print(mean_distance(b[[i]]))
        # print(V(b[[i]])$name)
        # cc<-cluster_edge_betweenness(b[[i]])
        # dendPlot(cc, mode="hclust")
        # plot(cc, b[[i]])
        # length(cc)
        # membership(cc)
      }
      
      res<-foreach (i = 1:length(clusters), .combine=rbind) %do%
      {
        spos<-median(indel.sub[id %in% clusters[[i]], spos])
        epos<-median(indel.sub[id %in% clusters[[i]], epos])
        mean_dist<-clusters.mean_dist[i]
        samples.num<-length(clusters[[i]])
        samples.id<-paste0(clusters[[i]], collapse = ",")
        groups.num<-length(unique(indel.sub[id %in% clusters[[i]], group]))
        groups.id<-paste0(unique(indel.sub[id %in% clusters[[i]], group]), collapse = ",")
        c(indel, chrom, spos, epos, mean_dist, samples.num, groups.num, samples.id, groups.id)
      }
      if (length(clusters)==1) {res<-rbind(res,c(indel, chrom, spos=NA, epos=NA, mean_dist=NA, samples.num=NA, groups.num=NA, samples.id=NA, groups.id=NA))}
      res<-data.table(res)
      setnames(res, colnames(res), c("indel","chr","spos","epos","mean_dist","samples.num","groups.num","samples.id","groups.id"))
      res
      
    }
  }
  
}


# get samples ids from the "unique" indel using the indel table all.indels.sub
get_samples_mean_dist<-function(all.indels.sub.row)
{
  samps<-strsplit(all.indels.sub.row[, samples.id], split=":")[[1]][rep(c(F,T),as.numeric(all.indels.sub.row[, samples.num])+1)]
  samps[1:length(samps)-1]
}

# inspect a chromosome for a specific sample
# prints the window based genome coverage for the respective chromosome and mark the identified indel regions
inspect_chr<- function(sample.c=c("CH32"), chrom=31, genome.cov=genome.cov, indels=indels)
{
  gg<-ggplot(data=genome.cov[chr==chrom & sample %in% sample.c], aes(x=win, y=cov.chrSomyNorm.diff, col=indel)) +
    geom_point(size=0.6)+
    labs(x=paste0("Window [",winsize/1000,"kb] on chromsome ",chrom), y="Copy number difference, normalised by somy") +
    theme_bw() +
    # geom_vline(xintercept = c(indels[sample==sample.c & chr==chrom, spos], indels[sample==sample.c & chr==chrom, epos])) +
    facet_wrap(~id,ncol=2)
  # theme(strip.background =element_rect(fill=c("red")))+
  # facet_grid(sample~chr)#, scales="free")
  # theme(legend.position="none") #+ 
  plot(gg)
}

# inspect a chromosome for a specific sample
# prints the window based genome coverage for the respective chromosome and mark the identified indel regions
# region vector of the window numbers
# cap specifies the maximum coverage to be plotted, larger values will be set to the specified cap value
inspect_chr_region <- function(sample.c=c("CH32"), chrom=31, genome.cov=genome.cov, region=c(), cap=40, ptsize=0.6, orient="v", ncol=2, LD1=F)
{
  genome.cov.sub<-genome.cov[chr %in% chrom & sample %in% sample.c]
  genome.cov.sub[cov.chrSomyNorm.diff>cap, cov.chrSomyNorm.diff:=cap]
  gg<-ggplot(data=genome.cov.sub, aes(x=win, y=cov.chrSomyNorm.diff, col=indel)) +
    geom_point(size=ptsize)+ 
    theme_bw() 
  if (length(chrom)==1)
  {
    gg<-gg+facet_wrap(~id,ncol=ncol) +
      labs(x=paste0("Window [",winsize/1000,"kb] on chromsome ",chrom), y="Copy number difference, normalised by somy")
  } else {
    gg<-gg+facet_grid(chr~id, scales="free") +
      labs(x=paste0("Window [",winsize/1000,"kb] on chromsome ",chrom[1],"-",chrom[36]), y="Copy number difference, normalised by somy")
  }
  if (length(unique(genome.cov.sub$indel))==3) {
    gg<-gg+scale_colour_manual(values=c("red","blue","gray"))
  }
  if (length(unique(genome.cov.sub$indel))==2 & "insertion" %in% unique(genome.cov.sub$indel)) {
    gg<-gg+scale_colour_manual(values=c("blue","gray"))
  }
  if (length(unique(genome.cov.sub$indel))==2 & "deletion" %in% unique(genome.cov.sub$indel)) {
    gg<-gg+scale_colour_manual(values=c("red","gray"))
  }
  if (!is.null(region)) {
    gg<-gg+geom_vline(xintercept = region, size=0.2, linetype="dotted") 
  }
  if (LD1==T)
  {
    # gg<-gg+geom_rect(aes(xmin=393, xmax=413, ymin=-3, ymax=6), fill="green", alpha=0.5, inherit.aes = FALSE)
    gg<-gg+annotate("rect", xmin=393, xmax=413, ymin=-3.3, ymax=6, alpha = 0.2, fill="green")
  }
  plot(gg)
}

#----------------------------------
# variables
#----------------------------------

# get chr sizes
chr.size<-data.table(read.table("../../09_group_popgen/TriTrypDB-38_LinfantumJPCM5_Genome.fasta.fai"))
chr.size[,V3:=NULL]; chr.size[,V4:=NULL]; chr.size[,V5:=NULL]
setnames(chr.size,colnames(chr.size),c("chrom","size"))
chr.size<-chr.size[order(chrom)]
chr.size<-chr.size[!1:nrow(chr.size) %in% grep("LinJ.00",chr.size$chrom)]
chr.size[,chrom:=as.numeric(unlist(lapply(strsplit(as.character(chr.size[,chrom]),split="J."), function(x) x[2])))]


#----------------------------------


# load somy and median chr coverage data
# this contains somy and median coverage data information by chromosome and sample and has to be generated beforehand using script: a10_raw_somies.r
# however, the results file is also provided
load(file="../somy_all.RData") 
setnames(all, c("chr","sample"), c("chrom","samp"))

#
all[,RDmedian_haploid:=round(RDmedian/Somy)]
all[,RDmedian_haploid_samp_med:=median(RDmedian_haploid), by=.(samp)]
haploid.med.cov<-unique(all[,.(samp,RDmedian_haploid_samp_med)])
summary(haploid.med.cov[,RDmedian_haploid_samp_med]) # all samples
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.0    20.0    26.0    29.0    33.5   144.0 
x_BD<-c("BD09", "BD12", "BD14", "BD15", "BD17", "BD21" , "BD24" ,"BD25" )
x_CUK<-haploid.med.cov[,samp][grep("CUK",haploid.med.cov[,samp])]
x_restImamura_R_P<-c("BHU1064A1","BHU1065A1","BHU1137A1","BHU1139A1","BHU220A1","BHU816A1","BHU824A1","BHU931A1","BPK029A1","BPK035A1","BPK067A1","BPK077A1","BPK156A1","BPK157A1","BPK164A1","BPK294A1","BPK406A1","BPK413A1","BPK471A1","BPK512A1","BPK562A1","BPK612A1","BPK623A1","BPK648A1","BPK649A1","LinJPCM5","LdonLV9")
x_Zack_ZHang<-c("X356WTV","X363SKWTI","X364SPWTII","X383WTI","AM560WTI","AM563WTI","CL.SL","OVN3")
x_all_other_studies<-c(x_BD,x_CUK,x_restImamura_R_P,x_Zack_ZHang)
length(x_all_other_studies)
# [1] 54
length(x_all_other_studies)+97
# [1] 151
haploid.med.cov<-haploid.med.cov[!samp %in% x_all_other_studies] #haploid median coverage of all samples new for this study
summary(haploid.med.cov[,RDmedian_haploid_samp_med])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   20.00   26.50   28.35   33.00   88.00
###<<<<<<



load("../../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")

winsize=5000

######################
#
# generation of the main data table "genome.cov" based on which all indels will be called
#
######################
#
# Note: this can only be run when the intermediate files as described in a14_gene_cov.sh have been generated before
# otherwise the results table that was generated in the following if statement is provided and can be loaded  with the load cmd after the if statement.
if (F)
{
  genome.cov<-foreach(sample.c=sample.info$sample, .combine=rbind) %do%
    # genome.cov<-foreach(sample.c=c("EP"), .combine=rbind) %do%
  {
    # sample.c="EP"
    
    somy<-all[samp==sample.c, .(chrom, RDmedian, Somy)]
    somy[, haploid.cov:=round(RDmedian/Somy)] # haploid coverage
    
    ss<-sample.c
    ss<-gsub("[.]", "-", ss)
    ss<-gsub("X", "", ss)
    
    dat<-data.table(read.table(paste0("../../data/CNV_cov_file_examples/linj.",sample.c,".sorted.markdup.realigned.PP.rmdup.genome.cov")))
    setnames(dat, colnames(dat),c("chr","pos","cov"))
    dat<-dat[chr !="LinJ.00"][order(chr)]
    
    dat[,win:=round(pos/winsize)]
    dat[,cov:=as.double(cov)]
    dat[,cov.median.win:=median(cov), by=.(chr,win)]
    dat[,cov.median.chr:=median(cov), by=chr]
    
    dat.win<-unique(dat[,.(chr,win,cov.median.win,cov.median.chr)])
    dat.win[,win.chr:=1:nrow(dat.win)]
    dat.win[, cov.median.win.somynorm:=cov.median.win/median(somy[,haploid.cov])]
    
    # pdf(paste0("genomecov_",sample.c,"_",winsize/1000,"kb.pdf"),width=10,height=2)
    # gg<-ggplot(dat.win, aes(win.chr, cov.median.win.somynorm, col=chr)) +
    #   geom_point(size=0.5) +
    #   labs(x=paste0("Window [",winsize/1000," kb]"), y="Median coverage") +
    #   scale_colour_manual(values=rep(c("red","blue"),18)) +
    #   theme(legend.position="none") +
    #   ylim(0.5,max(somy$Somy)+0.5)
    # print(gg)
    # dev.off()
    
    for (cc in c(paste0("LinJ.0",1:9),paste0("LinJ.",10:36)))
    {
      dat.win[chr==cc, ploidy:=somy[chrom==cc,Somy]]
    }
    dat.win[,ploidy:=as.factor(ploidy)]
    
    # pdf(paste0("genomecov_",sample.c,"_",winsize/1000,"kb_somy.pdf"),width=10,height=2)
    # gg<-ggplot(dat.win, aes(win.chr, cov.median.win.somynorm, col=ploidy)) +
    #   geom_point(size=0.5) +
    #   labs(x=paste0("Window [",winsize/1000," kb]"), y="Median coverage") +
    #   # scale_colour_manual(values=rep(c("red","blue"),18)) +
    #   ylim(0.5,max(somy$Somy)+0.5) +
    #   theme(legend.position = "none")
    # print(gg)
    # dev.off()
    
    dat.win[,sample:=sample.c]
    dat.win[,chr:=as.numeric(str_sub(as.character(dat.win[,chr]), -2,-1))]
    dat.win
  }
  
  for (sample.c in sample.info[,sample]) {genome.cov[sample==sample.c, group:=sample.info[sample==sample.c,groups]]}
  genome.cov<-genome.cov[order(group)]
  save(genome.cov,file=paste0("all.genome.cov.median_",winsize/1000,"kb.RData"))
} 
load(file=paste0("../all.genome.cov.median_",winsize/1000,"kb.RData"))
#
genome.cov[,cov.median.win.somynorm:=NULL]
genome.cov[,ploidy:=as.numeric(as.character(ploidy))]
# normalise the median window coverage by the median haploid coverage of the respective chromosome -> cov.chrSomyNorm
# note cov.median.win.somynorm was normalised by the haploid coverage of the respective chromosome using the haploid coverage calulated ACRoss chromosomes using the window based data from here
for(sample.c in sample.info$sample)
{
  # sample.c="EP"
  # sample.c="Peking"
  print(sample.c)
  # somy<-all[samp==sample.c, .(chrom, RDmedian, Somy)]
  # somy[, haploid.cov:=round(RDmedian/Somy)] # haploid coverage
  # somy[,chr:=as.numeric(unlist(lapply(strsplit(as.character(somy[,chrom]), split="J."), function(x) x[2])))]
  
  for(chromosome in 1:36)
  {
    medcov_across_win<-median(genome.cov[sample==sample.c & chr==chromosome, cov.median.win])
    genome.cov[sample==sample.c & chr==chromosome, cov.chrSomyNorm:=cov.median.win*ploidy/medcov_across_win]
  }
}
genome.cov[, cov.chrSomyNorm.diff:=cov.chrSomyNorm-ploidy]


######################
#
# looking at detailed coverage for the two examples from figure 4
#
######################
#
# Note: this can only be run when the intermediate files as described in a14_gene_cov.sh have been generated before
if (F)  
{
  ex<-"indel150"
  ex<-"indel215"
  ex<-"del234"
  ex<-"del220"
  
  if (ex=="indel150")
  {
    focal_samples=c("IECDR1","BD09","BD12","BD14","BD15","BD17","BD21","BD24","BD25","BHU41","Don201","Inf206","Don081","Don038","LRC.L51p","Inf007","Inf152","Inf001","Inf045","Inf004","Inf055","Cha001","BD27","LRC_L47")
    focal_chrom="LinJ.27"; reg_spos=190000; reg_epos=300000#; reg_spos=1800000; reg_epos=1850000
    # read in positional information on repeats
    repeats<-data.table(read.table("../../a24_repeats/nucmer_genomeTT38_RepeatsUbeda/genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ27.coords.bed"))
    # load(file=paste0("genome.cov.indel_",ex,".RData"))
    new_rep<-data.table(read.table("../../a24_repeats/nucmer_LinJ27regTT38_LinJ27regTT38/nucmer_self_LinJ27regTT38_LinJ27regTT38_mincl30_minm7_newrepeats_minlen200.txt", header = T))
    png(paste0("genome.cov.",ex,"_repeats.png"),width=2000,height=2000)
  } else if (ex=="indel215")
  {
    focal_samples=c("BHU41","LRC_L47")
    focal_chrom="LinJ.35"; reg_spos=1820000; reg_epos=1850000# reg_spos=1800000; reg_epos=1850000png(paste0("genome.cov.",ex,"_repeats.png"),width=1100,height=600)
    # read in positional information on repeats
    repeats<-data.table(read.table("../../a24_repeats/nucmer_genomeTT38_RepeatsUbeda/genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ35.coords.bed"))
    png(paste0("genome.cov.",ex,"_repeats.png"),width=1300,height=800)
  } else if (ex=="del234")
  {
    focal_samples=c("BUMM3")
    focal_chrom="LinJ.35"; reg_spos=490000; reg_epos=	530000 #515000
    # read in positional information on repeats
    repeats<-data.table(read.table("../../a24_repeats/nucmer_genomeTT38_RepeatsUbeda/genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ35.coords.bed"))
    png(paste0("genome.cov.",ex,"_repeats.png"),width=1300,height=450) 
  } else if (ex=="del220")
  {
    focal_samples=c("Malta33","SUKKAR2")
    focal_chrom="LinJ.35"; reg_spos=1925000; reg_epos=	1995000 #1965000
    # read in positional information on repeats
    repeats<-data.table(read.table("../../a24_repeats/nucmer_genomeTT38_RepeatsUbeda/genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ35.coords.bed"))
    png(paste0("genome.cov.",ex,"_repeats.png"),width=1300,height=800) 
  }
  setnames(repeats,colnames(repeats),c("chr","spos","epos","name","strand","sposx","eposx"))
  repeats[,spos:=spos+1]
  repeats[,epos:=epos+1]
  
  region<-repeats[chr==focal_chrom & reg_spos<=spos & reg_epos>=epos]
  genome.cov.indel<-foreach(sample.c=focal_samples, .combine=rbind) %do%
  {
    # sample.c="BHU41"
    # sample.c="BHU41"
    
    somy<-all[samp==sample.c, .(chrom, RDmedian, Somy)]
    somy[, haploid.cov:=round(RDmedian/Somy)] # haploid coverage
    
    # files in /lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38
    ss<-sample.c
    ss<-gsub("[.]", "-", ss)
    ss<-gsub("X", "", ss)
    
    if (server)
    {
      file=paste0("/lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/genomecov/linj.",ss,".sorted.markdup.realigned.PP.rmdup.genome.cov")
      print(file)
    } else {
      file=paste0("/lustre/scratch118/infgen/team133/sf18/leish_donovaniComplex/A14_gene_cov_annotv38/genomecov/linj.",ss,".sorted.markdup.realigned.PP.rmdup.genome.cov")
      print(file)
      file=paste0("/Users/sf18/work/tmp/linj.",ss,".sorted.markdup.realigned.PP.rmdup.genome.cov")
    }
    #
    if (!file.exists(file))
    {
      print(paste("file", file, " does not exist!"))
    } else {
      print(paste("sample", file))
    }
    
    dat<-data.table(read.table(file))
    setnames(dat, colnames(dat),c("chr","pos","cov"))
    dat<-dat[chr==focal_chrom]
    # 
    dat[, ploidy:=somy[chrom==focal_chrom,Somy]]
    dat[,ploidy:=as.factor(ploidy)]
    #
    dat[,win:=round(pos/winsize)+1]
    dat[,cov:=as.double(cov)]
    dat[,cov.median.win:=median(cov), by=.(chr,win)]
    dat[,cov.median.chr:=median(cov), by=chr]
    #
    dat[, cov.median.win.somynorm:=cov.median.win/median(somy[,haploid.cov])]
    dat[, cov.somynorm:=cov/median(somy[,haploid.cov])]
    
    dat[,sample:=sample.c]
    dat[,chr:=as.numeric(str_sub(as.character(dat[,chr]), -2,-1))]
    
    dat[,cov.somynorm.round:=round(cov.somynorm,0)]
    
    dat
  }
  
  for (samp in unique(genome.cov.indel$sample))
  {
    genome.cov.indel[sample==samp,group:=sample.info[sample==samp,groups]]
  }
  genome.cov.indel[,name:=paste0(group,"_",sample)]
  
  color.palette  <- c(colorRampPalette(c("yellow","orange","green4","blue","red"))(6) ,c("purple","pink"))
  #
  # pdf(paste0("genome.cov.",ex,"_repeats.pdf"),width=11,height=6)
  # png(paste0("genome.cov.",ex,"_repeats.png"),width=1100,height=600)
  # png(paste0("genome.cov.",ex,"_repeats.png"),width=1300,height=800)
  gg<-ggplot(genome.cov.indel[pos>=reg_spos & pos<=reg_epos], aes(pos/1000, cov.somynorm, col=as.factor(cov.somynorm.round))) +
    # geom_point(size=0.5) +
    geom_point(size=1) +
    labs(x=paste0("Chromosome position [k bp]"), y="Coverage by base, normalised by somy", size=2) +
    scale_colour_manual(values=c("gray",
                                 color.palette,
                                 rep(color.palette[8],sum(sort(unique(genome.cov.indel$cov.somynorm.round))>8)),
                                 "black"),
                        name="Somy\nEquivalent") +
    theme_bw() + theme(text = element_text(size=20)) +
    facet_wrap(.~name, ncol=1) 
  xx<-unique(region[,.(spos,epos)])
  gg<-gg+geom_segment(aes(x = xx$spos[1]/1000, y = 0, xend = xx$epos[1]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[2]/1000, y = 0, xend = xx$epos[2]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[3]/1000, y = 0, xend = xx$epos[3]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[4]/1000, y = 0, xend = xx$epos[4]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[5]/1000, y = 0, xend = xx$epos[5]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[6]/1000, y = 0, xend = xx$epos[6]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[7]/1000, y = 0, xend = xx$epos[7]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[8]/1000, y = 0, xend = xx$epos[8]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[9]/1000, y = 0, xend = xx$epos[9]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[10]/1000, y = 0, xend = xx$epos[10]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[11]/1000, y = 0, xend = xx$epos[11]/1000, yend = 0), size=2, color="black")
  gg<-gg+geom_segment(aes(x = xx$spos[12]/1000, y = 0, xend = xx$epos[12]/1000, yend = 0), size=2, color="black")
  # for (i in 1:nrow(xx))
  # {
  #   gg<-gg+geom_segment(aes(x = xx$spos[i]/1000, y = 0, xend = xx$epos[i]/1000, yend = 0))
  # }
  if (ex=="indel150") # ad newly found repeats
  {
    gg<-gg+geom_segment(aes(x = new_rep$spos1[1]/1000, y = 0, xend = new_rep$epos1[1]/1000, yend = 0), size=2, color="blue")
    gg<-gg+geom_segment(aes(x = new_rep$spos1[2]/1000, y = 0, xend = new_rep$epos1[2]/1000, yend = 0), size=2, color="blue")
    gg<-gg+geom_segment(aes(x = new_rep$spos1[3]/1000, y = 0, xend = new_rep$epos1[3]/1000, yend = 0), size=2, color="blue")
    gg<-gg+geom_segment(aes(x = new_rep$spos1[4]/1000, y = 0, xend = new_rep$epos1[4]/1000, yend = 0), size=2, color="blue")
  }
  if (ex %in% c("del220")) # ad newly found repeats
  {
    gg<-gg+annotate("rect", xmin=1967, xmax=2000, ymin=-0.2, ymax=11, alpha = 0.15, fill="green")
  }
  print(gg)
  dev.off()
  
  save(genome.cov.indel,file=paste0("genome.cov.indel_",ex,".RData"))
}




######################
#
# large indel identification
#
######################
#
thresh.win=0.5 # minimum  copy number difference of window (coverage normalised to haploid coverage units) to chromosome copy number
thresh.indel=0.9 # minimum median difference (coverage normalised to haploid coverage units) of called indel
#
indels<-data.table(matrix(rep(NA,10), ncol=10))
setnames(indels,colnames(indels),c("spos","epos","length","meanSomyDiff","medianSomyDiff","sdSomyDiff","lengthKB","sample","chr","ploidy"))
for (sample.c in sample.info$sample)
{
  print(sample.c)
  for (chrom in 1:36)
  {
    aa<-genome.cov[chr==chrom & sample==sample.c,]

    # find insertions
    if (sum(aa$cov.chrSomyNorm.diff>thresh.win) > 0) {
      ins<-get_spos_epos(aa$cov.chrSomyNorm.diff>thresh.win, aa$cov.chrSomyNorm.diff)
      ins[,lengthKB:=length*winsize/1000]
      ins[,sample:=sample.c]
      ins[,chr:=chrom]
      ins[,ploidy:=unique(aa$ploidy)]
      indels<-rbind(indels,ins)
    }
    # find deletions
    if (sum(aa$cov.chrSomyNorm.diff<(-thresh.win)) > 0) {
      dels<-get_spos_epos(aa$cov.chrSomyNorm.diff<(-thresh.win), aa$cov.chrSomyNorm.diff)
      dels[,lengthKB:=length*winsize/1000]
      dels[,sample:=sample.c]
      dels[,chr:=chrom]
      dels[,ploidy:=unique(aa$ploidy)]
      indels<-rbind(indels,dels)
    }
    indels
  }
}
indels<-indels[lengthKB>=25 & abs(medianSomyDiff)>thresh.indel]
for (gg in unique(sample.info[,groups]))
{
  indels[sample %in% sample.info[groups==gg, sample], group:=gg]
  indels[sample %in% sample.info[groups==gg, sample], colour:=unique(sample.info[groups==gg, group.col])]
}
indels<-indels[,.(group, sample,chr, ploidy, spos, epos, length, lengthKB, meanSomyDiff, medianSomyDiff, sdSomyDiff)]
indels[, lengthKB:=as.double(lengthKB)]
indels[,lengthKB.mean.chr:=mean(lengthKB), by=.(chr)] # length
indels[,insert.sum.samp.chr:=sum(medianSomyDiff>0), by=.(sample,chr)] # medianSomyDiff --> number of inserts per sample and chr
indels[,deletion.sum.samp.chr:=sum(medianSomyDiff<0), by=.(sample,chr)] # medianSomyDiff --> number of deletions per sample and chr
write.table(indels, file=paste0("indels_thresh.win",thresh.win,"_thresh.indel",thresh.indel,".txt"), quote = F, row.names = F)
indels<-data.table(read.table(file=paste0("indels_thresh.win",thresh.win,"_thresh.indel",thresh.indel,".txt"), header = T))
indels[, lengthKB:=as.double(lengthKB)]
indels<-indels[order(group)]
indels[,name:=paste0(group,":",sample)]

# basic stats
hist(indels$lengthKB, nclass=200)
summary(indels$lengthKB)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 25.00   25.00   30.00   40.74   40.00  675.00
quantile(indels$lengthKB, probs=c(0.25,0.5,0.75,0.8,0.9,0.95,0.99))
# 25% 50% 75% 80% 90% 95% 99% 
# 25  30  40  45  60  65 235  
nrow(indels[lengthKB>=45]) / nrow(indels)
# [1] 0.2255319
nrow(indels[lengthKB>=50]) / nrow(indels)
# [1] 0.1670213
nrow(indels[lengthKB>=100]) / nrow(indels)
# [1] 0.03297872
sort(table(indels[lengthKB>=100, chr]))
# 4 20 26 36 11 18 23 29 31 35 
# 1  1  1  1  2  2  2  2  4 15 
#
######################
# indels stats
#
######################
# typical lengths of indels
pdf("indels_meanLengthPerChr.pdf")
ggplot(indels, aes(chr, lengthKB, group=chr)) +
  geom_boxplot()
dev.off()
#
######################
# have some chr indels more often than others? (does not differentiate how many are on a specific chr)
#
######################
# across all samples
indels[, insert.pres:=medianSomyDiff>0] # insert presence info
indels[, deletion.pres:=medianSomyDiff<0] # deletion presence info
indels.mod<-unique(indels[,.(group,sample,chr,insert.pres, deletion.pres)])
indels.mod[, insert.pres.mean:=sum(insert.pres)/151, by=.(chr)] # insert presence frac across samples
indels.mod[, deletion.pres.mean:=sum(deletion.pres)/151, by=.(chr)] # deletion presence frac across samples
indels.red<-unique(indels.mod[,.(chr,insert.pres.mean, deletion.pres.mean)])[order(chr)]
indels.red.melt<-melt(indels.red, id.vars = c("chr")); setnames(indels.red.melt, c("variable","value"), c("Indeltype","presence_fraction"))
pdf("indels_indelpresfrac_Acrosssamples.pdf")
ggplot(indels.red.melt, aes(as.factor(chr),presence_fraction,fill=Indeltype)) +
  geom_bar(stat = "identity", position = "dodge") + theme_bw() + labs(x="Chromsome", y="Fraction of samples with Indels") +
  scale_fill_manual(values=c("darkgreen", "red"), 
                    labels=c("insertion", "deletion"))
dev.off()
#
######################
# adding type of insertion information to genome.cov
#
######################
# 
genome.cov[, indel:="none"]
for (ii in 1:nrow(indels))
{
  if (indels[ii, medianSomyDiff>0])
  {
    genome.cov[sample==indels[ii,sample] & chr==indels[ii,chr] & win>=indels[ii,spos] & win<=indels[ii,epos], indel:="insertion" ]
  } else if (indels[ii, medianSomyDiff]<0)
  {
    genome.cov[sample==indels[ii,sample] & chr==indels[ii,chr] & win>=indels[ii,spos] & win<=indels[ii,epos], indel:= "deletion"]
  }
}
genome.cov$sample <-factor(genome.cov$sample, levels=unique(genome.cov$sample[order(genome.cov$group)]))
genome.cov[, id:=paste0(group,"_",sample)]



######################
#
# plot chr 27 with many dels
#
######################
# --> 22 of all the sample share the same deletion acros both species + two additional samples that do not have the deletion
pdf("indelPlots_del_nodel_chr27.pdf",width=8, height=8)
sample.c<-c(c(as.character(indels[chr==27 & deletion.pres==T, sample])), "BD27","LRC_L47")
inspect_chr_region(sample.c=sample.c, chrom=27, genome.cov=genome.cov, cap=2)
dev.off()



######################
#
# identify all indel that are "identical" by position
#
######################
# 
# they have a quality value associated,
# the mean distance between all nodes in a cluster, best if ==1 all indels in a cluster have the same start and stop +-2 windows, 
# i.e. every 2 nodes have and edge between them
# if larger indels borders are less conserved and it could be potentially different clusters, i.e. indels
all.indels<-foreach (chri=1:36, .combine = rbind) %do%
{
  print(chri)
  # calculating "identical" indels conditionion on max 1 or 2 windows accuracy of start and endpoint of the indel
  foreach (win_dist=1:2, .combine = rbind) %do%
  {
    zz1<-get_unique_indels(chrom=chri, indel="insertion", indels, win_dist)
    zz2<-get_unique_indels(chrom=chri, indel="deletion", indels, win_dist)
    zz<-rbind(zz1,zz2)
    zz[, win.dist:=win_dist]
    zz
  }
}
all.indels<-all.indels[!is.na(spos)]
# mark "telomere amplifications" --> indels starrting within 15 kb (<=3 5kb wins) from the chr border
for (cc in 1:36) {all.indels[chr==cc, chrsize:=chr.size[chrom==cc,size]]}
all.indels[, telomere_ampl:="no"]
all.indels[spos<=3, telomere_ampl:="start"]
all.indels[chrsize/winsize-epos<=3, telomere_ampl:="end"]
all.indels[, length:=epos-spos+1]
#
all.indels[,Linf.gr:=grepl("Linf",all.indels$groups.id)]
all.indels[,Ldon.gr:=grepl("Ldon",all.indels$groups.id)]
all.indels<-all.indels[order(-win.dist)]
all.indels[,id:=1:nrow(all.indels)]
all.indels[, spos_genome:=ifelse(spos-1>=0,(spos-1)*winsize+1,0+1)]
all.indels[, epos_genome:=epos*winsize]
all.indels[, lengthKB:=epos_genome-spos_genome]
write.table(all.indels, file="all.indels.txt")
#
all.indels<-data.table(read.table(file="all.indels.txt"))
all.indels[,chr:=as.numeric(chr)]
all.indels[,spos:=as.numeric(spos)]
all.indels[,epos:=as.numeric(epos)]
all.indels[,mean_dist:=as.numeric(mean_dist)]
all.indels[,samples.num:=as.numeric(samples.num)]
all.indels[,groups.num:=as.numeric(groups.num)]
all.indels[,samples.id:=as.character(samples.id)]
all.indels[, lengthKB:=as.numeric(lengthKB)]
# use win_dist=2
all.indels.sub<-all.indels[win.dist==2]
write.table(all.indels.sub, file="all.indels.sub.txt", row.names = F, quote = F)
table(all.indels.sub[,.(mean_dist,groups.num)])
table(all.indels.sub[,.(groups.num, samples.num)])

######################
# summary of shared indels between samples
#
######################
# 
pdf("all.indels.sub_groups_chr.pdf", width=8,height=5)
heatmap.2(table(all.indels.sub[,.(groups.num, samples.num)]), col= c("lightgrey",colorRampPalette(c("yellow","red"))(200)),
          key.xlab = "INDEL presence", key.title ="",  keysize = 1.4,
          trace='none', Colv = F, 
          cellnote=table(all.indels.sub[,.(groups.num, samples.num)]), notecol="black",
          dendrogram="none", Rowv=FALSE,
          xlab = "Presence in number of samples") 
dev.off()

table(all.indels.sub[,groups.num]) # number of unique indels in few and several groups
# 1   2   3   4   5   6   7   8 
# 154  49  16   7   9   7   2   1
table(all.indels.sub[,indel])
# deletion insertion 
# 62       183
table(all.indels.sub[,.(indel,chr)])
# chr
# indel        1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36
# deletion   0  1  1  3  3  0  0  1  3  1  2  3  1  1  2  0  1  1  4  0  1  1  3  1  1  1  1  1  1  1  9  0  2  1  7  3
# insertion  2  5  2  3  2  5  5  5  4  1  4  8  6  2  2  4  3  7  6  6  6  3 10  2  4  5  1  3  6  5 16  5  7  5 17  6
sort(table(all.indels.sub[,.(chr)]))
# indels.somySd.unique<-cbind(indels.somySd, table(all.indels.sub[,.(chr)]))

table(all.indels.sub[Linf.gr==T & Ldon.gr==T,chr]) # chromosomes with indels that are present in Ldon and Linf samples (as called permissively need to be confirmed that indels are indeed the same across samples)
# 2  3  4  5  6  7  8  9 10 11 12 13 14 15 17 19 20 21 22 23 24 25 26 27 28 29 30 31 33 34 35 
# 2  1  2  2  1  1  2  2  1  2  5  2  1  1  2  4  1  1  1  2  2  1  2  1  1  2  1  9  3  2  9

# number of indels shared between groups
sum(table(all.indels.sub[Linf.gr==T & Ldon.gr==T,.(groups.num,samples.num)]))
# [1] 69 ->69/245=28%
# save regions of all indels shared between both species
yy<-all.indels.sub[Linf.gr==T & Ldon.gr==T,.(id,chr,spos_genome,epos_genome,indel, samples.num, groups.num,length)]
yy[,chr:=ifelse(chr<10,paste0("LinJ.0",chr),paste0("LinJ.",chr))]
yy[,spos.bed:=spos_genome-1]
yy[,epos.bed:=epos_genome-1]
yy[,info:=paste0("id",id,"_",indel,"_samp.num",samples.num,"_group.num",groups.num,"_lengthKB",length*winsize/1000)]
write.table(yy[,.(chr, spos.bed, epos.bed, info)], file="indels_in_both_species.bed", quote=F, row.names = F, sep = "\t")
#
######################
# summary of shared indels between samples that are present in both species
#
######################
# 
pdf("all.indels.sub_groups_chr_bothspecies.pdf", width=8,height=5)
heatmap.2(table(all.indels.sub[Linf.gr==T & Ldon.gr==T,.(groups.num,samples.num)]), col= c("lightgrey",colorRampPalette(c("yellow","red"))(200)),
          key.xlab = "INDEL presence", key.title ="",  keysize = 1.4,
          trace='none', Colv = F, 
          cellnote=table(all.indels.sub[Linf.gr==T & Ldon.gr==T,.(groups.num,samples.num)]), notecol="black",
          dendrogram="none", Rowv=FALSE,
          xlab = "Presence in number of samples") 
dev.off()



######################
#
# plot all indels that are "shared" between both species
#
######################
# 
dir.create("indels_both_species", showWarnings = F)
for (i in all.indels.sub[Linf.gr==T & Ldon.gr==T,][,id])
{
  print(i)
  row<-all.indels.sub[id==i]
  sample.c<-get_samples_mean_dist(row)
  pdf(paste0("indels_both_species/id_",as.character(row[,id]),"_chr",row[,chr],"_type",row[,indel],"_nGroups",row[,groups.num],".pdf"),
      width=8,height=10)
  inspect_chr_region(sample.c, chrom=row[,chr], genome.cov=genome.cov, region=c(row[,spos],row[,epos]), cap=8, ncol=2)#8)
  dev.off()
  
  if (row[,indel]=="insertion") {
    xxx<-indels[sample %in% sample.c & chr==row[,chr] & insert.pres==T & abs(spos-row[,spos])<=4 & abs(epos-row[,epos])<=4]
  } else {
    xxx<-indels[sample %in% sample.c & chr==row[,chr] & insert.pres==F & abs(spos-row[,spos])<=4 & abs(epos-row[,epos])<=4]
  }
  xxx[,abs_copy:=ploidy+medianSomyDiff]
  write.table(xxx, file=paste0("indels_both_species/id_",as.character(row[,id]),"_chr",row[,chr],"_type",row[,indel],"_nGroups",row[,groups.num],".txt"),
              quote=F, row.names = F)
}



######################
#
# test if the "indels" are at the chr bordres and therefore might represent rather telomeric amplifications
#
######################
# 
table(all.indels.sub$telomere_ampl)
# end    no start 
# 64   127    54 
round(table(all.indels.sub$telomere_ampl)/nrow(all.indels.sub)*100)
# end    no start 
# 26 51 22
table(all.indels.sub[,.(indel,telomere_ampl)])
# telomere_ampl
# indel       end no start
# deletion   14 34    14
# insertion  50 93    40
#
# number of samples with the most "telomeric indels"
teloindels.id<-all.indels.sub[telomere_ampl!="no", samples.id]
teloindels.samp<-unlist(lapply(strsplit(strsplit(paste(teloindels.id, collapse = ","), split=",")[[1]], split=":"),function(x) x[2]))
# number of samples with any end indeld amplifications
length(unique(teloindels.samp))
# [1] 100 ->100/151=66%
# samples with more or equal to 10 end indels and number
sort(table(teloindels.samp), decreasing = T)[sort(table(teloindels.samp), decreasing = T)>=10]
teloindels.samp
# HN167    HN336 BPK406A1 LRC.L51p BPK164A1 BPK029A1  LdonLV9     BD09     CH35 
# 54       40       33       27       14       12       12       10       10 
length(sort(table(teloindels.samp), decreasing = T)[sort(table(teloindels.samp), decreasing = T)>=10])
# [1] 9 ->9/151=6%
#
# samples with most telomeric amplifications
sort(table(teloindels.samp), decreasing = T)[1:10]
# teloindels.samp
# HN167    HN336 BPK406A1 LRC.L51p BPK164A1 BPK029A1  LdonLV9     BD09     CH35     CUK8 
# 54       40       33       27       14       12       12       10       10        9
# number of samples with only a single telomeric amplification
sum(sort(table(teloindels.samp), decreasing = T)==1)
# [1] 43 ->43/151=28%
pdf("TOPtelomereAmpl_samples.pdf", width=6, height=11)
inspect_chr_region(sample.c=c("HN167","HN336"), chrom=1:36, genome.cov=genome.cov, region=c(), cap=8, ptsize = 0.2)
dev.off()

inspect_chr_region(sample.c=c("BPK077A1"), chrom=1:36, genome.cov=genome.cov, region=c(), cap=8, ptsize = 0.2)




######################
#
# plot largest indels on chr 35
#
######################

all.indels.sub[chr==35 & epos-spos+1>=20]
pos<-c()
# for (i in 1:nrow(all.indels.sub[chr==35 & epos-spos+1>=20]))
# {
#   pdf(paste0("Chr35_largestIndels",i,".pdf"))
#   indel.id<-all.indels.sub[chr==35 & epos-spos+1>=20][i,samples.id]
#   spos<-all.indels.sub[chr==35 & epos-spos+1>=20][i,spos]
#   epos<-all.indels.sub[chr==35 & epos-spos+1>=20][i,epos]
#   if (grepl(",",indel.id)) {pos<-c(pos, spos,epos)}
#   indel.samp<-unlist(lapply(strsplit(strsplit(paste(indel.id, collapse = ","), split=",")[[1]], split=":"),function(x) x[2]))
#   inspect_chr_region(sample.c=indel.samp, chrom=35, genome.cov=genome.cov, region=c(spos,epos), cap=4, ptsize = 0.2)
#   dev.off()
# }
indel.samp<-unlist(lapply(strsplit(strsplit(paste(all.indels.sub[chr==35 & epos-spos+1>=20,samples.id],collapse = ","),split=",")[[1]],split=":"),function(x) x[2]))
pdf(paste0("Chr35_largestIndels_all.pdf"))
inspect_chr_region(sample.c=indel.samp, chrom=35, genome.cov=genome.cov, region=pos, cap=6, ptsize = 0.2, orient="h", LD1=T)
dev.off()
write.table(all.indels.sub[chr==35 & epos-spos+1>=20], file="Chr35_largestIndels_all.txt", quote=F, row.names = F)

# save(genome.cov, file=paste0("genome.cov_w",winsize,".RData"))


library(data.table)
library(foreach)
# library(matrixStats)
library(ggplot2)

setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/08_genetically_distinct_subclones")



######################
#
# get allele freqs

sample.all<-c("Peking","D_2","STRAIN_B","STRAIN_A","RACOON_DOG","SKIN","DOG_STRAIN","IECDR1","BD09","BD12","BD14",
              "BD15","BD17","BD21","BD22","BD24","BD25","BD27","Ldon282cl2","BPK029A1","BPK035A1","BPK067A1","BPK077A1",
              "BPK156A1","BPK157A1","BPK164A1","BPK282I9","BPK294A1","BPK406A1","BPK413A1","BPK471A1","BPK512A1","BPK562A1",
              "BPK612A1","BPK623A1","BPK648A1","BPK649A1","L60b","OVN3","CL-SL","LRC-L51p","Chowd5","STL2-78","STL2-79","DD8",
              "Nandi","AG83","BHU41","BHU220A1","BHU1062-4","Don201","BHU816A1","BHU824A1","BHU931A1","BHU1064A1","BHU1065A1",
              "BHU1137A1","BHU1139A1","LRC-L1311","LRC-L1313","Inf206","BUMM3","SUKKAR2","Don081","LdonLV9","GEBRE1","Don038",
              "356WTV","AM560WTI","363SKWTI","383WTI","364SPWTII","AM563WTI","LRC-L53","LRC-L57","MRC74","NLB-323","LRC-L445",
              "LRC-L61","1S","GILANI","GE","38-UMK","45-UMK","452BM","597-2","597LN","762L","LEM3472","1026-8","855-9","SUDAN1",
              "Inf007","LRC-L699","LRC-L740","LRC-L1275","TH4","TH5","TH6","NT10","NT16","LRC-L1312","LRC-L1296","LRC_L1303",
              "CH32","CH33","CH34","CH35","CH36","EP","Inf152","CUK2","CUK3","CUK4","CUK5","CUK6","CUK7","CUK8","CUK9","CUK10",
              "CUK11","CUK12","ISS174","ISS2420","ISS2426","ISS2429","ISS2508","Malta33","BUCK","Inf001","LRC_L47","RM1",
              "LEM1985","LPN114","LEM3278","Inf045","Inf004","Inf055","BCN83","BCN87","LinJPCM5","IMT260","IMT373cl1","ITMAP26",
              "Cha001","MAM","ARL","WC","WR285","HN167","HN336")
load("../05_heterozygosities/hets_sample_heterozygosities.RData")
sample.het<-names(sort(hets[hets>0.004], decreasing = T))
sample.het<-sub("X","", sample.het)
sample.het<-sub(".","-", sample.het, fixed=T)




############################################################

# boxplot for sample distr

sample_SNP_pos=T # if T only posiions that wre called as heterozygous for the respecitve sample are used otherwise all SNP across the entire data set
#
load("../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")
dir.create("Reffreq_plots", showWarnings = F)
min.SNP=100 # the one that was used as a threshold for the peak calculation before
#
# load(file=paste0("peakx.samples_sampleSNPpos",sample_SNP_pos,".RData"))
het.index<-data.table(cbind(sort(sample.het), 1:46)) # sample index to load the correct estimated peak data in (red lines)
het.index[,V2:=as.numeric(V2)]
#
# load("/Volumes/sf18/leish_donovaniComplex/A10_raw_somies/somy_all.RData") #all
# setnames(all, "sample","samp")
somy1<-data.table(read.table("../04_aneuploidy/somies_updated.txt", header = T)); somy1[,chr:=1:.N]
somy<-melt(somy1, measure.vars = sample.info$sample)
setnames(somy,colnames(somy),c("chr","samp","somy"))
somy[chr<10,chrom:=paste0("LinJ.0",chr)]
somy[chr>=10,chrom:=paste0("LinJ.",chr)]
#
somy.col<-colorRampPalette(c("yellow","orange","green4","blue","red"))(6)
#

# foreach (sample = sort(sample.het), .combine=rbind) %do%

# for (sample in c("BPK512A1","LRC-L1311","LRC-L1312","LRC-L1313")) # for these one no peak estimates exist as there are tooo few SNPs
# for (sample in c("MAM","EP","CH32","CH34","BPK157A1","GILANI","Malta33","364SPWTII", "SUKKAR2","BUMM3","GE",
#                  "LEM3472","LRC-L740","LRC-L53","ISS2429","ISS2426","ISS174","Inf055","Inf152", "BPK512A1"))
# for (sample in paste0("CUK",3:12))
# for (sample in c("BPK157A1","Inf152","ISS174","ISS2426","ISS2429","LRC-L53","MAM","BUMM3","LRC-L740","Malta33","SUKKAR2"))
# for (sample in c("LRC-L53","MAM","BUMM3","LRC-L740","Malta33","SUKKAR2"))
for (sample in c("EP","MAM"))
{
  # sample="MAM"
  # sample="EP"

  file=paste0("../data/read_coverage/linj.",sample,".sorted.markdup.realigned.sort.SNPs.pileup.sync")
  
  dat<-data.table(read.table(file))
  setnames(dat,colnames(dat),c("chr","pos","ref","sync")) # sync entry A:T:C:G:N:del
  sample<-gsub(sample, pattern = "364SPWTII",replacement = "X364SPWTII")
  sample<-gsub(sample, pattern = "LRC-L740",replacement = "LRC.L740")
  sample<-gsub(sample, pattern = "LRC-L53",replacement = "LRC.L53")
  
  ss<-strsplit(as.character(dat$sync), split=":")
  dat[,a:=as.integer(sapply(ss, function(x) x[1]))]
  dat[,t:=as.integer(sapply(ss, function(x) x[2]))]
  dat[,c:=as.integer(sapply(ss, function(x) x[3]))]
  dat[,g:=as.integer(sapply(ss, function(x) x[4]))]
  dat[,ref:=as.character(ref)]
  dat[,chr:=as.character(chr)]
  dat[,sync:=NULL]
  dat[,sum:=as.integer(sapply(ss, function(x) sum(as.integer(x[1:4]))))]
  # t2<-top2(dat[,4:7,with=F])
  dat[ref=="a",ref.count:=a]
  dat[ref=="t",ref.count:=t]
  dat[ref=="c",ref.count:=c]
  dat[ref=="g",ref.count:=g]
  dat[,ref.freq:=round(ref.count/sum,2)]
  dat[,ref.freq.raw:=ref.count/sum]
  
  dat<-dat[!ref.freq %in% c(0,1),]
  dat[,samp:=sample]
  dat[,chrom:=sapply( dat$chr, function(x) strsplit(x, split = "J.")[[1]][2])]
  dat[,cc:=1]
  dat[,nSNP.chr:=sum(cc), by=.(chr)]
  # dat[,highnSNP.chr:=nSNP.chr>=min.SNP]
  
  for (ch in c(paste0("LinJ.0",1:9),paste0("LinJ.",10:36)))
  {
    sss<-gsub(sample, pattern = "LRC-L131",replacement = "LRC.L131")
    sss<-gsub(sample, pattern = "CL-SL",replacement = "CL.SL")
    ss<-somy[chrom==ch & samp==sss, somy]
    dat[chr==ch,somy:=ss] 
  }
  levels(dat$somy) <-c(1,2,3,4,5,6,7,8)
  
  xx<-dat
  
  if (sample_SNP_pos)
  {
    sss<-gsub(sample, pattern = "LRC-L131",replacement = "LRC.L131")
    sss<-gsub(sample, pattern = "CL-SL",replacement = "CL.SL")
    SNP.pos<-data.table(read.table(paste0("../05_heterozygosities/het_pos_per-sample/SNP_pos_",sss,".txt"), header = T))
    setnames(SNP.pos, colnames(SNP.pos), c("chr","pos"))
    print(paste(sample,dim(xx)))
    xx<-merge(xx,SNP.pos, by = c("chr","pos")) # reduce data set to positions called as SNPs for the respective sample
    print(paste(sample,dim(xx)))
  } else {
    xx<-xx[ref.freq>0.05 & ref.freq<0.95] # remove extreme freq loci
  }
  
  pdf(paste0("Reffreq_plots/Reffreq_",sample,"_minSNP",min.SNP,"_SNPpos",sample_SNP_pos,".pdf"), width=2, height=9)
  # pdf(paste0("Reffreq_plots/Reffreq_",sample,"_minSNP",min.SNP,"_SNPpos",sample_SNP_pos,"_extra.pdf"), width=5, height=15)#width=2, height=9)
  # pdf(paste0("Reffreq_plots/Reffreq_",sample,"_minSNP",min.SNP,".pdf"), width=5, height=9)
  gg<-ggplot(data=xx, aes(x=ref.freq, fill=as.factor(somy))) +
    geom_histogram(binwidth = 0.025) + xlim(0,1) +
    facet_grid(chrom~samp, scales = "free_y") +
    scale_fill_manual(values=somy.col[sort(unique(dat$somy))]) +
    # scale_fill_manual(values=col.use) +
    guides(fill=FALSE) +
    xlab("Reference frequency") #+
    # geom_vline(xintercept = c(1/3,2/3), colour="blue", size=0.1) +#, linetype="dotted") +
    # geom_vline(xintercept = c(1/4,3/4), colour="black", size=0.1) + #, linetype="dashed")
    # geom_vline(aes(xintercept=vlines),  xx, col="darkred", size=0.4) #+ theme_bw()
    plot(gg)
  dev.off()

  # for samples that present clone mixtures only
  # get the histograms for high qual SNPs vs lower qual SNPs
  if (sample %in%  c("BPK157A1","Inf152","ISS174","ISS2426","ISS2429","LRC.L53","MAM","BUMM3","LRC.L740","Malta33","SUKKAR2"))
  {
    library(patchwork)
        
    xx[,chr_pos:=paste0(chr,"_",pos)]
    highqual_SNPs<-data.table(read.table(paste0("/Volumes/sf18/leish_donovaniComplex/A00_SNP_stats/SNP_qual_bias_check/highConf_hetSNPs",sample,".txt"), header = T))
    highqual_SNPs[,chr_pos:=paste0(chr,"_",pos)]
    xx[,highQual:=ifelse(chr_pos %in% highqual_SNPs$chr_pos, T, F)]
    
    pdf(paste0("Reffreq_plots/Reffreq_",sample,"_minSNP",min.SNP,"_SNPpos",sample_SNP_pos,"_qualComp.pdf"), 
        width=6, height=13)#width=2, height=9)
    g1<-ggplot(data=xx[highQual==T], aes(x=ref.freq, fill=as.factor(somy))) +
      geom_histogram(binwidth = 0.025) + xlim(0,1) +
      facet_grid(chrom~samp, scales = "free_y") + 
      facet_grid(chrom~samp, scales = "free_y") + 
      scale_fill_manual(values=somy.col[sort(unique(dat$somy))]) +
      guides(fill=FALSE) + 
      xlab("Reference frequency")  + ggtitle("High quality SNPs") + theme(axis.text = element_text(size = 6))
    g2<-ggplot(data=xx[highQual==F], aes(x=ref.freq, fill=as.factor(somy))) +
      geom_histogram(binwidth = 0.025) + xlim(0,1) +
      facet_grid(chrom~samp, scales = "free_y") + 
      facet_grid(chrom~samp, scales = "free_y") + 
      scale_fill_manual(values=somy.col[sort(unique(dat$somy))]) +
      guides(fill=FALSE) +
      xlab("Reference frequency") + ggtitle("Remaining SNPs") + theme(axis.text = element_text(size = 6))
    plot(g1 +g2)
    dev.off()
    
    sample
    # [1] "BPK157A1"
    table(xx[,highQual])
    # FALSE  TRUE 
    # 3139   216 
  }
  
}



##############################
#
# phasing based on unequal somies
#
############################################################

# phasing of subset of samples - get frequencies

# NA.frac<-0.2
load(file = paste0("../01_phylogenetic_reconstruction/l.gtStampp_cov0.2_poly.RData")) # l.gtStampp.freq
for (chrom in 1:36) 
{
  l.gtStampp.freq[[chrom]]<-data.table(l.gtStampp.freq[[chrom]])  
}

dir.create("phase", showWarnings = F)

# args <- commandArgs(TRUE)
# input.sample <- args[1] 

# for (samp in sample.het)
for (samp in c("MAM","ISS174","ISS2426","ISS2429","BPK157A1","LRC.L53","Inf152"))
{
  #   samp="MAM"
  print(samp)

  file=paste0("../data/read_coverage/linj.",sample,".sorted.markdup.realigned.sort.SNPs.pileup.sync") # note: bam files have beenquality filtered
  
  if (substr(samp,1,1) %in% c(1,3,4,5,7,8)) {sampR<-paste0("X",samp)} else {sampR<-samp}
  sampR<-sub("-",".", sampR, fixed=T)
  
  if (T)
  {
    dat<-data.table(read.table(file))
    setnames(dat,colnames(dat),c("chr","pos","ref","sync")) # sync entry A:T:C:G:N:del
    
    ss<-strsplit(as.character(dat$sync), split=":")
    dat[,a:=as.integer(sapply(ss, function(x) x[1]))]
    dat[,t:=as.integer(sapply(ss, function(x) x[2]))]
    dat[,c:=as.integer(sapply(ss, function(x) x[3]))]
    dat[,g:=as.integer(sapply(ss, function(x) x[4]))]
    dat[,ref:=as.character(ref)]
    dat[,chr:=as.character(chr)]
    dat[,sync:=NULL]
    dat[,sum:=as.integer(sapply(ss, function(x) sum(as.integer(x[1:4]))))]
    # t2<-top2(dat[,4:7,with=F])
    dat[ref=="a",ref.count:=a]
    dat[ref=="t",ref.count:=t]
    dat[ref=="c",ref.count:=c]
    dat[ref=="g",ref.count:=g]
    dat[,ref.freq:=round(ref.count/sum,2)]
    dat[,ref.freq.raw:=ref.count/sum]
    
    dat<-dat[!ref.freq %in% c(0,1),] # --> excluding homozygous loci
    dat[,chr:=as.integer(lapply(strsplit(chr, split="J."), function(x) x[2]))]
    
    all.het.i=list()
      
    # dat with only the called het genotypes for plotting the allele frequency profiles  
    dat.het.pos <-foreach (chrom = 1:36, .combine=rbind) %do% # get all truely heterozygous loci (based on the called genotype and not just if the frequency is different from 0, 1)
    {
      print(chrom)
      pos_list <- strsplit(colnames(l.gtStampp.freq[[chrom]])[6:ncol(l.gtStampp.freq[[chrom]])], split="_")
      posi <- as.integer(unlist(lapply(pos_list,function(x) x[2]))) # all positions in the data table
      print(length(posi))
      
      index<-(1:length(posi))[posi %in% dat[chr == chrom, pos]] # SNP columns
      index=index+5 # meta data columns
      het.i=c() # index of heterozygous loci
      hom.i=c() # index of loci that are called as homozygotes even though their frequency is not 0 or 1
      het.pos=c() # pos of heterozygous loci
      hom.pos=c() # pos of loci that are called as homozygotes even though their frequency is not 0 or 1
      vals<-c()
      for (cc in index)
      {
        val<-l.gtStampp.freq[[chrom]][Sample==sampR,cc,with=F]
        vals<-c(vals,val)
        #       print(val)
        if (!is.na(val))
        {
          if (val==1)
          {
            hom.i<-c(hom.i,cc)
            hom.pos<-c(hom.pos,posi[cc-5])
          } 
          else if (val >0 & val <1)
          {
            het.i<-c(het.i,cc)
            het.pos<-c(het.pos,posi[cc-5])
          }
        }
        #       l.gtStampp.freq[[chrom]][Sample=="MAM",het.i,with=F]
        #       l.gtStampp.freq[[chrom]][Sample=="MAM",hom.i,with=F]
      }
      all.het.i[[chrom]]<-het.i
      
      dat[chr==chrom & pos %in% het.pos,]
    }

    write.table(dat.het.pos, file=paste0("phase/dat.het.pos_",sampR,".txt"), quote=FALSE, row.names=F)
    save(all.het.i, file=paste0("all.het.i_",samp,".RData"))
  } else {
    dat.het.pos <-data.table(read.table(file=paste0("phase/dat.het.pos_",sampR,".txt"),header = T))
    load(file=paste0("all.het.i_",samp,".RData"))
  }

}


############################################################

# phasing of subset of samples - mixed samples

l.gtStampp.freq.haplo<-list()

# for (samp in c("MAM","ISS174","ISS2426","ISS2429","BPK157A1","LRC-L53","Inf152"))
for (samp in c("MAM","ISS174","ISS2426","ISS2429","BPK157A1","LRC-L53","Inf152"))
{
  #   samp="MAM"
  print(samp)
  
  file=paste0("../data/read_coverage/linj.",sample,".sorted.markdup.realigned.sort.SNPs.pileup.sync") # note: bam files have beenquality filtered
  if (substr(samp,1,1) %in% c(1,3,4,5,7,8)) {sampR<-paste0("X",samp)} else {sampR<-samp}
  sampR<-sub("-",".", sampR, fixed=T)
  
  dat.het.pos <-data.table(read.table(file=paste0("phase/dat.het.pos_",sampR,".txt"),header = T))
  load(file=paste0("all.het.i_",samp,".RData"))  
  
  for (chrom in 1:36)
  {
    print(paste(sampR, chrom))
    ref<-l.gtStampp.freq[[chrom]][Sample==sampR,]; ref[,Sample:=paste0(Sample,"_high")]
    alt<-l.gtStampp.freq[[chrom]][Sample==sampR,]; alt[,Sample:=paste0(Sample,"_low")]
    
    for (i in all.het.i[[chrom]] )
    {
#       print(paste(sampR,chrom,i))
      col.n<-colnames(l.gtStampp.freq[[chrom]])[i] # colname current
      pos.c<-as.integer(strsplit(col.n, split="_")[[1]][2])
      
      if (dat.het.pos[chr==chrom & pos==pos.c,ref.freq.raw]  > 0.5)
      {
        ref[,(col.n):=1]
        alt[,(col.n):=0]
      } else {
        ref[,(col.n):=0]
        alt[,(col.n):=1]        
      }
    }
    
    if (length(l.gtStampp.freq.haplo) <36) # in the first round (once for each chr), the tables foreach chr do not yet exist
    {
      l.gtStampp.freq.haplo[[chrom]] <-rbind( ref, alt)
    } else {
      l.gtStampp.freq.haplo[[chrom]] <-rbind( l.gtStampp.freq.haplo[[chrom]], ref, alt)
    }
    l.gtStampp.freq.haplo[[chrom]][,ploidy:=1]
    
#     dat.het.pos[chr==chrom,] 
  }
}


for (chrom in 1:36) # save new haplos to the genotype table
{
  l.gtStampp.freq[[chrom]] <-rbind( l.gtStampp.freq[[chrom]], l.gtStampp.freq.haplo[[chrom]])
}


dir.create("haplo", showWarnings = F)
save(l.gtStampp.freq, file = "haplo/l.gtStampp.freq0.2_poly_addhaplo.RData") # l.gtStampp.freq
load(file = "haplo/l.gtStampp.freq0.2_poly_addhaplo.RData") # l.gtStampp.freq with phased haplotypes





############################################################

# phasing of subset of samples - only Ldon complex data - for triploid chromosome samples
#AAAAAAAA

# load StAMPP formatted data with the added whole genome haplotypes in it
load(file = "haplo/l.gtStampp.freq0.2_poly_addhaplo.RData") # l.gtStampp.freq

# check with chr are triploid in het samples --> for phasing
load("../05_heterozygosities/hets_sample_heterozygosities.RData")
sample.het<-names(sort(hets[hets>0.004], decreasing = T))
ploidy<-data.table(read.table("../04_aneuploidy/somies_updated.txt", header = T))
ploidy.t<-data.table(melt(ploidy, measure.vars = colnames(ploidy))) #
setnames(ploidy.t,colnames(ploidy.t),c("sample","somy"))
ploidy.t[,chr:=rep(1:36,151)]
ploidy.t[,sample:=as.character(sample)]
triploid<-table(ploidy.t[sample %in% sample.het & somy==3, .(sample,chr)])
sort(colSums(triploid), decreasing = T)
sort(ploidy.t[sample %in% sample.het & chr==26 & somy==3]$sample)
sort(ploidy.t[sample %in% sample.het & chr==23 & somy==3]$sample)
sort(ploidy.t[sample %in% sample.het & chr==8 & somy==3]$sample)
sort(ploidy.t[sample %in% sample.het & chr==9 & somy==3]$sample)

# initialize empty list for new haplotypes to be added
l.gtStampp.freq.haplo<-list() # new haplotypes to be added
l.gtStampp.freq.haplo[[37]] <-NA
l.gtStampp.freq.haplo[[37]] <-NULL

# save both possible haplotype into the new list
triploid.chr<-as.integer(names(sort(colSums(triploid), decreasing = T)[1:4]))
for (chrom in triploid.chr) # the 4 most commonlt triploid chromosomes
{
  print(chrom) 
  samples<-sort(ploidy.t[sample %in% sample.het & chr==chrom & somy==3]$sample) # all samples that are triploid for this chromosome
  print(length(samples))  

  for (samp in samples)
  {
    #   samp="MAM"
    print(samp)
    
    file=paste0("../data/read_coverage/linj.",sample,".sorted.markdup.realigned.sort.SNPs.pileup.sync") # note: bam files have beenquality filtered
    if (substr(samp,1,1) %in% c("X")) {sampR<-substr(samp,2,10)} else {sampR<-samp}
    sampR<-sub(".","-", sampR, fixed=T)
    
    dat.het.pos <-data.table(read.table(file=paste0("phase/dat.het.pos_",samp,".txt"),header = T))
    all.het.i<-NULL
    load(file=paste0("all.het.i_",sampR,".RData"))  
    print(length(all.het.i))

    print(paste(sampR, chrom))
    ref<-l.gtStampp.freq[[chrom]][Sample==samp,]; ref[,Sample:=paste0(Sample,"_high")]
    alt<-l.gtStampp.freq[[chrom]][Sample==samp,]; alt[,Sample:=paste0(Sample,"_low")]
    
    for (i in all.het.i[[chrom]] )
    {
      #       print(paste(sampR,chrom,i))
      col.n<-colnames(l.gtStampp.freq[[chrom]])[i] # colname current, e.g. "LinJ.19_2947"
      pos.c<-as.integer(strsplit(col.n, split="_")[[1]][2])
      
      if (dat.het.pos[chr==chrom & pos==pos.c,ref.freq.raw]  > 0.5)
      {
        ref[,(col.n):=1]
        alt[,(col.n):=0]
      } else {
        ref[,(col.n):=0]
        alt[,(col.n):=1]        
      }
    }
    
    if (is.null(l.gtStampp.freq.haplo[[chrom]])) # in the first round (once for each chr), the tables foreach chr do not yet exist
    {
      l.gtStampp.freq.haplo[[chrom]] <-rbind( ref, alt)
    } else {
      l.gtStampp.freq.haplo[[chrom]] <-rbind( l.gtStampp.freq.haplo[[chrom]], ref, alt)
    }
    l.gtStampp.freq.haplo[[chrom]][,ploidy:=1]
      
  }
}

# add haplotypes from new list into the genotype list
for (chrom in 1:36) # save new haplos to the genotype table
{
  print(chrom)
  if(!is.null(l.gtStampp.freq.haplo[[chrom]]))
  {
    print("in")
    l.gtStampp.freq[[chrom]] <-rbind( l.gtStampp.freq[[chrom]], l.gtStampp.freq.haplo[[chrom]])
  }
}

save(l.gtStampp.freq, file = paste0("haplo/l.gtStampp.freq0.2_poly_addhaplo_phasetriploid.RData")) # l.gtStampp.freq



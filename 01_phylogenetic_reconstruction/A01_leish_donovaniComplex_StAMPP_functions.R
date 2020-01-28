
#### libraries #######################

library(ggplot2)



#### variables #######################




#### FUNCTIONS #######################

# format SNP data to StAMPP format
format_StAMPP<-function(dat, iloci=NULL)
{
  # first 4 columns information for StAMPP
  Sample=names(dat)[9:ncol(dat)]
  Pop=rep(NA,length(Sample))
  Format=rep("BiA")
  
  if (is.null(iloci)) iloci=c(1:nrow(dat))
  # processing per individual
  ## iloci=1:100
  ## ind=dat[iloci,9:ncol(dat),with=F]$CL13a
  # unique(str_count(as.character(dat[iloci,9:ncol(dat),with=F]$CL2), ":"))
  individuals=foreach (ind=dat[iloci,9:ncol(dat),with=F]) %do% 
  {
    #print(as.character(ind))
    #somy=unique(as.integer(matrix(unlist(strsplit(as.character(ind), split=":")), ncol=6, byrow=T)[,4]))
    somy=unique(as.integer(sapply(strsplit(as.character(ind), split=":"),function(x) x[4])))
    
    if (length(somy>1)) 
      #gt=matrix(unlist(strsplit(as.character(ind), split=":")), ncol=6, byrow=T)[,1]
      gt=sapply(strsplit(as.character(ind), split=":"),function(x) x[1])
    else (stop("Somy has to be identical within one individual!"))
    # encode alleles with AB
    gt=gsub("1","B", gsub("0","A",gsub("/","",gt)))
    gt=gsub( paste0(rep("[.]",somy), collapse=""), "-9",gt)
  }

  # getting somy for each individual
  somy=foreach (ind=dat[iloci,9:ncol(dat),with=F]) %do% 
  {
    #unique(as.integer(matrix(unlist(strsplit(as.character(ind), split=":")), ncol=6, byrow=T)[,4]))
    unique(as.integer(sapply(strsplit(as.character(ind), split=":"),function(x) x[4])))
  }
  somy=unlist(somy)
  Ploidy=somy
  
  # set locus matrix
  loci=matrix(unlist(individuals), nrow=length(Sample), byrow=T)
  lname=paste(dat$CHR,dat$POS.CHR, sep="_")
  colnames(loci)<-lname[iloci]
  
  data.table(cbind(Sample,Pop,Ploidy,Format,loci))
}



#' Reduces gtStampp formated matrix to sufficiently covered SNPs across all input samples
#'
#' @param gtStampp.freq Stampp formated frequency matrix
#' @param NA.frac the maximal fractions of Nas across samples allowed
#' @param col.offset number of colums with header / general information (default=5; as expected of gtStampp matrix)
#' @return gtStampp.freq gtStampp input matrix with columns with to many NAs across samples removed
reduce.cov<-function(gtStampp.freq, NA.frac, coll.offset=5)
{
  gtStampp.freq<-data.table(gtStampp.freq)
  n.samp=nrow(gtStampp.freq)
  n.SNPs=ncol(gtStampp.freq)-coll.offset
  take <-foreach(snp=1:n.SNPs, .combine=c) %do%
  {
    if(sum(is.na(gtStampp.freq[,snp+coll.offset]))/n.samp <=NA.frac) {T} else F
  }
  gtStampp.freq[,c(rep(T,coll.offset),take),with=F]
}

#' Reduces gtStampp formated matrix to polymorphic SNPs across all input samples
#'
#' @param gtStampp.freq Stampp formated frequency matrix
#' @param col.offset number of colums with header / general information (default=5; as expected of gtStampp matrix)
#' @return gtStampp.freq gtStampp input matrix with columns with to many NAs across samples removed
reduce.poly<-function(gtStampp.freq, coll.offset=5)
{
  gtStampp.freq<-data.table(gtStampp.freq)
  n.SNPs=ncol(gtStampp.freq)-coll.offset
  take <-foreach(snp=1:n.SNPs, .combine=c) %do%
  {
    if (nrow(unique(na.omit(gtStampp.freq[,snp+coll.offset,with=F]))) > 1) T else F
  }
  gtStampp.freq[,c(rep(T,coll.offset),take),with=F]
}





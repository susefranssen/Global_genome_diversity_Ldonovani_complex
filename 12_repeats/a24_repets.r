


library(data.table)
library(ggplot2)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/12_repeats/")

#--------------------------------
# functions
#--------------------------------

plot_nucmer_self <-function(dir, name, genomic_spos, lines, image)
{
  dat<-data.table(read.table(paste0(dir,name)))
  setnames(dat,colnames(dat),c("spos1","epos1","spos2","epos2","len1","len2","identity","seqlen1","seqlen2","chr1","chr2"))
  setnames(repeats, colnames(repeats)[1:3], c("chr","spos","epos"))
  #
  dat[,id:=paste0(spos1,"_",spos2,"_",identity)]
  xx<-unique(dat[,.(spos1,epos1,spos2,epos2,identity,id)])
  x1<-xx[,.(spos1,spos2,identity,id)]
  x2<-xx[,.(epos1,epos2,identity,id)]
  setnames(x1, colnames(x1), c("pos1","pos2","identity","id"))
  setnames(x2, colnames(x2), c("pos1","pos2","identity","id"))
  yy<-rbind(x1,x2)
  yy[,pos1:=pos1+genomic_spos]
  yy[,pos2:=pos2+genomic_spos]
  #
  rep_sub<-repeats[spos>=genomic_spos & epos<=max(yy$pos1,yy$pos2)]
  #
  if (image=="pdf")
  {
    pdf(paste0(dir,"nucmer_self_",gsub("_tab.coords","",name),".pdf"),width=6,height=5)
  } else {
    png(paste0(dir,"nucmer_self_",gsub("_tab.coords","",name),".png"),width=600,height=500)
  }
  
  gg<-ggplot(yy, aes(x=pos1, y=pos2, group=id, color=identity)) +
    geom_rect(aes(xmin=lines[1], xmax=lines[2], ymin=lines[1], ymax=lines[2]), color="lightgray", fill="lightgray", alpha=0.9) +
    geom_line() + 
    # geom_vline(xintercept = lines[1], color="black", lwd=0.2) +
    # geom_hline(yintercept = lines[1], color="black", lwd=0.2) +
    # geom_vline(xintercept = lines[2], color="black", lwd=0.2) +
    # geom_hline(yintercept = lines[2], color="black", lwd=0.2) +
    scale_color_gradient(low="red", high="blue") +
    theme_bw() +
    labs(x="Genomic position [bp]", y="Genomic position [bp]")
  if (nrow(rep_sub)==2)
  {
    gg<-gg +
      geom_segment(aes(x=rep_sub[1,spos],xend=rep_sub[1,epos],y=genomic_spos,yend=genomic_spos), color="chartreuse4", lwd=1.5) +
      geom_segment(aes(x=rep_sub[2,spos],xend=rep_sub[2,epos],y=genomic_spos,yend=genomic_spos), color="chartreuse4", lwd=1.5) +
      geom_vline(xintercept = rep_sub[1,spos], color="chartreuse4", lwd=0.2) +
      geom_vline(xintercept = rep_sub[2,spos], color="chartreuse4", lwd=0.2)
  } else {
    print(nrow(rep_sub))
  }
  if (name=="LinJ27regTT38_LinJ27regTT38_mincl30_minm7_tab.coords") 
  {
    new_rep<-dat[identity>=80 & len1<5000 & len1>=200]
    new_rep[,spos1:=spos1+genomic_spos]
    new_rep[,epos1:=epos1+genomic_spos]
    new_rep[,spos2:=spos2+genomic_spos]
    new_rep[,epos2:=epos2+genomic_spos]
    write.table(new_rep, paste0(dir,"nucmer_self_",gsub("_tab.coords","",name),"_newrepeats_minlen200.txt"),quote = F, row.names = F)
    new_rep<-dat[identity>=80 & len1<5000 & len1>=100]
    new_rep[,spos1:=spos1+genomic_spos]
    new_rep[,epos1:=epos1+genomic_spos]
    new_rep[,spos2:=spos2+genomic_spos]
    new_rep[,epos2:=epos2+genomic_spos]
    write.table(new_rep, paste0(dir,"nucmer_self_",gsub("_tab.coords","",name),"_newrepeats_minlen100.txt"),quote = F, row.names = F)
    gg<-gg +
      geom_segment(aes(x=new_rep[1,spos1],xend=new_rep[1,epos1],y=genomic_spos,yend=genomic_spos), color="chartreuse3", lwd=1.5) +
      geom_segment(aes(x=new_rep[2,spos1],xend=new_rep[2,epos1],y=genomic_spos,yend=genomic_spos), color="chartreuse3", lwd=1.5) +
      geom_segment(aes(x=new_rep[3,spos1],xend=new_rep[3,epos1],y=genomic_spos,yend=genomic_spos), color="chartreuse3", lwd=1.5) +
      geom_segment(aes(x=new_rep[4,spos1],xend=new_rep[4,epos1],y=genomic_spos,yend=genomic_spos), color="chartreuse3", lwd=1.5) +
      geom_vline(xintercept = new_rep[1,spos1], color="chartreuse3", lwd=0.2) +
      geom_vline(xintercept = new_rep[2,spos1], color="chartreuse3", lwd=0.2) +
      geom_vline(xintercept = new_rep[3,spos1], color="chartreuse3", lwd=0.2) +
      geom_vline(xintercept = new_rep[4,spos1], color="chartreuse3", lwd=0.2) +
      geom_vline(xintercept = new_rep[5,spos1], color="chartreuse3", lwd=0.2) +
      geom_vline(xintercept = new_rep[6,spos1], color="chartreuse3", lwd=0.2)
  }
  print(gg)
  dev.off()
  
  # dat
  rep_sub
}



#--------------------------------



#####################
# 
# 
# 
#####################
# 
name="LinJ27regTT38_LinJ27regTT38_mincl30_minm7_tab.coords"
dir=""
genomic_spos=190000
lines=c(206000,264000)
repeats<-data.table(read.table("genomeTT9.38_RepeatsUbeda_tab_TT38.LinJ27.coords.bed"))
#
plot_nucmer_self(dir, name, genomic_spos, lines, "pdf")
plot_nucmer_self(dir, name, genomic_spos, lines, "png")







datrep<-data.table(read.table("nucmer_LinJ27regTT38_LinJ27regTT38/LinJ27regTT38_LinJ27regTT38.repeats", header = T))
datrep[,Start2:=as.character(Start2)]
datrep[, rev:=F]
datrep[grepl("r", Start2), rev:=T]
datrep[rev==T, Start2:=gsub("r","",Start2)]
datrep[, Start2:=as.numeric(Start2)]
datrep[,gs1:=Start1+190000]
datrep[,gs2:=Start2+190000]
datrep[,xx:=gs2-gs1+1]


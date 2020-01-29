



library(data.table)
library(foreach)
library(ggplot2)
library(reshape2)
library(topGO)
library(Rgraphviz)




#############
#
# GSE - gene set enrichment analysis
# this script will focus on the result os weightFish using p=0.05
#
#############
# comp="Linf_vs_Ldon"
# comp="Linf1star_43"
# comp="CH.Linf.nohy_3"
# comp="CUK.Linf_11"
# comp="Ldon1star_44"
# comp="Ldon2_7"
# comp="Ldon3star_18"
# comp="Ldon4_4"
# comp="Ldon5star_7"



for (comp in c("Linf_vs_Ldon","Linf1star_43","CH.Linf.nohy_3","CUK.Linf_11","Ldon1star_44","Ldon2_7","Ldon3star_18","Ldon4_4","Ldon5star_7"))
{
  setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/Global_genome_diversity_Ldonovani_complex/11_selection_signatures/")
  dir.create(paste0("SNPeff_",comp), showWarnings = F)
  setwd(paste0("SNPeff_",comp))
  
  
  #----------------------------------------------
  # input data
  
  # 1 columns with genes of interest
  tag=paste0("SNPeff_",comp,"_stats.genes_high_mod")
  genes<-data.table(read.table(paste0("SNPeff_",comp,"_stats.genes_high_mod.id")))
  
  
  
  #----------------------------------------------
  # variables
  annV=38 # gff annotation version from TritryDB
  # function to discriminate gene set of interest from the background genes
  geneSelFunc <-function (score) { return(round(score) >0) }
  
  
  
  
  
  #----------------------------------------------
  # setup
  #
  # 1) load geneID to GO term annotation
  geneID2GO <-readMappings(file = paste0("../../10_CNVs/CNV_gene/TriTrypDB-",annV,"_LinfantumJPCM5_geneID_to_GOterm.txt"))
  # GO2geneID <-readMappings(file = paste0("~/work/references/GO_annotation/TriTrypDB-",annV,"_LinfantumJPCM5_GOterm_to_geneID.txt"))
  #
  # numbers tested for v38:
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
  
  # prepare geneList with genes set of interest: val=1 bg genes val=0
  genes[,val:=1]
  allGenes<-data.table(names(geneID2GO))
  bgGenes<-allGenes[!V1 %in% as.character(unlist(genes[,1]))]
  bgGenes[,val:=0]
  geneList<-rbind(genes,bgGenes)
  geneL<-unlist(geneList[,2])
  names(geneL)<-as.character(unlist(geneList[,1]))
  
  GO_analysis <-function(geneL, geneSelFunc)
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
  
  #
  # testing for enrichment of genes where the majority of samples has an increased copy number change >=4
  #
  res<-GO_analysis(geneL, geneSelFunc)
  res[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)]
  
  # add genes ids of genes responsible for the enrichment 
  for (goID in res[[2]][classicFish<0.05 | weightFish<0.05][,GO.ID])
  {
    as.character(unlist(genes[,1])) # all genes from set of interest
    genesInTerm(res[[1]], whichGO = goID)[[1]] # all genes from that GO term
    sigGenesOfGOid <- intersect( as.character(unlist(genes[,1])), genesInTerm(res[[1]], whichGO = goID)[[1]] ) 
    res[[2]][GO.ID==goID, sig_Genes:=paste0(sigGenesOfGOid, collapse = ",")]
  }
  
  write.table(res[[2]][classicFish<0.05 | weightFish<0.05][order(weightFish)], file=paste0(tag,"_topGO_classic.weight.0.05.txt"), 
              quote=T, row.names=F)
  
  #
  # pdf(paste0(tag,"_topGO_algos.classic.weight.weight01.elim.pdf"),width=24,height=20)
  # par(mfrow=c(2,2))
  # showSigOfNodes(res[[1]], score(res[[3]]), firstSigNodes = nrow(res[[2]][classicFish<0.05]), useInfo = "all")
  # showSigOfNodes(res[[1]], score(res[[4]]), firstSigNodes = nrow(res[[2]][weightFish<0.05]), useInfo = "all")
  # showSigOfNodes(res[[1]], score(res[[5]]), firstSigNodes = nrow(res[[2]][weight01Fish<0.05]), useInfo = "all")
  # showSigOfNodes(res[[1]], score(res[[6]]), firstSigNodes = nrow(res[[2]][elimFish<0.05]), useInfo = "all")
  # dev.off()
  # 
  # pdf(paste0(tag,"_topGO_algos.weight.pdf"),width=10,height=8)
  # par(mfrow=c(1,1))
  # showSigOfNodes(res[[1]], score(res[[4]]), firstSigNodes = nrow(res[[2]][weightFish<0.05]), useInfo = "all")
  # dev.off()
  
  # general stats for weightFish
  # geneData(res[[4]])
  
  stats<-rbind(c("Number genes of interest",paste0("Number bg genes in version v",annV)), c(nrow(genes), nrow(bgGenes)))
  write.table(stats, file=paste0(tag,"_topGO_stats.txt"), quote=T, row.names=F, col.names = F)
  
}




library(data.table)
library(ggplot2)
library(ade4)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/02_IBD_analysis/")

#-----------------------------
# functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
{
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

ReplaceLowerOrUpperTriangle <- function(m, triangle.to.replace)
{
  # If triangle.to.replace="lower", replaces the lower triangle of a square matrix with its upper triangle.
  # If triangle.to.replace="upper", replaces the upper triangle of a square matrix with its lower triangle.
  
  if (nrow(m) != ncol(m)) stop("Supplied matrix must be square.")
  if      (tolower(triangle.to.replace) == "lower") tri <- lower.tri(m)
  else if (tolower(triangle.to.replace) == "upper") tri <- upper.tri(m)
  else stop("triangle.to.replace must be set to 'lower' or 'upper'.")
  m[tri] <- t(m)[tri]
  return(m)
}

GeoDistanceInMetresMatrix <- function(df.geopoints)
{
  # Returns a matrix (M) of distances between geographic points.
  # M[i,j] = M[j,i] = Distance between (df.geopoints$lat[i], df.geopoints$lon[i]) and
  # (df.geopoints$lat[j], df.geopoints$lon[j]).
  # The row and column names are given by df.geopoints$name.
  
  GeoDistanceInMetres <- function(g1, g2){
    # Returns a vector of distances. (But if g1$index > g2$index, returns zero.)
    # The 1st value in the returned vector is the distance between g1[[1]] and g2[[1]].
    # The 2nd value in the returned vector is the distance between g1[[2]] and g2[[2]]. Etc.
    # Each g1[[x]] or g2[[x]] must be a list with named elements "index", "lat" and "lon".
    # E.g. g1 <- list(list("index"=1, "lat"=12.1, "lon"=10.1), list("index"=3, "lat"=12.1, "lon"=13.2))
    DistM <- function(g1, g2){
      require("Imap")
      return(ifelse(g1$index > g2$index, 0, gdist(lat.1=g1$lat, lon.1=g1$lon, lat.2=g2$lat, lon.2=g2$lon, units="km")))
    }
    return(mapply(DistM, g1, g2))
  }
  
  n.geopoints <- nrow(df.geopoints)
  
  # The index column is used to ensure we only do calculations for the upper triangle of points
  df.geopoints$index <- 1:n.geopoints
  
  # Create a list of lists
  list.geopoints <- by(df.geopoints[,c("index", "lat", "lon")], 1:n.geopoints, function(x){return(list(x))})
  
  # Get a matrix of distances (in metres)
  mat.distances <- ReplaceLowerOrUpperTriangle(outer(list.geopoints, list.geopoints, GeoDistanceInMetres), "lower")
  
  # Set the row and column names
  rownames(mat.distances) <- df.geopoints$name
  colnames(mat.distances) <- df.geopoints$name
  
  return(mat.distances)
}

adj_name_notation<-function(old_name)
{
  old_name<-as.character(old_name)
  old_name[startsWith(old_name, "X")]<-gsub("X", "", old_name[startsWith(old_name, "X")])
  old_name<-gsub("\\.", "-", old_name)
  old_name<-gsub("BHU1062-4", "BHU1062_4", old_name)
  old_name<-gsub("LRC_L1303", "LRC-L1303", old_name)
  old_name<-gsub("LRC_L47", "LRC-L47", old_name)
}

manteltest_for_samps<-function(all, samps, nperm)
{
  sub_Neisd<-Neisd_m[sample %in% samps, colnames(Neisd_m) %in% samps, with=F]
  sub_sampdist<-samp_dist_m[sample %in% samps, colnames(samp_dist_m) %in% samps, with=F]
  print(dim(sub_Neisd))
  # mantel.rtest(as.dist(sub_Neisd),as.dist(sub_sampdist),nrepet = 100)
  mantel.randtest(as.dist(sub_Neisd),as.dist(sub_sampdist),nrepet = nperm)
}

mantel_summary<-function(name, res)
{
  c(group=name, correlation=res$obs, pvalue=res$pvalue, res$expvar, nperm=res$rep)
}
#-----------------------------


# load geographical coordinates of the sampling locations
coord<-data.table(read.table("TableS1_coordinates.txt", header = T))
setnames(coord,colnames(coord)[c(1,5,6)],c("name","lat","lon"))
coord[,name:=as.character(name)]
coord<-coord[order(name)]
coord[, name:=gsub("Ldon282CL4", "Ldon282cl2", name)]
# calulate geographical distances between sampling locations
samp_dist_m<-GeoDistanceInMetresMatrix(coord[,c(1,5,6)])
samp_dist<-data.table(reshape2::melt(samp_dist_m))
setnames(samp_dist,colnames(samp_dist),c("samp1","samp2","dist"))

# load weighted NeisD that was also used to build the tree
Neisd_m<-data.table(read.table("../01_phylogenetic_reconstruction/NeisD/LmjF.allchr.weighted_NeisD.ind_NA.frac0.2_dphC.txt"))
setnames(Neisd_m,colnames(Neisd_m),c("sample",as.character(Neisd_m$V1)))
#
# adjust name notation to match distance matrix
Neisd_m$sample<-adj_name_notation(Neisd_m$sample)
colnames(Neisd_m)<-adj_name_notation(colnames(Neisd_m))
#
Neisd<-melt(Neisd_m, id.vars = "sample")
setnames(Neisd,colnames(Neisd),c("samp1","samp2","NeisD"))
Neisd<-Neisd[order(samp1,samp2)]

# check that sample names between genetic and geographic distances match
comp_names<-data.table(cbind(unique(Neisd$samp1),as.character(unique(samp_dist$samp1))))
comp_names[V1!=V2]

all<-merge(Neisd,samp_dist)
all[,samp1:=as.character(samp1)]
all[,samp2:=as.character(samp2)]
all[, samp1_species:=unlist(lapply(all$samp1, function(x) coord[name==x, Leish_species]))]
all[, samp2_species:=unlist(lapply(all$samp2, function(x) coord[name==x, Leish_species]))]
all[, samp1_group:=unlist(lapply(all$samp1, function(x) coord[name==x, Leish_group]))]
all[, samp2_group:=unlist(lapply(all$samp2, function(x) coord[name==x, Leish_group]))]
all[, samp1_group_colour:=unlist(lapply(all$samp1, function(x) coord[name==x, Leish_group_colour]))]
all[, samp2_group_colour:=unlist(lapply(all$samp2, function(x) coord[name==x, Leish_group_colour]))]
all[,comb_group:=paste(samp1_group,samp2_group, sep="_")]
all[samp1_species!=samp2_species,comp_species:="both"]
all[samp1_species==samp2_species & samp1_species=="L.infantum",comp_species:="Linf"]
all[samp1_species==samp2_species & samp1_species=="L.donovani",comp_species:="Ldon"]
all<-all[order(samp1_group, samp2_group)]
# Linf core definition
all[,Linf_Mon1:=F]
all[samp1_group %in% c("Linf1") & samp2_group %in% c("Linf1"),Linf_Mon1:=T]
all[samp1 %in% c("Inf152","WR285","Inf045","LRC-L47","Inf055","Inf004") | samp2 %in% c("Inf152","WR285","Inf045","LRC-L47","Inf055","Inf004"),Linf_Mon1:=F]
#


# prepare data for Mantel test
Neisd_m
samp_dist_m
#
samp_dist_m<-GeoDistanceInMetresMatrix(coord[,c(1,5,6)])
a<-colnames(samp_dist_m)
samp_dist_m<-data.table(samp_dist_m)
samp_dist_m[,sample:=a]
setcolorder(samp_dist_m, c("sample",sort(a)))
samp_dist_m<-samp_dist_m[order(sample)]
#
setcolorder(Neisd_m, c("sample",sort(a)))
Neisd_m<-Neisd_m[order(sample)]
# check sample names and orders are identical
sum(colnames(Neisd_m)!=colnames(samp_dist_m))
# [1] 0
sum(Neisd_m$sample!=samp_dist_m$sample)
# [1] 0



pdf("IBD_overview.pdf")
ggplot(all, aes(dist, NeisD)) + 
  geom_point() +
  facet_grid(samp1_species~samp2_species)
dev.off()






##########
# Ldon
Ldon_gg<-ggplot(all[comp_species=="Ldon"], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) + 
  labs(title=expression(paste(italic("L. donovani"), " all")), x="Distance [km]")# + theme_update(plot.title = element_text(hjust = 0.5))
Ldon_mantel<-manteltest_for_samps(all, unique(c(all[comp_species=="Ldon",samp1],all[comp_species=="Ldon",samp2])), 10000)

##########
# Linf
Linfall_gg<-ggplot(all[comp_species=="Linf"], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  labs(title=expression(paste(italic("L. infantum"), " all")), x="Distance [km]") 
Linfall_mantel<-manteltest_for_samps(all, unique(c(all[comp_species=="Linf",samp1],all[comp_species=="Linf",samp2])), 10000)
#
Linf1_gg<-ggplot(all[(samp1_group=="Linf1") & (samp2_group=="Linf1")], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  labs(title=expression(paste(italic("L. infantum"), " Linf1: MON1 & non-MON1")), x="Distance [km]") 
Linf1_samps<-unique(c(all[samp1_group=="Linf1",samp1],all[samp2_group=="Linf1",samp2]))
Linf1_mantel<-manteltest_for_samps(all, Linf1_samps, 10000)
#
american_samps<-c("MAM","WR285","Cha001","WC","ARL","HN167","HN336")
#
Linf1_noAmerica_gg<-ggplot(all[(samp1_group=="Linf1") & (samp2_group=="Linf1")
                               & !samp1 %in% american_samps & !samp2 %in% american_samps], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  labs(title=expression(paste(italic("L. infantum"), " Linf1: MON1 & non-MON1 no America")), x="Distance [km]") 
Linf1_noAmerica_samps<-unique(c(all[samp1_group=="Linf1" & !samp1 %in% american_samps,samp1],
                                all[samp2_group=="Linf1" & !samp2 %in% american_samps,samp2]))
Linf1_noAmerica_mantel<-manteltest_for_samps(all, Linf1_noAmerica_samps, 10000)
#
Linf_Mon1_gg<-ggplot(all[Linf_Mon1==T,], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  labs(title=expression(paste(italic("L. infantum"), " Linf1: MON1")), x="Distance [km]") 
Linf_Mon1_mantel<-manteltest_for_samps(all, unique(c(all[Linf_Mon1==T,samp1],all[Linf_Mon1==T,samp2])), 10000)
#
Linf_Mon1_noAmerica_gg<-ggplot(all[Linf_Mon1==T & !samp1 %in% american_samps & !samp2 %in% american_samps,], aes(dist, NeisD)) + 
  geom_point() + geom_smooth(method='lm',formula=y~x) +
  labs(title=expression(paste(italic("L. infantum"), " Linf1: MON1 no America")), x="Distance [km]") 
Linf_Mon1_noAmerica_samps<-unique(c(all[Linf_Mon1==T & !samp1 %in% american_samps,samp1],
                                    all[Linf_Mon1==T & !samp2 %in% american_samps,samp2]))
Linf_Mon1_noAmerica_mantel<-manteltest_for_samps(all, Linf_Mon1_noAmerica_samps, 10000)


cor.test(all[samp1 %in% Linf1_samps & samp2 %in% Linf1_samps,NeisD],
         all[samp1 %in% Linf1_samps & samp2 %in% Linf1_samps,dist])
cor.test(all[samp1 %in% Linf1_noAmerica_samps & samp2 %in% Linf1_noAmerica_samps,NeisD],
         all[samp1 %in% Linf1_noAmerica_samps & samp2 %in% Linf1_noAmerica_samps,dist])




res_mantel<-data.table(rbind(mantel_summary("L. donovani all",Ldon_mantel),
      mantel_summary("L. infantum all",Linfall_mantel),
      mantel_summary("L. infantum Linf1: MON1 & non-MON1",Linf1_mantel),
      mantel_summary("L. infantum Linf1: MON1 & non-MON1 - no America",Linf1_noAmerica_mantel),
      mantel_summary("L. infantum Linf1: MON1",Linf_Mon1_mantel),
      mantel_summary("L. infantum Linf1: MON1 - no America",Linf_Mon1_noAmerica_mantel)))
res_mantel[,correlation:=as.numeric(correlation)]
res_mantel[,pvalue:=as.numeric(pvalue)]
res_mantel[,Std.Obs:=as.numeric(Std.Obs)]
res_mantel[,Expectation:=as.numeric(Expectation)]
res_mantel[,Variance:=as.numeric(Variance)]
res_mantel
write.table(res_mantel, file="IBD_groups_manteltest.txt", quote = T, row.names = F, sep = ",")

# pdf("IBD_groups.pdf",width=8, height=12)
png("IBD_groups.png",width=800, height=1200)
multiplot(Ldon_gg, Linf1_gg, Linf1_noAmerica_gg, Linfall_gg, Linf_Mon1_gg, Linf_Mon1_noAmerica_gg, cols = 2)
dev.off()

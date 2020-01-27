
library(foreach)
library(data.table)
library(reshape)
library(ggplot2)
library(circlize)
library(ape)


setwd("~/work/Leish_donovaniComplex/MS01_globalDiversity_Ldonovani/github_scripts/05_heterozygosities/")

load("../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")



#-----------------------------------
# functions

# calculates heterozygosities
get_hets <- function(dat)
{
  i=1; gg=c()
  hets <- foreach (cc=dat[,9:ncol(dat)], .combine=c) %do%
  {
    # cc=unlist(dat[,15])
    het_vector <-strsplit(as.character(cc), split = c(":"))
    het_vector <-unlist(het_vector)[ c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE) ]
    
    print(paste("sample",i)); i=i+1
    #     print(unique(het_vector))
    gg <- unique(c(gg,het_vector))
    # [1] "1/1" "0/0" "./." "0/1"
    
    het_vector <- replace(het_vector, het_vector=="0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/1", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="1/1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="./.", NA)
    het_vector <- replace(het_vector, het_vector=="0/2", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="2/2", 0^2+1^2)
    
    het_vector <- replace(het_vector, het_vector=="0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/1", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="././.", NA)
    
    het_vector <- replace(het_vector, het_vector=="0/0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/1", 0.25^2+0.75^2)
    het_vector <- replace(het_vector, het_vector=="0/0/1/1", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/1", 0.25^2+0.75^2)
    het_vector <- replace(het_vector, het_vector=="./././.", NA)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/2", 0.25^2+0.75^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2/2", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/2", 0.25^2+0.75^2)
    
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1/1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/1/1", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="0/0/1/1/1", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/1/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="././././.", NA)
    het_vector <- replace(het_vector, het_vector=="2/2/2/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/2/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1/2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2/2/1", 0.4^2+0.4^2+0.2^2)
    
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0/0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1/1/1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="./././././.", NA)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0/1", (1/6)^2+(5/6)^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/1/1/1", (1/6)^2+(5/6)^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0/2", (1/6)^2+(5/6)^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/1/1", (2/6)^2+(4/6)^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/2/2/2", 0^2+1^2)
    
    het_vector <- replace(het_vector, het_vector=="1/1/2", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="0", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="1", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector==".", NA)
    het_vector <- replace(het_vector, het_vector=="2", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="2/1/1", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="1/2", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="2/1", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/2", 0.75^2+0.25^2)
    het_vector <- replace(het_vector, het_vector=="2/1/1/1", 0.75^2+0.25^2)
    het_vector <- replace(het_vector, het_vector=="1/2/2", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="2/1/1/1/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/2/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="0/2/1", (1/3)^2+(1/3)^2+(1/3)^2)
    het_vector <- replace(het_vector, het_vector=="2/2/1", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="1/1/2/2", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/1/2", (1/3)^2+(1/3)^2+(1/3)^2) 
    het_vector <- replace(het_vector, het_vector=="1/1/2/2/2", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="1/2/2/2", 0.25^2+0.75^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/2", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/0/1/2", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/2/1/1", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/1/2/2", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/1", 0.25^2+0.75^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2/1", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="2/2/1/1", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="2/1/1/1/1/1", (1/6)^2+(5/6)^2)
    het_vector <- replace(het_vector, het_vector=="2/2/1/1/1", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="2/1/1/1/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/2/1", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="1/2/2/2/2", 0.2^2+0.8^2)
    het_vector <- replace(het_vector, het_vector=="0/0/2/1/1", 0.4^2+0.4^2+0.2^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/2/2", 0.4^2+0.4^2+0.2^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/1/1", 0.4^2+0.4^2+0.2^2)
    het_vector <- replace(het_vector, het_vector=="0/1/1/1/2", 0.6^2+0.2^2+0.2^2)
    het_vector <- replace(het_vector, het_vector=="2/2/2/1/1", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/1/1", 0.4^2+0.4^2+0.2^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/1", 0.25^2+0.25^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1/1/2", (5/6)^2+(1/6)^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/2/2", 0.4^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="1/1/3", (1/3)^2+(2/3)^2)
    het_vector <- replace(het_vector, het_vector=="3/3", 0^2+1^2)
    het_vector <- replace(het_vector, het_vector=="0/2/2/2/1", 0.2^2+0.2^2+0.6^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/2/1", (4/6)^2+(1/6)^2+(1/6)^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0/0/0/0", 1^2)
    het_vector <- replace(het_vector, het_vector=="./././././././.", NA)
    het_vector <- replace(het_vector, het_vector=="0/0/1/1/1/1/1/1", (2/8)^2+(6/8)^2)
    het_vector <- replace(het_vector, het_vector=="1/1/1/1/1/1/1/1", 1^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/1/1/1/1", 0.5^2+0.5^2)
    het_vector <- replace(het_vector, het_vector=="0/0/0/0/0/0/1/1", (2/8)^2+(6/8)^2)
    
    #     print(unique(het_vector[grep("/",het_vector)]))
    
    het_vector <- as.double(het_vector)
    # print(table(het_vector))
    
    hh <- 1 - sum(het_vector, na.rm = T) / sum(!is.na(het_vector))
    hh
  }
  names(hets) <-colnames(dat[,9:ncol(dat)])
  hets
}

# differenciats between heterozygotes, homozygotes all polarized by the reference
# 0: hom ref
# 1: het with ref and alt alleles
# 2: het with 2 diferent alt alleles
# 3: hom alt
get_genos <- function(dat)
{
  i=1; gg=c()
  het.info <- foreach (cc=dat[,9:ncol(dat)], .combine=cbind) %do%
  {
    
    geno <-strsplit(as.character(cc), split = c(":"))
    geno <-unlist(geno)[ c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE) ]
    
    print(paste("sample",i)); i=i+1
    #     print(unique(geno))
    gg <- unique(c(gg,geno))
    # [1] "1/1" "0/0" "./." "0/1"
    geno <- replace(geno, geno=="0", 0)
    geno <- replace(geno, geno=="1", 3)
    geno <- replace(geno, geno==".", NA)
    geno <- replace(geno, geno=="2", 3)
    
    geno <- replace(geno, geno=="0/0", 0)
    geno <- replace(geno, geno=="0/1", 1)
    geno <- replace(geno, geno=="1/1", 3)
    geno <- replace(geno, geno=="./.", NA)
    geno <- replace(geno, geno=="0/2", 1)
    geno <- replace(geno, geno=="2/2", 3)
    
    geno <- replace(geno, geno=="0/0/0", 0)
    geno <- replace(geno, geno=="0/0/1", 1)
    geno <- replace(geno, geno=="0/1/1", 1)
    geno <- replace(geno, geno=="1/1/1", 3)
    geno <- replace(geno, geno=="0/0/0", 0)
    geno <- replace(geno, geno=="0/0/2", 1)
    geno <- replace(geno, geno=="0/2/2", 1)
    geno <- replace(geno, geno=="2/2/2", 3)
    geno <- replace(geno, geno=="././.", NA)
    
    geno <- replace(geno, geno=="0/0/0/0", 0)
    geno <- replace(geno, geno=="1/1/1/1", 3)
    geno <- replace(geno, geno=="0/0/0/1", 1)
    geno <- replace(geno, geno=="0/0/1/1", 1)
    geno <- replace(geno, geno=="0/1/1/1", 1)
    geno <- replace(geno, geno=="./././.", NA)
    geno <- replace(geno, geno=="0/0/0/0", 0)
    geno <- replace(geno, geno=="2/2/2/2", 3)
    geno <- replace(geno, geno=="0/0/0/2", 1)
    geno <- replace(geno, geno=="0/0/2/2", 1)
    geno <- replace(geno, geno=="0/2/2/2", 1)
    
    geno <- replace(geno, geno=="0/0/0/0/0", 0)
    geno <- replace(geno, geno=="1/1/1/1/1", 3)
    geno <- replace(geno, geno=="0/0/0/0/1", 1)
    geno <- replace(geno, geno=="0/0/0/1/1", 1)
    geno <- replace(geno, geno=="0/0/1/1/1", 1)
    geno <- replace(geno, geno=="0/1/1/1/1", 1)
    geno <- replace(geno, geno=="././././.", NA)
    geno <- replace(geno, geno=="2/2/2/2/2", 3)
    geno <- replace(geno, geno=="0/0/0/0/2", 1)
    geno <- replace(geno, geno=="0/0/2/2/2", 1)
    geno <- replace(geno, geno=="0/0/0/2/2", 1)
    geno <- replace(geno, geno=="0/2/2/2/2", 1)
    geno <- replace(geno, geno=="1/1/1/1/2", 2)
    geno <- replace(geno, geno=="0/0/2/2/1", 1)
    
    geno <- replace(geno, geno=="0/0/0/0/0/0", 0)
    geno <- replace(geno, geno=="1/1/1/1/1/1", 3)
    geno <- replace(geno, geno=="./././././.", NA)
    geno <- replace(geno, geno=="0/0/0/0/0/1", 1)
    geno <- replace(geno, geno=="0/1/1/1/1/1", 1)
    geno <- replace(geno, geno=="0/0/0/0/0/2", 1)
    geno <- replace(geno, geno=="0/0/0/0/1/1", 1)
    geno <- replace(geno, geno=="2/2/2/2/2/2", 3)
    
    geno <- replace(geno, geno=="1/1/2", 2)
    geno <- replace(geno, geno=="2/1/1", 2)
    geno <- replace(geno, geno=="1/2", 2)
    geno <- replace(geno, geno=="2/1", 2)
    geno <- replace(geno, geno=="1/1/1/2", 2)
    geno <- replace(geno, geno=="2/1/1/1", 2)
    geno <- replace(geno, geno=="1/2/2", 2)
    geno <- replace(geno, geno=="2/1/1/1/1", 2)
    geno <- replace(geno, geno=="2/2/2/2/1", 2)
    geno <- replace(geno, geno=="0/2/1", 1)
    geno <- replace(geno, geno=="2/2/1", 2)
    geno <- replace(geno, geno=="1/1/2/2", 2)
    geno <- replace(geno, geno=="0/1/2", 1) 
    geno <- replace(geno, geno=="1/1/2/2/2", 2)
    geno <- replace(geno, geno=="1/2/2/2", 2)
    geno <- replace(geno, geno=="0/1/1/2", 1)
    geno <- replace(geno, geno=="0/0/1/2", 1)
    geno <- replace(geno, geno=="0/2/1/1", 1)
    geno <- replace(geno, geno=="0/1/2/2", 1)
    geno <- replace(geno, geno=="2/2/2/1", 2)
    geno <- replace(geno, geno=="0/0/2/1", 1)
    geno <- replace(geno, geno=="2/2/1/1", 2)
    geno <- replace(geno, geno=="2/1/1/1/1/1", 2)
    geno <- replace(geno, geno=="2/2/1/1/1", 2)
    geno <- replace(geno, geno=="2/1/1/1/1", 2)
    geno <- replace(geno, geno=="2/2/2/2/1", 2)
    geno <- replace(geno, geno=="1/2/2/2/2", 2)
    geno <- replace(geno, geno=="0/0/2/1/1", 1)
    geno <- replace(geno, geno=="0/1/1/2/2", 1)
    geno <- replace(geno, geno=="0/2/2/1/1", 1)
    geno <- replace(geno, geno=="0/1/1/1/2", 1)
    geno <- replace(geno, geno=="2/2/2/1/1", 2)
    geno <- replace(geno, geno=="0/2/2/1/1", 1)
    geno <- replace(geno, geno=="0/2/2/1", 1)
    geno <- replace(geno, geno=="1/1/1/1/1/2", 2)
    geno <- replace(geno, geno=="1/1/1/2/2", 2)
    geno <- replace(geno, geno=="1/1/3", 2)
    geno <- replace(geno, geno=="3/3", 3)
    geno <- replace(geno, geno=="0/2/2/2/1", 1)
    geno <- replace(geno, geno=="0/0/0/0/2/1", 1)
    geno <- replace(geno, geno=="0/0/0/0/0/0/0/0", 0)
    geno <- replace(geno, geno=="./././././././.", NA)
    geno <- replace(geno, geno=="0/0/1/1/1/1/1/1", 1)
    geno <- replace(geno, geno=="1/1/1/1/1/1/1/1", 3)
    geno <- replace(geno, geno=="0/0/0/0/1/1/1/1", 1)
    geno <- replace(geno, geno=="0/0/0/0/0/0/1/1", 1)
    
    #     print(unique(geno[grep("/",geno)]))
    
    geno <- as.double(geno)
    
    geno
  }
  colnames(het.info) <-colnames(dat[9:ncol(dat)])
  het.info
}

# get the percentage of reference alleles
get_refallels <- function(dat)
{
  i=1; gg=c()
  ref <- foreach (cc=dat[,9:ncol(dat)], .combine=c) %do%
  {
    geno <-strsplit(as.character(cc), split = c(":"))
    geno <-unlist(geno)[ c(TRUE,FALSE,FALSE,FALSE,FALSE,FALSE) ]
    
    print(paste("sample",i)); i=i+1
    #     print(unique(geno))
    gg <- unique(c(gg,geno))
    # [1] "1/1" "0/0" "./." "0/1"
    
    geno <- replace(geno, geno=="0/0", 1)
    geno <- replace(geno, geno=="0/1", 0.5)
    geno <- replace(geno, geno=="1/1", 0)
    geno <- replace(geno, geno=="./.", NA)
    geno <- replace(geno, geno=="0/2", 0.5)
    geno <- replace(geno, geno=="2/2", 0)
    
    geno <- replace(geno, geno=="0/0/0", 1)
    geno <- replace(geno, geno=="0/0/1", (2/3))
    geno <- replace(geno, geno=="0/1/1", (1/3))
    geno <- replace(geno, geno=="1/1/1", 0)
    geno <- replace(geno, geno=="0/0/0", 1)
    geno <- replace(geno, geno=="0/0/2", (2/3))
    geno <- replace(geno, geno=="0/2/2", (1/3))
    geno <- replace(geno, geno=="2/2/2", 0)
    geno <- replace(geno, geno=="././.", NA)
    
    geno <- replace(geno, geno=="0/0/0/0", 1)
    geno <- replace(geno, geno=="1/1/1/1", 0)
    geno <- replace(geno, geno=="0/0/0/1", 0.75)
    geno <- replace(geno, geno=="0/0/1/1", 0.5)
    geno <- replace(geno, geno=="0/1/1/1", 0.25)
    geno <- replace(geno, geno=="./././.", NA)
    geno <- replace(geno, geno=="0/0/0/0", 1)
    geno <- replace(geno, geno=="2/2/2/2", 0)
    geno <- replace(geno, geno=="0/0/0/2", 0.75)
    geno <- replace(geno, geno=="0/0/2/2", 0.5)
    geno <- replace(geno, geno=="0/2/2/2", 0.25)
    
    geno <- replace(geno, geno=="0/0/0/0/0", 1)
    geno <- replace(geno, geno=="1/1/1/1/1", 0)
    geno <- replace(geno, geno=="0/0/0/0/1", 0.8)
    geno <- replace(geno, geno=="0/0/0/1/1", 0.6)
    geno <- replace(geno, geno=="0/0/1/1/1", 0.4)
    geno <- replace(geno, geno=="0/1/1/1/1", 0.2)
    geno <- replace(geno, geno=="././././.", NA)
    geno <- replace(geno, geno=="2/2/2/2/2", 0)
    geno <- replace(geno, geno=="2/2/2/2/2", 0)
    geno <- replace(geno, geno=="0/0/0/0/2", 0.8)
    geno <- replace(geno, geno=="0/0/2/2/2", 0.4)
    geno <- replace(geno, geno=="0/0/0/2/2", 0.6)
    geno <- replace(geno, geno=="0/2/2/2/2", 0.2)
    geno <- replace(geno, geno=="1/1/1/1/2", 0)
    geno <- replace(geno, geno=="0/0/2/2/1", 0.4)
    
    geno <- replace(geno, geno=="0/0/0/0/0/0", 1)
    geno <- replace(geno, geno=="1/1/1/1/1/1", 0)
    geno <- replace(geno, geno=="./././././.", NA)
    geno <- replace(geno, geno=="0/0/0/0/0/1", (5/6))
    geno <- replace(geno, geno=="0/1/1/1/1/1", (1/6))
    geno <- replace(geno, geno=="0/0/0/0/0/2", (5/6))
    geno <- replace(geno, geno=="0/0/0/0/1/1", (4/6))
    geno <- replace(geno, geno=="2/2/2/2/2/2", 0)
    
    
    geno <- replace(geno, geno=="1/1/2", 0)
#     geno <- replace(geno, geno=="0", 1)
#     geno <- replace(geno, geno=="1", 0)
    geno <- replace(geno, geno==".", NA)
    geno <- replace(geno, geno=="2", 0)
    geno <- replace(geno, geno=="2/1/1", 0)
    geno <- replace(geno, geno=="1/2", 0)
    geno <- replace(geno, geno=="2/1", 0)
    geno <- replace(geno, geno=="1/1/1/2", 0)
    geno <- replace(geno, geno=="2/1/1/1", 0)
    geno <- replace(geno, geno=="1/2/2", 0)
    geno <- replace(geno, geno=="2/1/1/1/1", 0)
    geno <- replace(geno, geno=="2/2/2/2/1", 0)
    geno <- replace(geno, geno=="0/2/1", 1/3)
    geno <- replace(geno, geno=="2/2/1", 0)
    geno <- replace(geno, geno=="1/1/2/2", 0)
    geno <- replace(geno, geno=="0/1/2", 1/3) 
    geno <- replace(geno, geno=="1/1/2/2/2", 0)

    geno <- replace(geno, geno=="1/2/2/2", 0)
    geno <- replace(geno, geno=="0/1/1/2", 0.25)
    geno <- replace(geno, geno=="0/0/1/2", 0.5)
    geno <- replace(geno, geno=="0/2/1/1", 0.25)
    geno <- replace(geno, geno=="0/1/2/2", 0.25)
    geno <- replace(geno, geno=="2/2/2/1", 0)
    geno <- replace(geno, geno=="0/0/2/1", 0.5)
    geno <- replace(geno, geno=="2/2/1/1", 0)
    geno <- replace(geno, geno=="2/1/1/1/1/1", 0)
    geno <- replace(geno, geno=="2/2/1/1/1", 0)
    geno <- replace(geno, geno=="2/1/1/1/1", 0)
    geno <- replace(geno, geno=="2/2/2/2/1", 0)
    geno <- replace(geno, geno=="1/2/2/2/2", 0)
    geno <- replace(geno, geno=="0/0/2/1/1", 0.4)
    geno <- replace(geno, geno=="0/1/1/2/2", 0.2)
    geno <- replace(geno, geno=="0/2/2/1/1", 0.2)
    geno <- replace(geno, geno=="0/1/1/1/2", 0.2)
    geno <- replace(geno, geno=="2/2/2/1/1", 0)
    geno <- replace(geno, geno=="0/2/2/1/1", 02)
    geno <- replace(geno, geno=="0/2/2/1", 0.25)
    geno <- replace(geno, geno=="1/1/1/1/1/2", 0)
    geno <- replace(geno, geno=="1/1/1/2/2", 0)
    geno <- replace(geno, geno=="1/1/3", 0)
    geno <- replace(geno, geno=="3/3", 0)
    geno <- replace(geno, geno=="0/2/2/2/1", 0.2)
    geno <- replace(geno, geno=="0/0/0/0/2/1", (4/6))
    geno <- replace(geno, geno=="0/0/0/0/0/0/0/0", 1)
    geno <- replace(geno, geno=="./././././././.", NA)
    geno <- replace(geno, geno=="0/0/1/1/1/1/1/1", 0.25)
    geno <- replace(geno, geno=="1/1/1/1/1/1/1/1", 1)
    geno <- replace(geno, geno=="0/0/0/0/1/1/1/1", 0.5)
    geno <- replace(geno, geno=="0/0/0/0/0/0/1/1", 2/3)
    
    geno <- as.double(geno)
    pcref <- sum(geno, na.rm = T) / sum(!is.na(geno))
    
  }
  names(ref) <-colnames(dat[9:ncol(dat)])
  ref
}

plot_mds <- function(mydata, title="")
{
  d <- dist(mydata) # euclidean distances between the rows
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  fit # view results
  
  # plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
#   plot(x, y, xlab="MDS, Coordinate 1", ylab="MDS, Coordinate 2", main=title, type="n")
#   text(x, y, labels = row.names(mydata), cex=.7, col=group.col)
  
  dat<-data.table(fit$points)
  dat[,samp:=rownames(fit$points)]
  
  gg <-ggplot(data=dat, aes(x=V1, y=V2)) +
    geom_point(colour=group.col) + 
    geom_text(aes(label=dat$samp , colour = group.col, hjust=0.5, vjust=-0.5))  + 
    scale_color_manual(values=c(phyloCol[1],phyloCol[3],phyloCol[2],phyloCol[5],phyloCol[4],"gray34"), guide=FALSE) +
    xlab("MDS, Coordinate 1") + ylab("MDS, Coordinate 2") + theme_bw()
#     ggtitle(title)                       
  print(gg)  

  fit
}

plot_pca <- function(mydata, title="")
{
  fit <- princomp(mydata, cor=TRUE)
  summary(fit) # print variance accounted for
  loadings(fit) # pc loadings
  fit$scores # the principal components
#   biplot(fit) 
  plot(fit$scores[,1],fit$scores[,2], col="white", xlab="PC 1", ylab="PC 2", main=title)
  text(fit$scores[,1],fit$scores[,2], labels=names(hets), col=group.col)

  plot(fit,type="lines") # scree plot
}

print_sample_hets <-function(samples=c("CUK10","CUK8","CUK3"), cols=rainbow(length(samples)), hh=0.15)
{
  if(het.info[1,chr]=="LinJ.01"){
    xx<-sapply(strsplit(as.character(het.info[,chr]), split="J."),function(x) x[2])
    het.info[,chr:=xx]
  }
  
  circos.clear()
  par(mar = c(3,2,3,2)) 
  circos.par("track.height" = hh, start.degree=90)
  
  dcir <- het.info[,.(chr,pos,get(samples[1]))]
  setnames(dcir, colnames(dcir)[3], c("y"))
  dcir <-dcir[is.na(y), y:=0.5]
  #   dcir <-dcir[chr %in% c("01","02","03")]
  cc<-dcir$y; cc<-replace(cc,cc %in% c(0,3),"blue"); cc<-replace(cc,cc %in% c(1,2),"red"); cc<-replace(cc,cc==0.5,"white")
  dcir[,cols:=cc]
  # Step 1: initialize
  circos.initialize( factors=as.character(dcir$chr), x=dcir$pos )
  # Step 2: Build the regions. 
  circos.trackPlotRegion(factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
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
    aa = c(1, 0.5)
    if(theta < 90 || theta > 270)  aa = c(0, 0.5)
    circos.text(x=mean(xlim), y=4, labels=name, facing = dd, cex=2,  adj = aa)
  })
  # Step 3: Add points
  yy<-dcir[,cols]
  yy<-replace(yy,yy=="blue",cols[1])
  circos.trackPoints(as.character(dcir$chr), dcir$pos, dcir$y, col=yy , pch = 20, cex = 0.4)   
  for (i in 2:length(samples))
  {
    dcir <- het.info[,.(chr,pos,get(samples[i]))]
    setnames(dcir, colnames(dcir)[3], c("y"))
    dcir <-dcir[is.na(y), y:=0.5]
    #     dcir <-dcir[chr %in% c("01","02","03")]
    cc<-dcir$y; cc<-replace(cc,cc %in% c(0,3),"blue"); cc<-replace(cc,cc %in% c(1,2),"red"); cc<-replace(cc,cc==0.5,"white")
    dcir[,cols:=cc]
    circos.trackPlotRegion(factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
    {
      circos.axis( labels=FALSE, major.tick=FALSE)
    })
    # Step 3: Add points
    yy<-dcir[,cols]
    yy<-replace(yy,yy=="blue",cols[i])
    circos.trackPoints(as.character(dcir$chr), dcir$pos, dcir$y, col=yy , pch = 20, cex = 0.4)
  }
}


print_sample_hets_win <-function(samples=c("CUK10","CUK8","CUK3"), color=rainbow(length(samples)), hh=0.15)
{
  circos.clear()
  par(mar = c(3,2,3,2)) 
  circos.par("track.height" = hh, start.degree=90)
  
  dcir <- het.wins[,.(chr,w10kb,get(samples[1]))]
  setnames(dcir, colnames(dcir)[3], c("y"))
  #   dcir <-dcir[chr %in% c("01","02","03")]
  dcir[,cols:=color[1]]
  # Step 1: initialize
  circos.initialize( factors=as.character(dcir$chr), x=dcir$w10kb )
  # Step 2: Build the regions. 
  circos.trackPlotRegion(factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
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
  # Step 3: Add points
  circos.trackLines(as.character(dcir$chr), dcir$w10kb, dcir$y, col=dcir$cols, lwd=2)
  
  if(length(samples)>=2)
  {
    for (i in 2:length(samples))
    {
      dcir <- het.wins[,.(chr,w10kb,get(samples[i]))]
      setnames(dcir, colnames(dcir)[3], c("y"))
      #     dcir <-dcir[chr %in% c("01","02","03")]
      dcir[,cols:=color[i]]
      circos.trackPlotRegion(factors = as.character(dcir$chr), y = dcir$y, panel.fun = function(x, y) 
      {
        circos.axis( labels=FALSE, major.tick=FALSE)
      })
      # Step 3: Add points
      circos.trackLines(as.character(dcir$chr), dcir$w10kb, dcir$y, col=dcir$cols, lwd=2)
    }
  }

}

#-----------------------------------
# 
  
#-----------------------------------
# variables
# note: the vcflike file is provided as gzipped files split for each chromosome in he respective folder
# in order to read it in as done in the cmd below files have to be unzipped, combined and header information within the file removed
data<-read.table("../data/snp_files/snps.filt.leish_global.linj.complex.vcflike.txt", comment.char="", header=T)
dat <-data.table(data)
load("../01_phylogenetic_reconstruction/sample.info_NEWsubgroups.RData")
#
# follow three assignments take several minutes to caluclate, if needed repeatedly should once be calculated and saved as .RData file
hets <-get_hets(data)
het.info <-data.table(get_genos(data)); het.info[,chr:=data[,1]]; het.info[,pos:=data[,2]]
refs <-get_refallels(data) # gets percentage of reference alleles per sample
#
# resulting datatables:
sample.info
hets
het.info
refs
# save(hets, het.info, refs, file="calculated_tables.RData")
# load(file="calculated_tables.RData")

# get heterozigosities per chr used for analysis in 04_aneuploidy
# run once and save the data (runs several minutes ...)
chr.hets <- foreach (chr = c(paste0("LinJ.0",1:9), paste0("LinJ.",10:36)), .combine=cbind) %do%
{
  print(chr)
  get_hets(dat[CHR==chr])
}
chr.hets<-data.table(chr.hets)
setnames(chr.hets,colnames(chr.hets),paste0("chr",1:36))
chr.hets[,sample:=colnames(dat)[9:ncol(dat)]]
chr.hets<-melt(chr.hets, id.vars = "sample")
#
save(chr.hets, file="chr.hets.RData")





res <-data.table(cbind(hets,refs))
res[,groups:=sample.info[,groups]]
res[,group.col:=sample.info[,group.col]]
res[,samp:=colnames(het.info)[1:151]]
write.table(res, "heterozygosities_pcrefalleles.txt", quote=FALSE, row.names=F)

sum(res$hets <0.004)
# [1] 105
sum(res$hets >0.004)
# [1] 46
46/151
# [1] 0.3046358


res[samp=="MAM"]
# hets      refs     groups group.col samp
# 1: 0.0650172 0.9405769 other_Linf     black  MAM

# pdf("heterozygosities_hist.pdf",width=5, height=5)
# ggplot(data=res, aes(hets, fill=groups)) + 
#   geom_histogram(breaks=seq(0, 0.07, by = 0.001)) + theme_bw() +
#   scale_fill_manual(values=unique(sample.info[order(groups),group.col]))#, guide=FALSE) +
# dev.off()


pdf("heterozygosities_pcrefalleles.pdf",width=8, height=8)
ggplot(data=res, aes(x=refs, y=hets)) +
  geom_point(aes(colour=groups)) + 
  geom_text(aes(label=res$samp , colour = groups, hjust=0.5, vjust=-0.5))  + 
  scale_color_manual(values=unique(sample.info[order(groups),group.col]), guide=F) +
  xlab("Fraction of reference alleles") + ylab("Heterozygosity") + 
  theme_bw() + geom_hline(yintercept = 0.004, linetype="dotted") 
dev.off()

# pdf("heterozygosities_pcrefalleles_ylog.pdf",width=8, height=8)
# ggplot(data=res, aes(x=refs, y=hets)) +
#   geom_point(aes(colour=groups)) + 
#   geom_text(aes(label=res$samp , colour = groups, hjust=0.5, vjust=-0.5))  + 
#   scale_color_manual(values=unique(sample.info[order(groups),group.col]), guide=F) +
#   xlab("Fraction of reference alleles") + ylab("Heterozygosity") + 
#   theme_bw() + geom_hline(yintercept = 0.004, linetype="dotted") + scale_y_continuous(trans='log')
# dev.off()



#-----------------------------------
# analysis per chromosome

dat <-data.table(data)

id=c(paste0("0",1:9),10:36)
#id=c("01","02")

hets.chr <- foreach (chrom=paste0("LinJ.",id), .combine=rbind ) %do%
{
  print(chrom)
  het <- get_hets(data.frame(dat[CHR==chrom]))
  het
}

refs.chr <- foreach (chrom=paste0("LinJ.",id), .combine=rbind ) %do%
{
  print(chrom)
  ref <- get_refallels(data.frame(dat[CHR==chrom]))
  ref
}

hets.chr <-data.table(hets.chr)
refs.chr <-data.table(refs.chr)

write.table(hets.chr, "heterozygosities_chr.txt", quote=FALSE, row.names=F)
write.table(refs.chr, "pcrefallels_chr.txt", quote=FALSE, row.names=F)





#------------------------------------
# heterozygosity statistics, window analysis
dir.create("hetWindows")

# initialize windows
het.info[,w10kb:=pos%/%10000+1, by=chr]
samples <-colnames(het.info)

het.wins <-foreach(i=1:151, .combine=cbind) %do%
{
  het.info[,w10samp:=0]
  het.info[get(samples[i]) %in% c(1,2),w10samp:=sum(get(samples[i])), by=.(chr,w10kb)]
  xx<-unique(het.info[,.(chr,w10kb,w10samp)])
  xx[,y:=sum(w10samp),by=.(chr,w10kb)]
  aa<-unique(xx[,.(chr,w10kb,y)])
  aa$y
}
het.wins <-data.table(het.wins)
setnames(het.wins, colnames(het.wins), colnames(het.info[,1:151,with=FALSE]))
het.wins[,chr:=aa$chr]
het.wins[,w10kb:=aa$w10kb]

save(het.wins,file="het.wins.RData") # information needed in 06_hybridisation_signatures

# for(i in 1:151)
# {
#   #i=146
#   pdf(paste0("hetWindows/hetwin_",samples[i],".pdf"),width=10,height=8)
#   gg <-ggplot(data=het.wins, aes(x=w10kb, y=get(samples[i]))) +
#     geom_line() + facet_wrap(~chr, ncol=6, scales="free_x")
#   print(gg)
#   dev.off()
# }

#--------------------------------------------
#
# correlate hets with somies
#
ploidy <-data.table(read.table("../04_aneuploidy/somies_updated.txt", head=T))

samples <-colnames(het.info)[1:151]
cor.het_somy <-foreach(samp=samples, .combine=rbind) %do%
{
  xx<-cor.test(ploidy[,get(samp)], hets.chr[,get(samp)], method="spearman")
  c(xx$estimate,xx$p.value)
}
cor.het_somy <-data.table(cor.het_somy)
cor.het_somy[,sample:=colnames(het.info)[1:151]]
cor.het_somy[,groups:=groups]
cor.het_somy[,group.col:=sample.info[order(groups),group.col]]
setnames(cor.het_somy, c("rho","V2"), c("R.cor","p.val"))
cor.het_somy<-cor.het_somy[order(R.cor)]
cor.het_somy[,ord:=1:151]
cor.het_somy[,FDR:=p.adjust(p.val)]
#
pdf("cor.heterozy_somy.pdf",width=6, height=8)
ggplot(data=cor.het_somy, aes(x=ord, y=R.cor)) +
  geom_point(colour=cor.het_somy$group.col) + 
  geom_text(aes(label=ifelse(cor.het_somy$FDR<= 0.001,cor.het_somy$sample,""),
                colour = cor.het_somy$group.col, hjust=-0.1, vjust=0.5))  + 
  scale_color_manual(values=unique(sample.info[order(groups),group.col]), guide=FALSE) +
  xlab("Ordered samples") + ylab("Correlation between chromosomal heterozygosity and somy") + theme_bw() + xlim(0,166)
dev.off()


###
# 
#
het.wins.gg <-melt(het.wins, id.vars = c("chr","w10kb"))
setnames(het.wins.gg, c("variable","value"), c("sample","het.count"))
for (gg in  unique(sample.info[,groups]))
{
  print(gg)
  het.wins.gg[sample %in% sample.info[groups==gg,sample], group:=gg]
  het.wins.gg[sample %in% sample.info[groups==gg,sample], group.col:=unique(sample.info[groups==gg,group.col])]
}
het.wins.gg <- het.wins.gg[order(group)]
het.wins.gg[group %in% c("Ldon1","Ldon2","Ldon3","Ldon4","Ldon5","other_Ldon"), species:="donovani"]
het.wins.gg[group %in% c("Linf1","other_Linf","CH_Linf","CUK_Linf"), species:="infantum"]
sample.info[groups %in% c("Ldon1","Ldon2","Ldon3","Ldon4","Ldon5","other_Ldon"), species:="donovani"]
sample.info[groups %in% c("Linf1","other_Linf","CH_Linf","CUK_Linf"), species:="infantum"]
sample.info[sample %in% names(hets[hets>0.004]), het.cat:="high"]
sample.info[sample %in% names(hets[hets<0.004]), het.cat:="low"]

het.wins.gg$sample <- factor(het.wins.gg$sample, levels=unique(het.wins.gg$sample))


pdf("Hetcount_win_across_samples.pdf",width=17,height=7)
ggplot(het.wins.gg, aes(x=factor(sample), y=het.count)) + 
  geom_boxplot(colour=sample.info[order(species,groups),group.col],
               outlier.colour="black", outlier.shape=16,
               outlier.size=0.6) + 
  facet_grid(.~species, scales = "free_x") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Heterozygous counts across 10kb windows") + ylim(0,109)
dev.off()
#
# pdf("Hetcount_win_across_samples_ylog.pdf",width=17,height=7)
# ggplot(het.wins.gg, aes(x=factor(sample), y=het.count)) + 
#   geom_boxplot(colour=sample.info[order(species,groups),group.col],
#                outlier.colour="black", outlier.shape=16,
#                outlier.size=0.6) + 
#   facet_grid(.~species, scales = "free_x") + theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   xlab("") + ylab("Heterozygous counts across 10kb windows") + ylim(0,109) + scale_y_continuous(trans='log')
# dev.off()
#

# one outlier for EP removed in plot
het.wins.gg[order(-het.count)][1:5,]
# chr w10kb  sample het.count      group group.col  species
# 1:  34    70      EP       146 other_Linf     black infantum
# 2:  31    99  X452BM       109      Ldon3   #8B1C62 donovani
# 3:  31   108 X383WTI       108      Ldon3   #8B1C62 donovani
# 4:  34    71      EP       108 other_Linf     black infantum
# 5:  08    42     MAM       108 other_Linf     black infantum


het.wins.gg[,med.hetcount.samp:=median(het.count), by=.(sample)]
het.wins.gg[,quant75.hetcount.samp:=quantile(het.count, probs = 0.75), by=.(sample)]
het.wins.gg[sample %in% names(hets[hets>0.004]), het.cat:="high"]
het.wins.gg[sample %in% names(hets[hets<0.004]), het.cat:="low"]
#
pdf("Hetcount_win_across_samples_cat.pdf",width=17,height=7)
ggplot(het.wins.gg, aes(x=factor(sample), y=het.count)) +
  geom_boxplot(colour=sample.info[order(species,het.cat,groups),group.col],
               outlier.colour="black", outlier.shape=16,
               outlier.size=0.6) +
  facet_grid(.~species*het.cat, scales = "free_x", space = "free") + theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("") + ylab("Heterozygous counts across 10kb windows") + ylim(0,109)
dev.off()


# median het.counts by sample
unique(het.wins.gg[,.(sample, med.hetcount.samp, het.cat)])[med.hetcount.samp>0]
# sample med.hetcount.samp het.cat
# 1:     CH34                 1    high
# 2:   GILANI                 1    high
# 3:  LRC.L53                 2    high
# 4:   Inf152                 2    high
# 5:   ISS174                 1    high
# 6:  ISS2426                 1    high
# 7:  ISS2429                 1    high
# 8:   Inf055                 1    high
# 9: LRC.L740                 1    high
# 10:       EP                 4    high
# 11:      MAM                14    high
table(unique(het.wins.gg[,.(sample, med.hetcount.samp, het.cat)])[med.hetcount.samp>0][,.(med.hetcount.samp, het.cat)])
het.cat
# med.hetcount.samp high
# 1     7
# 2     2
# 4     1
# 14    1

unique(het.wins.gg[,.(sample, group, quant75.hetcount.samp, het.cat),])[quant75.hetcount.samp>0]

table(unique(het.wins.gg[,.(sample, group, quant75.hetcount.samp, het.cat),])[quant75.hetcount.samp>0][,.(quant75.hetcount.samp, het.cat)])
het.cat
# quant75.hetcount.samp high low
# 1    14   2
# 2     8   0
# 3     3   0
# 4     1   0
# 7     1   0
# 12    1   0
# 13    1   0
# 15    1   0
# 20    1   0
#-------------
# sums  31  2


#----------------------------------------------
#
# get called SNP positions for each sample
# use data
data<-data.table(data)

dir.create("het_pos_per-sample", showWarnings = F)
# for heterozygosous sample where allelefreqs are printed
sub<-data[,.(MAM,EP,CH32,CH34,BPK157A1,GILANI,Malta33,X364SPWTII,SUKKAR2,BUMM3,GE,LEM3472,LRC.L740,LRC.L53,ISS2429,ISS2426,ISS174,Inf055,Inf152,BPK512A1)]
# sub<-data

for (col in 1:ncol(sub))
{
  print(col)
  # genotype info of the respective sample
  cc<-as.character(unlist(sub[,col, with=F]))
  # extract genotypes
  aa<-unlist(lapply(cc, function(x) strsplit(x,split=":")[[1]][1]))
  
  # data[(grepl("1",aa) | grepl("2",aa)) & grepl("0",aa), ISS2429]
  # aa[(grepl("1",aa) | grepl("2",aa)) & grepl("0",aa)]
  # unique(aa[(grepl("1",aa) | grepl("2",aa)) & grepl("0",aa)])
  
  write.table(data[(grepl("1",aa) | grepl("2",aa)) & grepl("0",aa), .(CHR, POS.CHR)], 
              file=paste0("het_pos_per-sample/SNP_pos_",colnames(sub)[col],".txt"), quote=FALSE, row.names=F)
  print(sum((grepl("1",aa) | grepl("2",aa)) & grepl("0",aa)))
}


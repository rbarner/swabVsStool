#### modified from script written by Anthony Fodor (https://github.com/afodor/metagenomicsTools/blob/master/src/scripts/vanderbilt/pcoa.txt) ######
rm(list=ls())
library("vegan")
taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("../microbialClassfications/")
  inFileName <- paste("rdpClassifications_",taxa, "Level.RData", sep ="")
  myT <-readRDS(inFileName)
  myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  myPCOA <- capscale(myTLogged~1,distance="bray")

  setwd("../mds/")
  saveRDS(myPCOA$CA$u,file=paste("mds_", taxa, "_logged.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("eigenValues_", taxa, "_logged.RData", sep=""))
}

taxaLevels <- c("otu")
for(taxa in taxaLevels )
{
  setwd("../microbialClassfications/")
  inFileName <- paste("qiime_",taxa, "s.RData", sep ="")
  myT <-readRDS(inFileName)
  myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  myPCOA <- capscale(myTLogged~1,distance="bray")
  
  setwd("../mds/")
  saveRDS(myPCOA$CA$u,file=paste("16S/mds_", taxa, "_logged.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("eigenValues_", taxa, "_logged.RData", sep=""))
}



wgsLevels <- c("Families","Pathways_level4",
               "Pathways_level3",
               "Pathways_level2",
               "Pathways_level1",
               "MetabolicPathways_level2",
               "MetabolicPathways_level3")
for(wgs in wgsLevels )
{
  setwd("../metagenomeFunctions/")
  inFileName <- paste("kegg",wgs, ".RData", sep ="")
  myT <-readRDS(inFileName)
  myT <- myT[-42,]
  myTLogged <- log10((myT*175000) +1)
  myPCOA <- capscale(myTLogged~1,distance="euclidean")
  
  setwd("../mds/")
  saveRDS(myPCOA$CA$u,file=paste("WGS/mds_", wgs, "_loggedFiltered.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("eigenValues_", wgs, "_loggedFiltered.RData", sep=""))
}
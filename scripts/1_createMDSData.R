rm(list=ls())
library("vegan")
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")

################################### RDP classifications/ Qiime ################################################################
taxaLevels <- c("phylum","class","order","family","genus")
tools <- c("rdpClassifications","qiime")
for(tool in tools)
{
  for(taxa in taxaLevels )
  {
    setwd("data/microbialClassifications/")
    inFileName <- paste(tool,"_",taxa, "Level.txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header = TRUE, row.names=1)
    myT <- t(myT)
    myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
    #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    myPCOA <- capscale(myTLogged~1,distance="bray")
    
    setwd("../../mds/")
    saveRDS(myPCOA$CA$u,file=paste(tool,"_mds_", taxa, "_loggedFiltered.RData",sep=""))
    saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste(tool,"_eigenValues_", taxa, "_loggedFiltered.RData", sep=""))
    setwd("..")
  }
}
###################################################### OTUs #############################################################33
setwd("data/microbialClassifications/")
inFileName <- paste("qiime_otuLevel.txt", sep ="")
print(inFileName)
myT <-read.delim(inFileName,header=TRUE,row.names=1)
myT <- t(myT)
myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
myPCOA <- capscale(myTLogged~1,distance="bray")
  
setwd("../../mds/")
saveRDS(myPCOA$CA$u,file=paste("qiime_mds_otu_loggedFiltered.RData",sep=""))
saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("qiime_eigenValues_otu_loggedFiltered.RData", sep=""))
setwd("..")

###################################################### OTUs Rarefied #############################################################33
setwd("data/microbialClassifications/")
inFileName <- paste("qiime_otuRarefiedLevel.txt", sep ="")
print(inFileName)
myT <-read.delim(inFileName,header=TRUE,row.names=1)
myT <- t(myT)
myTLogged <- log10(myT +1)
myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
myPCOA <- capscale(myTLogged~1,distance="bray")

setwd("../../mds/")
saveRDS(myPCOA$CA$u,file=paste("qiime_mds_otuRarefied_loggedFiltered.RData",sep=""))
saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("qiime_eigenValues_otuRarefied_loggedFiltered.RData", sep=""))
setwd("..")

taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied")
tools <- c("qiime")
for(tool in tools)
{
  for(taxa in taxaLevels )
  {
    setwd("data/microbialClassifications/")
    inFileName <- paste(tool,"_",taxa, "Level.txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header = TRUE, row.names=1)
    myT <- t(myT)
    myTLogged <- log10(myT +1)
    #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    myPCOA <- capscale(myTLogged~1,distance="bray")
    
    setwd("../../mds/")
    saveRDS(myPCOA$CA$u,file=paste(tool,"_mds_", taxa, "_loggedFiltered.RData",sep=""))
    saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste(tool,"_eigenValues_", taxa, "_loggedFiltered.RData", sep=""))
    setwd("..")
  }
}


############################################ Metaphlan ############################################
taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
  inFileName <- paste("metaphlan_",taxa, "Level.txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  myTLogged <- log10((myT*175000) +1)
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.5]
  myPCOA <- capscale(myTLogged~1,distance="euclidean")
  
  setwd("../../mds/")
  saveRDS(myPCOA$CA$u,file=paste("metaphlan_mds_", taxa, "_loggedFiltered.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("metaphlan_eigenValues_", taxa, "_loggedFiltered.RData", sep=""))
  setwd("..")
}

############################################ Kraken ############################################
taxaLevels <- c("phylum","class","order","family","genus","species")
tools <- c("kraken16S","krakenWGS")
for(tool in tools)
{
  for(taxa in taxaLevels )
  {
    setwd("data/microbialClassifications/")
    inFileName <- paste(tool,"_",taxa, "Level.txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header=TRUE, row.names=1)
    myT <- t(myT)
    myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
    myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.5]
    myPCOA <- capscale(myTLogged~1,distance="euclidean")
    
    setwd("../../mds/")
    saveRDS(myPCOA$CA$u,file=paste(tool,"_mds_", taxa, "_loggedFiltered.RData",sep=""))
    saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste(tool,"_eigenValues_", taxa, "_loggedFiltered.RData", sep=""))
    setwd("..")
  }
}
############################################ WGS/PICRUsT functions ############################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
tools <- c("wgs","picrust")
for(tool in tools)
{
  for(wgs in wgsLevels )
  {
    setwd("data/metagenomeFunctions/")
    inFileName <- paste(tool,"_",wgs, ".txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header=TRUE, row.names=1)
    myT <- t(myT)
    if(tool %in% "wgs")
    {
      myT <- myT[!row.names(myT) %in% "ST00046",]
    }
    myTLogged <- log10((myT*175000) +1)
    myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    myPCOA <- capscale(myTLogged~1,distance="bray")
    
    setwd("../../mds/")
    saveRDS(myPCOA$CA$u,file=paste(tool,"_mds_", wgs, "_loggedFiltered.RData",sep=""))
    saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste(tool,"_eigenValues_", wgs, "_loggedFiltered.RData", sep=""))
    setwd("..")
  }
}

################################ 1. No tissue ##############################
############################################ Kraken ############################################
taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
  inFileName <- paste("krakenWGS","_",taxa, "Level.txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  samplesToBeIncluded <- row.names(myT)[!startsWith(row.names(myT),"FT")]
  myT <- myT[row.names(myT) %in% samplesToBeIncluded,]
  myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.5]
  myPCOA <- capscale(myTLogged~1,distance="euclidean")
  
  setwd("../../mds/")
  saveRDS(myPCOA$CA$u,file=paste("krakenWGSNoTissue","_mds_", taxa, "_loggedFiltered.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("krakenWGSNoTissue","_eigenValues_", taxa, "_loggedFiltered.RData", sep=""))
  setwd("..")
}

############################################ WGS/PICRUsT functions ############################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
for(wgs in wgsLevels )
{
  setwd("data/metagenomeFunctions/")
  inFileName <- paste("wgs_",wgs, ".txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  samplesToBeIncluded <- row.names(myT)[!startsWith(row.names(myT),"FT")]
  myT <- myT[row.names(myT) %in% samplesToBeIncluded,]
  myT <- myT[!row.names(myT) %in% "ST00046",]
  myTLogged <- log10((myT*175000) +1)
  myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  myPCOA <- capscale(myTLogged~1,distance="bray")
  
  setwd("../../mds/")
  saveRDS(myPCOA$CA$u,file=paste("wgsNoTissue_mds_", wgs, "_loggedFiltered.RData",sep=""))
  saveRDS(myPCOA$CA$eig/sum(myPCOA$CA$eig),file=paste("wgs_eigenValues_", wgs, "NoTissue_loggedFiltered.RData", sep=""))
  setwd("..")
}

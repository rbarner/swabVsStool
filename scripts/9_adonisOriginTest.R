rm(list=ls())
library("vegan")

################################### RDP classifications/ Qiime ################################################################
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData$visit <- unlist(strsplit(as.character(sampleData$type),split = "_"))[c(FALSE,TRUE)]
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied","otuRarefied")
#tools <- c("rdpClassifications","qiime")
tools <- c("qiime")
for(tool in tools)
{
  pValOriginList <- numeric(0);
  pValParticipantList <- numeric(0);
  pValBetaDispList <- numeric(0);
  pValBetaDispParticipantList <- numeric(0);
  for(taxa in taxaLevels )
  {
    setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
    inFileName <- paste(tool,"_",taxa, "Level.txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header = TRUE, row.names=1)
    myT <- myT[,order(colnames(myT))]
    myT <- t(myT)
    myT <- myT[rownames(myT) %in% rownames(sampleData),]
    sampleData3 <- sampleData[rownames(sampleData) %in% rownames(myT),]
    #myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
    myTLogged <- log10((myT) +1)
    myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    
    distTest <- vegdist(myTLogged,method="bray")
    testAdon <- adonis(distTest~Origin,data=sampleData3)
    testAdonParticipant <- adonis(distTest~study_id,data=sampleData3)
    pValOrigin <- pf(testAdon$aov.tab$F.Model[1],testAdon$aov.tab$Df[1],testAdon$aov.tab$Df[2],lower.tail = FALSE);
    pValParticipant <- pf(testAdonParticipant$aov.tab$F.Model[1],testAdonParticipant$aov.tab$Df[1],testAdonParticipant$aov.tab$Df[2],lower.tail = FALSE);
    
    
    testBetadisp <- betadisper(distTest,sampleData3$Origin)
    anovaBetaDisp <- anova(testBetadisp)
    
    testBetadispParticipant <- betadisper(distTest,sampleData3$study_id)
    anovaBetaDispParticipant <- anova(testBetadispParticipant)
    
    #pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
    pValBetaDisp <- pf(anova(testBetadisp)$"F value"[1],anova(testBetadisp)$"Df"[1],anova(testBetadisp)$"Df"[2],lower.tail = FALSE);
    pValBetaDispParticipant <- pf(anova(testBetadispParticipant)$"F value"[1],anova(testBetadispParticipant)$"Df"[1],anova(testBetadispParticipant)$"Df"[2],lower.tail = FALSE);
    
    
    pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
    pValParticipantList[[length(pValParticipantList)+1]]<- pValParticipant;
    pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
    pValBetaDispParticipantList[[length(pValBetaDispParticipantList)+1]]<- pValBetaDispParticipant;
  }
  makeTable <- data.frame(taxaLevels,pValOriginList,pValParticipantList,pValBetaDispList,pValBetaDispParticipantList);
  write("TaxaLevel\tOrigin p-value\tParticipant p-value\tBeta Dispersion Origin\tBeta Dispersion Participant",paste("../../statisticalModels/adonis_",tool,".txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/adonis_",tool,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
}
###################################################### OTUs #############################################################
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
inFileName <- paste("qiime_otuLevel.txt", sep ="")
print(inFileName)
pValOriginList <- numeric(0);
pValBetaDispList <- numeric(0); 
myT <-read.delim(inFileName,header=TRUE,row.names=1)
myT <- t(myT)
myT <- myT[rownames(myT) %in% rownames(sampleData),]
myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)

distTest <- vegdist(myTLogged,method="bray")
testAdon <- adonis(distTest~Origin+visit,data=sampleData)
testBetadisp <- betadisper(distTest,sampleData$Origin)
anovaBetaDisp <- anova(testBetadisp)

pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]

pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;

#setwd("../../statisticalModels/")
makeTable <- data.frame(c("otu"),pValOriginList,pValBetaDispList);
write("TaxaLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_qiime_otu.txt",sep=""));
write.table(makeTable,paste("../../statisticalModels/adonis_qiime_otu.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
#setwd("..")


############################################ Metaphlan (no tissue) ############################################
taxaLevels <- c("phylum","class","order","family","genus")
pValOriginList <- numeric(0);
pValBetaDispList <- numeric(0); 
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
  inFileName <- paste("metaphlan_",taxa, "Level.txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  samplesToBeIncluded <- row.names(myT)[!startsWith(row.names(myT),"FT")]
  myT <- myT[row.names(myT) %in% samplesToBeIncluded,]
  sampleData3 <- sampleData2[row.names(sampleData2) %in% samplesToBeIncluded,]
  myT <- myT[rownames(myT) %in% rownames(sampleData3),]
  myTLogged <- log10((myT*175000)+1)
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.85]
  distTest <- vegdist(myTLogged,method="euclidean")
  testAdon <- adonis(distTest~Origin,data=sampleData3)
  testBetadisp <- betadisper(distTest,sampleData3$Origin)
  anovaBetaDisp <- anova(testBetadisp)
  
  pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
  pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
  
  pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
  pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
}
#setwd("../../statisticalModels/")
makeTable <- data.frame(taxaLevels,pValOriginList,pValBetaDispList);
write("TaxaLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_metaphlan_noTissue.txt",sep=""));
write.table(makeTable,paste("../../statisticalModels/adonis_metaphlan_noTissue.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
#setwd("..")

############################################ Metaphlan ############################################
taxaLevels <- c("phylum","class","order","family","genus")
pValOriginList <- numeric(0);
pValBetaDispList <- numeric(0); 
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
  inFileName <- paste("metaphlan_",taxa, "Level.txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  myT <- myT[rownames(myT) %in% rownames(sampleData2),]
  myTLogged <- log10((myT*175000)+1)
  myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.85]
  distTest <- vegdist(myTLogged,method="euclidean")
  testAdon <- adonis(distTest~Origin,data=sampleData2)
  testBetadisp <- betadisper(distTest,sampleData2$Origin)
  anovaBetaDisp <- anova(testBetadisp)
  
  pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
  pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
  
  pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
  pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
}
#setwd("../../statisticalModels/")
makeTable <- data.frame(taxaLevels,pValOriginList,pValBetaDispList);
write("TaxaLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_metaphlan.txt",sep=""));
write.table(makeTable,paste("../../statisticalModels/adonis_metaphlan.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
#setwd("..")

############################################ Kraken ############################################
taxaLevels <- c("phylum","class","order","family","genus","species")
tools <- c("kraken16S","krakenWGS")
for(tool in tools)
{
  pValOriginList <- numeric(0);
  pValBetaDispList <- numeric(0); 
  for(taxa in taxaLevels )
  {
    setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
    inFileName <- paste(tool,"_",taxa, "Level.txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header=TRUE, row.names=1)
    myT <- t(myT)
    if(tool %in% "kraken16S")
    {
      myT <- myT[rownames(myT) %in% rownames(sampleData),]
    }else{
      myT <- myT[rownames(myT) %in% rownames(sampleData2),]
    }
    myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
    myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    distTest <- vegdist(myTLogged,method="euclidean")
    if(tool %in% "kraken16S")
    {
      testAdon <- adonis(distTest~Origin, data=sampleData)
      testBetadisp <- betadisper(distTest,sampleData$Origin)
    }else{
      testAdon <- adonis(distTest~Origin, data=sampleData2)
      testBetadisp <- betadisper(distTest,sampleData2$Origin)
    }
    anovaBetaDisp <- anova(testBetadisp)
    
    pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
    pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
    
    pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
    pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
  }
  #setwd("../../statisticalModels/")
  makeTable <- data.frame(taxaLevels,pValOriginList,pValBetaDispList);
  write("TaxaLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_",tool,".txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/adonis_",tool,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  #setwd("..")
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
  pValOriginList <- numeric(0);
  pValBetaDispList <- numeric(0); 
  for(wgs in wgsLevels )
  {
    setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/metagenomeFunctions/")
    inFileName <- paste(tool,"_",wgs, ".txt", sep ="")
    print(inFileName)
    myT <-read.delim(inFileName,header=TRUE, row.names=1)
    myT <- t(myT)
    if(tool %in% "wgs")
    {
      myT <- myT[!row.names(myT) %in% "ST00046",]
      myT <- myT[rownames(myT) %in% rownames(sampleData2),]
      sampleData3 <- sampleData2[rownames(sampleData2) %in% rownames(myT),]
      myTLogged <- log10((myT*175000)+1)
      
    }else{
      myT <- myT[rownames(myT) %in% rownames(sampleData),]
      myTLogged <- log10((myT/rowSums(myT))*(sum(colSums(myT))/dim(myT)[1])+1)
      sampleDataFiltered <- sampleData[rownames(sampleData) %in% rownames(myT),]
    }
    myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    
    distTest <- vegdist(myTLogged,method="bray")
    if(tool %in% "wgs"){
      testAdon <- adonis(distTest~Origin, data=sampleData3)
      testBetadisp <- betadisper(distTest,sampleData3$Origin)
    }else{
      testAdon <- adonis(distTest~Origin,data=sampleDataFiltered)
      testBetadisp <- betadisper(distTest,sampleDataFiltered$Origin)
    }
    anovaBetaDisp <- anova(testBetadisp)
    
    pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
    pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
    
    pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
    pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
  }
  #setwd("../../statisticalModels/")
  makeTable <- data.frame(wgsLevels,pValOriginList,pValBetaDispList);
  write("wgsLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_",tool,".txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/adonis_",tool,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  #setwd("..")
}

################################ 1. No tissue ##############################
############################################ Kraken ############################################
taxaLevels <- c("phylum","class","order","family","genus","species")
pValOriginList <- numeric(0);
pValBetaDispList <- numeric(0); 
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/microbialClassifications/")
  inFileName <- paste("krakenWGS","_",taxa, "Level.txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  samplesToBeIncluded <- row.names(myT)[!startsWith(row.names(myT),"FT")]
  myT <- myT[row.names(myT) %in% samplesToBeIncluded,]
  sampleData3 <- sampleData2[row.names(sampleData2) %in% samplesToBeIncluded,]
  myT <- myT[rownames(myT) %in% rownames(sampleData3),]
  myTLogged <- log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  distTest <- vegdist(myTLogged,method="bray")
  testAdon <- adonis(distTest~Origin,data=sampleData3)
  testBetadisp <- betadisper(distTest,sampleData3$Origin)
  anovaBetaDisp <- anova(testBetadisp)
  
  pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
  pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
  
  pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
  pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
}
#setwd("../../statisticalModels/")
makeTable <- data.frame(taxaLevels,pValOriginList,pValBetaDispList);
write("TaxaLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_krakenWGS_noTissue.txt",sep=""));
write.table(makeTable,paste("../../statisticalModels/adonis_krakenWGS_noTissue.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
#setwd("..")

############################################ WGS/PICRUsT functions ############################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
pValOriginList <- numeric(0);
pValBetaDispList <- numeric(0); 
for(wgs in wgsLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/data/metagenomeFunctions/")
  inFileName <- paste("wgs_",wgs, ".txt", sep ="")
  print(inFileName)
  myT <-read.delim(inFileName,header=TRUE, row.names=1)
  myT <- t(myT)
  myT <- myT[!row.names(myT) %in% "ST00046",]
  samplesToBeIncluded <- row.names(myT)[!startsWith(row.names(myT),"FT")]
  sampleData3 <- sampleData2[row.names(sampleData2) %in% samplesToBeIncluded,]
  myT <- myT[row.names(myT) %in% samplesToBeIncluded,]
  myT <- myT[rownames(myT) %in% rownames(sampleData3),]
  myTLogged <- log10((myT*175000) +1)
  myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  
  distTest <- vegdist(myTLogged,method="bray")
  testAdon <- adonis(distTest~Origin,data=sampleData3)
  testBetadisp <- betadisper(distTest,sampleData3$Origin)
  anovaBetaDisp <- anova(testBetadisp)
  
  pValOrigin <- testAdon$aov.tab$`Pr(>F)`[1]
  pValBetaDisp <- anovaBetaDisp$`Pr(>F)`[1]
  
  pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
  pValBetaDispList[[length(pValBetaDispList)+1]]<- pValBetaDisp;
}
#setwd("../../statisticalModels/")
makeTable <- data.frame(wgsLevels,pValOriginList,pValBetaDispList);
write("wgsLevel\tOrigin p-value\tBeta Dispersion",paste("../../statisticalModels/adonis_wgs_noTissue.txt",sep=""));
write.table(makeTable,paste("../../statisticalModels/adonis_wgs_noTissue.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);


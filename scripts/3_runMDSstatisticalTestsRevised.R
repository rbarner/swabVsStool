library(nlme)
library(matrixStats)

####### Functions used #######
mdsStatisticalModels <- function(dataSet,eigen,tool,level)
{
  bacteriaSwab <- split(dataSet,dataSet$Origin)$SWAB
  bacteriaStool <- split(dataSet,dataSet$Origin)$STOOL
  
  stoolMeans <- round(colMeans(bacteriaStool[,37:51]),3); 
  swabMeans <- round(colMeans(bacteriaSwab[,37:51]),3); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,37:51])),3); 
  swabSDs <- round(colSds(as.matrix(bacteriaSwab[,37:51])),3); 
  
  pValOriginList=numeric(0);
  pValIndividualList=numeric(0);
  for(i in 1:15)
  {
    f <- as.formula(paste("MDS",i,"~","Origin","+","visit",sep=""));
    simpleMod <- gls(f,method="REML",data=dataSet);
    mixedMod <- lme(f,method="REML",random=~1|study_id,data=dataSet);
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValParticipant <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)

    pValOriginList[[length(pValOriginList)+1]]<- format(pValOrigin,digits=3);
    pValIndividualList[[length(pValIndividualList)+1]]<- format(pValParticipant,digits=3);
  }
  originAdj <- p.adjust(pValOriginList,method = "BH")
  individualAdj <- p.adjust(pValIndividualList,method = "BH")
  makeTable=data.frame(eigen[1:15]*100,stoolMeans,stoolSDs,swabMeans,swabSDs,originAdj,individualAdj);
  names(makeTable) <- cbind("stool mean","stool sd","swab mean","swab sd","Origin adj p-value","Participant adj p-value");
  write("MDS Axis\t% variation explained\tstoolMeans\tstoolSDs\tswabMeans\tswabSDs\tOrigin adj p-value\tParticipant adj p-value",paste("../statisticalModels/3_mds_",taxa,"_",classifier,"_individual_origin_pVal.txt",sep=""));
  write.table(makeTable,paste("../statisticalModels/3_mds_",level,"_",tool,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
}

############ MDS classifications #######################
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData$visit <- unlist(strsplit(as.character(sampleData$type),split = "_"))[c(FALSE,TRUE)]
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

classifierList <- c("krakenWGS","krakenWGSNoTissue","kraken16S","rdpClassifications", "qiime","metaphlan")
classifierList <- c("qiime","krakenWGSNoTissue","kraken16S","rdpClassifications")

for(classifier in classifierList)
{
  if(classifier %in% c("qiime"))
  {
    taxaLevels <- c("phylum","phylumRarefied","class","classRarefied","order","orderRarefied","family","familyRarefied","genus","genusRarefied","otu","otuRarefied")
    
  }else{
    taxaLevels <- c("phylum","class","order","family","genus")
  }
  for(taxa in taxaLevels )
  {
    setwd("mds")
    mdsFile <- paste(classifier,"_mds_", taxa, "_loggedFiltered.RData",sep="");
    print(mdsFile)
    eigenFile <- paste(classifier,"_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    if(classifier %in% c("krakenWGS","metaphlan"))
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile)
    
    mdsStatisticalModels(dataSet = mdsMeta,eigen = eigen,tool = classifier,level=taxa);
    setwd("..") 
  }
}

setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
functionList <- c("wgs","picrust")
for(funct in functionList)
{
  if(funct %in% c("wgs"))
  {
    wgsLevels <- c("keggFamilies",
                   "keggPathwaysLevel3",
                   "keggPathwaysLevel2",
                   "keggPathwaysLevel1",
                   "metabolickeggPathwaysLevel2",
                   "metabolickeggPathwaysLevel3",
                   "keggFamiliesNoTissue",
                   "keggPathwaysLevel3NoTissue",
                   "keggPathwaysLevel2NoTissue",
                   "keggPathwaysLevel1NoTissue",
                   "metabolickeggPathwaysLevel2NoTissue",
                   "metabolickeggPathwaysLevel3NoTissue")
  }else{
    wgsLevels <- c("keggFamilies",
                   "keggPathwaysLevel3",
                   "keggPathwaysLevel2",
                   "keggPathwaysLevel1",
                   "metabolickeggPathwaysLevel2",
                   "metabolickeggPathwaysLevel3")
  }
  for(wgs in wgsLevels )
  {
    setwd("mds")
    mdsFile <- paste(funct,"_mds_", wgs, "_loggedFiltered.RData",sep="");
    print(mdsFile)
    eigenFile <- paste(funct,"_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    #sampleData  <- sampleData[-26,]
    if(funct %in% "wgs")
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }
    else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    mdsStatisticalModels(mdsMeta,mds,eigen,funct,wgs);
    setwd("..")
  }
}
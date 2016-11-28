############ MDS classifications #######################
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

taxaLevels <- c("phylum","class","order","family","genus")
classifierList <- c("rdpClassifications", "qiime","metaphlan")
for(classifier in classifierList)
{
  for(taxa in taxaLevels )
  {
    setwd("mds")
    mdsFile <- paste(classifier,"_mds_", taxa, "_loggedFiltered.RData",sep="");
    eigenFile <- paste(classifier,"_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    if(classifier %in% "metaphlan")
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }
    else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    
    pValOriginList=numeric(0);
    pValIndividualList=numeric(0);
    for(i in 1:8)
    {
      f<-as.formula(paste("MDS",i,"~","Origin",sep=""));
      aovtest <- aov(f,mdsMeta);
      pValOrigin <- anova(aovtest)$"Pr(>F)"[1];
      pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
      
      g<-as.formula(paste("MDS",i,"~","study_id",sep = ""));
      aovtest <- aov(g,mdsMeta);
      pValIndividual <- anova(aovtest)$"Pr(>F)"[1];
      pValIndividualList[[length(pValIndividualList)+1]]<- pValIndividual;
    }
    originAdj <- p.adjust(pValOriginList,method = "BH")
    individualAdj <- p.adjust(pValIndividualList,method = "BH")
    makeTable=data.frame(eigen[1:8]*100,originAdj,individualAdj);
    names(makeTable)=cbind("Origin adj p-value","Individual adj p-value");
    write("MDS Axis\t% variation explained\tOrigin adj p-value\tIndividual adj p-value",paste("../statisticalModels/3_mds_",taxa,"_",classifier,"_individual_origin_pVal.txt",sep=""));
    write.table(makeTable,paste("../statisticalModels/3_mds_",taxa,"_",classifier,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
    setwd("..")
    }
}

##################################### Functions ################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
functionList <- c("wgs","picrust")

for(funct in functionList)
{
  for(wgs in wgsLevels )
  {
    setwd("mds")
    mdsFile <- paste(funct,"_mds_", wgs, "_loggedFiltered.RData",sep="");
    print(mdsFile)
    eigenFile <- paste(funct,"_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    if(funct %in% "wgs")
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }
    else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    
    pValOriginList=numeric(0);
    pValIndividualList=numeric(0);
    for(i in 1:8)
    {
      f<-as.formula(paste("MDS",i,"~","Origin",sep=""));
      aovtest <- aov(f,mdsMeta);
      pValOrigin <- anova(aovtest)$"Pr(>F)"[1];
      pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
      
      g<-as.formula(paste("MDS",i,"~","study_id",sep = ""));
      aovtest2 <- aov(g,mdsMeta);
      pValIndividual <- anova(aovtest2)$"Pr(>F)"[1];
      pValIndividualList[[length(pValIndividualList)+1]]<- pValIndividual;
    }
    originAdj <- p.adjust(pValOriginList,method = "BH")
    individualAdj <- p.adjust(pValIndividualList,method = "BH")
    makeTable=data.frame(eigen[1:8]*100,originAdj,individualAdj);
    names(makeTable)=cbind("Origin adj p-value","Individual adj p-value");
    write("MDS Axis\t% variation explained\tOrigin adj p-value\tIndividual adj p-value",paste("../statisticalModels/3_mds_",wgs,"_",funct,"_individual_origin_pVal.txt",sep=""));
    write.table(makeTable,paste("../statisticalModels/3_mds_",wgs,"_",funct,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
    setwd("..")
    }
}
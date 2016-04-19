############ functions used #######################
setwd("C://Users//Roshonda/PPCCT_swabstool/data/mds/")
sampleData=readRDS("../key/mapping_key_16S.RData");
taxaLevels <- c("phylum","class","order","family","genus","otu")

for(taxa in taxaLevels )
{
  mdsFile <- paste(  "mds_", taxa, "_logged.RData",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, "_logged.RData",sep="");
  
  mds <-readRDS(mdsFile);
  mdsMeta <- cbind(sampleData,mds);
  eigen <-readRDS(eigenFile);

  pValOriginList=numeric(0);
  pValIndividualList=numeric(0);
  for(i in (dim(sampleData)[2]+1): (dim(sampleData)[2]+15))
  {
    f<-as.formula(paste(names(mdsMeta)[i],"~","Origin"));
    ttest <- t.test(f,mdsMeta);
    pValOrigin <- ttest$p.value[[1]];
    pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
    
    g<-as.formula(paste(names(mdsMeta)[i],"~","study_id"));
    aovtest <- aov(g,mdsMeta);
    pValIndividual <- anova(aovtest)$"Pr(>F)"[1];
    pValIndividualList[[length(pValIndividualList)+1]]<- pValIndividual;
  }
  originAdj <- p.adjust(pValOriginList,method = "BH")
  individualAdj <- p.adjust(pValIndividualList,method = "BH")
  makeTable=data.frame(eigen[1:15]*100,originAdj,individualAdj);
  names(makeTable)=cbind("Origin adj p-value","Individual adj p-value");
  write("MDS Axis\t% variation explained\tOrigin adj p-value\tIndividual adj p-value",paste("../../statisticalModels/mds_",taxa,"_individual_origin_pVal.txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/mds_",taxa,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
}


sampleData=readRDS("../key/mapping_key_WGS.RData");
sampleData <- sampleData[-42,]
wgsLevels <- c("Families","Pathways_level4",
               "Pathways_level3",
               "Pathways_level2",
               "Pathways_level1",
               "MetabolicPathways_level2",
               "MetabolicPathways_level3")

for(wgs in wgsLevels )
{
  mdsFile <- paste(  "mds_", wgs, "_loggedFiltered.RData",sep="");
  eigenFile <- paste(  "eigenValues_", wgs, "_loggedFiltered.RData",sep="");
  
  mds <-readRDS(mdsFile);
  mdsMeta <- cbind(sampleData,mds);
  eigen <-readRDS(eigenFile);
  
  pValOriginList=numeric(0);
  pValIndividualList=numeric(0);
  for(i in (dim(sampleData)[2]+1): length(names(mdsMeta)))
  {
    f<-as.formula(paste(names(mdsMeta)[i],"~","type"));
    aovtest <- aov(f,mdsMeta);
    pValOrigin <- anova(aovtest)$"Pr(>F)"[1];
    pValOriginList[[length(pValOriginList)+1]]<- pValOrigin;
    
    g<-as.formula(paste(names(mdsMeta)[i],"~","study_id"));
    aovtest2 <- aov(g,mdsMeta);
    pValIndividual <- anova(aovtest2)$"Pr(>F)"[1];
    pValIndividualList[[length(pValIndividualList)+1]]<- pValIndividual;
  }
  originAdj <- p.adjust(pValOriginList,method = "BH")
  individualAdj <- p.adjust(pValIndividualList,method = "BH")
  makeTable=data.frame(eigen*100,originAdj,individualAdj);
  names(makeTable)=cbind("Origin adj p-value","Individual adj p-value");
  write("MDS Axis\t% variation explained\tOrigin adj p-value\tIndividual adj p-value",paste("../../statisticalModels/mds_",wgs,"_individual_origin_pVal.txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/mds_",wgs,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
}
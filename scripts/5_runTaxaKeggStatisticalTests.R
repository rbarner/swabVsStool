library(nlme)

############################################ Classifications ######################
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]

taxaLevels <- c("phylum","class","order","family","genus","otu")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassfications/")
  if(taxa=="otu")
  {
    inFileName <- paste("qiime_",taxa, "Level.txt", sep ="")
    print(inFileName)
    bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
    bacteria <- t(bacteria)
    bacteriaLogged <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/nrow(bacteria)) +1)
    bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
    colnames(bacteriaLogged) <- paste("OTU",colnames(bacteriaLogged),sep="")
  }
  else
  {
    inFileName <- paste("rdpClassifications_",taxa, "Level.txt", sep ="")
    print(inFileName)
    bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
    bacteria <- t(bacteria)
    bacteriaLogged <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/nrow(bacteria)) +1)
    bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
    colnames(bacteriaLogged) <- gsub(" ","_",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("/","_",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("-","_",colnames(bacteriaLogged))
  }
  
  bacteriaMeta <- merge(sampleData16S,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$SWAB
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$STOOL
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValTimeList <- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  
  for(i in (ncol(sampleData16S)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
 
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,236,lower.tail = FALSE);
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,236,lower.tail = FALSE);
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,swabMeans,pValOriginList,pValTimeList, pValIndividualList,originAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tswabMeans\tOrigin p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,".txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}

#################### Kraken 16S ######################
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]

taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassfications/")
    inFileName <- paste("kraken16S_",taxa, "Level.txt", sep ="")
    print(inFileName)
    bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
    bacteria <- t(bacteria)
    colnames(bacteria)  <- gsub("\\[","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("\\]","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("\\(","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("\\)","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("-","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("\\+","\\.",colnames(bacteria))
    colnames(bacteria)  <- gsub("\\/","\\.",colnames(bacteria))
    bacteriaLogged <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/dim(bacteria)[1]) +1)
    #myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
    #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
    bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
  
  bacteriaMeta <- merge(sampleData16S,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$SWAB
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$STOOL
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValTimeList <- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  
  for(i in (ncol(sampleData16S)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,236,lower.tail = FALSE);
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,236,lower.tail = FALSE);
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,swabMeans,pValOriginList,pValTimeList, pValIndividualList,originAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tswabMeans\tOrigin p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_kraken16S.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_kraken16S.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}


########################################### kraken classifications ##############################
sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"

taxaLevels <- c("phylum","class","order","family","genus","species")
#taxaLevels <- c("species")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassfications/")
  inFileName <- paste("krakenWGS_",taxa, "Level.txt", sep ="")
  print(inFileName)
  bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
  bacteria <- t(bacteria)
  colnames(bacteria)  <- gsub("\\[","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("\\]","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("\\(","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("\\)","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("-","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("\\+","\\.",colnames(bacteria))
  colnames(bacteria)  <- gsub("\\/","\\.",colnames(bacteria))
  bacteriaLogged <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/dim(bacteria)[1]) +1)
  #myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.95]
  
  bacteriaMeta <- merge(sampleDataWGS,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$swab
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$stool
  bacteriaTissue <- split(bacteriaMeta,bacteriaMeta$Origin)$tissue
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValOriginSwabStoolList <- numeric(0);
  pValOriginTissueStoolList <- numeric(0);
  pValOriginTissueSwabList <- numeric(0);
  pValTimeList <- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  tissueMeans <- colMeans(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)]); 
  
  for(i in (ncol(sampleDataWGS)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
    
    anovaModForTukey <- aov(model,data=bacteriaMeta)
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[3,4];
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=bacteriaMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- pValOriginSwabStool;
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- pValOriginTissueStool;
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- pValOriginTissueSwab;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,swabMeans,tissueMeans,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tswabMeans\ttissueMeans\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}



########################################### Metaphlan classifications ##############################
sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"

taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassfications/")
  inFileName <- paste("metaphlan_",taxa, "Level.txt", sep ="")
  print(inFileName)
  bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
  bacteria <- t(bacteria)
  bacteriaLogged <- log10((bacteria*175000) +1)
  #bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.95]
  
  bacteriaMeta <- merge(sampleDataWGS,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$swab
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$stool
  bacteriaTissue <- split(bacteriaMeta,bacteriaMeta$Origin)$tissue
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValOriginSwabStoolList <- numeric(0);
  pValOriginTissueStoolList <- numeric(0);
  pValOriginTissueSwabList <- numeric(0);
  pValTimeList <- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  tissueMeans <- colMeans(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)]); 
  
  for(i in (ncol(sampleDataWGS)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
    
    anovaModForTukey <- aov(model,data=bacteriaMeta)
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[3,4];
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=bacteriaMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- pValOriginSwabStool;
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- pValOriginTissueStool;
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- pValOriginTissueSwab;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,swabMeans,tissueMeans,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tswabMeans\ttissueMeans\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_metaphlan.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_metaphlan.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}


################################################ WGS functions##############################
sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"

wgsLevels <- c("keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")

wgsLevels <- c("keggFamilies") 
for(wgs in wgsLevels )
{
  setwd("data/metagenomeFunctions/")
  inFileName <- paste("wgs_",wgs, ".txt", sep ="")
  print(inFileName)
  kegg <-read.delim(inFileName,header = TRUE, row.names=1)
  kegg <- t(kegg)
  keggLogged <- log10((kegg*175000) +1)
  keggLogged <- keggLogged[,(colSums(keggLogged==0)/nrow(keggLogged))<=0.75]
  keggLogged2 <- keggLogged
  colnames(keggLogged2) <- substr(colnames(keggLogged2),1,6)
  
  keggMeta <- merge(sampleDataWGS,keggLogged2,by="row.names");
  keggSwab <- split(keggMeta,keggMeta$Origin)$swab
  keggStool <- split(keggMeta,keggMeta$Origin)$stool
  keggTissue <- split(keggMeta,keggMeta$Origin)$tissue
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValOriginSwabStoolList <- numeric(0);
  pValOriginTissueStoolList <- numeric(0);
  pValOriginTissueSwabList <- numeric(0);
  pValTimeList <- numeric(0);
  
  stoolMeans <- colMeans(keggStool[,names(keggStool) %in% colnames(keggLogged2)]); 
  swabMeans <- colMeans(keggSwab[,names(keggSwab) %in% colnames(keggLogged2)]); 
  tissueMeans <- colMeans(keggTissue[,names(keggTissue) %in% colnames(keggLogged2)]); 
  
  for(i in (ncol(sampleDataWGS)+2): ncol(keggMeta))
  {
    model <- as.formula(paste(names(keggMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=keggMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=keggMeta,control = lmeControl(singular.ok=TRUE, returnObject = TRUE));
    
    anovaModForTukey <- aov(model,data=keggMeta)
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[3,4];
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=keggMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- pValOriginSwabStool;
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- pValOriginTissueStool;
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- pValOriginTissueSwab;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(keggLogged),stoolMeans,swabMeans,tissueMeans,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("wgs\tstoolMeans\tswabMeans\ttissueMeans\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgs.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgs.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..") 
}


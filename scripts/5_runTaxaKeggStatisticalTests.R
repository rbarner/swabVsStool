library(nlme)
sampleData=readRDS("../key/mapping_key_16S.RData");
taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("../data/microbialClassfications/")
  if(taxa=="otu")
  {
    inFileName <- paste("qiime_",taxa, "s.RData", sep ="")
  }
  else
  {
    inFileName <- paste("rdpClassifications_",taxa, "Level.RData", sep ="")
  }
  bacteria <-readRDS(inFileName)
  bacteria <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/dim(bacteria)[1]) +1)
  colnames(bacteria) <- gsub(" ","_",colnames(bacteria))
  colnames(bacteria) <- gsub("/","_",colnames(bacteria))
  colnames(bacteria) <- gsub("-","_",colnames(bacteria))
  bacteriaMeta <- cbind(sampleData,bacteria);
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$SWAB
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$STOOL
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValTreatmentList <- numeric(0);
  pValTimeList <- numeric(0);
  stoolMeans=colMeans(bacteriaStool[,(dim(sampleData)[2]+1):dim(bacteriaStool)[2]]);
  swabMeans <- colMeans(bacteriaSwab[,(dim(sampleData)[2]+1):dim(bacteriaSwab)[2]]);
  
  for(i in (dim(sampleData)[2]+1):dim(bacteriaMeta)[2])
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","treatment","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
 
    pValOrigin <-summary(aov(mixedMod,data=bacteriaMeta))[[1]][["Pr(>F)"]][1];
    pValTreatment <- summary(aov(mixedMod,data=bacteriaMeta))[[1]][["Pr(>F)"]][2];
    pValTime <- summary(aov(mixedMod,data=bacteriaMeta))[[1]][["Pr(>F)"]][3];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1, lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTreatmentList[[length(pValTreatmentList)+1]] <- pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] <- pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  treatmentAdj <- p.adjust(pValTreatmentList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteria),stoolMeans,swabMeans,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList,originAdj,treatmentAdj,timeAdj,individualAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tswabMeans\tOrigin p-value\tTreatment p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tTreatment (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("marginalStatisticalModels_taxaByTaxa_",taxa,".txt",sep=""));
  write.table(makeTable,paste("marginalStatisticalModels_taxaByTaxa_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
}


sampleData=read.delim("/Users/Roshonda/Dropbox/FodorLab/dataForProposal/chapter4SwabStool/key/mapping_file_wgs_ordered.txt");
#sampleData=read.delim("/Users/Roshonda/Dropbox/FodorLab/dataForProposal/chapter4SwabStool/key/mapping_file_wgs_ordered.txt");
names(sampleData)[2] <- "Origin"
sampleData  <- sampleData[-42,]
#taxaLevels <- c("keggFamilies")
taxaLevels <- c("pathways_level1","pathways_level2","pathways_level3")
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/Dropbox/FodorLab/dataForProposal/chapter4SwabStool/tables/wgs/")
  fileName <- paste(  "shrubsole_WGS_",taxa, "LogFiltered_wTissue.txt", sep ="")
  kegg <-read.delim(fileName,header=TRUE,row.names=1)
  #kegg <-read.table(fileName,header=TRUE,sep="\t")
  kegg <- t(kegg)
  kegg <- kegg[-42,]
  #keggNames  <- kegg[1,]
  #kegg <- kegg[-1,]
  #colnames(kegg) <- substr(colnames(kegg),1,6)
  keggMeta <- cbind(sampleData,kegg);
  keggSwab <- split(keggMeta,keggMeta$Origin)$swab
  keggStool <- split(keggMeta,keggMeta$Origin)$stool
  keggTissue <- split(keggMeta,keggMeta$Origin)$tissue
  
  pValIndividualList=numeric(0);
  pValOriginList=numeric(0);
  pValTreatmentList = numeric(0);
  pValTimeList = numeric(0);
  stoolMeans=colMeans(keggStool[,24:dim(keggStool)[2]]);
  swabMeans <- colMeans(keggSwab[,24:dim(keggSwab)[2]]);
  tissueMeans <- colMeans(keggTissue[,24:dim(keggTissue)[2]]);
  
  for(i in 1:dim(kegg)[2])
  {
    model=as.formula(paste(names(keggMeta)[i+23],"~","Origin","+","treatment","+","visit"));
    print(model);
    simpleMod=lm(model,data=keggMeta);
    simpleMod2=gls(model,method="REML",data=keggMeta);
    mixedMod=lme(model,method="REML",random=~1|study_id,data=keggMeta);
    
    pValOrigin <-summary(aov(mixedMod,data=keggMeta))[[1]][["Pr(>F)"]][1];
    pValTreatment <- summary(aov(mixedMod,data=keggMeta))[[1]][["Pr(>F)"]][2];
    pValTime <- summary(aov(mixedMod,data=keggMeta))[[1]][["Pr(>F)"]][3];
    pValIndividual  <- pchisq(anova(simpleMod2,mixedMod)$"L.Ratio"[2],1, lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]]=pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]]=pValIndividual;
    pValTreatmentList[[length(pValTreatmentList)+1]] = pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] = pValTime;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  treatmentAdj <- p.adjust(pValTreatmentList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  makeTable=data.frame(colnames(kegg),stoolMeans,swabMeans,tissueMeans,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList,originAdj,treatmentAdj,timeAdj,individualAdj);
  #names(makeTable)=cbind("t score","p-value");
  setwd("C://Users//Roshonda//Dropbox//FodorLab//dataForProposal/chapter4SwabStool/statisticalModels/")
  write("KeggFunction\tStoolMeans\tSwabMeans\tTissueMeans\tOrigin p-value\tTreatment p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tTreatment (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("anova_functions_",taxa,"_wMeans.txt",sep=""));
  write.table(makeTable,paste("anova_functions_",taxa,"_wMeans.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
}



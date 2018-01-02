library(matrixStats)
library(nlme)
library(MuMIn)
library(ggplot2)
############################################ Classifications ######################
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]

taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied","otuRarefied")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
  if(taxa=="otuRarefied")
  {
    inFileName <- paste("qiime_",taxa, "Level.txt", sep ="")
    print(inFileName)
    bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
    bacteria <- t(bacteria)
    bacteriaLogged <- log10(bacteria +1)
    bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
    colnames(bacteriaLogged) <- paste("OTU",colnames(bacteriaLogged),sep="")
  }else
  {
    inFileName <- paste("qiime_",taxa, "Level.txt", sep ="")
    print(inFileName)
    bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
    bacteria <- t(bacteria)
    bacteriaLogged <- log10(bacteria +1)
    bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
    colnames(bacteriaLogged) <- gsub(" ","_",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("/","_",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("-","_",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("\\]","\\.",colnames(bacteriaLogged))
    colnames(bacteriaLogged) <- gsub("\\[","\\.",colnames(bacteriaLogged))
  }
  
  bacteriaMeta <- merge(sampleData16S,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$SWAB
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$STOOL
  
  pValttestList <- numeric(0);
  pValWilcoxonList<- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)])),3); 
  swabSDs <- round(colMeans(as.matrix(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)])),3); 
  
  for(i in (ncol(sampleData16S)+2): ncol(bacteriaMeta))
  {
      model<-as.formula(paste(names(bacteriaMeta)[i],"~","Origin"));
      print(model)
      wilcoxtest <- wilcox.test(model,bacteriaMeta);
      ttest <- t.test(model,bacteriaMeta);
      
      pValWilcoxon=wilcoxtest$p.value;
      pValttest=ttest$p.value;
      
      pValttestList[[length(pValttestList)+1]] <- format(pValttest,3);
      pValWilcoxonList[[length(pValWilcoxonList)+1]] <- format(pValWilcoxon,3);  
  }
  ttestAdj <- p.adjust(pValttestList, method = "fdr")
  wilcoxonAdj <- p.adjust(pValWilcoxonList, method = "fdr")
  
  makeTable2=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,-log10(as.numeric(as.character(pValttestList))),-log10(as.numeric(as.character(pValWilcoxonList))), ttestAdj, wilcoxonAdj);
  names(makeTable2) <- c("taxa","stoolMeans","stoolSDs","swabMeans","swabSDs","pValttest","pValWilcoxon","adjPvalttest","adjPvalWilcoxon")
  correlationPvals <- cor(as.numeric(as.character(pValWilcoxonList)),as.numeric(as.character(pValttestList)))
  tiff(paste("../../plots/10_nonparametricVsparametric_pValues_",taxa,"_qiime.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350);
  p <- ggplot(makeTable2)
  print(p + geom_point(aes(pValWilcoxon,pValttest),colour = "black",size = 6) +
          geom_hline(yintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
          geom_vline(xintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
          xlab("-Log10(Wilcoxon p-values)") + ylab("-Log10(T-Test p-values)") +
          ggtitle(expression(atop(taxa, atop(italic(paste("cor= ",correlationPvals,sep="")), ""))))+
          ggtitle(taxa) +
          theme_classic(base_size = 28)+
          theme(plot.title = element_text(hjust = 0.5),
                axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=24),
                text=element_text(face="bold",size=24),
                legend.position="bottom",
                legend.title=element_blank()
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
  graphics.off()
  
  setwd("../../statisticalModels/")
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,pValttestList,pValWilcoxonList,ttestAdj, wilcoxonAdj);
  write("Taxa\tstoolMeans\tstoolSDs\tswabMeans\tswabSDs\tttest p-value\twilcoxon p-value\tTtest (adj p-value)\tWilcoxon (adj p-value)",paste("10_nonparametricVsparametric_",taxa,"_qiimeRarefied.txt",sep=""));
  write.table(makeTable,paste("10_nonparametricVsparametric_",taxa,"_qiimeRarefied.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}

#################### Kraken 16S ######################
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]

taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
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
  bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
  
  bacteriaMeta <- merge(sampleData16S,bacteriaLogged,by="row.names");
  bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$SWAB
  bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$STOOL
  
  pValttestList <- numeric(0);
  pValWilcoxonList<- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)])),3); 
  swabSDs <- round(colMeans(as.matrix(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)])),3); 
  
  for(i in (ncol(sampleData16S)+2): ncol(bacteriaMeta))
  {
    model<-as.formula(paste(names(bacteriaMeta)[i],"~","Origin"));
    print(model)
    wilcoxtest <- wilcox.test(model,bacteriaMeta);
    ttest <- t.test(model,bacteriaMeta);
    
    pValWilcoxon=wilcoxtest$p.value;
    pValttest=ttest$p.value;
    
    pValttestList[[length(pValttestList)+1]] <- format(pValttest,3);
    pValWilcoxonList[[length(pValWilcoxonList)+1]] <- format(pValWilcoxon,3);  
  }
  ttestAdj <- p.adjust(pValttestList, method = "fdr")
  wilcoxonAdj <- p.adjust(pValWilcoxonList, method = "fdr")
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,pValttestList,pValWilcoxonList, ttestAdj, wilcoxonAdj);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tstoolSDs\tswabMeans\tswabSDs\tttest p-value\twilcoxon p-value\tTtest (adj p-value)\tWilcoxon (adj p-value)",paste("10_nonparametricVsparametric_taxaByTaxa_",taxa,"_kraken16S.txt",sep=""));
  write.table(makeTable,paste("10_nonparametricVsparametric_",taxa,"_kraken16S.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}


########################################### kraken WGS ##############################
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
  bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
  
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
  rSquaredMarginalList <- numeric(0);
  rSquaredConditionalList <- numeric(0);
  
  stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
  swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
  tissueMeans <- colMeans(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)]); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)])),3); 
  swabSDs <- round(colSds(as.matrix(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)])),3); 
  tissueSDs <- round(colSds(as.matrix(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)])),3); 
  
  for(i in (ncol(sampleDataWGS)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin"));
    print(model);
    
    anovaModForTukey <- aov(model,data=bacteriaMeta)
    
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValTime <- pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE);
    
    #pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[3,4];
    #pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=bacteriaMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- format.pval(pValOrigin,3);
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- format.pval(pValOriginSwabStool,3);
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- format.pval(pValOriginTissueStool,3);
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- format.pval(pValOriginTissueSwab,3);
    pValIndividualList[[length(pValIndividualList)+1]] <- format.pval(pValIndividual,3);
    pValTimeList[[length(pValTimeList)+1]] <- format.pval(pValTime,3);
    
    rsquaredMarginal <- round(r.squaredGLMM(mixedMod)[1],3)
    rsquaredConditional <- round(r.squaredGLMM(mixedMod)[2],3)
    rSquaredMarginalList[[length(rSquaredMarginalList)+1]] <- rsquaredMarginal;
    rSquaredConditionalList[[length(rSquaredConditionalList)+1]] <- rsquaredConditional;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,tissueMeans,tissueSDs,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj,rSquaredMarginalList, rSquaredConditionalList);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tstoolSD\tswabMeans\tswabSD\ttissueMeans\ttissueSD\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)\trsquared marginal\trsquared conditional",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS_mixedMod.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS_mixedMod.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
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
  rSquaredMarginalList <- numeric(0);
  rSquaredConditionalList <- numeric(0);
  
  stoolMeans <- round(colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]),3); 
  swabMeans <- round(colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]),3); 
  tissueMeans <- round(colMeans(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)]),3); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)])),3); 
  swabSDs <- round(colSds(as.matrix(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)])),3); 
  tissueSDs <- round(colSds(as.matrix(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)])),3);
  for(i in (ncol(sampleDataWGS)+2): ncol(bacteriaMeta))
  {
    model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=bacteriaMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
    
    anovaModForTukey <- aov(model,data=bacteriaMeta)
    
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValTime <- pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE);
    
    #pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[3,4];
    #pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=bacteriaMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- format.pval(pValOrigin,3);
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- format.pval(pValOriginSwabStool,3);
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- format.pval(pValOriginTissueStool,3);
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- format.pval(pValOriginTissueSwab,3);
    pValIndividualList[[length(pValIndividualList)+1]] <- format.pval(pValIndividual,3);
    pValTimeList[[length(pValTimeList)+1]] <- format.pval(pValTime,3);
    
    rsquaredMarginal <- round(r.squaredGLMM(mixedMod)[1],3)
    rsquaredConditional <- round(r.squaredGLMM(mixedMod)[2],3)
    rSquaredMarginalList[[length(rSquaredMarginalList)+1]] <- rsquaredMarginal;
    rSquaredConditionalList[[length(rSquaredConditionalList)+1]] <- rsquaredConditional;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,tissueMeans,tissueSDs,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj,rSquaredMarginalList,rSquaredConditionalList);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tstoolSD\tswabMeans\tswabSD\ttissueMeans\ttissueSD\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)\trsquared marginal\trsquared conditional",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_metaphlan_mixedMod.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_metaphlan_mixedMod.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
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
  rSquaredMarginalList <- numeric(0);
  rSquaredConditionalList <- numeric(0);
  
  stoolMeans <- round(colMeans(keggStool[,names(keggStool) %in% colnames(keggLogged2)]),3); 
  swabMeans <- round(colMeans(keggSwab[,names(keggSwab) %in% colnames(keggLogged2)]),3); 
  tissueMeans <- round(colMeans(keggTissue[,names(keggTissue) %in% colnames(keggLogged2)]),3); 
  
  stoolSDs <- round(colSds(as.matrix(keggStool[,names(keggStool) %in% colnames(keggLogged2)])),3); 
  swabSDs <- round(colSds(as.matrix(keggSwab[,names(keggSwab) %in% colnames(keggLogged2)])),3); 
  tissueSDs <- round(colSds(as.matrix(keggTissue[,names(keggTissue) %in% colnames(keggLogged2)])),3); 
  
  for(i in (ncol(sampleDataWGS)+2): ncol(keggMeta))
  {
    model <- as.formula(paste(names(keggMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=keggMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=keggMeta,control = lmeControl(singular.ok=TRUE, returnObject = TRUE));
    
    anovaModForTukey <- aov(model,data=keggMeta)
    
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValTime <- pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE);
    
    #pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValOriginSwabStool <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[1,4];
    pValOriginTissueStool <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[2,4];
    pValOriginTissueSwab <- TukeyHSD(aov(model,data=keggMeta),"Origin")$Origin[3,4];
    #pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    #pValTime <- TukeyHSD(aov(model,data=keggMeta),"visit")$visit[4];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- format.pval(pValOrigin,3);
    pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- format.pval(pValOriginSwabStool,3);
    pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- format.pval(pValOriginTissueStool,3);
    pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- format.pval(pValOriginTissueSwab,3);
    pValIndividualList[[length(pValIndividualList)+1]] <- format.pval(pValIndividual,3);
    pValTimeList[[length(pValTimeList)+1]] <- format.pval(pValTime,3);
    
    rsquaredMarginal <- round(r.squaredGLMM(mixedMod)[1],3)
    rsquaredConditional <- round(r.squaredGLMM(mixedMod)[2],3)
    rSquaredMarginalList[[length(rSquaredMarginalList)+1]] <- rsquaredMarginal;
    rSquaredConditionalList[[length(rSquaredConditionalList)+1]] <- rsquaredConditional;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
  originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
  originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  timeAdj <- p.adjust(pValTimeList, method = "fdr")
  
  
  makeTable=data.frame(colnames(keggLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,tissueMeans,tissueSDs,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj,rSquaredMarginalList,rSquaredConditionalList);
  setwd("../../statisticalModels/")
  write("wgs\tstoolMeans\tstoolSD\tswabMeans\tswabSD\ttissueMeans\ttissueSD\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)\trsquared marginal\trsquared conditional",paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgs_mixedMod.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgs_mixedMod.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..") 
}

################### WGS functions (No Tissue)  ########################
sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"

wgsLevels <- c("keggPathwaysLevel1",
               "keggPathwaysLevel2",
               "keggPathwaysLevel3",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")

#wgsLevels <- c("keggFamilies") 
for(wgs in wgsLevels )
{
  setwd("data/metagenomeFunctions/")
  inFileName <- paste("wgs_",wgs, ".txt", sep ="")
  print(inFileName)
  kegg <-read.delim(inFileName,header = TRUE, row.names=1)
  kegg <- kegg[,-c(1:16)]
  kegg <- t(kegg)
  #keggLogged <- kegg
  keggLogged <- log10((kegg*175000) +1)
  keggLogged <- keggLogged[,(colSums(keggLogged==0)/nrow(keggLogged))<=0.75]
  keggLogged2 <- keggLogged
  #colnames(keggLogged2) <- substr(colnames(keggLogged2),1,6)
  
  keggMeta <- merge(sampleDataWGS,keggLogged2,by="row.names");
  keggSwab <- split(keggMeta,keggMeta$Origin)$swab
  keggStool <- split(keggMeta,keggMeta$Origin)$stool
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  rSquaredMarginalList <- numeric(0);
  rSquaredConditionalList <- numeric(0);
  
  stoolMeans <- round(colMeans(keggStool[,names(keggStool) %in% colnames(keggLogged2)]),3); 
  swabMeans <- round(colMeans(keggSwab[,names(keggSwab) %in% colnames(keggLogged2)]),3); 
  
  stoolSDs <- round(colSds(as.matrix(keggStool[,names(keggStool) %in% colnames(keggLogged2)])),3); 
  swabSDs <- round(colSds(as.matrix(keggSwab[,names(keggSwab) %in% colnames(keggLogged2)])),3); 
  for(i in (ncol(sampleDataWGS)+2): ncol(keggMeta))
  {
    model <- as.formula(paste(names(keggMeta)[i],"~","Origin","+","visit"));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=keggMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=keggMeta,control = lmeControl(singular.ok=TRUE, returnObject = TRUE));
    #r.squaredGLMM(mixedMod)
    anovaModForTukey <- aov(model,data=keggMeta)
    
    #pValOrigin <- pf(anova(mixedMod)$"F-value"[2],1,75,lower.tail = FALSE);
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- format.pval(pValOrigin,3);
    pValIndividualList[[length(pValIndividualList)+1]] <- format.pval(pValIndividual,3);
    
    rsquaredMarginal <- round(r.squaredGLMM(mixedMod)[1],3)
    rsquaredConditional <- round(r.squaredGLMM(mixedMod)[2],3)
    rSquaredMarginalList[[length(rSquaredMarginalList)+1]] <- rsquaredMarginal;
    rSquaredConditionalList[[length(rSquaredConditionalList)+1]] <- rsquaredConditional;
  }
  originAdj <- p.adjust(pValOriginList, method = "fdr")
  individualAdj <- p.adjust(pValIndividualList, method = "fdr")
  
  
  makeTable=data.frame(colnames(keggLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,pValOriginList, pValIndividualList,originAdj,individualAdj,rSquaredMarginalList,rSquaredConditionalList);
  setwd("../../statisticalModels/")
  write("wgs\tstoolMeans\tstoolSD\tswabMeans\tswabSD\tOrigin p-value\tIndividual p-value\tOrigin (adj p-value)\tIndividual (adj p-value)\trsquared marginal\trsquared conditional",paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgsNoTissue_mixedMod.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_wgsBywgs_",wgs,"_wgsNoTissue_mixedMod.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..") 
}
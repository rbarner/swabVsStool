


sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"

#bacteriaMeta <- merge(sampleDataWGS,bacteriaLogged,by="row.names");
#bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$swab
#bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$stool
#bacteriaTissue <- split(bacteriaMeta,bacteriaMeta$Origin)$tissue

#keggMetaSplit=split(keggMeta, f=keggMeta$Origin);

taxaLevels <- c("phylum","class","order","family","genus","species")
#taxaLevels <- c("species")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
  inFileName <- paste("krakenWGS_",taxa, "Level2.txt", sep ="")
  print(inFileName)
  bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
  bacteria <- t(bacteria)

  bacteriaMeta <- merge(sampleDataWGS,bacteria,by="row.names");
  bacteriaMetaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$swab
  bacteriaMetaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$stool
  bacteriaMetaTissue <- split(bacteriaMeta,bacteriaMeta$Origin)$tissue
  
  bacteriaMetaTissue=bacteriaMetaTissue[order(bacteriaMetaTissue[,6]),];
  tissuePatients = bacteriaMetaTissue$study_id;
  bacteriaMetaTissueSplits = split(bacteriaMetaTissue,f=bacteriaMetaTissue$visit)
  bacteriaMetaTissuePre = bacteriaMetaTissueSplits$Pre;
  bacteriaMetaTissuePost = bacteriaMetaTissueSplits$Post;
  

  bacteriaMetaSwab=bacteriaMetaSwab[order(bacteriaMetaSwab[,6]),];
  bacteriaMetaSwabSub = bacteriaMetaSwab[bacteriaMetaSwab$study_id %in% tissuePatients,];
  bacteriaMetaSwabSub2 = bacteriaMetaSwab[!(bacteriaMetaSwab$study_id %in% tissuePatients),];
  swabPatients = bacteriaMetaSwab$study_id;
  bacteriaMetaSwabSplits = split(bacteriaMetaSwabSub,f=bacteriaMetaSwabSub$visit)
  bacteriaMetaSwabSplits2 = split(bacteriaMetaSwabSub2,f=bacteriaMetaSwabSub2$visit);
  bacteriaMetaSwabPre = bacteriaMetaSwabSplits$Pre;
  bacteriaMetaSwabPost = bacteriaMetaSwabSplits$Post;
  bacteriaMetaSwabPre2 = bacteriaMetaSwabSplits2$Pre;
  bacteriaMetaSwabPost2 = bacteriaMetaSwabSplits2$Post;
  
  bacteriaMetaStool=bacteriaMetaStool[order(bacteriaMetaStool[,6]),];
  bacteriaMetaStoolSub = bacteriaMetaStool[bacteriaMetaStool$study_id %in% tissuePatients,];
  bacteriaMetaStoolSub2 = bacteriaMetaStool[!(bacteriaMetaStool$study_id %in% swabPatients),];
  bacteriaMetaStoolSplits = split(bacteriaMetaStoolSub,f=bacteriaMetaStoolSub$visit);
  bacteriaMetaStoolSplits2 = split(bacteriaMetaStoolSub2,f=bacteriaMetaStoolSub2$visit);
  bacteriaMetaStoolPre = bacteriaMetaStoolSplits$Pre;
  bacteriaMetaStoolPost = bacteriaMetaStoolSplits$Post;
  bacteriaMetaStoolPre2 = bacteriaMetaStoolSplits2$Pre;
  bacteriaMetaStoolPost2 = bacteriaMetaStoolSplits2$Post;
  
  
  bacteriaMetaPre = rbind(bacteriaMetaStoolPre2,bacteriaMetaSwabPre2,bacteriaMetaTissuePre);
  bacteriaMetaPairedPre = rbind(bacteriaMetaStoolPre,bacteriaMetaSwabPre,bacteriaMetaTissuePre);
  bacteriaMetaPre2 = bacteriaMetaPre[-16,];
  bacteriaMetaPost = rbind(bacteriaMetaStoolPost2,bacteriaMetaSwabPost2,bacteriaMetaTissuePost);
  bacteriaMetaPairedPost = rbind(bacteriaMetaStoolPost,bacteriaMetaSwabPost,bacteriaMetaTissuePost);
  
  bacteriaMetaPaired=rbind(bacteriaMetaPairedPost,bacteriaMetaPairedPre)
  bacteriaMetaDistinct=rbind(bacteriaMetaPost,bacteriaMetaPre)
  
  bacteriaSwab <- split(bacteriaMetaPaired,bacteriaMetaPaired$Origin)$swab
  bacteriaStool <- split(bacteriaMetaPaired,bacteriaMetaPaired$Origin)$stool
  bacteriaTissue <- split(bacteriaMetaPaired,bacteriaMetaPaired$Origin)$tissue
  
  bacteriaSwabDistinct <- split(bacteriaMetaDistinct,bacteriaMetaDistinct$Origin)$swab
  bacteriaStoolDistinct <- split(bacteriaMetaDistinct,bacteriaMetaDistinct$Origin)$stool
  bacteriaTissueDistinct <- split(bacteriaMetaDistinct,bacteriaMetaDistinct$Origin)$tissue
  
  bacteriaSwabDistinct2 <- bacteriaSwabDistinct[,-c(1:23)]
  bacteriaStoolDistinct2 <- bacteriaStoolDistinct[,-c(1:23)]
  bacteriaTissueDistinct2 <- bacteriaTissueDistinct[,-c(1:23)]
  
  mean(rowSums(bacteriaTissueDistinct2))
  sd(rowSums(bacteriaTissueDistinct2))
  
  mean(rowSums(bacteriaSwabDistinct2))
  sd(rowSums(bacteriaSwabDistinct2))
  
  mean(rowSums(bacteriaStoolDistinct2))
  sd(rowSums(bacteriaStoolDistinct2))
  
  
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
  
  
  makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,stoolSDs,swabMeans,swabSDs,tissueMeans,tissueSDs,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj,rSquaredMarginalList, rSquaredConditionalList);
  setwd("../../statisticalModels/")
  write("Taxa\tstoolMeans\tstoolSD\tswabMeans\tswabSD\ttissueMeans\ttissueSD\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)\trsquared marginal\trsquared conditional",paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS_mixedMod.txt",sep=""));
  write.table(makeTable,paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS_mixedMod.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  setwd("..")
}


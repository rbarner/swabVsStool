library(nlme)
#####################################################33
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")

sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"
sampleDataWGS$Origin <- factor(toupper(sampleDataWGS$Origin))
sampleDataWGS$visit <- tolower(sampleDataWGS$visit)
sampleDataWGS$participantID <- paste(sampleDataWGS$study_id,sampleDataWGS$visit,sampleDataWGS$Origin,sep = "_")

sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]
sampleData16S$participantID <- paste(sampleData16S$study_id,sampleData16S$visit,sampleData16S$Origin,sep = "_")

copyCounts <- read.delim("data/16ScopyCount.txt",header = TRUE)
taxaLevels <- c("phylum","class","order","family","genus","species")
#copyCountColumnIndex <- 15
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
  setwd("data/microbialClassfications/")
  inFileName16S <- paste("kraken16S_",taxa, "Level.txt", sep ="")
  print(inFileName16S)
  bacteria16S <-read.delim(inFileName16S,header = TRUE, row.names=1)
  bacteria16S <- t(bacteria16S)
  if(taxa %in% "species")
  {
    testMatrix <- sapply(strsplit(colnames(bacteria16S),"_"), `[`, 1:2)
    colnames(bacteria16S) <- apply(testMatrix, 2, paste, collapse=" ")
  }else{
    testMatrix <- sapply(strsplit(colnames(bacteria16S),"_"), `[`, 1)
    colnames(bacteria16S) <- testMatrix;
  }
  bacteria16SLogged <- log10((bacteria16S/rowSums(bacteria16S))* (sum(colSums(bacteria16S))/dim(bacteria16S)[1]) +1)
  #bacteria16SLogged <- bacteria16SLogged[,(colSums(bacteria16SLogged==0)/nrow(bacteria16SLogged))<=0.75]
  bacteria16S <- bacteria16S[,(colSums(bacteria16S==0)/nrow(bacteria16S))<=0.6]
  bacteria16SMeta <- merge(sampleData16S,bacteria16SLogged,by="row.names");
  bacteria16SMeta2 <- merge(sampleDataWGS,bacteria16S,by="row.names");
  
  inFileNameWGS <- paste("krakenWGS_",taxa, "Level.txt", sep ="")
  print(inFileNameWGS)
  bacteriaWGS <-read.delim(inFileNameWGS,header = TRUE, row.names=1)
  bacteriaWGS <- t(bacteriaWGS)
  if(taxa %in% "species")
  {
    testMatrix <- sapply(strsplit(colnames(bacteriaWGS),"_"), `[`, 1:2);
    colnames(bacteriaWGS) <- apply(testMatrix, 2, paste, collapse=" ");
  }else{
    testMatrix <- sapply(strsplit(colnames(bacteriaWGS),"_"), `[`, 1);
    colnames(bacteriaWGS) <- testMatrix;
  }
  
  bacteriaWGSLogged <- log10((bacteriaWGS/rowSums(bacteriaWGS))* (sum(colSums(bacteriaWGS))/dim(bacteriaWGS)[1]) +1)
  #bacteriaWGSLogged <- bacteriaWGSLogged[,(colSums(bacteriaWGSLogged==0)/nrow(bacteriaWGSLogged))<=0.75]
  bacteriaWGS <- bacteriaWGS[,(colSums(bacteriaWGS==0)/nrow(bacteriaWGS))<=0.6]
  bacteriaWGSMeta <- merge(sampleDataWGS,bacteriaWGSLogged,by="row.names");
  bacteriaWGSMeta2 <- merge(sampleDataWGS,bacteriaWGS,by="row.names");
  
  
  #bacteriaWGSMetaCopycount <- merge(copyCounts,bacteriaWGSMeta,by.x="species",by.y="col.names")
  subsetList <- intersect(sampleData16S$participantID,sampleDataWGS$participantID)
  bacteriaWGSMeta2 <- bacteriaWGSMeta2[bacteriaWGSMeta2$participantID %in% subsetList,]
  bacteria16SMeta2 <- bacteria16SMeta2[bacteria16SMeta2$participantID %in% subsetList,]
  
  bacteriaNames <- intersect(colnames(bacteriaWGS),colnames(bacteria16S))
  #bacteriaNames <- bacteriaNames[bacteriaNames %in% copyCounts[,copyCountColumnIndex]]
  
  correlationList=numeric();
  pValOriginList <- numeric();
  pVal16SList <- numeric()
  pValInteractionList <- numeric();
  pValParticipantList <- numeric();
  
  pValOriginCopyCountList <- numeric();
  pVal16SCopyCountList <- numeric()
  pValInteractionCopyCountList <- numeric();
  pValParticipantCopyCountList <- numeric();
  
  pdf(paste("7_",taxa,"_level_16SvsWGS_comparisonPlots.pdf",sep=""))
  for(thisTaxa in bacteriaNames)
  {
    #copyCount <- names(sort(-table(copyCounts[copyCounts[,copyCountColumnIndex] %in% thisTaxa,]$X16S.count)))[1]
    #copyCount <- mean(copyCounts[copyCounts[,copyCountColumnIndex] %in% thisTaxa,]$X16S.count);
    bacteriaWGSMetaSub <- bacteriaWGSMeta2[,names(bacteriaWGSMeta2) %in% c(thisTaxa,"participantID","Origin","study_id")]
    names(bacteriaWGSMetaSub)[4] <- "taxaWGS"
    bacteria16SMetaSub <- bacteria16SMeta2[,names(bacteria16SMeta2) %in% c(thisTaxa,"participantID","Origin","study_id")]
    names(bacteria16SMetaSub)[4] <- "taxa16S"
    #if(ncol(bacteria16SMetaSub)==2)
    #{
     # bacteria16SMetaSub$taxa <- 0.0
    #}
    #else
    #{
     # bacteriaWGSMetaSub$taxa <- 0.0
    #}
    subTable <- merge(bacteria16SMetaSub,bacteriaWGSMetaSub,by="participantID",all=TRUE)
    #subTable$taxa16SCopycount <- subTable$taxa16S/as.numeric(copyCount)
    #subTable <- subTable[,-4]
    #names(subTable)[3:4] <- c("taxa16S","taxaWGS")
    correlation <- cor(subTable$taxa16S,subTable$taxaWGS)
    
    plot(subTable$taxa16S,subTable$taxaWGS,
         col=ifelse(subTable$Origin.x=="SWAB","red2","black"),
         pch=19,
         ylab="Log10 (WGS sequence read count)",
         xlab="Log10 (16S rRNA sequence read count)",
         cex=1.5,
         main = paste("Taxa:",thisTaxa))
    
    correlationList[[length(correlationList)+1]] <- correlation;
    
    statMod <- as.formula("taxaWGS~taxa16S*Origin.x");
    simpleMod <- gls(statMod,method="REML",data=subTable,control=lmeControl(singular.ok=TRUE, returnObject = TRUE));
    mixedMod <- lme(statMod,method="REML",random=~1|study_id.x,data=subTable,control=lmeControl(singular.ok=TRUE, returnObject = TRUE));
    pVal16S <- format(pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE),digit=3);
    pValOrigin <- format(pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE),digit=3);
    pValInteraction <- format(pf(anova(mixedMod)$"F-value"[4],anova(mixedMod)$"numDF"[4],anova(mixedMod)$"denDF"[4],lower.tail = FALSE),digit=3);
    pValParticipant <- format(pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE),digit=3)
    
    
    #statMod <- lm(taxaWGS~taxa16S+Origin.x+Origin.x*taxa16S,data=subTable)
    #pValOrigin <- anova(statMod)$"Pr(>F)"[1]
    #pVal16S <- anova(statMod)$"Pr(>F)"[2]
    #pValInteraction <- anova(statMod)$"Pr(>F)"[3]
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pVal16SList[[length(pVal16SList)+1]] <- pVal16S;
    pValInteractionList[[length(pValInteractionList)+1]] <- pValInteraction;
    pValParticipantList[[length(pValParticipantList)+1]] <- pValParticipant;
    
    #statModCopyCount <- lm(taxaWGS~taxa16SCopycount+Origin.x+Origin.x*taxa16SCopycount,data=subTable)
    
    # statModCopyCount <- as.formula("taxaWGS~taxa16SCopycount*Origin.x");
    # simpleMod <- gls(statModCopyCount,method="REML",data=subTable,control=lmeControl(singular.ok=TRUE, returnObject = TRUE));
    # mixedMod <- lme(statModCopyCount,method="REML",random=~1|study_id.x,data=subTable,control=lmeControl(singular.ok=TRUE, returnObject = TRUE));
    # pVal16SCopyCount <- format(pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE),digit=3);
    # pValOriginCopyCount <- format(pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE),digit=3);
    # pValInteractionCopyCount <- format(pf(anova(mixedMod)$"F-value"[4],anova(mixedMod)$"numDF"[4],anova(mixedMod)$"denDF"[4],lower.tail = FALSE),digit=3);
    # pValParticipantCopyCount <- format(pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE),digit=3)
    # 
    
    #pValOriginCopyCount  <- anova(statModCopyCount)$"Pr(>F)"[1]
    #pVal16SCopyCount  <- anova(statModCopyCount)$"Pr(>F)"[2]
    #pValInteractionCopyCount  <- anova(statModCopyCount)$"Pr(>F)"[3]
    
    #pValOriginCopyCountList[[length(pValOriginCopyCountList)+1]] <- pValOriginCopyCount;
    #pVal16SCopyCountList[[length(pVal16SCopyCountList)+1]] <- pVal16SCopyCount;
    #pValInteractionCopyCountList[[length(pValInteractionCopyCountList)+1]] <- pValInteractionCopyCount;
    #pValParticipantCopyCountList[[length(pValParticipantCopyCountList)+1]] <- pValParticipantCopyCount;
    
  }
  dev.off()
  
  pdf(paste("7_",taxa,"_correlation_histogram.pdf",sep=""))
  hist(correlationList, col="grey76", main=paste(taxa, "level 16S vs WGS correlations"),xlab = "correlations")
  dev.off()
  
  
  makeTable=data.frame(bacteriaNames,correlationList);
  setwd("../../statisticalModels/")
  write("Taxa\tcorrelation",paste("7_correlationKraken16SvsWGS_",taxa,"2.txt",sep=""));
  write.table(makeTable,paste("7_correlationKraken16SvsWGS_",taxa,"2.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);

  pValOriginListAdj <- p.adjust(pValOriginList,method = "fdr")
  pVal16SListAdj <- p.adjust(pVal16SList,method = "fdr")
  pValInteractionListAdj <- p.adjust(pValInteractionList,method = "fdr")
  pValParticipantListAdj <- p.adjust(pValParticipantList,method = "fdr")
  
  makeTable=data.frame(bacteriaNames,pValOriginList,pVal16SList,pValInteractionList,pValParticipantList,pValOriginListAdj,pVal16SListAdj,pValInteractionListAdj,pValParticipantListAdj);
  write("Taxa\tOrigin P-value\t16S P-value\tInteraction P-value\tParticipant P-value\tOrigin p-val adj\t16S p-val adj\tInteraction p-val adj\tParticipantp-val adj",paste("7_pValues16SvsWGS_",taxa,".txt",sep=""));
  write.table(makeTable,paste("7_pValues16SvsWGS_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  
  # pValOriginCopyCountListAdj <- p.adjust(pValOriginCopyCountList,method = "fdr")
  # pVal16SCopyCountListAdj <- p.adjust(pVal16SCopyCountList,method = "fdr")
  # pValInteractionCopyCountListAdj <- p.adjust(pValInteractionCopyCountList,method = "fdr")
  # pValParticipantCopyCountListAdj <- p.adjust(pValParticipantCopyCountList,method = "fdr")
  # 
  # 
  # makeTable=data.frame(bacteriaNames,pValOriginCopyCountList,pVal16SCopyCountList,pValInteractionCopyCountList,pValParticipantCopyCountList,pValOriginCopyCountListAdj,pVal16SCopyCountListAdj,pValInteractionCopyCountListAdj,pValParticipantCopyCountListAdj);
  # write("Taxa\tOrigin P-value\t16S P-value\tInteraction P-value\tParticipant P-value\tOrigin p-val adj\t16S p-val adj\tInteraction p-val adj\tParticipant p-val adj",paste("7_pValuesCopyCount16SvsWGS_",taxa,".txt",sep=""));
  # write.table(makeTable,paste("7_pValuesCopyCount16SvsWGS_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
  # copyCountColumnIndex <- copyCountColumnIndex-1;
}

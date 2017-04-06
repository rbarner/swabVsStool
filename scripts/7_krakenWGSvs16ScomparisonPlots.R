setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")

sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"
sampleDataWGS$Origin <- factor(toupper(sampleDataWGS$Origin))
sampleDataWGS$visit <- tolower(sampleDataWGS$visit)
sampleDataWGS$participantID <- paste(sampleDataWGS$study_id,sampleDataWGS$visit,sampleDataWGS$Origin,sep = "_")

sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]
sampleData16S$participantID <- paste(sampleData16S$study_id,sampleData16S$visit,sampleData16S$Origin,sep = "_")

taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
  setwd("data/microbialClassfications/")
  inFileName16S <- paste("kraken16S_",taxa, "Level.txt", sep ="")
  print(inFileName16S)
  bacteria16S <-read.delim(inFileName16S,header = TRUE, row.names=1)
  bacteria16S <- t(bacteria16S)
  #colnames(bacteria16S)  <- gsub("\\[","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("\\]","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("\\(","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("\\)","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("-","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("\\+","\\.",colnames(bacteria16S))
  #colnames(bacteria16S)  <- gsub("\\/","\\.",colnames(bacteria16S))
  bacteria16SLogged <- log10((bacteria16S/rowSums(bacteria16S))* (sum(colSums(bacteria16S))/dim(bacteria16S)[1]) +1)
  #myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  bacteria16SLogged <- bacteria16SLogged[,(colSums(bacteria16SLogged==0)/nrow(bacteria16SLogged))<=0.75]
  bacteria16SMeta <- merge(sampleData16S,bacteria16SLogged,by="row.names");
  
  inFileNameWGS <- paste("krakenWGS_",taxa, "Level.txt", sep ="")
  print(inFileNameWGS)
  bacteriaWGS <-read.delim(inFileNameWGS,header = TRUE, row.names=1)
  bacteriaWGS <- t(bacteriaWGS)
  #colnames(bacteriaWGS)  <- gsub("\\[","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("\\]","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("\\(","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("\\)","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("-","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("\\+","\\.",colnames(bacteriaWGS))
  #colnames(bacteriaWGS)  <- gsub("\\/","\\.",colnames(bacteriaWGS))
  bacteriaWGSLogged <- log10((bacteriaWGS/rowSums(bacteriaWGS))* (sum(colSums(bacteriaWGS))/dim(bacteriaWGS)[1]) +1)
  #myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
  #myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
  bacteriaWGSLogged <- bacteriaWGSLogged[,(colSums(bacteriaWGSLogged==0)/nrow(bacteriaWGSLogged))<=0.75]
  bacteriaWGSMeta <- merge(sampleDataWGS,bacteriaWGSLogged,by="row.names");
  
  subsetList <- intersect(sampleData16S$participantID,sampleDataWGS$participantID)
  bacteriaWGSMeta <- bacteriaWGSMeta[bacteriaWGSMeta$participantID %in% subsetList,]
  bacteria16SMeta <- bacteria16SMeta[bacteria16SMeta$participantID %in% subsetList,]
  
  bacteriaNames <- intersect(colnames(bacteriaWGSLogged),colnames(bacteria16SLogged))
  correlationList=numeric();
  pdf(paste(taxa,"_level_16SvsWGS_comparisonPlots.pdf",sep=""))
  for(thisTaxa in bacteriaNames)
  {
    bacteriaWGSMetaSub <- bacteriaWGSMeta[,names(bacteriaWGSMeta) %in% c(thisTaxa,"participantID","Origin")]
    bacteria16SMetaSub <- bacteria16SMeta[,names(bacteria16SMeta) %in% c(thisTaxa,"participantID","Origin")]
    #if(ncol(bacteria16SMetaSub)==2)
    #{
     # bacteria16SMetaSub$taxa <- 0.0
    #}
    #else
    #{
     # bacteriaWGSMetaSub$taxa <- 0.0
    #}
    subTable <- merge(bacteria16SMetaSub,bacteriaWGSMetaSub,by="participantID",all=TRUE)
    subTable <- subTable[,-4]
    names(subTable)[3:4] <- c("taxa16S","taxaWGS")
    correlation <- cor(subTable$taxa16S,subTable$taxaWGS)
    
    plot(subTable$taxa16S,subTable$taxaWGS,
         col=ifelse(subTable$Origin.x=="SWAB","red2","black"),
         pch=19,
         ylab="Log10 (WGS sequence read count)",
         xlab="Log10 (16S rRNA sequence read count)",
         cex=1.5,
         main = paste("Taxa:",thisTaxa))
    correlationList[[length(correlationList)+1]] <- correlation;
  }
  dev.off()
  
  pdf(paste(taxa,"_correlation_histogram.pdf",sep=""))
  hist(correlationList, col="grey76", main=paste(taxa, "level 16S vs WGS correlations"),xlab = "correlations")
  dev.off()
  makeTable=data.frame(bacteriaNames,correlationList);
  setwd("../../statisticalModels/")
  write("Taxa\tcorrelation",paste("correlationKraken16SvsWGS_",taxa,".txt",sep=""));
  write.table(makeTable,paste("correlationKraken16SvsWGS_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
}

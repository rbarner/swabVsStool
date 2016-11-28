source("scripts/iwanthue.R")
#sampleData <- readRDS("../key/mapping_key_16S.RData");
#sampleData2 <- readRDS("../key/mapping_key_WGS.RData");
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);

##################### RDP Classifier and QIIME ###############################################
colorsPlot16S <- c("#cb5469","#60a75c","#8065cb","#ad963e","#c361ab","#cc693c","#6596cc")
toolList <- c("rdpClassifications", "qiime")

phylumNames <- c("Firmicutes","Bacteroidetes","Proteobacteria", "Actinobacteria","Fusobacteria","Verrucomicrobia")
for(tool in toolList)
{
  taxaLevels <- c("phylum")
  for(taxa in taxaLevels )
  {
    inFileName <- paste("data/microbialClassfications/",tool,"_",taxa,"Level.txt",sep="")
    classificationLevel <- read.delim(inFileName,header=TRUE,row.names = 1)
    classificationLevel <- t(classificationLevel)
    classificationPhylum <- classificationLevel[,colnames(classificationLevel) %in% phylumNames]
    classificationOthers <- classificationLevel[,!colnames(classificationLevel) %in% phylumNames]
    classificationPhylum <- cbind(classificationPhylum,as.vector(rowSums(classificationOthers)))
    colnames(classificationPhylum)[7] <- "Others"
    
    classificationMeta <- merge(sampleData,classificationPhylum, by = "row.names")
    classificationLevel <- classificationMeta[,-c(1:4,7:35)];
    row.names(classificationLevel) <- classificationMeta$Row.names;
    #classificationLevel <- rbind(as.character(sampleData$Origin),as.character(sampleData$study_id),classificationLevel);
    orderedTable <- classificationLevel[order(classificationLevel$Origin,classificationLevel$study_id),];
    normClassificationLevel <- t(prop.table(as.table(as.matrix(orderedTable[,3:ncol(orderedTable)])),1))
    #normClassificationLevel <- t(prop.table(t(as.table(as.matrix(orderedTable[3:(dim(orderedTable)[1]),]))),1))
    normClassificationLevel2 <- normClassificationLevel[order(-rowSums(normClassificationLevel)),];
    
    # meltedData <- melt(normClassificationLevel2, id.var=row.names(normClassificationsLevel2)) 
    # names(meltedData) <- c("Bacteria","Sample","value")
    # 
    # p <- ggplot(meltedData,aes(x=Sample,y=value,fill = Bacteria))
    # p + geom_bar(stat="identity")
    # ggplot(DF1, aes(x = Rank, y = value, fill = variable)) + 
    #   geom_bar(stat = "identity")
    
    
    tiff(paste("plots/",tool,"_",taxa,"distributions.tiff",sep=""),width=4800,height=2200,compression="lzw",res=300)
    layout(matrix(c(0,0,0,0,0,0,0,2,2,
                    0,1,1,1,1,1,1,2,2,
                    0,1,1,1,1,1,1,2,2,
                    0,0,0,0,0,0,0,2,2),
                  4,9,byrow=TRUE));
    barplot(as.matrix(normClassificationLevel2),col=colorsPlot,border="black",
            axisnames=FALSE, main=paste("Microbial Distributions Per Sample (",taxa," level)",sep=""),
            xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=119),0.75,rep(0,times=119)))
    axis(1,at = c(60,180),cex.axis=1.25, labels = c("Stool Samples", "Swab Samples"),tick=FALSE,padj = 2)
    plot.new();
    legend("center",box.col="white",pch=c(rep(15,times = dim(normClassificationLevel2)[1])),
           cex=1.25,pt.cex=2.5,col=colorsPlot,legend=rownames(normClassificationLevel2));
    dev.off()
  }
}

##################### Metaphlan ###############################################
phylumNames <- c("Firmicutes","Bacteroidetes","Proteobacteria", "Actinobacteria","Fusobacteria","Verrucomicrobia")

inFileName <- paste("data/microbialClassfications/metaphlan_phylumLevel.txt",sep="")
functionWGS <- read.delim(inFileName,header=TRUE,row.names = 1)
functionWGS <- t(functionWGS)
functionMeta <- merge(sampleData2,functionWGS, by = "row.names")
functionWGS <- functionMeta[,-c(1,4:23)];
row.names(functionWGS) <- functionMeta$Row.names;

orderedTable <- functionWGS[order(functionWGS$type,functionWGS$study_id),];
normWGSFunctionLevel <- t(prop.table(as.table(as.matrix(orderedTable[,3:ncol(orderedTable)])),1))
normWGSFunctionLevel2 <- normWGSFunctionLevel[order(-rowSums(normWGSFunctionLevel)),];

tiff("plots/metaphlan_phylum_distributions.tiff",width=4800,height=2200,compression="lzw",res=300)
layout(matrix(c(0,0,0,0,0,0,0,2,2,2,
                0,1,1,1,1,1,1,2,2,2,
                0,1,1,1,1,1,1,2,2,2,
                0,0,0,0,0,0,0,2,2,2),
              4,10,byrow=TRUE));
barplot(as.matrix(normWGSFunctionLevel2),col=iwanthue(dim(normWGSFunctionLevel2)[1]),border="black",
        axisnames=FALSE, main="Distributions of MetaPhlAn Phylum",
        xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=100),0.75,rep(0,times=26),0.75,rep(0,times=15)))
axis(1,at = c(50,113,135),cex.axis=1.25, labels = c("Stool", "Swab","Tissue"),tick=FALSE,padj = 2)
plot.new();
legend("center",box.col="white",pch=c(rep(15,times = dim(normWGSFunctionLevel2)[1])),
       cex=1.25,pt.cex=2.5,col=rev(iwanthue(dim(normWGSFunctionLevel2)[1])),legend=rev(rownames(normWGSFunctionLevel2)));
dev.off()



####################### WGS functions ##########################################
wgsLevels <- c("keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2")
for(wgs in wgsLevels )
{
  inFileName <- paste("data/metagenomeFunctions/wgs_",wgs,".txt",sep="")
  functionWGS <- read.delim(inFileName,header=TRUE,row.names = 1)
  functionWGS <- t(functionWGS)
  functionMeta <- merge(sampleData2,functionWGS, by = "row.names")
  functionWGS <- functionMeta[,-c(1,4:23)];
  row.names(functionWGS) <- functionMeta$Row.names;
  
  orderedTable <- functionWGS[order(functionWGS$type,functionWGS$study_id),];
  normWGSFunctionLevel <- t(prop.table(as.table(as.matrix(orderedTable[,3:ncol(orderedTable)])),1))
  normWGSFunctionLevel2 <- normWGSFunctionLevel[order(-rowSums(normWGSFunctionLevel)),];
  
  tiff(paste("plots/",wgs,"_distributions.tiff",sep=""),width=4800,height=2200,compression="lzw",res=300)
  layout(matrix(c(0,0,0,0,0,0,0,2,2,2,
                  0,1,1,1,1,1,1,2,2,2,
                  0,1,1,1,1,1,1,2,2,2,
                  0,0,0,0,0,0,0,2,2,2),
                4,10,byrow=TRUE));
  barplot(as.matrix(normWGSFunctionLevel2),col=iwanthue(dim(normWGSFunctionLevel2)[1]),border="black",
          axisnames=FALSE, main=paste("Distributions of ",wgs,sep=""),
          xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=100),0.75,rep(0,times=26),0.75,rep(0,times=15)))
  axis(1,at = c(50,113,135),cex.axis=1.25, labels = c("Stool", "Swab","Tissue"),tick=FALSE,padj = 2)
  plot.new();
  legend("center",box.col="white",pch=c(rep(15,times = dim(normWGSFunctionLevel2)[1])),
         cex=1.25,pt.cex=2.5,col=rev(iwanthue(dim(normWGSFunctionLevel2)[1])),legend=rev(rownames(normWGSFunctionLevel2)));
  dev.off()
}

############################## PICRUSt functions ################################################
wgsLevels <- c("keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2")
for(wgs in wgsLevels )
{
    inFileName <- paste("data/metagenomeFunctions/picrust_",wgs,".txt",sep="")
    classificationLevel <- read.delim(inFileName,header=TRUE,row.names = 1)
    classificationLevel <- t(classificationLevel)
    
    classificationMeta <- merge(sampleData,classificationLevel, by = "row.names")
    classificationLevel <- classificationMeta[,-c(1:4,7:35)];
    row.names(classificationLevel) <- classificationMeta$Row.names;
    #classificationLevel <- rbind(as.character(sampleData$Origin),as.character(sampleData$study_id),classificationLevel);
    orderedTable <- classificationLevel[order(classificationLevel$Origin,classificationLevel$study_id),];
    normClassificationLevel <- t(prop.table(as.table(as.matrix(orderedTable[,3:ncol(orderedTable)])),1))
    #normClassificationLevel <- t(prop.table(t(as.table(as.matrix(orderedTable[3:(dim(orderedTable)[1]),]))),1))
    normClassificationLevel2 <- normClassificationLevel[order(-rowSums(normClassificationLevel)),];
    
    # meltedData <- melt(normClassificationLevel2, id.var=row.names(normClassificationsLevel2)) 
    # names(meltedData) <- c("Bacteria","Sample","value")
    # 
    # p <- ggplot(meltedData,aes(x=Sample,y=value,fill = Bacteria))
    # p + geom_bar(stat="identity")
    # ggplot(DF1, aes(x = Rank, y = value, fill = variable)) + 
    #   geom_bar(stat = "identity")
    
    
    tiff(paste("plots/picrust_",wgs,"_distributions.tiff",sep=""),width=4800,height=2200,compression="lzw",res=300)
    layout(matrix(c(0,0,0,0,0,0,0,2,2,
                    0,1,1,1,1,1,1,2,2,
                    0,1,1,1,1,1,1,2,2,
                    0,0,0,0,0,0,0,2,2),
                  4,9,byrow=TRUE));
    barplot(as.matrix(normClassificationLevel2),col=iwanthue(dim(normClassificationLevel2)[1]),border="black",
            axisnames=FALSE, main=paste("Gene Pathway Distributions Per Sample (",wgs," level)",sep=""),
            xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=119),0.75,rep(0,times=119)))
    axis(1,at = c(60,180),cex.axis=1.25, labels = c("Stool Samples", "Swab Samples"),tick=FALSE,padj = 2)
    plot.new();
    legend("center",box.col="white",pch=c(rep(15,times = dim(normClassificationLevel2)[1])),
           cex=1.25,pt.cex=2.5,col=rev(iwanthue(dim(normClassificationLevel2)[1])),legend=rev(rownames(normClassificationLevel2)));
    dev.off()
}
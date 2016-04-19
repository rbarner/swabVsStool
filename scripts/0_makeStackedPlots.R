source("../../scripts/iwanthue.R")
sampleData <- readRDS("../key/mapping_key_16S.RData");
sampleData2 <- readRDS("../key/mapping_key_WGS.RData");


taxaLevels <- c("phylum","class","order")
for(taxa in taxaLevels )
{
  inFileName <- paste("../microbialClassfications/rdpClassifications_",taxa,"Level.RData",sep="")
  classificationLevel <- readRDS(inFileName)
  classificationLevel <- t(classificationLevel)
  classificationLevel <- rbind(sampleData$Origin,sampleData$study_id,classificationLevel);
  orderedTable <- classificationLevel[,order(classificationLevel[1,],classificationLevel[2,])];
  normClassificationLevel <- t(prop.table(t(as.table(as.matrix(orderedTable[3:(dim(orderedTable)[1]),]))),1))
  normClassificationLevel2 <- normClassificationLevel[order(-rowSums(normClassificationLevel)),];
  
  tiff(paste("../../plots/",taxa,"distributions.tiff"),width=4800,height=2200,compression="lzw",res=300)
  layout(matrix(c(0,0,0,0,0,0,0,2,2,
                  0,1,1,1,1,1,1,2,2,
                  0,1,1,1,1,1,1,2,2,
                  0,0,0,0,0,0,0,2,2),
                4,9,byrow=TRUE));
  barplot(as.matrix(normClassificationLevel2),col=iwanthue(dim(normClassificationLevel2)[1]),border="black",
          axisnames=FALSE, main=paste("Microbial Distributions Per Sample (",taxa," level)",sep=""),
          xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=119),0.5,rep(0,times=119)))
  axis(1,at = c(60,180),cex.axis=1.25, labels = c("Stool Samples", "Swab Samples"),tick=FALSE,padj = 2)
  plot.new();
  legend("center",box.col="white",pch=c(rep(15,times = dim(normClassificationLevel2)[1])),
         cex=1.25,pt.cex=2.5,col=rev(iwanthue(dim(normClassificationLevel2)[1])),legend=rev(rownames(normClassificationLevel2)));
  dev.off()
}


wgsLevels <- c("Pathways_level2",
               "Pathways_level1",
               "MetabolicPathways_level2")
for(wgs in wgsLevels )
{
  inFileName <- paste("../metagenomeFunctions/kegg",wgs,".RData",sep="")
  functionWGS <- readRDS(inFileName)
  functionWGS <- t(functionWGS)
  functionWGS <- rbind(sampleData2$type,sampleData2$study_id,functionWGS);
  orderedTable <- functionWGS[,order(functionWGS[1,],functionWGS[2,])];
  orderedTable  <- orderedTable[order(rownames(orderedTable)),];
  normWGSFunctionLevel <- t(prop.table(t(as.table(as.matrix(orderedTable[3:(dim(orderedTable)[1]),]))),1))
  normWGSFunctionLevel2 <- normWGSFunctionLevel[order(-rowSums(normWGSFunctionLevel)),];
  
  tiff(paste("../../plots/kegg",wgs,"distributions.tiff"),width=4800,height=2200,compression="lzw",res=300)
  layout(matrix(c(0,0,0,0,0,0,0,2,2,2,
                  0,1,1,1,1,1,1,2,2,2,
                  0,1,1,1,1,1,1,2,2,2,
                  0,0,0,0,0,0,0,2,2,2),
                4,10,byrow=TRUE));
  barplot(as.matrix(normWGSFunctionLevel2),col=iwanthue(dim(normWGSFunctionLevel2)[1]),border="black",
          axisnames=FALSE, main=paste("Distributions of Kegg ",wgs,sep=""),
          xaxt="n",yaxt="n",cex.main=3,space=c(rep(0,times=100),0.6,rep(0,times=26),0.6,rep(0,times=15)))
  axis(1,at = c(50,113,135),cex.axis=1.25, labels = c("Stool", "Swab","Tissue"),tick=FALSE,padj = 2)
  plot.new();
  legend("center",box.col="white",pch=c(rep(15,times = dim(normWGSFunctionLevel2)[1])),
         cex=1.25,pt.cex=2.5,col=rev(iwanthue(dim(normWGSFunctionLevel2)[1])),legend=rev(rownames(normWGSFunctionLevel2)));
  dev.off()
}
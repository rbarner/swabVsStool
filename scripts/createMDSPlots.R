rm(list=ls())
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
  title <- paste("MDS plot (",taxa," level)",sep="")
  comp1<-as.character(paste("MDS1"," ",(round(eigen[1],3))*100,"%",sep=""));
  comp2<-as.character(paste("MDS2"," ",(round(eigen[2],3))*100,"%",sep=""));
  
  tiff(paste("../../plots/mdsPlot_Axes12_",taxa,"_coloredByOrigin.tiff",sep=""),width=3600,height=3600,compression="lzw",res=350)
  layout(matrix(c(1,1,1,1,
                  1,1,1,1,
                  1,1,1,1,
                  2,2,2,2),
                ncol=4,nrow=4,byrow=TRUE));
  plot(mdsMeta$MDS1,mdsMeta$MDS2,
       cex=3.5,pch=19,col=ifelse(mdsMeta$Origin=='STOOL','red','blue'),
       xlab=comp1,ylab=comp2,main=title,
       las=1,cex.main=1.5,cex.lab=1.5,cex.axis=1.5); 
  plot.new()
  legend("center",c("Stool","Swab"),col=c("red","blue"),pch=c(19,19),cex=2,pt.cex=2,box.col="white",horiz=TRUE)
  dev.off()
  
  for(axisNum in seq(1:4))
  {
    comp1<-as.character(paste("MDS",axisNum," ",(round(eigen[axisNum],3))*100,"%"));
    tiff(paste("../../plots/mdsSeparation_",taxa,"_level_axis",axisNum,".tiff",sep=""),width=7000,height=3600,compression="lzw",res=400)
    layout(matrix(c(1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,
                    0,0,2,2,2,2,0,0),
                  ncol=8,nrow=4,byrow=TRUE));
    boxplot(mdsMeta[,35+axisNum]~mdsMeta$study_id,ylab=comp1,xlab="Participants", main=title,xaxt="n", cex=1.5, cex.main=1.75,cex.lab=1.5)
    points(mdsMeta[,35+axisNum]~mdsMeta$study_id,pch=19,cex=2,col=ifelse(mdsMeta$Origin=='STOOL','red','blue'))
    plot.new()
    legend("center",c("Stool","Swab"),col=c("red","blue"),pch=c(19,19),cex=2,pt.cex=2,box.col="white",horiz=TRUE)
    dev.off()
  } 
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
  title <- paste("MDS plot (",wgs," level)",sep="")
  comp1<-as.character(paste("MDS1"," ",(round(eigen[1],3))*100,"%",sep=""));
  comp2<-as.character(paste("MDS2"," ",(round(eigen[2],3))*100,"%",sep=""));
  
  tiff(paste("../../plots/mdsPlot_Axes12_",wgs,"_coloredByOrigin.tiff",sep=""),width=3600,height=3600,compression="lzw",res=350)
  layout(matrix(c(1,1,1,1,
                  1,1,1,1,
                  1,1,1,1,
                  2,2,2,2),
                ncol=4,nrow=4,byrow=TRUE));
  plot(mdsMeta$MDS1,mdsMeta$MDS2,
       cex=3.5,pch=19,col=ifelse(mdsMeta$type=='stool','red',ifelse(mdsMeta$type=='swab','cyan3','darkmagenta')),
       xlab=comp1,ylab=comp2,main=title,
       las=1,cex.main=1.5,cex.lab=1.5,cex.axis=1.5); 
  plot.new()
  legend("center",c("Stool","Swab","Tissue"),col=c("red","cyan3","darkmagenta"),pch=c(19,19,19),cex=2,pt.cex=2,box.col="white",horiz=TRUE)
  dev.off()
  
  for(axisNum in seq(1:4))
  {
    comp1<-as.character(paste("MDS",axisNum," ",(round(eigen[axisNum],3))*100,"%"));
    tiff(paste("../../plots/mdsSeparation_",wgs,"_level_axis",axisNum,".tiff",sep=""),width=7000,height=3600,compression="lzw",res=400)
    layout(matrix(c(1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,1,
                    0,0,2,2,2,2,0,0),
                  ncol=8,nrow=4,byrow=TRUE));
    boxplot(mdsMeta[,22+axisNum]~mdsMeta$study_id,ylab=comp1,xlab="Participants", main=title,xaxt="n", cex=1.5, cex.main=1.75,cex.lab=1.5)
    points(mdsMeta[,22+axisNum]~mdsMeta$study_id,pch=19,cex=2,col=ifelse(mdsMeta$type=='stool','red',ifelse(mdsMeta$type=='swab','cyan3','darkmagenta')))
    plot.new()
    legend("center",c("Stool","Swab","Tissue"),col=c("red","cyan3","darkmagenta"),pch=c(19,19,19),cex=2,pt.cex=2,box.col="white",horiz=TRUE)
    dev.off()
  }
}





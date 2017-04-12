#library(sfsmisc)
#rm(list=ls())
library(ggplot2)
library(nlme)

################### Functions Used #######################
originMDSPlots <- function(dataSet,MDSX,MDSY,eigen,tool,level)
{
  MDSX_string <- paste("MDS",MDSX,sep="")
  MDSY_string <- paste("MDS",MDSY,sep="")
  print("Stops at two MDS plots")
  statMod1 <- as.formula(paste("MDS",MDSX,"~","Origin","+","visit",sep=""));
  simpleModX <- gls(statMod1,method="REML",data=dataSet);
  mixedModX <- lme(statMod1,method="REML",random=~1|study_id,data=dataSet);
  pVal_originX <- pf(anova(mixedModX)$"F-value"[2],anova(mixedModX)$"numDF"[2],anova(mixedModX)$"denDF"[2],lower.tail = FALSE);

  statMod2 <- as.formula(paste("MDS",MDSY,"~","Origin","+","visit",sep=""));
  simpleModY <- gls(statMod2,method="REML",data=dataSet);
  mixedModY <- lme(statMod2,method="REML",random=~1|study_id,data=dataSet);
  pVal_originY <- pf(anova(mixedModY)$"F-value"[2],anova(mixedModY)$"numDF"[2],anova(mixedModY)$"denDF"[2],lower.tail = FALSE);
  
  title <- paste("MDS plot (",taxa," level)",sep="")
  comp1<-as.character(paste("MDS",MDSX," ", (round(eigen[MDSX],3))*100,"%, p-value = ",format(pVal_originX,digit=3),sep=""));
  comp2<-as.character(paste("MDS",MDSY," ", (round(eigen[MDSY],3))*100,"%, p-value = ",format(pVal_originY,digit=3),sep=""));

  p <- ggplot(dataSet,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes",MDSX,"_",MDSY,"_",tool,"_",level,"_coloredByOrigin.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes_string(MDSX_string,MDSY_string),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1) + ylab(comp2) +
          ggtitle(title) +
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=20),
                legend.position="bottom",
                legend.title=element_blank()
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
  graphics.off()
}

participantMDSBoxPlots <- function(dataSet,eigen,MDSX,tool,level)
{
  MDSX_string <- paste("MDS",MDSX,sep="")
  print("Stops at participant plots")
  statMod1 <- as.formula(paste("MDS",MDSX,"~","Origin","+","visit",sep=""));
  simpleModX <- gls(statMod1,method="REML",data=dataSet);
  mixedModX <- lme(statMod1,method="REML",random=~1|study_id,data=dataSet);
  pVal_origin <- pf(anova(mixedModX)$"F-value"[2],anova(mixedModX)$"numDF"[2],anova(mixedModX)$"denDF"[2],lower.tail = FALSE);
  pVal_participant <- pchisq(anova(simpleModX,mixedModX)$"L.Ratio"[2],1,lower.tail = FALSE)
  
  print("Stops after participant models")
  title <- paste("MDS plot (",level," level)",sep="")
  comp1Participant<-as.character(paste("MDS",MDSX ," ", (round(eigen[MDSX],3))*100,"%",sep=""));
  comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format(pVal_participant,digit=3),", Origin p-value = ",format(pVal_origin,digit=3),sep=""));
  
  print("Stops after making labels")
  p <- ggplot(mdsMeta,aes_string(x="study_id",y=MDSX_string))
  tiff(paste("2_mdsPlot_Axes",MDSX,"_",tool,"_",level,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1ParticipantMDS) +
          ylab(comp1Participant) +
          ggtitle(title) +
          theme_classic(base_size = 20)+
          theme(axis.line=element_line(size=1),
                axis.ticks.y=element_line(size=1),
                axis.ticks.x=element_blank(),
                axis.text.y=element_text(face="bold",size=24),
                axis.text.x=element_blank(),
                text=element_text(face="bold",size=28),
                legend.position="bottom",
                legend.title=element_blank()
          )+
          theme(axis.line.x = element_line(color="black", size = 2),
                axis.line.y = element_line(color="black", size = 2)
          )
  )
  graphics.off()
}

####################### Creating Plots #########################
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData$visit <- unlist(strsplit(as.character(sampleData$type),split = "_"))[c(FALSE,TRUE)]
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

classifierList <- c("krakenWGS","krakenWGSNoTissue","kraken16S","rdpClassifications", "qiime","metaphlan")
for(classifier in classifierList)
{
  if(classifier %in% c("krakenWGS","kraken16S"))
  {
    taxaLevels <- c("phylum","class","order","family","genus","species")
  }else{
    taxaLevels <- c("phylum","class","order","family","genus")
  }
  for(taxa in taxaLevels )
  {
    setwd("mds")
    mdsFile <- paste(classifier,"_mds_", taxa, "_loggedFiltered.RData",sep="");
    eigenFile <- paste(classifier,"_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    #sampleData  <- sampleData[-26,]
    if(classifier %in% c("krakenWGS","metaphlan"))
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    
    setwd("../plots")
    
    originMDSPlots(mdsMeta,1,2,eigen,classifier,taxa);
    originMDSPlots(mdsMeta,3,4,eigen,classifier,taxa);
    
    participantMDSBoxPlots(dataSet=mdsMeta,eigen=eigen,MDSX=1,tool=classifier,level=taxa)
    participantMDSBoxPlots(dataSet=mdsMeta,eigen=eigen,MDSX=2,tool=classifier,level=taxa)
    participantMDSBoxPlots(dataSet=mdsMeta,eigen=eigen,MDSX=3,tool=classifier,level=taxa)
    participantMDSBoxPlots(dataSet=mdsMeta,eigen=eigen,MDSX=4,tool=classifier,level=taxa)
    
    setwd("..") 
  }
}

setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData$visit <- unlist(strsplit(as.character(sampleData$type),split = "_"))[c(FALSE,TRUE)]
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

functionList <- c("wgs","picrust")
for(funct in functionList)
{
  if(funct %in% c("wgs"))
  {
    wgsLevels <- c("keggFamilies",
                   "keggPathwaysLevel3",
                   "keggPathwaysLevel2",
                   "keggPathwaysLevel1",
                   "metabolickeggPathwaysLevel2",
                   "metabolickeggPathwaysLevel3",
                   "keggFamiliesNoTissue",
                   "keggPathwaysLevel3NoTissue",
                   "keggPathwaysLevel2NoTissue",
                   "keggPathwaysLevel1NoTissue",
                   "metabolickeggPathwaysLevel2NoTissue",
                   "metabolickeggPathwaysLevel3NoTissue")
  }else{
    wgsLevels <- c("keggFamilies",
                   "keggPathwaysLevel3",
                   "keggPathwaysLevel2",
                   "keggPathwaysLevel1",
                   "metabolickeggPathwaysLevel2",
                   "metabolickeggPathwaysLevel3")
  }
  for(wgs in wgsLevels )
  {
    setwd("mds")
    mdsFile <- paste(funct,"_mds_", wgs, "_loggedFiltered.RData",sep="");
    print(mdsFile)
    eigenFile <- paste(funct,"_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
    
    mds <-readRDS(mdsFile);
    #sampleData  <- sampleData[-26,]
    if(funct %in% "wgs")
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }
    else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    
    setwd("../plots")
    
    originMDSPlots(mdsMeta,1,2,eigen,funct,wgs);
    originMDSPlots(mdsMeta,3,4,eigen,funct,wgs);
    
    statMod1 <- formula(paste("MDS1","~","Origin","+","visit",sep=""))
    statMod2 <- formula(paste("MDS2","~","Origin","+","visit",sep=""))
    statMod3 <- formula(paste("MDS3","~","Origin","+","visit",sep=""))
    statMod4 <- formula(paste("MDS4","~","Origin","+","visit",sep=""))
    
    participantMDSBoxPlots(mdsMeta,eigen,1,funct,wgs)
    participantMDSBoxPlots(mdsMeta,eigen,2,funct,wgs)
    participantMDSBoxPlots(mdsMeta,eigen,3,funct,wgs)
    participantMDSBoxPlots(mdsMeta,eigen,4,funct,wgs)
    
    setwd("..")
  }
}
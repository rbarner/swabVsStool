library(ggplot2)
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
######################## Classifications ################################
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

taxaLevels <- c("phylum","class","order","family","genus","species")
#classifierList <- c("rdpClassifications", "qiime","metaphlan")
classifierList <- c("krakenWGS","kraken16S")
for(classifier in classifierList)
{
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
    }
    else
    {
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <-readRDS(eigenFile);
    
    ################ Generate p-values for plot ##############  
    test_origin_mds1 <- aov(MDS1~Origin,mdsMeta);
    pVal_origin_mds1 <- anova(test_origin_mds1)$"Pr(>F)"[1];
    test_origin_mds2 <- aov(MDS2~Origin,mdsMeta);
    pVal_origin_mds2 <- anova(test_origin_mds2)$"Pr(>F)"[1];
    test_origin_mds3 <- aov(MDS3~Origin,mdsMeta);
    pVal_origin_mds3 <- anova(test_origin_mds3)$"Pr(>F)"[1];
    test_origin_mds4 <- aov(MDS4~Origin,mdsMeta);
    pVal_origin_mds4 <- anova(test_origin_mds4)$"Pr(>F)"[1];
    
    test_participant_mds1 <- aov(MDS1~study_id,mdsMeta);
    pVal_participant_mds1 <- anova(test_participant_mds1)$"Pr(>F)"[1];
    test_participant_mds3 <- aov(MDS3~study_id,mdsMeta);
    pVal_participant_mds3 <- anova(test_participant_mds3)$"Pr(>F)"[1];
    test_participant_mds4 <- aov(MDS4~study_id,mdsMeta);
    pVal_participant_mds4 <- anova(test_participant_mds4)$"Pr(>F)"[1];

    
    title <- paste("MDS plot (",taxa," level)",sep="")
    comp1<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
    comp2<-as.character(paste("MDS2 ", (round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_origin_mds2,3),sep=""));
    comp3<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%, p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
    comp4<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%, p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
    
    setwd("../plots")
    
    p <- ggplot(mdsMeta,aes(colour = Origin,shape=Origin))
    tiff(paste("2_mdsPlot_Axes12_",classifier,"_",taxa,"_coloredByOrigin.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
    print(p + geom_point(aes(MDS1,MDS2),size = 8) +
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
    
    p <- ggplot(mdsMeta,aes(x=MDS3,y=MDS4,colour = Origin,shape=Origin))
    tiff(paste("2_mdsPlot_Axes34_",classifier,"_",taxa,"_coloredByOrigin.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
    print(p + geom_point(size = 8) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp3) + ylab(comp4) +
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
    
    comp1Participant<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%",sep=""));
    comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds1,3),", Origin p-value = ",format.pval(pVal_origin_mds1,3),sep=""));

    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS1))
    tiff(paste("2_mdsPlot_Axes1_",classifier,"_",taxa,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
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
    
    comp3Participant<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%",sep=""));
    comp3ParticipantMDS <-as.character(paste("Participant\nParticipant p-value = ",format.pval(pVal_participant_mds3,3),", Origin p-value = ",format.pval(pVal_origin_mds3,3),sep=""));

    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS3))
    tiff(paste("2_mdsPlot_Axes3_",classifier,"_",taxa,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
    print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp3ParticipantMDS)+
            ylab(comp3Participant) +
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
    
    comp4Participant<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%",sep=""));
    comp4ParticipantMDS <-as.character(paste("Participant\nParticipant p-value = ",format.pval(pVal_participant_mds4,3),", Origin p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
    
    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS4))
    tiff(paste("2_mdsPlot_Axes4_",classifier,"_",taxa,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
    print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp4ParticipantMDS) + 
            ylab(comp4Participant) +
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
   setwd("..") 
  }
}

######################## Gene pathways ################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
functionList <- c("wgs","picrust")
for(funct in functionList)
{
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
    
    ################ Generate p-values for plot ##############  
    test_origin_mds1 <- aov(MDS1~Origin,mdsMeta);
    pVal_origin_mds1 <- anova(test_origin_mds1)$"Pr(>F)"[1];
    test_origin_mds2 <- aov(MDS2~Origin,mdsMeta);
    pVal_origin_mds2 <- anova(test_origin_mds2)$"Pr(>F)"[1];
    test_origin_mds3 <- aov(MDS3~Origin,mdsMeta);
    pVal_origin_mds3 <- anova(test_origin_mds3)$"Pr(>F)"[1];
    test_origin_mds4 <- aov(MDS4~Origin,mdsMeta);
    pVal_origin_mds4 <- anova(test_origin_mds4)$"Pr(>F)"[1];
    
    test_participant_mds1 <- aov(MDS1~study_id,mdsMeta);
    pVal_participant_mds1 <- anova(test_participant_mds1)$"Pr(>F)"[1];
    test_participant_mds3 <- aov(MDS3~study_id,mdsMeta);
    pVal_participant_mds3 <- anova(test_participant_mds3)$"Pr(>F)"[1];
    test_participant_mds4 <- aov(MDS4~study_id,mdsMeta);
    pVal_participant_mds4 <- anova(test_participant_mds4)$"Pr(>F)"[1];
    
    
    title <- paste("MDS plot (",wgs," level)",sep="")
    comp1<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
    comp2<-as.character(paste("MDS2 ", (round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_origin_mds2,3),sep=""));
    comp3<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%, p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
    comp4<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%, p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
    
    setwd("../plots")
    
    p <- ggplot(mdsMeta,aes(colour = Origin,shape=Origin))
    tiff(paste("2_mdsPlot_Axes12_",funct,"_",wgs,"_coloredByOrigin.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
    print(p + geom_point(aes(MDS1,MDS2),size = 8) +
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
    
    p <- ggplot(mdsMeta,aes(x=MDS3,y=MDS4,colour = Origin,shape=Origin))
    tiff(paste("2_mdsPlot_Axes34_",funct,"_",wgs,"_coloredByOrigin.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
    print(p + geom_point(size = 8) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp3) + ylab(comp4) +
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
    
    comp1Participant<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%",sep=""));
    comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds1,3),", Origin p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS1))
    tiff(paste("2_mdsPlot_Axes1_",funct,"_",wgs,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
    print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp1ParticipantMDS) + ylab(comp1Participant) +
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

    comp3Participant<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%",sep=""));
    comp3ParticipantMDS <-as.character(paste("Participant\nParticipant p-value = ",format.pval(pVal_participant_mds3,3),", Origin p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS3))
    tiff(paste("2_mdsPlot_Axes3_",funct,"_",wgs,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
    print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
            #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp3ParticipantMDS) + ylab(comp3Participant) +
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
    
    comp4Participant<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%",sep=""));
    comp4ParticipantMDS <-as.character(paste("Participant\nParticipant p-value = ",format.pval(pVal_participant_mds4,3),", Origin p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
    p <- ggplot(mdsMeta,aes(x=study_id,y=MDS4))
    tiff(paste("2_mdsPlot_Axes4_",funct,"_",wgs,"_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
    print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
            #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
            scale_colour_manual(values=c("red2","blue3","magenta4")) +
            xlab(comp4ParticipantMDS) + ylab(comp4Participant) +
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
    setwd("..") 
  }
}


################################# Make pathway MDS plots using distinct and matched samples ########################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
funct <- "wgs"
for(wgs in wgsLevels )
{
  setwd("mds")
  mdsFile <- paste(funct,"_mds_", wgs, "_loggedFiltered.RData",sep="");
  print(mdsFile)
  eigenFile <- paste(funct,"_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
  
  mds <-readRDS(mdsFile);
  #sampleData  <- sampleData[-26,]
  mdsMeta <- merge(sampleData2,mds, by = "row.names")
  
  allStoolPatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "stool"]))
  allSwabPatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "swab"]))
  allTissuePatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "tissue"]))
  matchedTissueSwabPatients <- intersect(allSwabPatients,allTissuePatients)
  matchedPatients <- intersect(matchedTissueSwabPatients,allStoolPatients)
  
  distinctSwabPatients <- allSwabPatients[!allSwabPatients %in% allTissuePatients]
  distinctStoolPatients <- allStoolPatients[!allStoolPatients %in% c(distinctSwabPatients,allTissuePatients)]
  
  mdsMetaDistinct <- mdsMeta[(mdsMeta$study_id %in% allTissuePatients & mdsMeta$Origin %in% "tissue") |(mdsMeta$study_id %in% distinctSwabPatients & mdsMeta$Origin %in% "swab") |(mdsMeta$study_id %in% distinctStoolPatients & mdsMeta$Origin %in% "stool"),];
  mdsMetaMatched <- mdsMeta[(mdsMeta$study_id %in% matchedPatients),];
  
  mdsMetaDistinctPre <- split(mdsMetaDistinct,f=mdsMetaDistinct$visit)$Pre
  mdsMetaDistinctPost <- split(mdsMetaDistinct,f=mdsMetaDistinct$visit)$Post
  
  mdsMetaMatchedPre <- split(mdsMetaMatched,f=mdsMetaMatched$visit)$Pre
  mdsMetaMatchedPost <- split(mdsMetaMatched,f=mdsMetaMatched$visit)$Post
  
  ################ Generate p-values for plot ##############  
  test_distinctPre_origin_mds1 <- aov(MDS1~Origin,mdsMetaDistinct);
  pVal_distinctPre_origin_mds1 <- anova(test_distinctPre_origin_mds1)$"Pr(>F)"[1];
  test_distinctPre_origin_mds2 <- aov(MDS2~Origin,mdsMetaDistinct);
  pVal_distinctPre_origin_mds2 <- anova(test_distinctPre_origin_mds2)$"Pr(>F)"[1];
  test_matchedPre_origin_mds1 <- aov(MDS1~Origin,mdsMetaMatched);
  pVal_matchedPre_origin_mds1 <- anova(test_matchedPre_origin_mds1)$"Pr(>F)"[1];
  test_matchedPre_origin_mds2 <- aov(MDS2~Origin,mdsMetaMatched);
  pVal_matchedPre_origin_mds2 <- anova(test_matchedPre_origin_mds2)$"Pr(>F)"[1];
  
  test_distinctPost_origin_mds1 <- aov(MDS1~Origin,mdsMetaDistinct);
  pVal_distinctPost_origin_mds1 <- anova(test_distinctPost_origin_mds1)$"Pr(>F)"[1];
  test_distinctPost_origin_mds2 <- aov(MDS2~Origin,mdsMetaDistinct);
  pVal_distinctPost_origin_mds2 <- anova(test_distinctPost_origin_mds2)$"Pr(>F)"[1];
  test_matchedPost_origin_mds1 <- aov(MDS1~Origin,mdsMetaMatched);
  pVal_matchedPost_origin_mds1 <- anova(test_matchedPost_origin_mds1)$"Pr(>F)"[1];
  test_matchedPost_origin_mds2 <- aov(MDS2~Origin,mdsMetaMatched);
  pVal_matchedPost_origin_mds2 <- anova(test_matchedPost_origin_mds2)$"Pr(>F)"[1];
  
 
  eigen <-readRDS(eigenFile);
  

  comp1matchedPre <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_matchedPre_origin_mds1,3),sep=""));
  comp2matchedPre <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_matchedPre_origin_mds2,3),sep=""));
  comp1distinctPre <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_distinctPre_origin_mds1,3),sep=""));
  comp2distinctPre <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_distinctPre_origin_mds2,3),sep=""));
  
  comp1matchedPost <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_matchedPost_origin_mds1,3),sep=""));
  comp2matchedPost <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_matchedPost_origin_mds2,3),sep=""));
  comp1distinctPost <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_distinctPost_origin_mds1,3),sep=""));
  comp2distinctPost <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_distinctPost_origin_mds2,3),sep=""));
  
  
  setwd("../plots")
  p <- ggplot(mdsMetaDistinctPre,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_",funct,"_",wgs,"_coloredByOrigin_distinctPre.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1distinctPre) + ylab(comp2distinctPre) +
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

  p <- ggplot(mdsMetaMatchedPre,aes(x=MDS1,y=MDS2,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_",funct,"_",wgs,"_coloredByOrigin_matchedPre.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1matchedPre) + ylab(comp2matchedPre) +
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
  
  p <- ggplot(mdsMetaDistinctPost,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_",funct,"_",wgs,"_coloredByOrigin_distinctPost.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1distinctPost) + ylab(comp2distinctPost) +
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
  
  p <- ggplot(mdsMetaMatchedPost,aes(x=MDS1,y=MDS2,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_",funct,"_",wgs,"_coloredByOrigin_matchedPost.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1matchedPost) + ylab(comp2matchedPost) +
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
  
  setwd("..") 
}

########################## Make Metaphlan 
taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("mds")
  mdsFile <- paste("metaphlan_mds_", taxa, "_loggedFiltered.RData",sep="");
  print(mdsFile)
  eigenFile <- paste("metaphlan_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
  
  mds <-readRDS(mdsFile);
  #sampleData  <- sampleData[-26,]
  mdsMeta <- merge(sampleData2,mds, by = "row.names")
  
  allStoolPatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "stool"]))
  allSwabPatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "swab"]))
  allTissuePatients <- levels(factor(mdsMeta$study_id[mdsMeta$Origin %in% "tissue"]))
  matchedTissueSwabPatients <- intersect(allSwabPatients,allTissuePatients)
  matchedPatients <- intersect(matchedTissueSwabPatients,allStoolPatients)
  
  distinctSwabPatients <- allSwabPatients[!allSwabPatients %in% allTissuePatients]
  distinctStoolPatients <- allStoolPatients[!allStoolPatients %in% c(distinctSwabPatients,allTissuePatients)]
  #distinctPatients <- c(distinctStoolPatients,distinctTissueSwabPatients)
  
  mdsMetaDistinct <- mdsMeta[(mdsMeta$study_id %in% allTissuePatients & mdsMeta$Origin %in% "tissue") |(mdsMeta$study_id %in% distinctSwabPatients & mdsMeta$Origin %in% "swab") |(mdsMeta$study_id %in% distinctStoolPatients & mdsMeta$Origin %in% "stool"),];
  mdsMetaMatched <- mdsMeta[(mdsMeta$study_id %in% matchedPatients),];
  
  mdsMetaDistinctPre <- split(mdsMetaDistinct,f=mdsMetaDistinct$visit)$Pre
  mdsMetaDistinctPost <- split(mdsMetaDistinct,f=mdsMetaDistinct$visit)$Post
  
  mdsMetaMatchedPre <- split(mdsMetaMatched,f=mdsMetaMatched$visit)$Pre
  mdsMetaMatchedPost <- split(mdsMetaMatched,f=mdsMetaMatched$visit)$Post
  
  ################ Generate p-values for plot ##############  
  test_distinctPre_origin_mds1 <- aov(MDS1~Origin,mdsMetaDistinctPre);
  pVal_distinctPre_origin_mds1 <- anova(test_distinctPre_origin_mds1)$"Pr(>F)"[1];
  test_distinctPre_origin_mds2 <- aov(MDS2~Origin,mdsMetaDistinctPre);
  pVal_distinctPre_origin_mds2 <- anova(test_distinctPre_origin_mds2)$"Pr(>F)"[1];
  
  test_matchedPre_origin_mds1 <- aov(MDS1~Origin,mdsMetaMatchedPre);
  pVal_matchedPre_origin_mds1 <- anova(test_matchedPre_origin_mds1)$"Pr(>F)"[1];
  test_matchedPre_origin_mds2 <- aov(MDS2~Origin,mdsMetaMatchedPre);
  pVal_matchedPre_origin_mds2 <- anova(test_matchedPre_origin_mds2)$"Pr(>F)"[1];
  
  test_distinctPost_origin_mds1 <- aov(MDS1~Origin,mdsMetaDistinctPost);
  pVal_distinctPost_origin_mds1 <- anova(test_distinctPost_origin_mds1)$"Pr(>F)"[1];
  test_distinctPost_origin_mds2 <- aov(MDS2~Origin,mdsMetaDistinctPost);
  pVal_distinctPost_origin_mds2 <- anova(test_distinctPost_origin_mds2)$"Pr(>F)"[1];
  
  test_matchedPost_origin_mds1 <- aov(MDS1~Origin,mdsMetaMatchedPost);
  pVal_matchedPost_origin_mds1 <- summary(test_matchedPost_origin_mds1)[[1]][["Pr(>F)"]];
  test_matchedPost_origin_mds2 <- aov(MDS2~Origin,mdsMetaMatchedPost);
  pVal_matchedPost_origin_mds2 <- summary(test_matchedPost_origin_mds2)[[1]][["Pr(>F)"]];
  
  eigen <-readRDS(eigenFile);
  
  comp1matchedPre <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_matchedPre_origin_mds1,3),sep=""));
  comp2matchedPre <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_matchedPre_origin_mds2,3),sep=""));
  comp1distinctPre <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_distinctPre_origin_mds1,3),sep=""));
  comp2distinctPre <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_distinctPre_origin_mds2,3),sep=""));
  
  comp1matchedPost <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_matchedPost_origin_mds1,3),sep=""));
  comp2matchedPost <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_matchedPost_origin_mds2,3),sep=""));
  comp1distinctPost <-as.character(paste("MDS1 ",(round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_distinctPost_origin_mds1,3),sep=""));
  comp2distinctPost <-as.character(paste("MDS2 ",(round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_distinctPost_origin_mds2,3),sep=""));
  
  
  setwd("../plots")
  p <- ggplot(mdsMetaDistinctPre,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_metaphlan_",taxa,"_coloredByOrigin_distinctPre.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1distinctPre) + ylab(comp2distinctPre) +
          ggtitle(taxa) +
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
  
  p <- ggplot(mdsMetaMatchedPre,aes(x=MDS1,y=MDS2,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_metaphlan_",taxa,"_coloredByOrigin_matchedPre.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1matchedPre) + ylab(comp2matchedPre) +
          ggtitle(taxa) +
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
  
  p <- ggplot(mdsMetaDistinctPost,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_metaphlan_",taxa,"_coloredByOrigin_distinctPost.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1distinctPost) + ylab(comp2distinctPost) +
          ggtitle(taxa) +
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
  
  p <- ggplot(mdsMetaMatchedPost,aes(x=MDS1,y=MDS2,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_metaphlan_",taxa,"_coloredByOrigin_matchedPost.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1matchedPost) + ylab(comp2matchedPost) +
          ggtitle(taxa) +
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
  
  setwd("..") 
}

############################### 2. No Tissue ###################################
######################## Classifications ################################
sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("mds")
  mdsFile <- paste("krakenWGS_mds_", taxa, "_loggedFiltered_NoTissue.RData",sep="");
  eigenFile <- paste("krakenWGS_eigenValues_", taxa, "_loggedFiltered_NoTissue.RData",sep="");
  
  mds <-readRDS(mdsFile);
  mdsMeta <- merge(sampleData2,mds, by = "row.names")
  eigen <-readRDS(eigenFile);
  
  ################ Generate p-values for plot ##############  
  test_origin_mds1 <- aov(MDS1~Origin,mdsMeta);
  pVal_origin_mds1 <- anova(test_origin_mds1)$"Pr(>F)"[1];
  test_origin_mds2 <- aov(MDS2~Origin,mdsMeta);
  pVal_origin_mds2 <- anova(test_origin_mds2)$"Pr(>F)"[1];
  test_origin_mds3 <- aov(MDS3~Origin,mdsMeta);
  pVal_origin_mds3 <- anova(test_origin_mds3)$"Pr(>F)"[1];
  test_origin_mds4 <- aov(MDS4~Origin,mdsMeta);
  pVal_origin_mds4 <- anova(test_origin_mds4)$"Pr(>F)"[1];
  
  test_participant_mds1 <- aov(MDS1~study_id,mdsMeta);
  pVal_participant_mds1 <- anova(test_participant_mds1)$"Pr(>F)"[1];
  test_participant_mds3 <- aov(MDS3~study_id,mdsMeta);
  pVal_participant_mds3 <- anova(test_participant_mds3)$"Pr(>F)"[1];
  test_participant_mds4 <- aov(MDS4~study_id,mdsMeta);
  pVal_participant_mds4 <- anova(test_participant_mds4)$"Pr(>F)"[1];
  
  
  title <- paste("MDS plot (",taxa," level)",sep="")
  comp1<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
  comp2<-as.character(paste("MDS2 ", (round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_origin_mds2,3),sep=""));
  comp3<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%, p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
  comp4<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%, p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
  
  setwd("../plots")
  
  p <- ggplot(mdsMeta,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_krakenWGS_",taxa,"_coloredByOrigin_NoTissue.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
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
  
  p <- ggplot(mdsMeta,aes(x=MDS3,y=MDS4,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes34_krakenWGS_",taxa,"_coloredByOrigin_NoTissue.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp3) + ylab(comp4) +
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
  
  comp1Participant<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%",sep=""));
  comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds1,3),", Origin p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS1))
  tiff(paste("2_mdsPlot_Axes1_krakenWGS_",taxa,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1ParticipantMDS) + ylab(comp1Participant) +
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
  
  comp3Participant<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%",sep=""));
  comp3ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds3,3),", Origin p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS3))
  tiff(paste("2_mdsPlot_Axes3_krakenWGS_",taxa,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp3ParticipantMDS) + ylab(comp3Participant) +
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
  
  comp4Participant<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%",sep=""));
  comp4ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds4,3),", Origin p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS4))
  tiff(paste("2_mdsPlot_Axes4_krakenWGS_",taxa,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp4ParticipantMDS) + ylab(comp4Participant) +
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
  setwd("..") 
}

######################## Gene pathways ################################
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")
for(wgs in wgsLevels )
{
  setwd("mds")
  mdsFile <- paste("wgs_mds_", wgs, "_loggedFiltered_NoTissue.RData",sep="");
  print(mdsFile)
  eigenFile <- paste("wgs_eigenValues_", wgs, "_loggedFiltered_NoTissue.RData",sep="");
  
  mds <-readRDS(mdsFile);
  #sampleData  <- sampleData[-26,]
  
  mdsMeta <- merge(sampleData2,mds, by = "row.names")
  
  eigen <-readRDS(eigenFile);
  
  ################ Generate p-values for plot ##############  
  test_origin_mds1 <- aov(MDS1~Origin,mdsMeta);
  pVal_origin_mds1 <- anova(test_origin_mds1)$"Pr(>F)"[1];
  test_origin_mds2 <- aov(MDS2~Origin,mdsMeta);
  pVal_origin_mds2 <- anova(test_origin_mds2)$"Pr(>F)"[1];
  test_origin_mds3 <- aov(MDS3~Origin,mdsMeta);
  pVal_origin_mds3 <- anova(test_origin_mds3)$"Pr(>F)"[1];
  test_origin_mds4 <- aov(MDS4~Origin,mdsMeta);
  pVal_origin_mds4 <- anova(test_origin_mds4)$"Pr(>F)"[1];
  
  test_participant_mds1 <- aov(MDS1~study_id,mdsMeta);
  pVal_participant_mds1 <- anova(test_participant_mds1)$"Pr(>F)"[1];
  test_participant_mds4 <- aov(MDS4~study_id,mdsMeta);
  pVal_participant_mds4 <- anova(test_participant_mds4)$"Pr(>F)"[1];
  test_participant_mds3 <- aov(MDS3~study_id,mdsMeta);
  pVal_participant_mds3 <- anova(test_participant_mds3)$"Pr(>F)"[1];
  
  title <- paste("MDS plot (",wgs," level)",sep="")
  comp1<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%, p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
  comp2<-as.character(paste("MDS2 ", (round(eigen[2],3))*100,"%, p-value = ",format.pval(pVal_origin_mds2,3),sep=""));
  comp3<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%, p-value = ",format.pval(pVal_origin_mds3,3),sep=""));
  comp4<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%, p-value = ",format.pval(pVal_origin_mds4,3),sep=""));
  
  setwd("../plots")
  
  p <- ggplot(mdsMeta,aes(colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes12_wgs_",wgs,"_coloredByOrigin_NoTissue.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(aes(MDS1,MDS2),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
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
  
  p <- ggplot(mdsMeta,aes(x=MDS3,y=MDS4,colour = Origin,shape=Origin))
  tiff(paste("2_mdsPlot_Axes34_wgs_",wgs,"_coloredByOrigin_NoTissue.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350)
  print(p + geom_point(size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp3) + ylab(comp4) +
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
  
  comp1Participant<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%",sep=""));
  comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds1,3),", Origin p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS1))
  tiff(paste("2_mdsPlot_Axes1_wgs_",wgs,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1ParticipantMDS) + ylab(comp1Participant) +
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
  
  comp1Participant<-as.character(paste("MDS1 ", (round(eigen[1],3))*100,"%",sep=""));
  comp1ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds1,3),", Origin p-value = ",format.pval(pVal_origin_mds1,3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS1))
  tiff(paste("2_mdsPlot_Axes1_wgs_",wgs,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp1ParticipantMDS) + ylab(comp1Participant) +
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
  
  comp3Participant<-as.character(paste("MDS3 ", (round(eigen[3],3))*100,"%",sep=""));
  comp3ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds3,3),", Origin p-value = ",format.pval(pVal_origin_mds3),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS3))
  tiff(paste("2_mdsPlot_Axes3_wgs_",wgs,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp3ParticipantMDS) + ylab(comp3Participant) +
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
  
  comp4Participant<-as.character(paste("MDS4 ", (round(eigen[4],3))*100,"%",sep=""));
  comp4ParticipantMDS <-as.character(paste("Participants\nParticipant p-value = ",format.pval(pVal_participant_mds4,3),", Origin p-value = ",format.pval(pVal_origin_mds4),sep=""));
  p <- ggplot(mdsMeta,aes(x=study_id,y=MDS4))
  tiff(paste("2_mdsPlot_Axes4_wgs_",wgs,"_coloredByOrigin_barchart_NoTissue.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
  print(p +geom_boxplot()+ geom_point(aes(colour = Origin,shape=Origin),size = 8) +
          #scale_colour_manual(values=c("#00728F","#DE3A6E")) +
          scale_colour_manual(values=c("red2","blue3","magenta4")) +
          xlab(comp4ParticipantMDS) + ylab(comp4Participant) +
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
  setwd("..") 
}





library(RColorBrewer);
require(ggplot2);
require(nlme);

sampleData <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData$visit <- unlist(strsplit(as.character(sampleData$type),split = "_"))[c(FALSE,TRUE)]
sampleData2 <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleData2)[1] <- "Origin"

taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied","otuRarefied")
#classifierList <- c("qiime","rdpClassifications")
classifierList <- c("qiime")
for(classifier in classifierList)
{
  for(taxa in taxaLevels )
  {
    setwd("mds")
    mdsFile <- paste(classifier,"_mds_", taxa, "_loggedFiltered.RData",sep="");
    eigenFile <- paste(classifier,"_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
    
    mds <- data.frame(readRDS(mdsFile));
    if(classifier %in% "metaphlan")
    {
      mdsMeta <- merge(sampleData2,mds, by = "row.names")
    }
    else{
      mdsMeta <- merge(sampleData,mds, by = "row.names")
    }
    eigen <- readRDS(eigenFile);
    
    pValIndividualList <- numeric(0);
    pValOriginList <- numeric(0);
    pValTimeList <- numeric(0);
    for(i in 1:16)
    {     
      model=as.formula(paste("MDS",i,"~","Origin","+","visit",sep=""));
      print(model);
      
      simpleMod <- gls(model,method="REML",data=mdsMeta);
      mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
      
      #pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,236,lower.tail = FALSE);
      pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
      
      #pValTime <- pf(anova(simpleMod)$"F-value"[3],1,236,lower.tail = FALSE);
      pValTime <- pf(anova(mixedMod)$"F-value"[3],anova(mixedMod)$"numDF"[3],anova(mixedMod)$"denDF"[3],lower.tail = FALSE);
      
      pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
      
      pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
      pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
      pValTimeList[[length(pValTimeList)+1]] <-  pValTime;
    }
    makeTable=data.frame((eigen[1:16])*100,pValOriginList,pValTimeList, pValIndividualList);
    setwd("../statisticalModels/")
    write("Axis\t% explained\tOrigin p-value\tTime p-value\tIndividual p-value",paste("4_anovaMDSindividual_",taxa,"_",classifier,".txt",sep=""));
    write.table(makeTable,paste("4_anovaMDSindividual_",taxa,"_",classifier,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
    
    setwd("../plots")
    
    tiff(paste("4_mdsPvalues_LineGraph_",taxa,"_",classifier,".tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350);
    maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
    minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;

    p <- ggplot(makeTable)
    print(p + geom_point(aes(1:16,-log10(pValOriginList)),colour = "black",size = 6) +
            geom_line(aes(1:16,-log10(pValOriginList)),colour = "black",size = 3) +
            ylim(minLimit,maxLimit)+
            geom_point(aes(1:16,-log10(pValIndividualList)),colour = "red2",size = 6) +
            geom_line(aes(1:16,-log10(pValIndividualList)),colour = "red2",size = 3) +
            geom_point(aes(1:16,-log10(pValTimeList)),colour = "blue3",size = 6) +
            geom_line(aes(1:16,-log10(pValTimeList)),colour = "blue3",size = 3) +
            geom_hline(yintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
            xlab("MDS axes") + ylab("-Log10(pValue)") +
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
    setwd("..")
  } 
}   

######################################## Metaphlan/KrakenWGS ##########################
#classifierList <- c("krakenWGS","metaphlan")
classifierList <- c("krakenWGS")
for(classifier in classifierList)
{
  for(taxa in taxaLevels )
  {
    setwd("mds")
    mdsFile <- paste(classifier,"_mds_", taxa, "_loggedFiltered.RData",sep="");
    eigenFile <- paste(classifier,"_eigenValues_", taxa, "_loggedFiltered.RData",sep="");
    
    mds <- data.frame(readRDS(mdsFile));
    mdsMeta <- merge(sampleData2,mds, by = "row.names")
    eigen <- readRDS(eigenFile);
    
    pValIndividualList <- numeric(0);
    pValOriginList <- numeric(0);
    pValTimeList <- numeric(0);
    for(i in 1:8)
    {     
      model=as.formula(paste("MDS",i,"~","Origin","+","visit",sep=""));
      print(model);
      
      simpleMod <- gls(model,method="REML",data=mdsMeta);
      mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
      
      pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
      pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
      pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
      
      pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
      pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
      pValTimeList[[length(pValTimeList)+1]] <-  pValTime;
    }
    makeTable=data.frame((eigen[1:8])*100,pValOriginList,pValTimeList, pValIndividualList);
    setwd("../statisticalModels/")
    write("Axis\t% explained\tOrigin p-value\tTime p-value\tIndividual p-value",paste("4_anovaMDSindividual_",taxa,"_",classifier,".txt",sep=""));
    write.table(makeTable,paste("4_anovaMDSindividual_",taxa,"_",classifier,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
    
    setwd("../plots")
    
    tiff(paste("4_mdsPvalues_LineGraph_",taxa,"_",classifier,".tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350);
    maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
    minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;
    
    p <- ggplot(makeTable)
    print(p + geom_point(aes(1:8,-log10(pValOriginList)),colour = "black",size = 6) +
            geom_line(aes(1:8,-log10(pValOriginList)),colour = "black",size = 3) +
            #ylim(minLimit,maxLimit)+
            ylim(minLimit,30)+
            geom_point(aes(1:8,-log10(pValIndividualList)),colour = "red2",size = 6) +
            geom_line(aes(1:8,-log10(pValIndividualList)),colour = "red2",size = 3) +
            geom_point(aes(1:8,-log10(pValTimeList)),colour = "blue3",size = 6) +
            geom_line(aes(1:8,-log10(pValTimeList)),colour = "blue3",size = 3) +
            geom_hline(yintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
            xlab("MDS axes") + ylab("-Log10(pValue)") +
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
    setwd("..")
  } 
}
################################### WGS functions ##################################  
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")

for(wgs in wgsLevels )
{
  setwd("mds")
  mdsFile <- paste("wgs_mds_", wgs, "_loggedFiltered.RData",sep="");
  eigenFile <- paste("wgs_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
  
  mds <- data.frame(readRDS(mdsFile));
  mdsMeta <- merge(sampleData2,mds, by = "row.names")
  
  mdsMeta2 <- mdsMeta[!mdsMeta$Origin %in% "tissue",]
  eigen <- readRDS(eigenFile);
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValTimeList <- numeric(0);
  for(i in 1:8)
  {     
    model=as.formula(paste("MDS",i,"~","Origin","+","visit",sep=""));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=mdsMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <-  pValTime;
  }
  makeTable=data.frame((eigen[1:8])*100,pValOriginList,pValTimeList, pValIndividualList);
  setwd("../statisticalModels/")
  write("Axis\t% explained\tOrigin p-value\tTime p-value\tIndividual p-value",paste("4_anovaMDSindividual_",wgs,"_wgs.txt",sep=""));
  write.table(makeTable,paste("4_anovaMDSindividual_",wgs,"_wgs.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
  
  setwd("../plots")
  
  tiff(paste("4_mdsPvalues_LineGraph_",wgs,"_wgs.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350);
  maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
  minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;
  
  p <- ggplot(makeTable)
  print(p + geom_point(aes(1:8,-log10(pValOriginList)),colour = "black",size = 6) +
          geom_line(aes(1:8,-log10(pValOriginList)),colour = "black",size = 3) +
          ylim(minLimit,maxLimit)+
          geom_point(aes(1:8,-log10(pValIndividualList)),colour = "red2",size = 6) +
          geom_line(aes(1:8,-log10(pValIndividualList)),colour = "red2",size = 3) +
          geom_point(aes(1:8,-log10(pValTimeList)),colour = "blue3",size = 6) +
          geom_line(aes(1:8,-log10(pValTimeList)),colour = "blue3",size = 3) +
          geom_hline(yintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
          xlab("MDS axes") + ylab("-Log10(pValue)") +
          ggtitle(wgs) +
          theme_classic(base_size = 28)+
          theme(axis.line=element_line(size=1),
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
  setwd("..")
} 

################################### picrust functions ##################################  
wgsLevels <- c("keggFamilies",
               "keggPathwaysLevel3",
               "keggPathwaysLevel2",
               "keggPathwaysLevel1",
               "metabolickeggPathwaysLevel2",
               "metabolickeggPathwaysLevel3")

for(wgs in wgsLevels )
{
  setwd("mds")
  mdsFile <- paste("picrust_mds_", wgs, "_loggedFiltered.RData",sep="");
  eigenFile <- paste("picrust_eigenValues_", wgs, "_loggedFiltered.RData",sep="");
  
  mds <- data.frame(readRDS(mdsFile));
  mdsMeta <- merge(sampleData,mds, by = "row.names")
  
  #mdsMeta2 <- mdsMeta[!mdsMeta$Origin %in% "tissue",]
  eigen <- readRDS(eigenFile);
  
  pValIndividualList <- numeric(0);
  pValOriginList <- numeric(0);
  pValTimeList <- numeric(0);
  for(i in 1:16)
  {     
    model=as.formula(paste("MDS",i,"~","Origin","+","visit",sep=""));
    print(model);
    
    simpleMod <- gls(model,method="REML",data=mdsMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
    
    pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,236,lower.tail = FALSE);
    pValTime <- pf(anova(simpleMod)$"F-value"[3],1,236,lower.tail = FALSE);
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
    pValTimeList[[length(pValTimeList)+1]] <-  pValTime;
  }
  makeTable=data.frame((eigen[1:16])*100,pValOriginList,pValTimeList, pValIndividualList);
  setwd("../statisticalModels/")
  write("Axis\t% explained\tOrigin p-value\tTime p-value\tIndividual p-value",paste("4_anovaMDSindividual_",wgs,"_picrust.txt",sep=""));
  write.table(makeTable,paste("4_anovaMDSindividual_",wgs,"_picrust.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
  
  setwd("../plots")
  
  tiff(paste("4_mdsPvalues_LineGraph_",wgs,"_picrust.tiff",sep=""),width=200,height=200,units="mm",compression="lzw",res=350);
  maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
  minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;
  
  p <- ggplot(makeTable)
  print(p + geom_point(aes(1:16,-log10(pValOriginList)),colour = "black",size = 6) +
          geom_line(aes(1:16,-log10(pValOriginList)),colour = "black",size = 3) +
          ylim(minLimit,maxLimit)+
          geom_point(aes(1:16,-log10(pValIndividualList)),colour = "red2",size = 6) +
          geom_line(aes(1:16,-log10(pValIndividualList)),colour = "red2",size = 3) +
          geom_point(aes(1:16,-log10(pValTimeList)),colour = "blue3",size = 6) +
          geom_line(aes(1:16,-log10(pValTimeList)),colour = "blue3",size = 3) +
          geom_hline(yintercept=-log10(0.05),colour="grey67",linetype="dashed",size=2) +
          xlab("MDS axes") + ylab("-Log10(pValue)") +
          ggtitle(wgs) +
          theme_classic(base_size = 28)+
          theme(axis.line=element_line(size=1),
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
  setwd("..")
} 




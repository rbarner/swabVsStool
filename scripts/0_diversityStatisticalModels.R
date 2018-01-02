library(vegan)
library(ggplot2)

diversityStatisticalModels <- function(dataSet,level)
{
  bacteriaSwab <- split(dataSet,dataSet$Origin)$SWAB
  bacteriaStool <- split(dataSet,dataSet$Origin)$STOOL
  
  stoolMeans <- round(colMeans(bacteriaStool[,37:41]),3); 
  swabMeans <- round(colMeans(bacteriaSwab[,37:41]),3); 
  
  stoolSDs <- round(colSds(as.matrix(bacteriaStool[,37:41])),3); 
  swabSDs <- round(colSds(as.matrix(bacteriaSwab[,37:41])),3); 
  
  pValOriginList=numeric(0);
  pValIndividualList=numeric(0);
  for(i in 1:5)
  {
    f <- as.formula(paste(names(dataSet)[i+36],"~","Origin","+","visit",sep=""));
    simpleMod <- gls(f,method="REML",data=dataSet);
    mixedMod <- lme(f,method="REML",random=~1|study_id,data=dataSet);
    pValOrigin <- pf(anova(mixedMod)$"F-value"[2],anova(mixedMod)$"numDF"[2],anova(mixedMod)$"denDF"[2],lower.tail = FALSE);
    pValParticipant <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
    
    pValOriginList[[length(pValOriginList)+1]]<- format(pValOrigin,digits=3);
    pValIndividualList[[length(pValIndividualList)+1]]<- format(pValParticipant,digits=3);
  }
  originAdj <- p.adjust(pValOriginList,method = "BH")
  individualAdj <- p.adjust(pValIndividualList,method = "BH")
  makeTable=data.frame(stoolMeans,stoolSDs,swabMeans,swabSDs,pValOriginList,pValIndividualList);
  names(makeTable) <- cbind("stool mean","stool sd","swab mean","swab sd","Origin adj p-value","Participant adj p-value");
  write("DiversityIndice\tstoolMeans\tstoolSDs\tswabMeans\tswabSDs\tOrigin adj p-value\tParticipant adj p-value",paste("../../statisticalModels/0_diversityIndex_",taxa,"_individual_origin_pVal.txt",sep=""));
  write.table(makeTable,paste("../../statisticalModels/0_diversityIndex_",taxa,"_individual_origin_pVal.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE); 
}
##############################################
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]

taxaLevels <- c("phylum","class","order","family","genus","otu")
taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied","otuRarefied")

taxaLevels <- c("otu","phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied","otuRarefied")
for(taxa in taxaLevels )
{
  setwd("data/microbialClassifications/")
  inFileName <- paste("qiime_",taxa, "Level.txt", sep ="")
  print(inFileName)
    
  bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
  bacteria <- t(bacteria)

    # Assign diversity indices
  shannonDiversity <- diversity(bacteria);
  simpsonDiversity <- diversity(bacteria,index = "simpson");
  invSimpsonDiversity <- diversity(bacteria,index= "invsimpson");
  richness <- specnumber(bacteria);
  evenness <- shannonDiversity/log(richness);
  diversityTable <- data.frame(shannonDiversity, simpsonDiversity, invSimpsonDiversity, richness, evenness)
  #write.table(diversityTable,"diversityIndices.txt",row.names = TRUE,col.names = TRUE,sep = "\t")
  write.table(diversityTable,"../../statisticalModels/0_diversityIndices.txt",row.names = TRUE,col.names = TRUE,sep = "\t")

  diversityMeta <- merge(sampleData16S,diversityTable, by = "row.names")
  diversityStatisticalModels(diversityMeta,taxa);
  
  setwd("../..")
  } 

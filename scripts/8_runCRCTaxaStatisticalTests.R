library(nlme)

sampleDataWGS <- read.delim("data/key/mapping_key_WGS.txt",header = TRUE, row.names=1);
names(sampleDataWGS)[1] <- "Origin"



setwd("data/microbialClassfications/")
inFileName <- paste("krakenWGS_speciesLevel.txt", sep ="")
print(inFileName)
bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
bacteria <- t(bacteria)
colnames(bacteria)  <- gsub("\\[","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("\\]","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("\\(","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("\\)","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("-","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("\\+","\\.",colnames(bacteria))
colnames(bacteria)  <- gsub("\\/","\\.",colnames(bacteria))

crcBacteria <- c(".Eubacterium._eligens",".Ruminococcus._torques","Akkermansia_muciniphila","Bacteroides_fragilis","Bacteroides_vulgatus",
                 "Bifidobacterium_longum","Butyvibrio_fibrisolvens","Enterococcus_faecalis","Escherichia_coli","Eubacterium_eligens",
                 "Eubacterium_rectale","Faecalibacterium_prausnitzii","Fusobacterium_nucleatum","Ruminococcus_torques","Streptococcus_thermophilus",
                 "Streptococcus_galloyticus")

bacteriaLogged <- log10((bacteria/rowSums(bacteria))* (sum(colSums(bacteria))/dim(bacteria)[1]) +1)
#myTLogged <-log10((myT/rowSums(myT))* (sum(colSums(myT))/dim(myT)[1]) +1)
#myTLogged <- myTLogged[,(colSums(myTLogged==0)/dim(myTLogged)[1])<=0.75]
bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
bacteriaLogged <- bacteriaLogged[,colnames(bacteriaLogged) %in% crcBacteria]

bacteriaMeta <- merge(sampleDataWGS,bacteriaLogged,by="row.names");
bacteriaSwab <- split(bacteriaMeta,bacteriaMeta$Origin)$swab
bacteriaStool <- split(bacteriaMeta,bacteriaMeta$Origin)$stool
bacteriaTissue <- split(bacteriaMeta,bacteriaMeta$Origin)$tissue

pValIndividualList <- numeric(0);
pValOriginList <- numeric(0);
pValOriginSwabStoolList <- numeric(0);
pValOriginTissueStoolList <- numeric(0);
pValOriginTissueSwabList <- numeric(0);
pValTimeList <- numeric(0);

stoolMeans <- colMeans(bacteriaStool[,names(bacteriaStool) %in% colnames(bacteriaLogged)]); 
swabMeans <- colMeans(bacteriaSwab[,names(bacteriaSwab) %in% colnames(bacteriaLogged)]); 
tissueMeans <- colMeans(bacteriaTissue[,names(bacteriaTissue) %in% colnames(bacteriaLogged)]); 

for(i in (ncol(sampleDataWGS)+2): ncol(bacteriaMeta))
{
  model <- as.formula(paste(names(bacteriaMeta)[i],"~","Origin","+","visit"));
  print(model);
  
  simpleMod <- gls(model,method="REML",data=bacteriaMeta);
  mixedMod <- lme(model,method="REML",random=~1|study_id,data=bacteriaMeta);
  
  anovaModForTukey <- aov(model,data=bacteriaMeta)
  
  pValOrigin <- pf(anova(simpleMod)$"F-value"[2],2,139,lower.tail = FALSE);
  pValOriginSwabStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[1,4];
  pValOriginTissueStool <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[2,4];
  pValOriginTissueSwab <- TukeyHSD(aov(model,data=bacteriaMeta),"Origin")$Origin[3,4];
  pValTime <- pf(anova(simpleMod)$"F-value"[3],1,139,lower.tail = FALSE);
  #pValTime <- TukeyHSD(aov(model,data=bacteriaMeta),"visit")$visit[4];
  pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1,lower.tail = FALSE)
  
  pValOriginList[[length(pValOriginList)+1]] <- pValOrigin;
  pValOriginSwabStoolList[[length(pValOriginSwabStoolList)+1]] <- pValOriginSwabStool;
  pValOriginTissueStoolList[[length(pValOriginTissueStoolList)+1]] <- pValOriginTissueStool;
  pValOriginTissueSwabList[[length(pValOriginTissueSwabList)+1]] <- pValOriginTissueSwab;
  pValIndividualList[[length(pValIndividualList)+1]] <- pValIndividual;
  pValTimeList[[length(pValTimeList)+1]] <- pValTime;
}
originAdj <- p.adjust(pValOriginList, method = "fdr")
originSwabStoolAdj <- p.adjust(pValOriginSwabStoolList, method = "fdr")
originTissueStoolAdj <- p.adjust(pValOriginTissueStoolList, method = "fdr")
originTissueSwabAdj <- p.adjust(pValOriginTissueSwabList, method = "fdr")
individualAdj <- p.adjust(pValIndividualList, method = "fdr")
timeAdj <- p.adjust(pValTimeList, method = "fdr")


makeTable=data.frame(colnames(bacteriaLogged),stoolMeans,swabMeans,tissueMeans,pValOriginList,pValOriginSwabStoolList,pValOriginTissueStoolList,pValOriginTissueSwabList,pValTimeList, pValIndividualList,originAdj,originSwabStoolAdj,originTissueStoolAdj,originTissueSwabAdj,timeAdj,individualAdj);
setwd("../../statisticalModels/")
write("Taxa\tstoolMeans\tswabMeans\ttissueMeans\tOrigin p-value\tOriginSwabStool p-value\tOriginTissueStool p-value\tOriginTissueSwab p-value\tTime p-value\tIndividual p-value\tOrigin (adj p-value)\tOriginSwabStool (adj p-value)\tOriginTissueStool (adj p-value)\tOriginTissueSwab (adj p-value)\tTime (adj p-value)\tIndividual (adj p-value)",paste("8_marginalStatisticalModels_CRCTaxa_species_krakenWGS.txt",sep=""));
write.table(makeTable,paste("8_marginalStatisticalModels_CRCTaxa_species_krakenWGS.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE, row.names = FALSE);
setwd("..")



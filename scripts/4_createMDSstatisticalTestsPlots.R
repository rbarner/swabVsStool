library(RColorBrewer)
require(nlme);

tiff(paste("mdsVariables_allNoTreatment3.tiff",sep=""),width=9600,height=4000,compression="lzw",res=400);
layout(matrix(c(0,1,1,2,2,3,3,4,4,5,5,6,6,
                0,1,1,2,2,3,3,4,4,5,5,6,6,
                7,7,7,7,7,7,7,7,7,7,7,7,7),
              ncol=13,nrow=3,byrow=TRUE));
for(taxa in taxaLevels )
{
  
  mdsFile <- paste(  "mds_", taxa, ".txt",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, ".txt",sep="");
  
  mds <-read.table(mdsFile,header=TRUE,sep="\t",row.names=1);
  mdsMeta <- cbind(sampleData,mds);
  eigen <-read.table(eigenFile,header=TRUE,sep="\t");
  
  pValIndividualList=numeric(0);
  pValOriginList=numeric(0);
  pValTreatmentList = numeric(0);
  pValTimeList = numeric(0);
  
  for(i in 37:52)
  {     
    #model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","treatment","+","visit","+","study_id"));
    #model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","treatment","+","visit"));
    model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","visit"));
    
    print(model);
    #simpleMod=lm(model,data=mdsMeta)
    
    simpleMod <- gls(model,method="REML",data=mdsMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
    
    pValOrigin <- summary(mixedMod)$tTable[2,5];
    #pValTreatment <- summary(mixedMod)$tTable[3,5];
    pValTime <- summary(mixedMod)$tTable[3,5];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1, lower.tail = FALSE)
    
    #pValOrigin=anova(simpleMod)$"Pr(>F)"[1];
    #pValTreatment=anova(simpleMod)$"Pr(>F)"[2];
    #pValTime=anova(simpleMod)$"Pr(>F)"[3];
    #pValIndividual=anova(simpleMod)$"Pr(>F)"[4];
    
    pValOriginList[[length(pValOriginList)+1]]=pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]]=pValIndividual;
    #pValTreatmentList[[length(pValTreatmentList)+1]] = pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] = pValTime;
  }
  #makeTable=data.frame((eigen$x[1:16])*100,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList);
  makeTable=data.frame((eigen$x[1:16])*100,pValOriginList,pValTimeList, pValIndividualList);
  maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
  minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;
  
  write("Axis\t% explained\tOrigin p-value\tTreatment p-value\tTime p-value\tIndividual p-value",paste("anovaMDSindividual2_",taxa,".txt",sep=""));
  write.table(makeTable,paste("anovaMDSindividualNoTreatment_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
  
  plot(1:16,-log10(pValOriginList),xlab="MDS Axes",ylab="-Log10(p-value)",col="black",ylim=c(minLimit,maxLimit),xlim=c(0,16),pch=19,cex=2,cex.lab=1.75,cex.axis = 2.25,cex.main=3, main = paste("(",taxa,")",sep=""))
  lines(1:16,-log10(pValOriginList),col="black",lwd=4)
  
  #points(1:16,-log10(pValTreatmentList),col="green3",pch=19,cex=2)
  #lines(1:16,-log10(pValTreatmentList),col="green3",lwd=4)
  
  points(1:16,-log10(pValIndividualList),col="red",pch=19,cex=2)
  lines(1:16,-log10(pValIndividualList),col="red",lwd=4)
  
  points(1:16,-log10(pValTimeList),col="blue",pch=19,cex=2)
  lines(1:16,-log10(pValTimeList), col="blue",lwd=4)
  
  abline(h=-log10(0.05),col="grey67",lwd=4,lty="dashed")
}
plot.new()
#legend("center",c("Sample Origin","Magnesium Treatment","Patient ID","Sample Timepoint","0.05 p-value"),lty=c(1,1,1,1,2),lwd=c(3,3,3,3,3),col=c("black","green3","red","blue","black"),cex=3,pt.cex=3,box.col="white",horiz=TRUE)
legend("center",c("Sample Origin","Participant","Sample Timepoint","0.05 p-value"),lty=c(1,1,1,2),lwd=c(3,3,3,3),col=c("black","red","blue","grey67"),cex=3,pt.cex=3,box.col="white",horiz=TRUE)
dev.off()


sampleData=read.delim("/Users/Roshonda/Dropbox/FodorLab/dataForProposal/chapter4SwabStool/key/mapping_file_wgs_ordered.txt");
names(sampleData)[2] <- "Origin"
sampleData  <- sampleData[-26,]
#taxaLevels <- c("modules","pathways")
taxaLevels <- c("keggFamilies","pathways_level1",
                "pathways_level2",
                "pathways_level3",
                "pathways_metabolicLevel2",
                "pathways_metabolicLevel3")


taxaLevels <- c("pathways_level1LogFiltered",
                "pathways_level2LogFiltered",
                "pathways_level3LogFiltered",
                "pathways_metabolicLevel2LogFiltered",
                "pathways_metabolicLevel3LogFiltered")

#taxaLevels <- c("keggFamilies_logfiltered")
for(taxa in taxaLevels )
{
  mdsFile <- paste(  "mds_", taxa, "_wTissue.txt",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, "_wTissue.txt",sep="");
  
  mds <-read.table(mdsFile,header=TRUE,sep="\t",row.names=1);
  #mds <- mds[-26,]
  mdsMeta <- cbind(sampleData,mds);
  eigen <-read.table(eigenFile,header=TRUE,sep="\t");
  
  mdsMetaSplit=split(mdsMeta, f=mdsMeta$Origin);
  
  mdsMetaSwab = mdsMetaSplit$"swab";
  mdsMetaSwab=mdsMetaSwab[order(mdsMetaSwab[,6]),];
  swabPatients = mdsMetaSwab$study_id;
  mdsMetaSwabSplits = split(mdsMetaSwab,f=mdsMetaSwab$visit)
  mdsMetaSwabPre = mdsMetaSwabSplits$Pre;
  mdsMetaSwabPost = mdsMetaSwabSplits$Post;
  
  mdsMetaStool = mdsMetaSplit$"stool";
  mdsMetaStool=mdsMetaStool[order(mdsMetaStool[,6]),];
  mdsMetaStoolSub = mdsMetaStool[mdsMetaStool$study_id %in% swabPatients,];
  mdsMetaStoolSub2 = mdsMetaStool[!(mdsMetaStool$study_id %in% swabPatients),];
  mdsMetaStoolSplits = split(mdsMetaStoolSub,f=mdsMetaStoolSub$visit);
  mdsMetaStoolSplits2 = split(mdsMetaStoolSub2,f=mdsMetaStoolSub2$visit);
  mdsMetaStoolPre = mdsMetaStoolSplits$Pre;
  mdsMetaStoolPost = mdsMetaStoolSplits$Post;
  mdsMetaStoolPre2 = mdsMetaStoolSplits2$Pre;
  mdsMetaStoolPost2 = mdsMetaStoolSplits2$Post;
  
  
  mdsMetaPre = rbind(mdsMetaStoolPre2,mdsMetaSwabPre);
  mdsMetaPairedPre = rbind(mdsMetaStoolPre,mdsMetaSwabPre);
  mdsMetaPre2 = mdsMetaPre[-16,];
  mdsMetaPost = rbind(mdsMetaStoolPost2,mdsMetaSwabPost);
  mdsMetaPairedPost = rbind(mdsMetaStoolPost,mdsMetaSwabPost);
  
  
  tiff(paste(taxa,"_5plotMDS.tiff",sep=""),width=4000,height=6300,compression="lzw",res=400);
  layout(matrix(c(0,1,1,0,
                  0,1,1,0,
                  2,2,3,3,
                  2,2,3,3,
                  4,4,5,5,
                  4,4,5,5,
                  6,6,6,6),
                ncol=4,nrow=7,byrow=TRUE));
  originMDSPlots3(mdsMeta, eigen,title="");
  originMDSPlots3(mdsMetaPairedPre, eigen,title="");
  originMDSPlots3(mdsMetaPairedPost, eigen,title="");
  originMDSPlots3(mdsMetaPre, eigen,title="");
  originMDSPlots3(mdsMetaPost, eigen,title="");
  
  plot.new();
  legend("center",pch=c(19,19),col=c('red','blue'),
         legend=c("STOOL","SWAB"),horiz=TRUE,cex=3,pt.cex=3,
         border=NA, box.col="white");
  dev.off();
  
  #originMDStTestTable(mdsMeta,eigen,newfile=paste("ttestMDSOrigin_",taxa,".txt",sep=""));
}


for(taxa in taxaLevels )
{
  mdsFile <- paste(  "mds_", taxa, ".txt",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, ".txt",sep="");
  
  mds <-read.table(mdsFile,header=TRUE,sep="\t",row.names=1);
  mds <- mds[-26,]
  mdsMeta <- cbind(sampleData,mds);
  eigen <-read.table(eigenFile,header=TRUE,sep="\t");
  
  originMDSPlots(mdsMeta, eigen,newfile=paste("originMDS_",taxa,".tiff",sep=""),title=paste(taxa," level PCoA plot",sep=""), location="bottomleft");
  originMDStTestTable(mdsMeta,eigen,newfile=paste("ttestMDSOrigin_",taxa,".txt",sep=""));
  
  pValIndividualList=numeric(0);
  pValOriginList=numeric(0);
  pValTreatmentList = numeric(0);
  pValTimeList = numeric(0);
  
  for(i in 1:15)
  {
    model=as.formula(paste(names(mdsMeta)[i+23],"~","Origin","+","treatment","+","visit"));
    print(model);
    simpleMod=gls(model,method="REML",data=mdsMeta);
    mixedMod=lme(model,method="REML",random=~1|study_id,data=mdsMeta);
    
    pValOrigin=summary(mixedMod)$tTable[2,5];
    pValTreatment=summary(mixedMod)$tTable[3,5];
    pValTime=summary(mixedMod)$tTable[4,5];
    #pValIndividual=anova(simpleMod,mixedMod)$"p-value"[2];
    pValIndividual = 0.5*(1-pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1))
    
    pValOriginList[[length(pValOriginList)+1]]=pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]]=pValIndividual;
    pValTreatmentList[[length(pValTreatmentList)+1]] = pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] = pValTime;
  }
  makeTable=data.frame((eigen$loading[1:15])*100,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList);
  maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
  minLimit  <- summary(unlist(-log10(makeTable)))[1]-0.5;
  #names(makeTable)=cbind("t score","p-value");
  write("Axis\t% explained\tOrigin p-value\tTreatment p-value\tTime p-value\tIndividual p-value",paste("anovaMDSindividual_",taxa,".txt",sep=""));
  write.table(makeTable,paste("anovaMDSindividual_",taxa,".txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
  
  tiff(paste("MDS_pVals_plot_",taxa,".tiff",sep=""),width=2200,height=2200,compression="lzw",res=300);
  
  plot(1:15,-log10(pValOriginList),xlab="MDS Axes",ylab="-Log10(p-value)",col="black",ylim=c(minLimit,maxLimit),xlim=c(0,16),pch=19,cex=2)
  lines(1:15,-log10(pValOriginList),col="black",lwd=4)
  
  points(1:15,-log10(pValTreatmentList),col="green3",pch=19,cex=2)
  lines(1:15,-log10(pValTreatmentList),col="green3",lwd=4)
  
  points(1:15,-log10(pValIndividualList),col="red",pch=19,cex=2)
  lines(1:15,-log10(pValIndividualList),col="red",lwd=4)
  
  points(1:15,-log10(pValTimeList),col="blue",pch=19,cex=2)
  lines(1:15,-log10(pValTimeList), col="blue",lwd=4)
  
  abline(h=-log10(0.05),col="black",lwd=4,lty="dashed")
  
  legend("topright",c("Sample Origin","Magnesium Treatment","Patient ID","Sample Time","0.05 p-value"),lty=c(1,1,1,1,2),lwd=c(3,3,3,3,3),col=c("black","green3","red","blue","black"),cex=1.25)
  dev.off()
  
}

taxaLevels <- c("pathways_level1LogFiltered",
                "pathways_level2LogFiltered",
                "pathways_level3LogFiltered",
                "pathways_metabolicLevel2LogFiltered",
                "pathways_metabolicLevel3LogFiltered")

taxaLevels <- c("keggFamilies_logfiltered")
for(taxa in taxaLevels )
{
  mdsFile <- paste(  "mds_", taxa, "_wTissue.txt",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, "_wTissue.txt",sep="");
  
  mds <-read.table(mdsFile,header=TRUE,sep="\t",row.names=1);
  #mds <- mds[-26,]
  mdsMeta <- cbind(sampleData,mds[,1:10]);
  eigen <-read.table(eigenFile,header=TRUE,sep="\t");
  
  mdsMetaSplit=split(mdsMeta, f=mdsMeta$Origin);
  
  mdsMetaTissue = mdsMetaSplit$"tissue";
  mdsMetaTissue=mdsMetaTissue[order(mdsMetaTissue[,6]),];
  tissuePatients = mdsMetaTissue$study_id;
  mdsMetaTissueSplits = split(mdsMetaTissue,f=mdsMetaTissue$visit)
  mdsMetaTissuePre = mdsMetaTissueSplits$Pre;
  mdsMetaTissuePost = mdsMetaTissueSplits$Post;
  
  mdsMetaSwab = mdsMetaSplit$"swab";
  mdsMetaSwab=mdsMetaSwab[order(mdsMetaSwab[,6]),];
  mdsMetaSwabSub = mdsMetaSwab[mdsMetaSwab$study_id %in% tissuePatients,];
  mdsMetaSwabSub2 = mdsMetaSwab[!(mdsMetaSwab$study_id %in% tissuePatients),];
  swabPatients = mdsMetaSwab$study_id;
  mdsMetaSwabSplits = split(mdsMetaSwabSub,f=mdsMetaSwabSub$visit)
  mdsMetaSwabSplits2 = split(mdsMetaSwabSub2,f=mdsMetaSwabSub2$visit);
  mdsMetaSwabPre = mdsMetaSwabSplits$Pre;
  mdsMetaSwabPost = mdsMetaSwabSplits$Post;
  mdsMetaSwabPre2 = mdsMetaSwabSplits2$Pre;
  mdsMetaSwabPost2 = mdsMetaSwabSplits2$Post;
  
  mdsMetaStool = mdsMetaSplit$"stool";
  mdsMetaStool=mdsMetaStool[order(mdsMetaStool[,6]),];
  mdsMetaStoolSub = mdsMetaStool[mdsMetaStool$study_id %in% tissuePatients,];
  mdsMetaStoolSub2 = mdsMetaStool[!(mdsMetaStool$study_id %in% swabPatients),];
  mdsMetaStoolSplits = split(mdsMetaStoolSub,f=mdsMetaStoolSub$visit);
  mdsMetaStoolSplits2 = split(mdsMetaStoolSub2,f=mdsMetaStoolSub2$visit);
  mdsMetaStoolPre = mdsMetaStoolSplits$Pre;
  mdsMetaStoolPost = mdsMetaStoolSplits$Post;
  mdsMetaStoolPre2 = mdsMetaStoolSplits2$Pre;
  mdsMetaStoolPost2 = mdsMetaStoolSplits2$Post;
  
  
  mdsMetaPre = rbind(mdsMetaStoolPre2,mdsMetaSwabPre2,mdsMetaTissuePre);
  mdsMetaPairedPre = rbind(mdsMetaStoolPre,mdsMetaSwabPre,mdsMetaTissuePre);
  mdsMetaPre2 = mdsMetaPre[-16,];
  mdsMetaPost = rbind(mdsMetaStoolPost2,mdsMetaSwabPost2,mdsMetaTissuePost);
  mdsMetaPairedPost = rbind(mdsMetaStoolPost,mdsMetaSwabPost,mdsMetaTissuePost);
  
  
  tiff(paste(taxa,"_withTissue_5plotMDS2.tiff",sep=""),width=4000,height=6300,compression="lzw",res=400);
  layout(matrix(c(0,1,1,0,
                  0,1,1,0,
                  2,2,3,3,
                  2,2,3,3,
                  4,4,5,5,
                  4,4,5,5,
                  6,6,6,6),
                ncol=4,nrow=7,byrow=TRUE));
  originMDSPlots4(mdsMeta, eigen,title="");
  legend("top",box.col=NA,cex=4,legend=c("A"));
  originMDSPlots4(mdsMetaPairedPre, eigen,title="");
  legend("top",box.col=NA,cex=4,legend=c("B"));
  originMDSPlots4(mdsMetaPairedPost, eigen,title="");
  legend("top",box.col=NA,cex=4,legend=c("C"));
  originMDSPlots4(mdsMetaPre, eigen,title="");
  legend("top",box.col=NA,cex=4,legend=c("D"));
  originMDSPlots4(mdsMetaPost, eigen,title="");
  legend("top",box.col=NA,cex=4,legend=c("E"));
  
  plot.new();
  legend("center",pch=c(19,19,19),col=c('red','cyan3','darkmagenta'),
         legend=c("STOOL","SWAB","TISSUE"),horiz=TRUE,cex=3,pt.cex=3,
         border=NA, box.col="white");
  dev.off();
  
  #originMDStTestTable(mdsMeta,eigen,newfile=paste("ttestMDSOrigin_",taxa,".txt",sep=""));
}



sampleData=read.delim("/Users/Roshonda/Dropbox/FodorLab/dataForProposal/chapter4SwabStool/key/mapping_file_wgs_orderedNoTissue.txt");
names(sampleData)[2] <- "Origin"
#sampleData  <- sampleData[-42,]
sampleData  <- sampleData[-26,]
taxaLevels <- c("pathways_level1LogFiltered",
                "pathways_level2LogFiltered",
                "pathways_level3LogFiltered",
                "pathways_logFiltered","keggFamilies_logFiltered")

mainTitles <- c("pathways_level1",
                "pathways_level2",
                "pathways_level3",
                "pathways_level4","keggFamilies")
k <- 1

tiff(paste("mdsVariables_functions_noTissue_all3.tiff",sep=""),width=9600,height=4000,compression="lzw",res=400);
layout(matrix(c(0,1,1,2,2,3,3,4,4,5,5,
                0,1,1,2,2,3,3,4,4,5,5,
                6,6,6,6,6,6,6,6,6,6,6),
              ncol=11,nrow=3,byrow=TRUE));
for(taxa in taxaLevels )
{
  
  mdsFile <- paste(  "mds_", taxa, "_noTissue.txt",sep="");
  eigenFile <- paste(  "eigenValues_", taxa, "_noTissue.txt",sep="");
  
  mds <-read.table(mdsFile,header=TRUE,sep="\t",row.names=1);
  if(!dim(sampleData)[1] %in% dim(mds)[1])
  {
    #mds <- mds[-27,]
    mds <- mds[-26,]
  }
  mdsMeta <- cbind(sampleData,mds);
  eigen <-read.table(eigenFile,header=TRUE,sep="\t");
  
  pValIndividualList=numeric(0);
  pValOriginList=numeric(0);
  pValTreatmentList = numeric(0);
  pValTimeList = numeric(0);
  
  for(i in 24:39)
  {     
    #model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","treatment","+","visit","+","study_id"));
    #model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","treatment","+","visit"));
    model=as.formula(paste(names(mdsMeta)[i],"~","Origin","+","visit"));
    
    print(model);
    #simpleMod=lm(model,data=mdsMeta)
    
    simpleMod <- gls(model,method="REML",data=mdsMeta);
    simpleMod2 <- lm(model,data=mdsMeta);
    mixedMod <- lme(model,method="REML",random=~1|study_id,data=mdsMeta);
    
    pValOrigin <- anova(simpleMod2)$"Pr(>F)"[1];
    #pValTreatment <- summary(mixedMod)$tTable[3,5];
    pValTime <- anova(simpleMod2)$"Pr(>F)"[2];
    pValIndividual <- pchisq(anova(simpleMod,mixedMod)$"L.Ratio"[2],1, lower.tail = FALSE)
    
    #pValOrigin=anova(simpleMod)$"Pr(>F)"[1];
    #pValTreatment=anova(simpleMod)$"Pr(>F)"[2];
    #pValTime=anova(simpleMod)$"Pr(>F)"[3];
    #pValIndividual=anova(simpleMod)$"Pr(>F)"[4];
    
    pValOriginList[[length(pValOriginList)+1]]=pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]]=pValIndividual;
    #pValTreatmentList[[length(pValTreatmentList)+1]] = pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] = pValTime;
  }
  #makeTable=data.frame((eigen$loading[1:16])*100,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList);
  makeTable=data.frame((eigen$loading[1:16])*100,pValOriginList,pValTimeList, pValIndividualList);
  maxLimit  <- summary(unlist(-log10(makeTable)))[6]+3;
  minLimit  <- summary(unlist(-log10(makeTable)))[1]-2;
  
  #write("Axis\t% explained\tOrigin p-value\tTreatment p-value\tTime p-value\tIndividual p-value",paste("anovaMDSindividual5_",taxa,".txt",sep=""));
  write("Axis\t% explained\tOrigin p-value\tTime p-value\tIndividual p-value",paste("anovaMDSindividualFunctions_",taxa,"_wTissue.txt",sep=""));
  write.table(makeTable,paste("anovaMDSindividualFunctions_",taxa,"_wTissue.txt",sep=""),quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
  
  plot(1:16,-log10(pValOriginList),xlab="MDS Axes",ylab="-Log10(p-value)",col="black",ylim=c(minLimit,maxLimit),xlim=c(0,16),pch=19,cex=2,cex.lab=1.75,cex.axis = 2.25,cex.main=3, main = paste("(",mainTitles[k],")",sep=""))
  lines(1:16,-log10(pValOriginList),col="black",lwd=4)
  
  #points(1:16,-log10(pValTreatmentList),col="green3",pch=19,cex=2)
  #lines(1:16,-log10(pValTreatmentList),col="green3",lwd=4)
  
  points(1:16,-log10(pValIndividualList),col="red",pch=19,cex=2)
  lines(1:16,-log10(pValIndividualList),col="red",lwd=4)
  
  points(1:16,-log10(pValTimeList),col="blue",pch=19,cex=2)
  lines(1:16,-log10(pValTimeList), col="blue",lwd=4)
  
  abline(h=-log10(0.05),col="grey67",lwd=4,lty="dashed")
  
  k <- k+1
}
plot.new()
#legend("center",c("Sample Origin","Magnesium Treatment","Patient ID","Sample Timepoint","0.05 p-value"),lty=c(1,1,1,1,2),lwd=c(3,3,3,3,3),col=c("black","green3","red","blue","black"),cex=3,pt.cex=3,box.col="white",horiz=TRUE)
legend("center",c("Sample Origin","Participant","Sample Timepoint","0.05 p-value"),lty=c(1,1,1,2),lwd=c(3,3,3,3),col=c("black","red","blue","grey67"),cex=3,pt.cex=3,box.col="white",horiz=TRUE)
dev.off()

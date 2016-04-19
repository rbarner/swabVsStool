############ functions used #######################
originMDStTestTable = function(dataFrame, componentPercentage, newfile)
{
  tValList=numeric(0);
  pValList=numeric(0);
  for(i in (dim(sampleData)[2]+1): length(names(dataFrame)))
  {
    f<-as.formula(paste(names(dataFrame)[i],"~","Origin"));
    ttest=t.test(f,dataFrame);
    tVal=ttest$statistic[[1]];
    pVal=ttest$p.value[[1]];
    tValList[[length(tValList)+1]]= tVal;
    pValList[[length(pValList)+1]] = pVal;  
  }
  makeTable=data.frame(componentPercentage$loading*100,tValList,pValList);
  names(makeTable)=cbind("t score","p-value");
  write("Component\tCumulative % explained\tt stat\tp-value",newfile);
  write.table(makeTable,newfile,quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
}

individualMDStTestTable = function(dataFrame, componentPercentage, newfile)
{
  pValIndividualList=numeric(0);
  pValOriginList=numeric(0);
  pValTreatmentList = numeric(0);
  pValTimeList = numeric(0);
  
  for(i in (dim(sampleData)[2]+1):dim(sampleData)[2]+14)
  {
    model=as.formula(paste(names(dataFrame)[i],"~","Origin","+","treatment","+","visit","+","study_id"));
    print(model);
    #simpleMod=gls(as.formula(paste(names(dataFrame)[i],"~","Origin","+","treatment")),method="REML",data=dataFrame);
    #mixedMod=lme(as.formula(paste(names(dataFrame)[i],"~","Origin","+","treatment")),method="REML",random=~1|study_id,data=dataFrame);
    simpleMod=lm(model,data=dataFrame)
    
    #pValOrigin=anova(mixedMod)$"p-value"[2];
    #pValTreatment=anova(mixedMod)$"p-value"[3];
    #pValIndividual=anova(simpleMod,mixedMod)$"p-value"[2];
    
    pValOrigin=anova(simpleMod)$"Pr(>F)"[1];
    pValTreatment=anova(simpleMod)$"Pr(>F)"[2];
    pValTime=anova(simpleMod)$"Pr(>F)"[3];
    pValIndividual=anova(simpleMod)$"Pr(>F)"[4];
    
    pValOriginList[[length(pValOriginList)+1]]=pValOrigin;
    pValIndividualList[[length(pValIndividualList)+1]]=pValIndividual;
    pValTreatmentList[[length(pValTreatmentList)+1]] = pValTreatment;
    pValTimeList[[length(pValTimeList)+1]] = pValTime;
  }
  makeTable=data.frame((componentPercentage$x[1:15])*100,pValOriginList,pValTreatmentList,pValTimeList, pValIndividualList);
  #names(makeTable)=cbind("t score","p-value");
  write("Axis\t% explained\tOrigin p-value\tTreatment p-value\tVisit\tIndividual p-value",newfile);
  write.table(makeTable,newfile,quote=FALSE, sep="\t",append=TRUE, col.names=FALSE);
}


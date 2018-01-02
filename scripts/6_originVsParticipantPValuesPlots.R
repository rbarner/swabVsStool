rm(list=ls())
library(ggplot2)
library(ggrepel)
library(dplyr)

################################# 16S classifications ########################################################
taxaLevels <- c("phylumRarefied","classRarefied","orderRarefied","familyRarefied","genusRarefied")
for(taxa in taxaLevels )
{
  pValuesTable <- read.delim(paste("statisticalModels/5_marginalStatisticalModels_taxaByTaxa_",taxa,"_mixedMod.txt",sep=""))
  
  pValuesTable$Origin_c <- -log10(pValuesTable$Origin.p.value)
  pValuesTable$Origin_c  <- ifelse(pValuesTable$stoolMeans > pValuesTable$swabMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)

  p <- ggplot(pValuesTable)
  
  tiff(paste("plots/6_",taxa,"_OriginVsParticipant.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       colour = ifelse(pValuesTable$Origin..adj.p.value. < 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. < 1e-20,"red2","gray53")),
                       size = ifelse(pValuesTable$Origin..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 1e-20,5,3))) +
          xlab("-Log10(Stool Vs Swab p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,Origin..adj.p.value. < 0.1 | Individual..adj.p.value. < 1e-20),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()
}

############################################# Metaphlan ##########################################################
taxaLevels <- c("phylum","class","order","family","genus")
for(taxa in taxaLevels )
{
  setwd("statisticalModels/")
  pValuesTable <- read.delim(paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_metaphlan.txt",sep=""))
  
  pValuesTable$Origin_c <- -log10(pValuesTable$Origin.p.value)
  pValuesTable$Origin_c  <- ifelse(pValuesTable$stoolMeans > pValuesTable$swabMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)
  
  setwd("../plots/")
  p <- ggplot(pValuesTable)
  
  tiff(paste("6_",taxa,"_OriginVsParticipant_metaphlan.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       colour = ifelse(pValuesTable$Origin..adj.p.value. < 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. < 0.1,"red2","gray53")),
                       size = ifelse(pValuesTable$Origin..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 0.1,5,3))) +
          xlab("-Log10(Stool Vs Swab p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,Origin..adj.p.value. < 0.1 | Individual..adj.p.value. < 0.1),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()
  setwd("..")
}

############################################# kraken WGS ##########################################################
taxaLevels <- c("phylum","class","order","family","genus","species")
for(taxa in taxaLevels )
{
  setwd("statisticalModels/")
  pValuesTable <- read.delim(paste("5_marginalStatisticalModels_taxaByTaxa_",taxa,"_krakenWGS.txt",sep=""))
  
  pValuesTable$Origin_c <- -log10(pValuesTable$Origin.p.value)
  #pValuesTable$Origin_c  <- ifelse(pValuesTable$stoolMeans > pValuesTable$swabMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)
  
  setwd("../plots/")
  p <- ggplot(pValuesTable)
  
  tiff(paste("6_",taxa,"_OriginVsParticipant_krakenWGS_allSampleOrigins.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       colour = ifelse(pValuesTable$Origin..adj.p.value. <= 0.1 & pValuesTable$Individual..adj.p.value. > 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$Origin..adj.p.value. > 0.1,"blue3",
                                              ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$Origin..adj.p.value. <= 0.1,"magenta4",
                                                     "gray53"))),
                       size = ifelse(pValuesTable$Origin..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 0.1,5,3))) +
          xlab("-Log10(Sample Origin p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,Origin..adj.p.value. < 1e-7 | Individual..adj.p.value. < 1e-7),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()

  
  pValuesTable$Origin_c <- -log10(pValuesTable$OriginSwabStool.p.value)
  pValuesTable$Origin_c  <- ifelse(pValuesTable$stoolMeans > pValuesTable$swabMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)
  
  p <- ggplot(pValuesTable)
  
  tiff(paste("6_",taxa,"_OriginVsParticipant_krakenWGS_swabStool.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       #colour = ifelse(pValuesTable$OriginSwabStool..adj.p.value. < 0.1,"red2",
                                       #ifelse(pValuesTable$Individual..adj.p.value. < 0.1,"red2","gray53")),
                       colour = ifelse(pValuesTable$OriginSwabStool..adj.p.value. <= 0.1 & pValuesTable$Individual..adj.p.value. > 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginSwabStool..adj.p.value. > 0.1,"blue3",
                                              ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginSwabStool..adj.p.value. <= 0.1,"magenta4",
                                                     "gray53"))),
                       size = ifelse(pValuesTable$OriginSwabStool..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 0.1,5,3))) +
          xlab("-Log10(Stool vs Swab p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,OriginSwabStool..adj.p.value. < 1e-7 | Individual..adj.p.value. < 1e-7),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()
  
  pValuesTable$Origin_c <- -log10(pValuesTable$OriginTissueSwab.p.value)
  pValuesTable$Origin_c  <- ifelse(pValuesTable$swabMeans > pValuesTable$tissueMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)
  
  p <- ggplot(pValuesTable)
  
  tiff(paste("6_",taxa,"_OriginVsParticipant_krakenWGS_swabTissue.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       #colour = ifelse(pValuesTable$OriginTissueSwab..adj.p.value. < 0.1,"red2",
                                      # ifelse(pValuesTable$Individual..adj.p.value. < 0.1,"red2","gray53")),
                       colour = ifelse(pValuesTable$OriginTissueSwab..adj.p.value. <= 0.1 & pValuesTable$Individual..adj.p.value. > 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginTissueSwab..adj.p.value. > 0.1,"blue3",
                                              ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginTissueSwab..adj.p.value. <= 0.1,"magenta4",
                                                     "gray53"))),
                       size = ifelse(pValuesTable$OriginTissueSwab..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 0.1,5,3))) +
          xlab("-Log10(Swab vs Tissue p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,OriginTissueSwab..adj.p.value. < 1e-7 | Individual..adj.p.value. < 1e-7),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()
  
  pValuesTable$Origin_c <- -log10(pValuesTable$OriginTissueStool.p.value)
  pValuesTable$Origin_c  <- ifelse(pValuesTable$stoolMeans > pValuesTable$tissueMeans,pValuesTable$Origin_c,-1*pValuesTable$Origin_c)
  
  pValuesTable$Individual_c <- -log10(pValuesTable$Individual.p.value)
  
  p <- ggplot(pValuesTable)
  
  tiff(paste("6_",taxa,"_OriginVsParticipant_krakenWGS_stoolTissue.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300);
  print(p + geom_point(aes(Origin_c,Individual_c),
                       #colour = ifelse(pValuesTable$OriginTissueStool..adj.p.value. < 0.1,"red2",
                                     #  ifelse(pValuesTable$Individual..adj.p.value. < 0.1,"red2","gray53")),
                       colour = ifelse(pValuesTable$OriginTissueStool..adj.p.value. <= 0.1 & pValuesTable$Individual..adj.p.value. > 0.1,"red2",
                                       ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginTissueStool..adj.p.value. > 0.1,"blue3",
                                              ifelse(pValuesTable$Individual..adj.p.value. <= 0.1 & pValuesTable$OriginTissueStool..adj.p.value. <= 0.1,"magenta4",
                                                     "gray53"))),
                       size = ifelse(pValuesTable$OriginTissueStool..adj.p.value. < 0.1,5,
                                     ifelse(pValuesTable$Individual..adj.p.value. < 0.1,5,3))) +
          xlab("-Log10(Stool vs Tissue p-value)")+ ylab("-Log10(Participant p-value)") +
          geom_hline(yintercept=0)+
          geom_vline(xintercept=0)+
          geom_text_repel(data = filter(pValuesTable,OriginTissueStool..adj.p.value. < 1e-7| Individual..adj.p.value. < 1e-7),
                          aes(Origin_c,Individual_c,label = Taxa),size=5,force=8,box.padding= unit(0.5,"lines"),fontface='bold') +
          theme_bw(base_size = 24)+
          theme(axis.line=element_line(size=1),
                axis.ticks=element_line(size=1),
                axis.text=element_text(face="bold",size=16),
                text=element_text(face="bold",size=24)
          )
  )
  graphics.off()
  
  setwd("..")
}

##################################################
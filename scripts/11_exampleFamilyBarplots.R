rm(list=ls())

library(ggplot2)
library(ggrepel)
library(dplyr)

setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")
sampleData16S <- read.delim("data/key/mapping_key_16S.txt",header = TRUE, row.names=1);
sampleData16S$visit <- unlist(strsplit(as.character(sampleData16S$type),split = "_"))[c(FALSE,TRUE)]


setwd("data/microbialClassifications/")

inFileName <- paste("qiime_familyRarefiedLevel.txt", sep ="")
print(inFileName)
bacteria <-read.delim(inFileName,header = TRUE, row.names=1)
bacteria <- t(bacteria)
bacteriaLogged <- log10(bacteria +1)
bacteriaLogged <- bacteriaLogged[,(colSums(bacteriaLogged==0)/nrow(bacteriaLogged))<=0.75]
colnames(bacteriaLogged) <- gsub(" ","_",colnames(bacteriaLogged))
colnames(bacteriaLogged) <- gsub("/","_",colnames(bacteriaLogged))
colnames(bacteriaLogged) <- gsub("-","_",colnames(bacteriaLogged))
colnames(bacteriaLogged) <- gsub("\\]","\\.",colnames(bacteriaLogged))
colnames(bacteriaLogged) <- gsub("\\[","\\.",colnames(bacteriaLogged))

bacteriaMeta <- merge(sampleData16S,bacteriaLogged,by="row.names");
p <- ggplot(bacteriaMeta,aes(x=study_id,y=Thermaceae))
tiff(paste("plots/11_Thermaceae_coloredByOrigin_barchart.tiff",sep=""),width=4200,height=3200,compression="lzw",res=300)
print(p +geom_boxplot(size=1.5)+ geom_point(aes(colour = Origin,shape=Origin),size = 7) +
        scale_colour_manual(values=c("red2","blue3")) +
        xlab("Participant") + ylab("Log-Normalized Abundance") +
        ggtitle("Thermaceae") +
        theme_classic(base_size = 28)+
        theme(axis.line=element_line(size=1),
              axis.ticks.y=element_line(size=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_text(face="bold",size=24),
              axis.text.x=element_blank(),
              text=element_text(face="bold",size=32),
              legend.text=element_text(face="bold",size=24),
              legend.position="bottom",
              legend.title=element_blank()
        )+
        theme(axis.line.x = element_line(color="black", size = 2),
              axis.line.y = element_line(color="black", size = 2)
        )
)
graphics.off()

p <- ggplot(bacteriaMeta,aes(x=study_id,y=Desulfovibrionaceae))
tiff(paste("plots/Desulfovibrionaceae_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
print(p +geom_boxplot(size=1.5)+ geom_point(aes(colour = Origin,shape=Origin),size = 7) +
        scale_colour_manual(values=c("red2","blue3")) +
        xlab("Participant") + ylab("Log-Normalized Abundance") +
        ggtitle("Desulfovibrionaceae") +
        theme_classic(base_size = 28)+
        theme(axis.line=element_line(size=1),
              axis.ticks.y=element_line(size=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_text(face="bold",size=24),
              axis.text.x=element_blank(),
              text=element_text(face="bold",size=32),
              legend.text=element_text(face="bold",size=24),
              legend.position="bottom",
              legend.title=element_blank()
        )+
        theme(axis.line.x = element_line(color="black", size = 2),
              axis.line.y = element_line(color="black", size = 2)
        )
)
graphics.off()

p <- ggplot(bacteriaMeta,aes(x=study_id,y=Enterobacteriaceae))
tiff(paste("plots/Enterobacteriaceae_coloredByOrigin_barchart.tiff",sep=""),width=400,height=200,units="mm",compression="lzw",res=350)
print(p +geom_boxplot(size=1.5)+ geom_point(aes(colour = Origin,shape=Origin),size = 7) +
        scale_colour_manual(values=c("red2","blue3")) +
        xlab("Participant") + ylab("Log-Normalized Abundance") +
        ggtitle("Enterobacteriaceae") +
        theme_classic(base_size = 28)+
        theme(axis.line=element_line(size=1),
              axis.ticks.y=element_line(size=1),
              axis.ticks.x=element_blank(),
              axis.text.y=element_text(face="bold",size=24),
              axis.text.x=element_blank(),
              text=element_text(face="bold",size=32),
              legend.text=element_text(face="bold",size=24),
              legend.position="bottom",
              legend.title=element_blank()
        )+
        theme(axis.line.x = element_line(color="black", size = 2),
              axis.line.y = element_line(color="black", size = 2)
        )
)
graphics.off()

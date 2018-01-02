library(phyloseq)
library(ggplot2)
library(plyr)
setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")

otuData <- read.table("data/otu_table_10000.txt",row.names = 1, header=TRUE)
otuData2 <- as.matrix(otuData)
otuTable <- otu_table(otuData2,taxa_are_rows = TRUE)

taxaData <- read.table("data/taxonomyTableCorrected.txt",row.names = 1, header=TRUE)
taxaData2 <- as.matrix(taxaData)
taxaTable <- tax_table(taxaData2)

phyloseqData <- phyloseq(otuTable,taxaTable)

phylumData <- tax_glom(phyloseqData, taxrank = "Phylum")
phylumTaxa <- as.data.frame(tax_table(phylumData))
phylumCounts <- as.data.frame(otu_table(phylumData))
phylumCounts$taxaNames <- phylumTaxa$Phylum
phylumCounts2 <- ddply(phylumCounts,"taxaNames",numcolwise(sum))
row.names(phylumCounts2) <- phylumCounts2$taxaNames
phylumCounts <- phylumCounts2[,-1]
write.table(phylumCounts,paste("data/microbialClassifications/qiime_phylumRarefiedLevel.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);



classData <- tax_glom(phyloseqData, taxrank = "Class")
classTaxa <- as.data.frame(tax_table(classData))
classCounts <- as.data.frame(otu_table(classData))
classCounts$taxaNames <- classTaxa$Class
classCounts2 <- ddply(classCounts,"taxaNames",numcolwise(sum))
row.names(classCounts2) <- classCounts2$taxaNames
classCounts <- classCounts2[,-1]
write.table(classCounts,paste("data/microbialClassifications/qiime_classRarefiedLevel.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


orderData <- tax_glom(phyloseqData, taxrank = "Order")
orderTaxa <- as.data.frame(tax_table(orderData))
orderCounts <- as.data.frame(otu_table(orderData))
orderCounts$taxaNames <- orderTaxa$Order
orderCounts2 <- ddply(orderCounts,"taxaNames",numcolwise(sum))
row.names(orderCounts2) <- orderCounts2$taxaNames
orderCounts <- orderCounts2[,-1]
write.table(orderCounts,paste("data/microbialClassifications/qiime_orderRarefiedLevel.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


familyData <- tax_glom(phyloseqData, taxrank = "Family")
familyTaxa <- as.data.frame(tax_table(familyData))
familyCounts <- as.data.frame(otu_table(familyData))
familyCounts$taxaNames <- familyTaxa$Family
familyCounts2 <- ddply(familyCounts,"taxaNames",numcolwise(sum))
row.names(familyCounts2) <- familyCounts2$taxaNames
familyCounts <- familyCounts2[,-1]
write.table(familyCounts,paste("data/microbialClassifications/qiime_familyRarefiedLevel.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


genusData <- tax_glom(phyloseqData, taxrank = "Genus")
genusTaxa <- as.data.frame(tax_table(genusData))
genusCounts <- as.data.frame(otu_table(genusData))
genusCounts$taxaNames <- genusTaxa$Genus
genusCounts2 <- ddply(genusCounts,"taxaNames",numcolwise(sum))
row.names(genusCounts2) <- genusCounts2$taxaNames
genusCounts <- genusCounts2[,-1]
write.table(genusCounts,paste("data/microbialClassifications/qiime_genusRarefiedLevel.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);

############################################### Not Rarefied #####################

setwd("C://Users/Roshonda/swabVsStoolMicrobiome/")

otuData <- read.table("data/otu_table_withouttax.txt",row.names = 1, header=TRUE)
otuData2 <- as.matrix(otuData)
otuTable <- otu_table(otuData2,taxa_are_rows = TRUE)

taxaData <- read.table("data/taxonomyTableCorrected.txt",row.names = 1, header=TRUE)
taxaData2 <- as.matrix(taxaData)
taxaTable <- tax_table(taxaData2)

phyloseqData <- phyloseq(otuTable,taxaTable)

phylumData <- tax_glom(phyloseqData, taxrank = "Phylum")
phylumTaxa <- as.data.frame(tax_table(phylumData))
phylumCounts <- as.data.frame(otu_table(phylumData))
phylumCounts$taxaNames <- phylumTaxa$Phylum
phylumCounts2 <- ddply(phylumCounts,"taxaNames",numcolwise(sum))
row.names(phylumCounts2) <- phylumCounts2$taxaNames
phylumCounts <- phylumCounts2[,-1]
write.table(phylumCounts,paste("data/microbialClassifications/qiime_phylum2Level.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);



classData <- tax_glom(phyloseqData, taxrank = "Class")
classTaxa <- as.data.frame(tax_table(classData))
classCounts <- as.data.frame(otu_table(classData))
classCounts$taxaNames <- classTaxa$Class
classCounts2 <- ddply(classCounts,"taxaNames",numcolwise(sum))
row.names(classCounts2) <- classCounts2$taxaNames
classCounts <- classCounts2[,-1]
write.table(classCounts,paste("data/microbialClassifications/qiime_class2Level.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


orderData <- tax_glom(phyloseqData, taxrank = "Order")
orderTaxa <- as.data.frame(tax_table(orderData))
orderCounts <- as.data.frame(otu_table(orderData))
orderCounts$taxaNames <- orderTaxa$Order
orderCounts2 <- ddply(orderCounts,"taxaNames",numcolwise(sum))
row.names(orderCounts2) <- orderCounts2$taxaNames
orderCounts <- orderCounts2[,-1]
write.table(orderCounts,paste("data/microbialClassifications/qiime_order2Level.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


familyData <- tax_glom(phyloseqData, taxrank = "Family")
familyTaxa <- as.data.frame(tax_table(familyData))
familyCounts <- as.data.frame(otu_table(familyData))
familyCounts$taxaNames <- familyTaxa$Family
familyCounts2 <- ddply(familyCounts,"taxaNames",numcolwise(sum))
row.names(familyCounts2) <- familyCounts2$taxaNames
familyCounts <- familyCounts2[,-1]
write.table(familyCounts,paste("data/microbialClassifications/qiime_family2Level.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);


genusData <- tax_glom(phyloseqData, taxrank = "Genus")
genusTaxa <- as.data.frame(tax_table(genusData))
genusCounts <- as.data.frame(otu_table(genusData))
genusCounts$taxaNames <- genusTaxa$Genus
genusCounts2 <- ddply(genusCounts,"taxaNames",numcolwise(sum))
row.names(genusCounts2) <- genusCounts2$taxaNames
genusCounts <- genusCounts2[,-1]
write.table(genusCounts,paste("data/microbialClassifications/qiime_genus2Level.txt",sep=""),quote=FALSE, sep="\t",col.names=TRUE, row.names = TRUE);

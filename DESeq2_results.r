library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(ellipse)
library(ggplot2)
library(RColorBrewer)
library(gplots)

######### DESeq2 differential expression for only coral sequences for all comparsion groups #########

sample_table=read_csv("./Input_files/Metadata_CBASS_HH.csv")
samplenames<-pull(sample_table, Sample)
filenames<-paste0(samplenames, '/quant.sf')

#provides column headers of the sample names
names(filenames) = pull(sample_table, Sample)
names(filenames)

#split sample table for desired comparison group e.g. temp 30 C vs 33 C samples n=8 

#splitsample_table = subset(sample_table, sample_table$Temperature == "30") #change value accordingly
splitsample_table = subset(sample_table, sample_table$Temperature == "30" & sample_table$Treatment == "Normoxia"| sample_table$Temperature == "33" & sample_table$Treatment == "Normoxia")

splitfilenames = subset(filenames, names(filenames) %in% splitsample_table$Sample)

setwd("./Input_files/quants_500_all")
split_salmon_data = tximport(files = splitfilenames, 
                          type="salmon", 
                          txOut = TRUE, 
                          dropInfReps = TRUE)

subsetlist <- read.table ("./Input_files/scler_500_contigs_uniq.txt", header=F)

#tibble to dataframe
splitsample_table = as.data.frame(splitsample_table)

#make another column with temp and treatment combined
splitsample_table$Condition=paste(splitsample_table$Treatment, splitsample_table$Temperature, sep="_")
splitsample_table$Condition<- as.factor(splitsample_table$Condition)

#insert split count data into deseq
dds=DESeqDataSetFromTximport(txi=split_salmon_data,
                               colData=splitsample_table, 
                               design= ~Condition)

#make rownames match colnames
rownames(splitsample_table)<-colnames(split_salmon_data$counts)

#remove non-coral transcripts
genesTokeep <- which(rownames(dds) %in% subsetlist$V1)
ddsC <- dds[genesTokeep, ]
rownames(ddsC)

ddsC=estimateSizeFactors(ddsC)
ddsC=DESeq(ddsC) 
res=results(ddsC, contrast=c("Condition", "Normoxia_33","Normoxia_30"), lfcThreshold = 0.0, alpha=0.05) #resultsNames(ddsC) to check contrast terms
message(summary(res, alpha =0.05))
results=as.data.frame(subset(res, res$padj<0.05))
write.table(res, "DE_output_33C_30C_comp.txt", sep = "\t", quote = FALSE, row.names= TRUE) #change accordingly

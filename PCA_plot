library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(ggplot2)

######### pca plot using all samples count data #########

#tximport data for deseq2

#not data.frame but tibble
sample_table=read_csv("./Input_files/Metadata_CBASS_HH.csv")
samplenames<-pull(sample_table, Sample)
filenames<-paste0(samplenames, '/quant.sf')

#provides column headers of the sample names
names(filenames) = pull(sample_table, Sample)
names(filenames)

setwd("./Input_files/quants_500_all")
#imports quant.sf data including numReads to use for count data
salmon_data = tximport(files = filenames, 
                       type="salmon", 
                       txOut = TRUE, 
                       dropInfReps = TRUE) #use dropInfReps to use older version than 0.8 salmon

#turn tibble into dataframe to use in DESeq2
sample_table = as.data.frame(sample_table)

#make another column with temp and treatment combined
sample_table$Condition=paste(sample_table$Treatment, sample_table$Temperature, sep="_")
sample_table$Condition<- as.factor(sample_table$Condition)
sample_table$Treatment<- as.factor(sample_table$Treatment)
sample_table$Temperature<- as.factor(sample_table$Temperature)

#put data into deseq
#need to use all salmon data first then deseq chooses to use counts and average transcript lengths for analysis

dds=DESeqDataSetFromTximport(txi=salmon_data,
                             colData=sample_table, 
                             design= ~Condition)

#estimate size factors e.g. each library (sample) depth
dds=estimateSizeFactors(dds)
vst=varianceStabilizingTransformation(dds)

#return data from plotPCA (DESeq package for simple pca plots) to customise
pca_data<-plotPCA(vst, intgroup= 'Condition', returnData = TRUE)

COL2=c("lightblue", "yellow","peachpuff","pink2", "blue", "gold2", "darkorange", "red")
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$Treatment<- sample_table$Treatment[match(pca_data$name, sample_table$Sample)]

#change info for key whether referring to deoxygenation or hypoxia treatment
pca_data$Treatment<- gsub("Deoxygenation", "Hypoxia", pca_data$Treatment)
colnames(pca_data)<-gsub("group", "Group", colnames(pca_data))
pca_data$Group<-gsub("Deoxygenation", "Hypoxia", pca_data$Group)
pca_data$Group<-gsub("Hypoxia_30", "HH_30", pca_data$Group)
pca_data$Group<-gsub("Hypoxia_33", "HH_33", pca_data$Group)
pca_data$Group<-gsub("Hypoxia_36", "HH_36", pca_data$Group)
pca_data$Group<-gsub("Hypoxia_39", "HH_39", pca_data$Group)
pca_data$Group<-gsub("Normoxia_30", "H0_30", pca_data$Group)
pca_data$Group<-gsub("Normoxia_33", "H0_33", pca_data$Group)
pca_data$Group<-gsub("Normoxia_36", "H0_36", pca_data$Group)
pca_data$Group<-gsub("Normoxia_39", "H0_39", pca_data$Group)

ggplot(pca_data,(aes(x=PC1,y=PC2, group= Condition)))+ 
  geom_point(aes(shape=Group, fill=Group), size=7)+ #geom_point(aes(fill= Group, shape=Treatment), size=5, color="black")+
  theme_classic()+
  labs(col="Condition")+
  scale_shape_manual(values=c(21,21, 21, 21, 24, 24, 24, 24)) +
  scale_fill_manual(values=COL2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  theme(axis.text = element_text(size=26))+
  theme(axis.title = element_text(size=26))+
  theme(legend.title = element_text(size=22))+
  theme(text = element_text(size=22))+
  theme(aspect.ratio = 1) #theme(aspect.ratio = 3/4)

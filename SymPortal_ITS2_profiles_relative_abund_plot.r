library(tidyr)
library(ggplot2)
library(reshape2)
library(dbplyr)

data<- read.delim("183_20211007_03_DBV_20211008T013033.profiles.relative.abund_and_meta.txt", sep="\t", header = FALSE) #183_20211007_03_DBV_20211008T013033.profiles.absolute.abund_and_meta.txt
data2<- data[7:15, 2:5] #subset dataframe for plotting info
colnames(data2)<- data2[1,] #make colnames the first row
colnames(data2)[1]<- "Samples" #add missing colanme
data2<- data2[-c(1),] #remove redundant first row of labels
data3<-pivot_longer(data2, cols=2:4, names_to = "Profile", values_to = "Abundance") #use tidyr to pivot wide table to long format for plotting

data3$Samples <- factor(data3$Samples, levels=c("H0_T1_30_R1", "H0_T1_30_R2", "H0_T1_30_R3", "H0_T1_30_R4", "HH_T1_30_R1", "HH_T1_30_R2", "HH_T1_30_R3", "HH_T1_30_R4")) #reorder samples
data3$Profile <- factor(data3$Profile, levels=unique(data3$Profile))
data3$Abundance <- as.numeric(data3$Abundance) #not integer as they are normalised to 1

colours<- c("#CCFF66", "#66CCFF", "#FF66FF")

 ggplot(data3, aes(x = Samples, fill = Profile, y = Abundance)) + 
  geom_bar(stat = "identity", colour = "black", position="fill") + 
  scale_y_continuous(labels = scales::percent) + #change to % scale
  theme(axis.text.x = element_text(angle = 90, size = 14, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 16), legend.title = element_text(size = 16), 
        axis.title.x = element_text(size = 16), 
        legend.text = element_text(size = 12, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 12)) + 
  labs(x = "Samples", y = "Relative Abundance (%)", fill = "ITS2 type profiles") +
  scale_fill_manual(values = colours)
 
 ggsave("FigureS5_ITS2typeprofile.pdf", width=9, height=6, dpi= 300)

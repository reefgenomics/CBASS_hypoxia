library("ggplot2")
library(grid)

######### Selecting genes of interest fpkms and average/SE for groups #########

cnt=read.table("./Input_files/fpkm_matrix_CBASS_HH.txt", header=TRUE, row.names = 1, sep= "\t")
cnts.sor=cnt[ , order(names(cnt))]
map.pre=read.delim("./Input_files/Metadata_CBASS_HH.csv", row.names = 1, header = TRUE, sep=",")
map=map.pre[order(rownames(map.pre)), ]

map.A=subset(map, map$Temperature == "30" & map$Treatment == "Normoxia")
map.B=subset(map, map$Temperature == "30" & map$Treatment == "Deoxygenation")
map.C=subset(map, map$Temperature == "33" & map$Treatment == "Normoxia")
map.D=subset(map, map$Temperature == "33" & map$Treatment == "Deoxygenation")
map.E=subset(map, map$Temperature == "36" & map$Treatment == "Normoxia")
map.F=subset(map, map$Temperature == "36" & map$Treatment == "Deoxygenation")
map.G=subset(map, map$Temperature == "39" & map$Treatment == "Normoxia")
map.H=subset(map, map$Temperature == "39" & map$Treatment == "Deoxygenation")

#Eggnog emapper annotation output for filtered transcriptome
gene=read.delim("./Input_files/MM_m5iwotms.emapper.annotations.csv", header=TRUE, row.names= 1, sep=",")
rownames(gene)=gsub("\\.p1", "", rownames(gene))
gene$geneid<-rownames(gene)

#use any of these depending on which type of annotation the GOI has
genes=subset(gene, gene$Preferred_name == "SCARB2")
#genes=subset(gene, gene == "HIF1A")
#genes=subset(gene, gene$Description == "PFAM FAD dependent oxidoreductase")

#30 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.A))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("Normoxia",length(nrow(out)))
C_30 <- out

#30 deoxy
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.B))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("30",length(nrow(out)))
out$Treatment = rep("Deoxygenation",length(nrow(out)))
T_30 <- out

#33 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.C))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("33",length(nrow(out)))
out$Treatment = rep("Normoxia",length(nrow(out)))
C_33 <- out

#33 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.D))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("33",length(nrow(out)))
out$Treatment = rep("Deoxygenation",length(nrow(out)))
T_33 <- out

#36 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.E))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("36",length(nrow(out)))
out$Treatment = rep("Normoxia",length(nrow(out)))
C_36 <- out

#36 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.F))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("36",length(nrow(out)))
out$Treatment = rep("Deoxygenation",length(nrow(out)))
T_36 <- out

#39 control
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.G))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("39",length(nrow(out)))
out$Treatment = rep("Normoxia",length(nrow(out)))
C_39 <- out

#39 treatment
cnts=cnts.sor[,(names(cnts.sor)[(names(cnts.sor) %in% row.names(map.H))])]
cnts.a=subset(cnts, row.names(cnts) %in% genes$geneid)
cnts.colsum <- apply(cnts.a, 2, sum)
Mean <- mean(cnts.colsum)
x <- sd(cnts.colsum)
SE <- x/sqrt(ncol(cnts.a)) 
out=data.frame(cbind(Mean, SE))
out$Temperature = rep("39",length(nrow(out)))
out$Treatment = rep("Deoxygenation",length(nrow(out)))
T_39 <- out

out=rbind(C_30, T_30, C_33, T_33, C_36, T_36, C_39, T_39)
out$Type = rep("SCARB2",length(nrow(out))) #change accordingly
write.table(out, "SCARB2_fpkm_all_temps.txt", sep = "\t", quote = FALSE, row.names= FALSE, col.names = TRUE) #change accordingly

######### scatterplots of fpkm expression patterns over timepoints #########

data <- read.delim("SCARB2_fpkm_all_temps.txt") #change file accordingly for different genes
data$Temperature <- as.character(data$Temperature) #as.character() identifies as character then can be made into a factor
data$Temperature <- factor(data$Temperature, levels=unique(data$Temperature))

#change scale/colour accordingly 
a<-ggplot(data, aes(x = Temperature, y = Mean, color = Treatment, group = Treatment)) + ylim(0, 60) + geom_point(aes(fill = Treatment), shape=21, colour= "black", size=7) + geom_line(size=2) + ylab(label="FPKMs") + xlab("Temperature") + theme_classic(base_size = 29) + theme(aspect.ratio = 1) + geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), width=0.2, color= "black") + scale_color_manual(values=c("#9933CC", "#CCCCFF")) + scale_fill_manual(values=c("#9933CC", "#CCCCFF")) 

# add gene label, SCARB2 is also known as CD36
grob <- grobTree(textGrob("CD36", x=0.1,  y=0.95, hjust=0,
                          gp=gpar(col="black", fontsize=22, fontface="italic")))
a + annotation_custom(grob)

#ggsave("CD36_fpkms_plot.png", width=9, height=10)

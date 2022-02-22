library(ggplot2)
library(RColorBrewer)
library(gplots)

#subset of interesting GO_enriched terms between heat+hypoxia vs heat only
GO_list=read.delim("./Input_files/GO_enriched_subset_list_HHvH02.txt", sep="\t")

#manually curated catergoes of GO enriched terms
a<-c("Iron binding", "Light response", "Protein synthesis", "Immunity", "Energy", "Structural", "ROS", "Stress signalling", "Epigenetics", "Cell death")
GO_list$type<-factor(GO_list$type, levels=a)
GO_list$temp= gsub(" HHvH0", "", GO_list$temp)
b<-c("30", "33", "36", "39")
GO_list$temp<-factor(GO_list$temp, levels=b)

ggplot(GO_list, aes(x=temp, y=Term, fill=weightFisher)) + geom_tile(color="gray90", size=0.05)+
  guides(fill=guide_legend(title="P-value"))+
  scale_fill_gradientn(colours=c("turquoise4", "turquoise1", "paleturquoise1"))+
  theme_bw()+ 
  scale_y_discrete(position = 'right')+
  theme(axis.text.x = element_blank())+
  theme(legend.title = element_text(size=10), axis.title=element_blank())+
  facet_grid(type~temp, scales="free",space="free", switch="y", labeller=label_wrap_gen(width=10, multi_line = TRUE))+
  theme(strip.text.y = element_text(angle = 0), strip.background.x = element_rect(color="gray20"))+
  theme(strip.text.y.left = element_text(angle = 0))+
  theme(strip.text.x = element_text(angle = 0))+
  theme(axis.text.y=element_text(size=10, hjust=0))+
  theme(panel.spacing.x=unit(0, "lines"))+
  theme(legend.text.align = 0, legend.text = element_text(size=8))

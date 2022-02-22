library(data.table)

######### Selecting annotated genes of associated with heat stress from DE comparisons #########

annot=read.csv("MM_m5iwotms.emapper.annotations.csv")

#calcium signalling
ca=annot[annot$Preferred_name %like% "ATP2A2|TMEM110|SLC24A2", ]
ca$type = rep("calcium_signalling",length(nrow(ca)))

#heat shock proteins
hsp=annot[annot$Preferred_name %like% "HYOU1|HSPA5|HSPA4L|HSP90B1|DNAJA4|CRYAB", ]
hsp$type = rep("heat_shock",length(nrow(hsp)))
hsp$Preferred_name<-gsub("CRYAB", "CRYAB/HSP20", hsp$Preferred_name)

#mitochondria-associated
mito=annot[annot$Preferred_name %like% "UQCRH|UQCRFS1|Ugcrb|HIGD1A", ]
mito$type = rep("mitochondria",length(nrow(mito)))
mito2=annot[annot$Description %like% "Thromboxane A synthase 1", ]
mito2$type = rep("mitochondria",length(nrow(mito2)))

#cytoskeletal
cyto=annot[annot$Preferred_name %like% "EPB41L2|CCDC155|MSN", ]
cyto$type = rep("cytoskeleton",length(nrow(cyto)))

#necrosis
nec=annot[annot$Preferred_name %like% "TRAF4|TNFAIP3|TNFAIP8", ]
nec$type = rep("necrosis",length(nrow(nec)))

#glycolysis
gly=annot[annot$Preferred_name %like% "SDHB|PDK2|GAPDH", ]
gly$type = rep("glycolysis",length(nrow(gly)))
gly2=annot[annot$Description %like% "fructose-bisphosphate aldolase", ]
gly2$type = rep("glycolysis",length(nrow(gly2)))
gly2$Preferred_name<-gsub("-", "ALDO", gly2$Preferred_name)

#ROS
ros=annot[annot$Preferred_name %like% "HMOX1|CCS", ]
ros$type = rep("ROS",length(nrow(ros)))
ros3=annot[annot$Description %like% "serves to protect cells from the toxic effects of hydrogen peroxide", ]
ros3$type = rep("ROS",length(nrow(ros3)))

out=rbind(ca, hsp, mito, mito2, cyto, nec, gly, gly2, ros, ros3)

#add geneidentifer to distinguish between genes with same annotation
out$geneidentifer= paste(out$Preferred_name, out$query, sep = "_")
#to make with DE geneIDs but geneidentifier keeps transcript level
out$query=gsub("\\.p1", "", out$query)
out$query=gsub("\\.p2", "", out$query)

DE5=read.table("DE_output_33C_30C_comp.txt", sep="\t", header=T)
DE5$comp = rep("33v30 H0",length(nrow(DE5)))
DE5$geneid = rownames(DE5)
DE5$descrip<- out$Description[match(DE5$geneid, out$query)] #adds to those DE with annotation, to filter out those with none
DE5=subset(DE5, !(DE5$descrip == ""))
DE5$type<- out$type[match(DE5$geneid, out$query)]
DE5$gene<- out$Preferred_name[match(DE5$geneid, out$query)]
DE5$ko<- out$KEGG_ko[match(DE5$geneid, out$query)]
DE5$geneidentifer<- out$geneidentifer[match(DE5$geneid, out$query)]

DE6=read.table("DE_output_36C_30C_comp.txt", sep="\t", header=T)
DE6$comp = rep("36v30 H0",length(nrow(DE6)))
DE6$geneid = rownames(DE6)
DE6$descrip<- out$Description[match(DE6$geneid, out$query)]
DE6=subset(DE6, !(DE6$descrip == ""))
DE6$type<- out$type[match(DE6$geneid, out$query)]
DE6$gene<- out$Preferred_name[match(DE6$geneid, out$query)]
DE6$ko<- out$KEGG_ko[match(DE6$geneid, out$query)]
DE6$geneidentifer<- out$geneidentifer[match(DE6$geneid, out$query)]

DE7=read.table("DE_output_39C_30C_comp.txt", sep="\t", header=T)
DE7$comp = rep("39v30 H0",length(nrow(DE7)))
DE7$geneid = rownames(DE7)
DE7$descrip<- out$Description[match(DE7$geneid, out$query)]
DE7=subset(DE7, !(DE7$descrip == ""))
DE7$type<- out$type[match(DE7$geneid, out$query)]
DE7$gene<- out$Preferred_name[match(DE7$geneid, out$query)]
DE7$ko<- out$KEGG_ko[match(DE7$geneid, out$query)]
DE7$geneidentifer<- out$geneidentifer[match(DE7$geneid, out$query)]

DE8=read.table("DE_output_33T_30T_comp.txt", sep="\t", header=T)
DE8$comp = rep("33v30 HH",length(nrow(DE8)))
DE8$geneid = rownames(DE8)
DE8$descrip<- out$Description[match(DE8$geneid, out$query)] 
DE8=subset(DE8, !(DE8$descrip == ""))
DE8$type<- out$type[match(DE8$geneid, out$query)]
DE8$gene<- out$Preferred_name[match(DE8$geneid, out$query)]
DE8$ko<- out$KEGG_ko[match(DE8$geneid, out$query)]
DE8$geneidentifer<- out$geneidentifer[match(DE8$geneid, out$query)]

DE9=read.table("DE_output_36T_30T_comp.txt", sep="\t", header=T)
DE9$comp = rep("36v30 HH",length(nrow(DE9)))
DE9$geneid = rownames(DE9)
DE9$descrip<- out$Description[match(DE9$geneid, out$query)] 
DE9=subset(DE9, !(DE9$descrip == ""))
DE9$type<- out$type[match(DE9$geneid, out$query)]
DE9$gene<- out$Preferred_name[match(DE9$geneid, out$query)]
DE9$ko<- out$KEGG_ko[match(DE9$geneid, out$query)]
DE9$geneidentifer<- out$geneidentifer[match(DE9$geneid, out$query)]

DE10=read.table("DE_output_39T_30T_comp.txt", sep="\t", header=T)
DE10$comp = rep("39v30 HH",length(nrow(DE10)))
DE10$geneid = rownames(DE10)
DE10$descrip<- out$Description[match(DE10$geneid, out$query)] 
DE10=subset(DE10, !(DE10$descrip == ""))
DE10$type<- out$type[match(DE10$geneid, out$query)]
DE10$gene<- out$Preferred_name[match(DE10$geneid, out$query)]
DE10$ko<- out$KEGG_ko[match(DE10$geneid, out$query)]
DE10$geneidentifer<- out$geneidentifer[match(DE10$geneid, out$query)]

output=rbind(DE5, DE6, DE7, DE8, DE9, DE10)
#write.table(output, "DE_all_comp_heatstress_heatgrid_input_SELECTED_genes_NEW3.txt", row.names = F, sep="\t")

######### heatgrid #########

a<-c("calcium_signalling", "heat_shock", "mitochondria", "cytoskeleton", "necrosis", "glycolysis", "ROS")
output$type<-factor(output$type, levels=a)
b<-c("33v30 H0", "33v30 HH", "36v30 H0", "36v30 HH", "39v30 H0", "39v30 HH")
output$comp<-factor(output$comp, levels=b)

ggplot(output, aes(x=comp, y=geneidentifer, fill=log2FoldChange)) + geom_tile(color="black", size=0.05)+
  scale_fill_gradientn(colours=c("lightyellow1", "grey90", "grey60", "turquoise4"))+
  guides(fill=guide_legend(title="log2fc", reverse = T))+
  theme_bw()+
  scale_y_discrete(position = 'right')+
  theme(axis.text.x = element_blank())+
  theme(legend.title = element_text(size=10), axis.title=element_blank())+
  facet_grid(type~comp, scales="free",space="free", switch="y", labeller=label_wrap_gen(width=10, multi_line = TRUE))+
  theme(strip.text.y = element_text(angle = 0), strip.background.x = element_rect(color="black"))+
  theme(strip.text.y.left = element_text(angle = 0, size = 8))+
  theme(strip.text.x = element_text(angle = 0))+
  theme(axis.text.y=element_text(size=8, hjust=0))+
  theme(panel.spacing.x=unit(0, "lines"))+
  theme(legend.text.align = 0, legend.text = element_text(size=8))

#ggsave("heatgrid_heat_hypoxia_SELECTED_NEW3.png", width=9, height=6)

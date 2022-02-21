library(topGO)
library(data.table)

######### GO enrichment #########

#GO annotations for all transcripts
annot=read.csv("./Input_files/MM_m5iwotms.emapper.annotations.csv") #eggnog mapper results that include GO annotations
annot$query=gsub("\\.p1", "", annot$query)
annot$query=gsub("\\.p2", "", annot$query)
Allgenes=annot[c("query", "GOs")]
Allgenes=subset(Allgenes, !(Allgenes$GOs == "-"))
write.table(Allgenes, "Allgenes_GO_annot_new.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

geneID2GO = readMappings(file = "./Input_files/Allgenes_GO_annot_new.txt") #gene names and GO list in 2 columns of all samples
#all possible gene names with GO annotation
geneUniverse = names(geneID2GO)
message("Number of genes with GO annotation: ", length(geneUniverse))

DE_output=read.delim("./Input_files/DE_output_30T_30C_comp.txt", header = TRUE, sep="\t") #GOI, change accordingly
DE_genes=DE_output$gene_id

#remove any GOI without GO annotation
keep=DE_genes %in% geneUniverse
keep=which(keep==TRUE)
DE_genes=DE_genes[keep]
#make named vector list of factors showing which are GOI
geneList=factor(as.integer(geneUniverse %in% DE_genes))
names(geneList)=geneUniverse
message("Number of DE genes with GO annotation: ", length(intersect(geneUniverse,DE_genes)))

#BP
#create topgo data object
myGOdata= new("topGOdata", description="My project", ontology="BP", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
#test for significance
#run weighted algorithm as classic doesnt take into consideration GO hierarchy so could overrepresent enrichment
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
#generate a table of results using Gentable function
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultsFisher", ranksOf = "weightFisher", topNodes = length(allGO))
#change format to non-scientific digits
options(scipen = 999)
#correct for multiple testing e.g. a p-adjusted value
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
BP=allRes
BP$ontology="BP"

#MF
myGOdata= new("topGOdata", description="My project", ontology="MF", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
options(scipen = 999)
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
MF=allRes
MF$ontology="MF"

#CC
myGOdata= new("topGOdata", description="My project", ontology="CC", allGenes=geneList,  annot = annFUN.gene2GO, gene2GO=geneID2GO)
resultFisher = runTest(myGOdata, algorithm="weight01", statistic="fisher")
allGO = usedGO(object = myGOdata)
allRes = GenTable(myGOdata, weightFisher = resultFisher, orderBy = "resultFisher", ranksOf = "weightFisher", topNodes = length(allGO))
options(scipen = 999)
allRes$adjusted.p= p.adjust(allRes$weightFisher, method = "bonferroni", n = length(allRes$weightFisher))
allRes$q.value= p.adjust(allRes$weightFisher, method = "fdr", n = length(allRes$weightFisher))
CC=allRes
CC$ontology="CC"

out=rbind(BP, CC, MF)
out.s=subset(out, out$weightFisher <0.001)
write.table(out.s, "GO_results_De_30T_30C_fisher.txt", row.names = FALSE, sep = "\t", quote = FALSE) #change accordingly

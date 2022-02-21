
######### filtering for only contigs with sequence length =/> 500 #########

library(dplyr)
library(seqinr)

blast_pre=read.delim("./Input_files/out_blast_tax_CONTROLS.txt", header = FALSE, sep = "\t")
#needed to split up sequences matching to multiple tax ids with ; in tax id file to run ete3 ncbiquery.
#so formed duplicate sequences which need removed from blast output

blastuniq=blast_pre[!duplicated(blast_pre$V1),] #remove duplicate sequences, leaves the first occurring one
write.csv(blastuniq, "out_blast_tax_CONTROLS_UNIQUE.csv", row.names = FALSE, quote = FALSE)
blastuniq2<-read.csv("./Input_files/out_blast_tax_CONTROLS_UNIQUE.csv", row.names = 1) #now seq names = rownames so only unqiue seq
blast=subset(blastuniq2, !(grepl("scaffold", rownames(blastuniq2)))) #removes scaffolds so only unique contigs

#subset blast by seq length
blast500=subset(blast, blast$V2 > 499)
write.table(rownames(blast500), "blast_500_contigs_uniq.txt", sep= "\t", quote=FALSE, col.names = FALSE, row.names = FALSE)

#use following bash command lines to get list of blast500 contig headers with length values
#rg -f blast_500_contigs_uniq.txt CBASS_assembly_CONTROLS.scafSeq.fasta > all_500_whole_fa_headers_CONTROLS.txt
#cat all_500_whole_fa_headers_CONTROLS.txt | tr -d '>,' > all_500_CONTROLS_fa_headers.txt

myfasta <- read.fasta (file = "./Input_files/CBASS_assembly_CONTROLS.scafSeq.fasta", seqtype = "DNA", as.string = TRUE, set.attributes = FALSE, forceDNAtolower = FALSE, whole.header=TRUE)
blast500names<- read.table("./Input_files/all_500_CONTROLS_fa_headers.txt", header=F, sep="\t")
my_fasta_sub <- myfasta [names(myfasta) %in% blast500names$V1]

#subset assembly fasta for subsequent mapping
write.fasta (sequences = my_fasta_sub, names = names (my_fasta_sub), file.out = "CBASS_assembly_CONTROLS_500.fasta")

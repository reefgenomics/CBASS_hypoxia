######### checking number of contigs per coral holobiont taxa group #########

#cnidarian contigs only
cnidaria=read.delim("./Input_files/checking_cnidaria_taxa.txt", header = FALSE) 
blastcnidaria=subset(blast, blast$V3 %in% cnidaria$V1) 
blast500cni=subset(blast500, rownames(blast500) %in% rownames(blastcnidaria)) 

#dinophyceae contigs only
dino=read.delim("./Input_files/checking_dino_taxa.txt", header=FALSE)
blastdino=subset(blast, blast$V3 %in% dino$V1)
blast500dino=subset(blast500, blast500$V3 %in% dino$V1)

#bacteria contigs only
bacteria=read.delim("./Input_files/checking_bacteria_taxa.txt", header=FALSE)
blastbacteria=subset(blast, blast$V3 %in% bacteria$V1)
blast500bacteria=subset(blast500, blast500$V3 %in% bacteria$V1)

#fungi contigs only
fungi=read.delim("./Input_files/checking_fungi_taxa.txt", header=FALSE)
blastfungi=subset(blast, blast$V3 %in% fungi$V1)
blast500fungi=subset(blast500, blast500$V3 %in% fungi$V1)

#virus contigs only
virus=read.delim("./Input_files/checking_virus_taxa.txt", header=FALSE)
blastvirus=subset(blast, blast$V3 %in% virus$V1)
blast500virus=subset(blast500, blast500$V3 %in% virus$V1)

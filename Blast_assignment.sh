#!/bin/bash

######### Blastn query sequences to assign taxa IDs #########

wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xf taxdb.tar.gz
blastn -query CBASS_assembly_CONTROLS.scafSeq.fasta -db /share/databases/nt/nt -evalue 1e-3 -num_threads 64 -max_target_seqs 1 -out out_blast_tax_CONTROLS -outfmt "6 qseqid qlen staxids sscinames evalue bitscore"

######### Assign taxa IDs to the main taxa groups of coral holobiont #########

cut -f2 out_blast_tax_CONTROLS | sort | uniq > taxIDs_all #get a list of unique taxIDs 
cat taxIDs_all | tr ";" "\n" > taxIDs_all2 #to split multiple taxIDs into several lines
cat taxIDs_all2 | while read line ; do ete3 ncbiquery --search $line --info >> tax_info_all_final ; done

grep "Cnidaria" tax_info_all_final > checking_cnidaria_taxa.txt
grep "Dinophyceae"  tax_info_all_final > checking_dino_taxa.txt
grep "Bacteria"  tax_info_all_final > checking_bacteria_taxa.txt
grep "Virus"  tax_info_all_final > checking_virus_taxa.txt
grep "Ascomycota\|Basidiomycota"  tax_info_all_final > checking_fungi_taxa.txt

#use these lists for checking number of contigs across coral holobiont taxa groups in Rscript Contigs_per_taxa_group

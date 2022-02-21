#!/bin/bash

######### Trimming adaptors ##########

for i in ./*_R1.fastq.gz ; do java -jar /home/cardena/software/Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads 16 $i ${i//_R1./_R2.} -baseout ${i//_R1.fastq.gz/.fastq} ILLUMINACLIP:/home/alderdicer/adapters/NovaseqPE150.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50; done

######### Assemble de novo transcriptome for non-bleached samples only ##########

cat *1P.fastq >ALL_1P.fastq
cat *2P.fastq >ALL_2P.fastq

cat > config_file
#maximal read length
max_rd_len=150
[LIB]
#maximal read length in this lib
rd_len_cutof=150
#average insert size
avg_ins=300
#if sequence needs to be reversed
reverse_seq=0
#in which part(s) the reads are used
asm_flags=3
#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)
map_len=32
q1=/home/alderdicer/CBASS_HH/X201SC21061862-Z01-F001/soap/assembly_CONTROLS/ALL_CONTROLS_1P.fastq
q2=/home/alderdicer/CBASS_HH/X201SC21061862-Z01-F001/soap/assembly_CONTROLS/ALL_CONTROLS_2P.fastq

/home/alderdicer/miniconda3/envs/soap/bin/SOAPdenovo-Trans-31mer all -s config_file -o CBASS_assembly_CONTROLS

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

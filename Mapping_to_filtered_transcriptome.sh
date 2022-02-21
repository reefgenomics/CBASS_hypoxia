#!/bin/bash

######### Map reads to transcriptome of only non-bleached coral samples >500 bp #########

bowtie2-build CBASS_assembly_CONTROLS_500.fasta  CBASS_assembly_CONTROLS_500_DB

for i in /home/alderdicer/CBASS_HH/X201SC21061862-Z01-F001/trimmed/*_1P.fastq.gz ; do bowtie2 -p 20 -x CBASS_assembly_CONTROLS_500_DB -1 $i -2 ${i//1P.fastq.gz/2P.fastq.gz} -S ${i//1P.fastq.gz/500_all.sam} ; done

######### Estimate read counts #########

for i in *500_all.sam; do /usr/bin/salmon quant -t CBASS_assembly_CONTROLS_500.fasta -l A -a $i -o quants_500_all/$i ; done

#edit folder names from salmon output so can use as input for tximport in R
for i in *.sam ; do mv -v "$i" "${i/_500_all.sam//}"; done

# CBASS_hypoxia

Data and code for RNASeq analysis and determining ED50 of FvFm over temperature.

Project: Short-term thermal stress assays using CBASS under normal and reduced oxygen conditions.

### Workflow

1.	De novo transcriptome was assembled using Denovo_transcriptome_assembly.sh
2.	Then filtered to only sequences =>500 bp using Filtering_transcriptome.r
3.	Reads were then mapped to transcriptome using Mapping_to_filtered_transcriptome.sh
4.	Principal components analysis including all samples plotted using PCA_plot.r
5.	Transcripts were annotated via BLAST query & assigned to taxa groups using Blast_assignment.sh
6.	Number of transcripts per taxa group of coral holobiont estimated using Contigs_per_taxa_group.r
7.	Differential expression of transcripts performed using DESeq2_results.r
8.	GO enrichment analysised performed using GO_enrichment_topGO.r
9.	Heatgrid of GO enriched terms using Heatgrid_GO_enrichment.r
10.	Heatgrid of DE transcripts using Heatgrid_DETS_heat_stress.r
11.	Scatterplot of FPKM across heating temperatures using FPKM_scatter_plots.r
12. Dose-response curves to determine ED50 of FvFm over temperature using DRC_PAM_ED50.r

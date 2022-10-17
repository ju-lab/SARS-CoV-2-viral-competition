# SARS-CoV-2-viral-competition
Scripts used for the paper titled "Relative infectivity of the SARS-CoV-2 Omicron variant in human alveolar cells"

### Consensus clonal mutations of viral variants 
Total 118 clonal mutations from four viral variants were written in unique_mutations.tsv

### VCF files of GISAID
Each viral variant's mutations were acquired through the GISAID fasta file aligned at NC_045512v2 
File name is (viral_variant).vcf

### Fraction of each viral variant by two different methods (RNA sequencing and Plaque assay)
Composition of each viral variant was calculated by two methods. Information is stored in live_virus_assay.xlsx
RNA sequencing counts both infectious and defective viruses. But Plaque assay only counts infectious viruses. 

### Expression of human type-2 alveolar cell genes 
Expression level of human type-2 alveolar organoid is stored in 20221004_fig_dat_per10k.tsv

### GTF file used in this study 
GTF file used for gene annotation of SARS-CoV-2 was uploaded as ncbiGenes_non_nestedexon.gtf 

### Depth coverage of SARS-CoV-2 transcriptome by two different methods (SMART-seq3 and 10X) 
Average value of normalized depth of SARS-CoV-2 infected alveolar cells was stored in samtools_depth_final.csv 
10X data from the paper "Three-Dimensional Human Alveolar Stem Cell Culture Models Reveal Infection Response to SARS-CoV-2" was used. (MOI 1, 60 hours post infection cells were used).

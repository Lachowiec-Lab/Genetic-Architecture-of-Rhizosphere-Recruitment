# Code and datasets supporting "Conditional genetic architecture of bacterial and fungal rhizosphere recruitment for crop performance"


## Preparing phyloseq objects
Microbial abundance based on 16S and ITS amplicon sequencing was processed to create phyloseq objects in the phyloseq folder.

### 16S
For 16S data, in the ~phyloseq/16S folder
1. 2021.R and 2022.R clean the raw reads with dada2 to create seqtabAll_F files for each location-year.

2. merging.R combines the seqtabAll_F files into a single seqtab called seqtab.merged. Chimeras are removed from this file, and then taxonomic assignment is completed using the silva nr99 v138.1 with species database to create the file taxTab. 

3. Processing.Rmd combines seqtab.merged.NoC.RDS and taxTab.RDS into a phyloseq object. It also uses the FieldData to tie into all the plot information.
	a. ps.RDS is complete at this point and saved. 
	
	b. Three-step filtering to create ps.filt.RDS:
		1. Remove reads that are not classified to Phylum.
		2. Prevalence threshold of 2 is used to remove single-instance ASVs
		3. Remove samples with fewer than 10,000 reads.

### ITS
For ITS data, in the ~phyloseq/ITS folder

## Preparing data for GWAS
Scripts for processing phyloseq objects for GWAS analyses are found in the ~/prepForGWAS/ folder

## Analyses
Scripts for examining microbial abundances, heritability calculations, GWAS, Wilcoxon tests are found in the ~/analyses/ folder

## Citations
This work supports the manuscript currently under review at Theoretical and Applied Genetics.

Raw 16S and ITS reads as fastq files are available in the NCBI Sequence Read Archive database under BioProject PRJNA909304, https://www.ncbi.nlm.nih.gov/bioproject/PRJNA909304.


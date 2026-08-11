library(BiocManager)
library(knitr)
library(gridExtra)
library(dada2)
library(ggplot2)
library(ggpubr)
library(lme4)
library(plyr)
library(dplyr)
library(vegan)





doitagain = function(miseq_path = "", 
                        filt_path = "",ErrorRate) {
  
  fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))

  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)
  
  # Specify the full path to the fnFs and fnRs
  fnFs <- file.path(miseq_path, fnFs)

  filt_path <- file.path(filt_path) # Place filtered files in filtered/ subdirectory
  if(!file_test("-d", filt_path)) dir.create(filt_path)
  filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))

  out <- filterAndTrim(fnFs,filtFs, maxN = 0, maxEE = ErrorRate, truncQ = 20, truncLen = 200, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
  
  saveRDS(out, file.path(filt_path,"out.RDS"))
  
  derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
  saveRDS(derepFs, file.path(filt_path,"derepFs.RDS"))
  
  # Name the derep-class objects by the sample names
  names(derepFs) <- sampleNames
  
  errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)
  
  saveRDS(errF, file.path(filt_path,"errF.RDS"))
  
  dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
  saveRDS(dadaFs, file.path(filt_path,"dadaFs.RDS"))
  
  # Create seqtable on forward reads
  seqtabAll <- makeSequenceTable(dadaFs)
  saveRDS(seqtabAll,file.path(filt_path,"seqtabAll_F.RDS"))
  
  rm(list = ls())
  
}

# 2021
doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2021/BZ/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/BZ21", ErrorRate = 6)

doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2021/CARC/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/MO21", ErrorRate = 6)

doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2021/SD/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/SD21", ErrorRate = 6)

doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2021/HI/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/HI21", ErrorRate = 6)

# 2022
doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2022/BZ/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/BZ22", ErrorRate = 6)

doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2022/CARC/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/MO22", ErrorRate = 6)

doitagain(miseq_path = "/home/j22f487/Desktop/ITS/Reads/2022/SD/cutadapt",
          filt_path = "/home/j22f487/Desktop/Phyloseq2.0/ITS/Quality_20_Len200/SD22", ErrorRate = 6)



# MERGE dada2 seqtabs

bz21 <- readRDS("BZ21/seqtabAll_F.RDS")
bz22 <- readRDS("BZ22/seqtabAll_F.RDS")
mo21 <- readRDS("MO21/seqtabAll_F.RDS")
mo22 <- readRDS("MO22/seqtabAll_F.RDS")

hi21 <- readRDS("HI21/seqtabAll_F.RDS")

sd21 <- readRDS("SD21/seqtabAll_F.RDS")
row.names(sd21) = substr(row.names(sd21),6,length(row.names(sd21))) # Fixed naming
row.names(sd21) = gsub("-RX","",row.names(sd21))

sd22 <- readRDS("SD22/seqtabAll_F.RDS")
row.names(sd22) = substr(row.names(sd22),6,length(row.names(sd22))) # Fixed naming
row.names(sd22) = gsub("-RX","",row.names(sd22))

# Add locyr tag to row names
row.names(bz21) = paste(row.names(bz21),"BZ21",sep="_")
row.names(bz22) = paste(row.names(bz22),"BZ22",sep="_")
row.names(mo21) = paste(row.names(mo21),"MO21",sep="_")
row.names(mo22) = paste(row.names(mo22),"MO22",sep="_")
row.names(hi21) = paste(row.names(hi21),"hi21",sep="_")
row.names(sd21) = paste(row.names(sd21),"SD21",sep="_")
row.names(sd22) = paste(row.names(sd22),"SD22",sep="_")


seqtab.merged = mergeSequenceTables(bz21,bz22,mo21,mo22,sd21,sd22,hi21)
saveRDS(seqtab.merged, "seqtab.merged.RDS")

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.merged, method="consensus", multithread=FALSE)
saveRDS(seqtab, "seqtab.merged.NoC.RDS")

# Assign taxonomy
library(dada2)
seqtab = readRDS("seqtab.merged.NoC.RDS")
tax <- assignTaxonomy(seqtab, "sh_general_release_dynamic_all_04.04.2024.fasta", multithread=TRUE)
saveRDS(tax, "taxTab.RDS")


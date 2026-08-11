library(BiocManager)
library(knitr)
library(gridExtra)
library(dada2)
library(phyloseq)
library(ggplot2)
library(ggpubr)
library(lme4)
library(plyr)
library(dplyr)
library(vegan)


#########
# BOZEMAN
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2022/Bozeman/KillianBozemanV6V8"
dir.exists(miseq_path)

fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)

plotQualityProfile(fnFs[1:2])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))

# Set truncLen for F and R based on plotQualityProfile
out.BZ <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     truncLen=200, maxEE=2, truncQ=2, maxN=0, trimLeft=5, rm.phix=TRUE,
                     compress=TRUE, verbose=FALSE, multithread=FALSE)

saveRDS(out.BZ, "BZ2022/out_BZ.RDS")
saveRDS(fnFs, "BZ2022/fnFs.RDS")
saveRDS(filtFs, "BZ2022/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "BZ2022/derepFs_BZ.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "BZ2022/errF_BZ.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "BZ2022/dadaFs_BZ.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"BZ2022/seqtabAll_F.RDS")

rm(list = ls())

#########
# CARC
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2022/CARC/KillianCARCV6V8"
dir.exists(miseq_path)

fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)

plotQualityProfile(fnFs[1:2])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))

# Set truncLen for F and R based on plotQualityProfile
out.MO <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     truncLen=200, maxEE=2, truncQ=2, maxN=0, trimLeft=5, rm.phix=TRUE,
                     compress=TRUE, verbose=FALSE, multithread=FALSE)

saveRDS(out.MO, "MO2022/out_MO.RDS")
saveRDS(fnFs, "MO2022/fnFs.RDS")
saveRDS(filtFs, "MO2022/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "MO2022/derepFs_MO.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "MO2022/errF_MO.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "MO2022/dadaFs_MO.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"MO2022/seqtabAll_F.RDS")

rm(list = ls())

#########
# South Dakota
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2022/SD/Ewing2022ALL"
dir.exists(miseq_path)

fnFs <- sort(list.files(miseq_path, pattern="_R1_001.fastq"))


# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sampleNames <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# Specify the full path to the fnFs and fnRs
fnFs <- file.path(miseq_path, fnFs)

plotQualityProfile(fnFs[1:2])

filt_path <- file.path(miseq_path, "filtered") # Place filtered files in filtered/ subdirectory
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sampleNames, "_F_filt.fastq.gz"))

# Set truncLen for F and R based on plotQualityProfile
out.SD <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     truncLen=200, maxEE=2, truncQ=2, maxN=0, trimLeft=5, rm.phix=TRUE,
                     compress=TRUE, verbose=FALSE, multithread=FALSE)

saveRDS(out.SD, "SD2022/out_SD.RDS")
saveRDS(fnFs, "SD2022/fnFs.RDS")
saveRDS(filtFs, "SD2022/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "SD2022/derepFs_SD.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "SD2022/errF_SD.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "SD2022/dadaFs_SD.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"SD2022/seqtabAll_F.RDS")

rm(list = ls())

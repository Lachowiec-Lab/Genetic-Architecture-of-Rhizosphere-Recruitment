library(BiocManager)
library(knitr)
library(gridExtra)
library(dada2)
library(ggplot2)
library(plyr)
library(dplyr)
library(vegan)


#########
# BOZEMAN
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2021/BZ2021/BZ_16S"
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

saveRDS(out.BZ, "BZ2021/out_BZ.RDS")
saveRDS(fnFs, "BZ2021/fnFs.RDS")
saveRDS(filtFs, "BZ2021/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "BZ2021/derepFs_BZ.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "BZ2021/errF_BZ.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "BZ2021/dadaFs_BZ.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"BZ2021/seqtabAll_F.RDS")

rm(list = ls())

#########
# CARC
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2021/CARC2021/CARC_16S"
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

saveRDS(out.MO, "MO2021/out_MO.RDS")
saveRDS(fnFs, "MO2021/fnFs.RDS")
saveRDS(filtFs, "MO2021/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "MO2021/derepFs_MO.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "MO2021/errF_MO.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "MO2021/dadaFs_MO.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"MO2021/seqtabAll_F.RDS")

rm(list = ls())

#########
# South Dakota
#########
setwd("~/Desktop/Phyloseq2.0/16S")
miseq_path <- "/home/j22f487/Desktop/2021/SD2021/16S/2021"
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

saveRDS(out.SD, "SD2021/out_SD.RDS")
saveRDS(fnFs, "SD2021/fnFs.RDS")
saveRDS(filtFs, "SD2021/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE) #n for low RAM PC's - try it out
saveRDS(derepFs, "SD2021/derepFs_SD.RDS")

# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames

errF <- learnErrors(filtFs, multithread=FALSE, verbose = TRUE)

saveRDS(errF, "SD2021/errF_SD.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=FALSE, verbose = TRUE)
saveRDS(dadaFs, "SD2021/dadaFs_SD.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"SD2021/seqtabAll_F.RDS")

rm(list = ls())

#########
# Hawaii
#########
setwd("~/Desktop/Phyloseq2.0/16S/all_maxEE_2")
miseq_path <- "/home/j22f487/Desktop/2021/HI2021_V2025/HI_16S"
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
out.HI <- filterAndTrim(fwd=file.path(fnFs), filt=file.path(filtFs),
                     truncLen=200, maxEE=2, truncQ=2, maxN=0, trimLeft=5, rm.phix=TRUE,
                     compress=TRUE, verbose=TRUE, multithread=FALSE)

saveRDS(out.HI, "HI2021/out_HI.RDS")
saveRDS(fnFs, "HI2021/fnFs.RDS")
saveRDS(filtFs, "HI2021/filtFs.RDS")

derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sampleNames
saveRDS(derepFs, "HI2021/derepFs_HI.RDS")


errF <- learnErrors(filtFs, multithread=15, verbose = TRUE, MAX_CONSIST = 15)
plotErrors(errF, nominalQ=F)

saveRDS(errF, "HI2021/errF_HI.RDS")

dadaFs <- dada(derepFs, err=errF, multithread=15, verbose = TRUE)
saveRDS(dadaFs, "HI2021/dadaFs_HI.RDS")

# Create seqtable on forward reads
seqtabAll <- makeSequenceTable(dadaFs)
saveRDS(seqtabAll,"HI2021/seqtabAll_F.RDS")

rm(list = ls())

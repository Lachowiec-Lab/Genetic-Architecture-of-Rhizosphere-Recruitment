library(dada2)
library(phyloseq)


bz21 <- readRDS("BZ2021/seqtabAll_F.RDS")
bz22 <- readRDS("BZ2022/seqtabAll_F.RDS")
mo21 <- readRDS("MO2021/seqtabAll_F.RDS")
mo22 <- readRDS("MO2022/seqtabAll_F.RDS")

hi21 <- readRDS("HI2021/seqtabAll_F.RDS")

sd21 <- readRDS("SD2021/seqtabAll_F.RDS")
row.names(sd21) = substr(row.names(sd21),6,length(row.names(sd21))) # Fixed naming
row.names(sd21) = gsub("-RX","",row.names(sd21))

sd22 <- readRDS("SD2022/seqtabAll_F.RDS")
row.names(sd22) = substr(row.names(sd22),6,length(row.names(sd22))) # Fixed naming
row.names(sd22) = gsub("-RX","",row.names(sd22))

# Add locyr tag to row names
row.names(bz21) = paste(row.names(bz21),"BZ21",sep="_")
row.names(bz22) = paste(row.names(bz22),"BZ22",sep="_")
row.names(mo21) = paste(row.names(mo21),"MO21",sep="_")
row.names(mo22) = paste(row.names(mo22),"MO22",sep="_")
row.names(hi21) = paste(row.names(hi21),"HI21",sep="_")
row.names(sd21) = paste(row.names(sd21),"SD21",sep="_")
row.names(sd22) = paste(row.names(sd22),"SD22",sep="_")


seqtab.merged = mergeSequenceTables(bz21,bz22,mo21,mo22,sd21,sd22,hi21)
saveRDS(seqtab.merged, "seqtab.merged.RDS")

# Remove chimeras
seqtab <- removeBimeraDenovo(seqtab.merged, method="consensus", multithread=FALSE)
saveRDS(seqtab, "seqtab.merged.NoC.RDS")
# Assign taxonomy
tax <- assignTaxonomy(seqtab, "silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=FALSE)
saveRDS(tax, "taxTab.RDS")

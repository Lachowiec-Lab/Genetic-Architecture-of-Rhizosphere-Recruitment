
### The following script contains code for generating lists of prevalent bacterial and fungal taxa for each location-year 16S or ITS dataset from the phyloseq objects,
# to be used in microbial phenotypic data preparation for GWAS.

# Here the swMT location is referred to as BZ and the cntrMT location is referred to as MO.
# Code for finding prevalent taxa is written out for each taxonomic rank/level, location-year (LocYr), and data type (16S or ITS) separately.


library(phyloseq)
library(tidyverse)



##### 16S, Bacteria #####



# Reading in the PS object for 16S from all LocYrs, this also includes all check and bulk samples

PS_16S_All = readRDS("ps.filt.RDS")


# Removing bulk and check samples

AllSamples <- sample_data(PS_16S_All) %>% data.frame

PS_16S_Exp = subset_samples(PS_16S_All, Line != "Bulk" & Line != "Hockett" & Line != "Lavina" & Line != "Merit57" & Line != "Odyssey")

ExpSamples <- sample_data(PS_16S_Exp) %>% data.frame




# BZ21

# Creating separate ps object for BZ21

BZ21 = subset_samples(PS_16S_Exp, Year == "2021" & Location == "BZ")

BZ21 = prune_taxa(taxa_sums(BZ21)>0, BZ21)      # so that OTU and tax tables only include samples from the given LocYr

# ntaxa(PSobj)    # will give the number of taxa in the OTU to check that they were indeed pruned

BZ21samples <- sample_data(BZ21) %>% data.frame


# Prevalent Phyla in BZ21

BZ21phy <- tax_glom(BZ21, "Phylum")   # agglomerate by Phylum

BZ21prevPhyPS <- filter_taxa(BZ21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated BZ21_16S ps object to phyla with a count of at least 1 in at least 80% of samples

BZ21prevPhyTax <- tax_table(BZ21prevPhyPS) %>% data.frame   # creating a data frame of that tax table


# Creating a data frame with LocYr and taxonomic level columns that additional sets of prevalent taxa can be added to

PrevBac <- BZ21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1)



# Now BZ21 by Class

BZ21cla <- tax_glom(BZ21, "Class")

BZ21prevClaPS <- filter_taxa(BZ21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ21prevClaTax <- tax_table(BZ21prevClaPS) %>% data.frame

BZ21prevClaTax <- BZ21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ21prevClaTax)



# BZ21 by Order

BZ21ord <- tax_glom(BZ21, "Order")

BZ21prevOrdPS <- filter_taxa(BZ21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ21prevOrdTax <- tax_table(BZ21prevOrdPS) %>% data.frame

BZ21prevOrdTax <- BZ21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ21prevOrdTax)



# BZ21 by Family

BZ21fam <- tax_glom(BZ21, "Family")

BZ21prevFamPS <- filter_taxa(BZ21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ21prevFamTax <- tax_table(BZ21prevFamPS) %>% data.frame

BZ21prevFamTax <- BZ21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ21prevFamTax)



# BZ21 by Genus

BZ21gen <- tax_glom(BZ21, "Genus")

BZ21prevGenPS <- filter_taxa(BZ21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ21prevGenTax <- tax_table(BZ21prevGenPS) %>% data.frame

BZ21prevGenTax <- BZ21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ21prevGenTax)




# BZ21 by ASV

# (no agglomeration for ASV level)

BZ21prevAsvPS <- filter_taxa(BZ21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ21prevAsvTax <- tax_table(BZ21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "BZ21") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ21prevAsvTax)









# BZ22

# Creating separate ps object for BZ22

BZ22 = subset_samples(PS_16S_Exp, Year == "2022" & Location == "BZ")

BZ22 = prune_taxa(taxa_sums(BZ22)>0, BZ22)   

BZ22samples <- sample_data(BZ22) %>% data.frame


# Prevalent Phyla in BZ22

BZ22phy <- tax_glom(BZ22, "Phylum")   # agglomerate by Phylum

BZ22prevPhyPS <- filter_taxa(BZ22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated BZ22_16S ps object to phyla with a count of at least 1 in at least 80% of samples

BZ22prevPhyTax <- tax_table(BZ22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

BZ22prevPhyTax <- BZ22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevPhyTax)



# Now BZ22 by Class

BZ22cla <- tax_glom(BZ22, "Class")

BZ22prevClaPS <- filter_taxa(BZ22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevClaTax <- tax_table(BZ22prevClaPS) %>% data.frame

BZ22prevClaTax <- BZ22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevClaTax)



# BZ22 by Order

BZ22ord <- tax_glom(BZ22, "Order")

BZ22prevOrdPS <- filter_taxa(BZ22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevOrdTax <- tax_table(BZ22prevOrdPS) %>% data.frame

BZ22prevOrdTax <- BZ22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevOrdTax)



# BZ22 by Family

BZ22fam <- tax_glom(BZ22, "Family")

BZ22prevFamPS <- filter_taxa(BZ22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevFamTax <- tax_table(BZ22prevFamPS) %>% data.frame

BZ22prevFamTax <- BZ22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevFamTax)



# BZ22 by Genus

BZ22gen <- tax_glom(BZ22, "Genus")

BZ22prevGenPS <- filter_taxa(BZ22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevGenTax <- tax_table(BZ22prevGenPS) %>% data.frame

BZ22prevGenTax <- BZ22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevGenTax)




# BZ22 by ASV

# (no agglomeration for ASV level)

BZ22prevAsvPS <- filter_taxa(BZ22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevAsvTax <- tax_table(BZ22prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, BZ22prevAsvTax)









# MO21

# Creating separate ps object for MO21

MO21 = subset_samples(PS_16S_Exp, Year == "2021" & Location == "MO")

MO21 = prune_taxa(taxa_sums(MO21)>0, MO21)   

MO21samples <- sample_data(MO21) %>% data.frame


# Prevalent Phyla in MO21

MO21phy <- tax_glom(MO21, "Phylum")   # agglomerate by Phylum

MO21prevPhyPS <- filter_taxa(MO21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated MO21_16S ps object to phyla with a count of at least 1 in at least 80% of samples

MO21prevPhyTax <- tax_table(MO21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

MO21prevPhyTax <- MO21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevPhyTax)



# Now MO21 by Class

MO21cla <- tax_glom(MO21, "Class")

MO21prevClaPS <- filter_taxa(MO21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevClaTax <- tax_table(MO21prevClaPS) %>% data.frame

MO21prevClaTax <- MO21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevClaTax)



# MO21 by Order

MO21ord <- tax_glom(MO21, "Order")

MO21prevOrdPS <- filter_taxa(MO21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevOrdTax <- tax_table(MO21prevOrdPS) %>% data.frame

MO21prevOrdTax <- MO21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevOrdTax)



# MO21 by Family

MO21fam <- tax_glom(MO21, "Family")

MO21prevFamPS <- filter_taxa(MO21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevFamTax <- tax_table(MO21prevFamPS) %>% data.frame

MO21prevFamTax <- MO21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevFamTax)



# MO21 by Genus

MO21gen <- tax_glom(MO21, "Genus")

MO21prevGenPS <- filter_taxa(MO21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevGenTax <- tax_table(MO21prevGenPS) %>% data.frame

MO21prevGenTax <- MO21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevGenTax)




# MO21 by ASV

# (no agglomeration for ASV level)

MO21prevAsvPS <- filter_taxa(MO21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevAsvTax <- tax_table(MO21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, MO21prevAsvTax)









# MO22

# Creating separate ps object for MO22

MO22 = subset_samples(PS_16S_Exp, Year == "2022" & Location == "MO")

MO22 = prune_taxa(taxa_sums(MO22)>0, MO22)   

MO22samples <- sample_data(MO22) %>% data.frame


# Prevalent Phyla in MO22

MO22phy <- tax_glom(MO22, "Phylum")   # agglomerate by Phylum

MO22prevPhyPS <- filter_taxa(MO22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated MO22_16S ps object to phyla with a count of at least 1 in at least 80% of samples

MO22prevPhyTax <- tax_table(MO22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

MO22prevPhyTax <- MO22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevPhyTax)



# Now MO22 by Class

MO22cla <- tax_glom(MO22, "Class")

MO22prevClaPS <- filter_taxa(MO22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevClaTax <- tax_table(MO22prevClaPS) %>% data.frame

MO22prevClaTax <- MO22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevClaTax)



# MO22 by Order

MO22ord <- tax_glom(MO22, "Order")

MO22prevOrdPS <- filter_taxa(MO22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevOrdTax <- tax_table(MO22prevOrdPS) %>% data.frame

MO22prevOrdTax <- MO22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevOrdTax)



# MO22 by Family

MO22fam <- tax_glom(MO22, "Family")

MO22prevFamPS <- filter_taxa(MO22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevFamTax <- tax_table(MO22prevFamPS) %>% data.frame

MO22prevFamTax <- MO22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevFamTax)



# MO22 by Genus

MO22gen <- tax_glom(MO22, "Genus")

MO22prevGenPS <- filter_taxa(MO22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevGenTax <- tax_table(MO22prevGenPS) %>% data.frame

MO22prevGenTax <- MO22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevGenTax)




# MO22 by ASV

# (no agglomeration for ASV level)

MO22prevAsvPS <- filter_taxa(MO22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevAsvTax <- tax_table(MO22prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, MO22prevAsvTax)









# SD21

# Creating separate ps object for SD21

SD21 = subset_samples(PS_16S_Exp, Year == "2021" & Location == "SD")

SD21 = prune_taxa(taxa_sums(SD21)>0, SD21)   

SD21samples <- sample_data(SD21) %>% data.frame


# Prevalent Phyla in SD21

SD21phy <- tax_glom(SD21, "Phylum")   # agglomerate by Phylum

SD21prevPhyPS <- filter_taxa(SD21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated SD21_16S ps object to phyla with a count of at least 1 in at least 80% of samples

SD21prevPhyTax <- tax_table(SD21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

SD21prevPhyTax <- SD21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevPhyTax)



# Now SD21 by Class

SD21cla <- tax_glom(SD21, "Class")

SD21prevClaPS <- filter_taxa(SD21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevClaTax <- tax_table(SD21prevClaPS) %>% data.frame

SD21prevClaTax <- SD21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevClaTax)



# SD21 by Order

SD21ord <- tax_glom(SD21, "Order")

SD21prevOrdPS <- filter_taxa(SD21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevOrdTax <- tax_table(SD21prevOrdPS) %>% data.frame

SD21prevOrdTax <- SD21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevOrdTax)



# SD21 by Family

SD21fam <- tax_glom(SD21, "Family")

SD21prevFamPS <- filter_taxa(SD21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevFamTax <- tax_table(SD21prevFamPS) %>% data.frame

SD21prevFamTax <- SD21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevFamTax)



# SD21 by Genus

SD21gen <- tax_glom(SD21, "Genus")

SD21prevGenPS <- filter_taxa(SD21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevGenTax <- tax_table(SD21prevGenPS) %>% data.frame

SD21prevGenTax <- SD21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevGenTax)




# SD21 by ASV

# (no agglomeration for ASV level)

SD21prevAsvPS <- filter_taxa(SD21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevAsvTax <- tax_table(SD21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, SD21prevAsvTax)









# SD22

# Creating separate ps object for SD22

SD22 = subset_samples(PS_16S_Exp, Year == "2022" & Location == "SD")

SD22 = prune_taxa(taxa_sums(SD22)>0, SD22)   

SD22samples <- sample_data(SD22) %>% data.frame


# Prevalent Phyla in SD22

SD22phy <- tax_glom(SD22, "Phylum")   # agglomerate by Phylum

SD22prevPhyPS <- filter_taxa(SD22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated SD22_16S ps object to phyla with a count of at least 1 in at least 80% of samples

SD22prevPhyTax <- tax_table(SD22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

SD22prevPhyTax <- SD22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevPhyTax)



# Now SD22 by Class

SD22cla <- tax_glom(SD22, "Class")

SD22prevClaPS <- filter_taxa(SD22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevClaTax <- tax_table(SD22prevClaPS) %>% data.frame

SD22prevClaTax <- SD22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevClaTax)



# SD22 by Order

SD22ord <- tax_glom(SD22, "Order")

SD22prevOrdPS <- filter_taxa(SD22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevOrdTax <- tax_table(SD22prevOrdPS) %>% data.frame

SD22prevOrdTax <- SD22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevOrdTax)



# SD22 by Family

SD22fam <- tax_glom(SD22, "Family")

SD22prevFamPS <- filter_taxa(SD22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevFamTax <- tax_table(SD22prevFamPS) %>% data.frame

SD22prevFamTax <- SD22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevFamTax)



# SD22 by Genus

SD22gen <- tax_glom(SD22, "Genus")

SD22prevGenPS <- filter_taxa(SD22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevGenTax <- tax_table(SD22prevGenPS) %>% data.frame

SD22prevGenTax <- SD22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevGenTax)




# SD22 by ASV

# (no agglomeration for ASV level)

SD22prevAsvPS <- filter_taxa(SD22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevAsvTax <- tax_table(SD22prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, SD22prevAsvTax)









# HI21

# Creating separate ps object for HI21

HI21 = subset_samples(PS_16S_Exp, Year == "2021" & Location == "HI")

HI21 = prune_taxa(taxa_sums(HI21)>0, HI21)   

HI21samples <- sample_data(HI21) %>% data.frame


# Prevalent Phyla in HI21

HI21phy <- tax_glom(HI21, "Phylum")   # agglomerate by Phylum

HI21prevPhyPS <- filter_taxa(HI21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated HI21_16S ps object to phyla with a count of at least 1 in at least 80% of samples

HI21prevPhyTax <- tax_table(HI21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

HI21prevPhyTax <- HI21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevPhyTax)



# Now HI21 by Class

HI21cla <- tax_glom(HI21, "Class")

HI21prevClaPS <- filter_taxa(HI21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevClaTax <- tax_table(HI21prevClaPS) %>% data.frame

HI21prevClaTax <- HI21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevClaTax)



# HI21 by Order

HI21ord <- tax_glom(HI21, "Order")

HI21prevOrdPS <- filter_taxa(HI21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevOrdTax <- tax_table(HI21prevOrdPS) %>% data.frame

HI21prevOrdTax <- HI21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevOrdTax)



# HI21 by Family

HI21fam <- tax_glom(HI21, "Family")

HI21prevFamPS <- filter_taxa(HI21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevFamTax <- tax_table(HI21prevFamPS) %>% data.frame

HI21prevFamTax <- HI21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevFamTax)



# HI21 by Genus

HI21gen <- tax_glom(HI21, "Genus")

HI21prevGenPS <- filter_taxa(HI21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevGenTax <- tax_table(HI21prevGenPS) %>% data.frame

HI21prevGenTax <- HI21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevGenTax)




# HI21 by ASV

# (no agglomeration for ASV level)

HI21prevAsvPS <- filter_taxa(HI21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevAsvTax <- tax_table(HI21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  relocate(ASV, .after = 8) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1) %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1)

PrevBac <- bind_rows(PrevBac, HI21prevAsvTax)






write_csv(PrevBac, "PrevalentBacteria.csv")









##### ITS, Fungi #####


# Reading in the PS object for ITS from all location-years, this also includes all check and bulk samples

PS_ITS_All = readRDS("ps.q20.trunc200.RDS")


# Removing bulk and check samples

AllSamples <- sample_data(PS_ITS_All) %>% data.frame

PS_ITS_Exp = subset_samples(PS_ITS_All, Line != "Bulk" & Line != "Hockett" & Line != "Lavina" & Line != "Merit57" & Line != "Odyssey")

ExpSamples <- sample_data(PS_ITS_Exp) %>% data.frame




# BZ21 ITS reads were too low quality for analysis





# BZ22

# Creating separate ps object for BZ22

BZ22 = subset_samples(PS_ITS_Exp, Year == "2022" & Location == "BZ")

BZ22samples <- sample_data(BZ22) %>% data.frame


# Prevalent Phyla in BZ22

BZ22phy <- tax_glom(BZ22, "Phylum")   # agglomerate by Phylum

BZ22prevPhyPS <- filter_taxa(BZ22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated BZ22_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

BZ22prevPhyTax <- tax_table(BZ22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

BZ22prevPhyTax <- BZ22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- BZ22prevPhyTax




# Now BZ22 by Class

BZ22cla <- tax_glom(BZ22, "Class")

BZ22prevClaPS <- filter_taxa(BZ22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevClaTax <- tax_table(BZ22prevClaPS) %>% data.frame

BZ22prevClaTax <- BZ22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, BZ22prevClaTax)



# BZ22 by Order

BZ22ord <- tax_glom(BZ22, "Order")

BZ22prevOrdPS <- filter_taxa(BZ22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevOrdTax <- tax_table(BZ22prevOrdPS) %>% data.frame

BZ22prevOrdTax <- BZ22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, BZ22prevOrdTax)



# BZ22 by Family

BZ22fam <- tax_glom(BZ22, "Family")

BZ22prevFamPS <- filter_taxa(BZ22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevFamTax <- tax_table(BZ22prevFamPS) %>% data.frame

BZ22prevFamTax <- BZ22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, BZ22prevFamTax)



# BZ22 by Genus

BZ22gen <- tax_glom(BZ22, "Genus")

BZ22prevGenPS <- filter_taxa(BZ22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

BZ22prevGenTax <- tax_table(BZ22prevGenPS) %>% data.frame

BZ22prevGenTax <- BZ22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, BZ22prevGenTax)



# BZ22 by ASV

BZ22prevAsvPS <- filter_taxa(BZ22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the BZ22 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

BZ22prevAsvTax <- tax_table(BZ22prevAsvPS) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "BZ22") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column

PrevFun <- bind_rows(PrevFun, BZ22prevAsvTax)









# MO21

# Creating separate ps object for MO21

MO21 = subset_samples(PS_ITS_Exp, Year == "2021" & Location == "MO")

MO21samples <- sample_data(MO21) %>% data.frame


# Prevalent Phyla in MO21

MO21phy <- tax_glom(MO21, "Phylum")   # agglomerate by Phylum

MO21prevPhyPS <- filter_taxa(MO21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated MO21_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

MO21prevPhyTax <- tax_table(MO21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

MO21prevPhyTax <- MO21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO21prevPhyTax)




# Now MO21 by Class

MO21cla <- tax_glom(MO21, "Class")

MO21prevClaPS <- filter_taxa(MO21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevClaTax <- tax_table(MO21prevClaPS) %>% data.frame

MO21prevClaTax <- MO21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO21prevClaTax)



# MO21 by Order

MO21ord <- tax_glom(MO21, "Order")

MO21prevOrdPS <- filter_taxa(MO21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevOrdTax <- tax_table(MO21prevOrdPS) %>% data.frame

MO21prevOrdTax <- MO21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO21prevOrdTax)



# MO21 by Family

MO21fam <- tax_glom(MO21, "Family")

MO21prevFamPS <- filter_taxa(MO21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevFamTax <- tax_table(MO21prevFamPS) %>% data.frame

MO21prevFamTax <- MO21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO21prevFamTax)



# MO21 by Genus

MO21gen <- tax_glom(MO21, "Genus")

MO21prevGenPS <- filter_taxa(MO21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO21prevGenTax <- tax_table(MO21prevGenPS) %>% data.frame

MO21prevGenTax <- MO21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO21prevGenTax)



# MO21 by ASV

MO21prevAsvPS <- filter_taxa(MO21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the MO21 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

MO21prevAsvTax <- tax_table(MO21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO21") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column


PrevFun <- bind_rows(PrevFun, MO21prevAsvTax)









# MO22

# Creating separate ps object for MO22

MO22 = subset_samples(PS_ITS_Exp, Year == "2022" & Location == "MO")

MO22samples <- sample_data(MO22) %>% data.frame


# Prevalent Phyla in MO22

MO22phy <- tax_glom(MO22, "Phylum")   # agglomerate by Phylum

MO22prevPhyPS <- filter_taxa(MO22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated MO22_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

MO22prevPhyTax <- tax_table(MO22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

MO22prevPhyTax <- MO22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO22prevPhyTax)




# Now MO22 by Class

MO22cla <- tax_glom(MO22, "Class")

MO22prevClaPS <- filter_taxa(MO22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevClaTax <- tax_table(MO22prevClaPS) %>% data.frame

MO22prevClaTax <- MO22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO22prevClaTax)



# MO22 by Order

MO22ord <- tax_glom(MO22, "Order")

MO22prevOrdPS <- filter_taxa(MO22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevOrdTax <- tax_table(MO22prevOrdPS) %>% data.frame

MO22prevOrdTax <- MO22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO22prevOrdTax)



# MO22 by Family

MO22fam <- tax_glom(MO22, "Family")

MO22prevFamPS <- filter_taxa(MO22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevFamTax <- tax_table(MO22prevFamPS) %>% data.frame

MO22prevFamTax <- MO22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO22prevFamTax)



# MO22 by Genus

MO22gen <- tax_glom(MO22, "Genus")

MO22prevGenPS <- filter_taxa(MO22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

MO22prevGenTax <- tax_table(MO22prevGenPS) %>% data.frame

MO22prevGenTax <- MO22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, MO22prevGenTax)



# MO22 by ASV

MO22prevAsvPS <- filter_taxa(MO22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the MO22 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

MO22prevAsvTax <- tax_table(MO22prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "MO22") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column

PrevFun <- bind_rows(PrevFun, MO22prevAsvTax)









# SD21

# Creating separate ps object for SD21

SD21 = subset_samples(PS_ITS_Exp, Year == "2021" & Location == "SD")

SD21samples <- sample_data(SD21) %>% data.frame


# Prevalent Phyla in SD21

SD21phy <- tax_glom(SD21, "Phylum")   # agglomerate by Phylum

SD21prevPhyPS <- filter_taxa(SD21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated SD21_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

SD21prevPhyTax <- tax_table(SD21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

SD21prevPhyTax <- SD21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD21prevPhyTax)




# Now SD21 by Class

SD21cla <- tax_glom(SD21, "Class")

SD21prevClaPS <- filter_taxa(SD21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevClaTax <- tax_table(SD21prevClaPS) %>% data.frame

SD21prevClaTax <- SD21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD21prevClaTax)



# SD21 by Order

SD21ord <- tax_glom(SD21, "Order")

SD21prevOrdPS <- filter_taxa(SD21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevOrdTax <- tax_table(SD21prevOrdPS) %>% data.frame

SD21prevOrdTax <- SD21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD21prevOrdTax)



# SD21 by Family

SD21fam <- tax_glom(SD21, "Family")

SD21prevFamPS <- filter_taxa(SD21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevFamTax <- tax_table(SD21prevFamPS) %>% data.frame

SD21prevFamTax <- SD21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD21prevFamTax)



# SD21 by Genus

SD21gen <- tax_glom(SD21, "Genus")

SD21prevGenPS <- filter_taxa(SD21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD21prevGenTax <- tax_table(SD21prevGenPS) %>% data.frame

SD21prevGenTax <- SD21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD21prevGenTax)



# SD21 by ASV

SD21prevAsvPS <- filter_taxa(SD21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the SD21 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

SD21prevAsvTax <- tax_table(SD21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD21") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column

PrevFun <- bind_rows(PrevFun, SD21prevAsvTax)









# SD22

# Creating separate ps object for SD22

SD22 = subset_samples(PS_ITS_Exp, Year == "2022" & Location == "SD")

SD22samples <- sample_data(SD22) %>% data.frame


# Prevalent Phyla in SD22

SD22phy <- tax_glom(SD22, "Phylum")   # agglomerate by Phylum

SD22prevPhyPS <- filter_taxa(SD22phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated SD22_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

SD22prevPhyTax <- tax_table(SD22prevPhyPS) %>% data.frame   # creating a data frame of that tax table

SD22prevPhyTax <- SD22prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD22prevPhyTax)




# Now SD22 by Class

SD22cla <- tax_glom(SD22, "Class")

SD22prevClaPS <- filter_taxa(SD22cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevClaTax <- tax_table(SD22prevClaPS) %>% data.frame

SD22prevClaTax <- SD22prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD22prevClaTax)



# SD22 by Order

SD22ord <- tax_glom(SD22, "Order")

SD22prevOrdPS <- filter_taxa(SD22ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevOrdTax <- tax_table(SD22prevOrdPS) %>% data.frame

SD22prevOrdTax <- SD22prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD22prevOrdTax)



# SD22 by Family

SD22fam <- tax_glom(SD22, "Family")

SD22prevFamPS <- filter_taxa(SD22fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevFamTax <- tax_table(SD22prevFamPS) %>% data.frame

SD22prevFamTax <- SD22prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD22prevFamTax)



# SD22 by Genus

SD22gen <- tax_glom(SD22, "Genus")

SD22prevGenPS <- filter_taxa(SD22gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

SD22prevGenTax <- tax_table(SD22prevGenPS) %>% data.frame

SD22prevGenTax <- SD22prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, SD22prevGenTax)



# SD22 by ASV

SD22prevAsvPS <- filter_taxa(SD22, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the SD22 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

SD22prevAsvTax <- tax_table(SD22prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "SD22") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column

PrevFun <- bind_rows(PrevFun, SD22prevAsvTax)









# HI21

# Creating separate ps object for HI21

HI21 = subset_samples(PS_ITS_Exp, Year == "2021" & Location == "HI")

HI21samples <- sample_data(HI21) %>% data.frame


# Prevalent Phyla in HI21

HI21phy <- tax_glom(HI21, "Phylum")   # agglomerate by Phylum

HI21prevPhyPS <- filter_taxa(HI21phy, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the Phylum agglomerated HI21_ITS ps object to phyla with a count of at least 1 in at least 80% of samples

HI21prevPhyTax <- tax_table(HI21prevPhyPS) %>% data.frame   # creating a data frame of that tax table

HI21prevPhyTax <- HI21prevPhyTax %>% 
  mutate(TaxLevel = "Phylum") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, HI21prevPhyTax)




# Now HI21 by Class

HI21cla <- tax_glom(HI21, "Class")

HI21prevClaPS <- filter_taxa(HI21cla, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevClaTax <- tax_table(HI21prevClaPS) %>% data.frame

HI21prevClaTax <- HI21prevClaTax %>% 
  mutate(TaxLevel = "Class") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, HI21prevClaTax)



# HI21 by Order

HI21ord <- tax_glom(HI21, "Order")

HI21prevOrdPS <- filter_taxa(HI21ord, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevOrdTax <- tax_table(HI21prevOrdPS) %>% data.frame

HI21prevOrdTax <- HI21prevOrdTax %>% 
  mutate(TaxLevel = "Order") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, HI21prevOrdTax)



# HI21 by Family

HI21fam <- tax_glom(HI21, "Family")

HI21prevFamPS <- filter_taxa(HI21fam, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevFamTax <- tax_table(HI21prevFamPS) %>% data.frame

HI21prevFamTax <- HI21prevFamTax %>% 
  mutate(TaxLevel = "Family") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, HI21prevFamTax)



# HI21 by Genus

HI21gen <- tax_glom(HI21, "Genus")

HI21prevGenPS <- filter_taxa(HI21gen, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)

HI21prevGenTax <- tax_table(HI21prevGenPS) %>% data.frame

HI21prevGenTax <- HI21prevGenTax %>% 
  mutate(TaxLevel = "Genus") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1)

PrevFun <- bind_rows(PrevFun, HI21prevGenTax)



# HI21 by ASV

HI21prevAsvPS <- filter_taxa(HI21, function(x) sum(x > 0) >= (0.80*length(x)), TRUE)   # filtering the HI21 ps object to ASVs with a count of at least 1 in at least 80% of samples (no agglomerating needed)

HI21prevAsvTax <- tax_table(HI21prevAsvPS) %>%
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  mutate(TaxLevel = "ASV") %>% 
  relocate(TaxLevel, .before = 1) %>% 
  mutate(LocYr = "HI21") %>% 
  relocate(LocYr, .before = 1) %>% 
  select("TaxLevel", "LocYr", "ASV")                          # creating a data frame of that tax table, with ASV names in a column

PrevFun <- bind_rows(PrevFun, HI21prevAsvTax)









write_csv(PrevFun, "PrevalentFungi.csv")


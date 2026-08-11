
### The following script contains code to get from phyloseq objects prepared from 16S and ITS sequencing data to phenotype files ready for GWAS.
# It utilizes lists of prevalent taxa generated in the "PrevalentTaxa.R" script.

# Here the swMT location is referred to as BZ and the cntrMT location is referred to as MO.


library(tidyverse)
library(phyloseq)
library(zCompositions)



# Reading in the PS object for 16S from all location-years, this also includes all check and bulk samples

PS_16S_All = readRDS("ps.filt.RDS")

# Removing bulk samples

PS_16S_Exp = subset_samples(PS_16S_All, Line != "Bulk")



# Reading in the PS object for ITS from all location-years, this also includes all check and bulk samples

PS_ITS_All = readRDS("ps.q20.trunc200.RDS")

# Removing bulk samples

PS_ITS_Exp = subset_samples(PS_ITS_All, Line != "Bulk")

# Removing an outlier sample from SD21_ITS that had far far fewer zeros than any other sample, creating issues during zero imputation and transformation

PS_ITS_Exp = subset_samples(PS_ITS_Exp, !(Plot == "208" & Year == "2021" & Location == "SD"))



##### PREP SECTION #####

# Set the phyloseq object (PSobj) to a data type (16S or ITS) and enter a year (2021 or 2022) and a location (BZ, MO, SD, or HI).
# Then run through the prep code (cycling through taxonomic ranks below), save the resulting dataframe with a name indicating the data type and location-year (LocYr), then repeat for the other type/LocYr datasets.

PSobj = subset_samples(PS_16S_Exp, Year == "2021" & Location == "BZ")     # set location, year, and 16S or ITS in this line

PSobj = prune_taxa(taxa_sums(PSobj)>0, PSobj)   # so that OTU and tax tables only include samples from the given LocYr

# ntaxa(PSobj)    # will give the number of taxa in the OTU to check that they were indeed pruned


# Determining which highly sparse ASVs to remove by filtering the PS object, then saving that list of taxa.
# The requirements for retaining an ASV are at least 0.7% prevalence (present in at least 0.7% of samples which works out to be 2 samples)
# and at least 0.1% maximum relative abundance (RA) (based on within-sample RA rather than counts or across-dataset RA so as not to be influenced by sample differences in read-depth).
# Since one of the filters is based on RA, it is important to convert OTU table counts to RA before prevalence filtering so that these RAs are not influenced by that filtering.

# Transforming data in the phyloseq object into relative abundances

PSobj_filt = transform_sample_counts(PSobj, function(x) x / sum(x))


# Filter to ASVs with at least 0.7% prevalence (RA greater than 0 in at least 0.7% of samples)

PSobj_filt <- filter_taxa(PSobj_filt, function(x) sum(x > 0) > 0.007 * length(x), prune = TRUE)   # prune = TRUE makes sure to create the filtered PS object, not just return a vector of the taxa fitting the function


# Further filter to taxa with a relative abundance of at least 0.1% in at least one sample

PSobj_filt = filter_taxa(PSobj_filt, function(x) max(x) >= 0.001, TRUE)

# ntaxa(PSobj_filt)   checking effect of filtering on ASV number


# Subsetting the original PS object (with the raw counts) to only contain the ASVs retained after filtering.

filtASVs <- taxa_names(PSobj_filt)   # making a vector of ASV names/numbers (these stay the same through processing steps) to keep

PSobj <- prune_taxa(taxa_names(PSobj) %in% filtASVs, PSobj)   # pruning the PS object to only include these ASVs



### TAXA RANKS ###

# Agglomerating the filtered PSobj of counts to a taxonomic rank.
# Set this, run through the rest of the prep code, save the dataframe for the rank, repeat this section for the other ranks, then join these into a final dataframe for the LocYr/type dataset.

PSobj_glom <- tax_glom(PSobj, "Phylum")   # agglomerate PS object by Phylum/Class/Order/Family/Genus


# Pull out otu table, impute zeros, clr transform, and reassemble PSobj with this edited otu table.

PSobj_otu <- otu_table(PSobj_glom) %>% data.frame   # pull out otu table as df


PSobj_otuGBM <- cmultRepl(PSobj_otu, label = 0, output = "prop", method = "GBM", z.warning = 1, adjust = TRUE)   # impute 0s
# label = 0 means that items labeled "0" in the count table are the things to be imputed
# output = "prop" means the output will be proportions.
# z.warning indicates what percentage of zeros is unacceptable for a column(taxon) to have. Setting to 1 will not remove any.
# adjust = TRUE indicates that non-zero values in the row should be adjusted in a way that preserves the ratio between them while also ensuring that the sum of proportions is 1 for each row


# Converting post-imputation proportions to pseudo-counts (by multiplying by pre-filtering sample count totals).
# This way (post ASV filtering) sample read depths stay the same and non-zero count ratios stay the same while zeros are imputed.

depths <- rowSums(PSobj_otu)

PSobj_otuGBM <- PSobj_otuGBM * depths


# Performing clr transformations

PSobj_otuGBM <- clr(PSobj_otuGBM) %>% data.frame


# Creating a new PSobj with the transformed otu table

PSobj_GBM <- phyloseq(otu_table(PSobj_otuGBM, taxa_are_rows = FALSE),
                      sample_data(sample_data(PSobj_glom)),
                      tax_table(tax_table(PSobj_glom)))


# Rename PSobj for the data type, LocYr, and taxonomic rank.
# Then repeat the "TAXA RANKS" section for each taxonomic rank in the type/LocYr dataset, and the "PREP SECTION" for each type/LocYr dataset.

PS_16S_BZ21_gen <- PSobj_GBM   # phy/cla/ord/fam/gen for rank



# This is a good place to save an .RData of transformed PS objects at each rank for each dataset.

save(PS_16S_BZ21_phy, PS_16S_BZ21_cla, PS_16S_BZ21_ord, PS_16S_BZ21_fam, PS_16S_BZ21_gen, 
     PS_16S_BZ22_phy, PS_16S_BZ22_cla, PS_16S_BZ22_ord, PS_16S_BZ22_fam, PS_16S_BZ22_gen,
     PS_16S_MO21_phy, PS_16S_MO21_cla, PS_16S_MO21_ord, PS_16S_MO21_fam, PS_16S_MO21_gen, 
     PS_16S_MO22_phy, PS_16S_MO22_cla, PS_16S_MO22_ord, PS_16S_MO22_fam, PS_16S_MO22_gen, 
     PS_16S_SD21_phy, PS_16S_SD21_cla, PS_16S_SD21_ord, PS_16S_SD21_fam, PS_16S_SD21_gen, 
     PS_16S_SD22_phy, PS_16S_SD22_cla, PS_16S_SD22_ord, PS_16S_SD22_fam, PS_16S_SD22_gen, 
     PS_16S_HI21_phy, PS_16S_HI21_cla, PS_16S_HI21_ord, PS_16S_HI21_fam, PS_16S_HI21_gen,
     PS_ITS_BZ22_phy, PS_ITS_BZ22_cla, PS_ITS_BZ22_ord, PS_ITS_BZ22_fam, PS_ITS_BZ22_gen, 
     PS_ITS_MO21_phy, PS_ITS_MO21_cla, PS_ITS_MO21_ord, PS_ITS_MO21_fam, PS_ITS_MO21_gen,  
     PS_ITS_MO22_phy, PS_ITS_MO22_cla, PS_ITS_MO22_ord, PS_ITS_MO22_fam, PS_ITS_MO22_gen,  
     PS_ITS_SD21_phy, PS_ITS_SD21_cla, PS_ITS_SD21_ord, PS_ITS_SD21_fam, PS_ITS_SD21_gen,  
     PS_ITS_SD22_phy, PS_ITS_SD22_cla, PS_ITS_SD22_ord, PS_ITS_SD22_fam, PS_ITS_SD22_gen,  
     PS_ITS_HI21_phy, PS_ITS_HI21_cla, PS_ITS_HI21_ord, PS_ITS_HI21_fam, PS_ITS_HI21_gen, 
     file = "FilteredTransformed_PSobjects.RData")







##### FILTERING TAXA FOR PREVALENCE #####

# Taxa within each datatype (16S or ITS) by location-year (LocYr) dataset/phyloseq object (PSobj) need to be filtered to only include prevalent taxa at each taxonomic rank.
# Lists of prevalent taxa for each dataset are derived in the "PrevalentTaxaListing.R" script.


# Read-in csv files with prevalent taxa lists.
# These are taxa that have a count of at least 1 in at least 80% of samples within a data type / location-year / taxonomic rank.

PrevBac <- read.csv("PrevalentBacteria.csv")

PrevFun <- read.csv("PrevalentFungi.csv")



# Load saved filtered, imputed, transformed PS objects if not loaded already.

load("FilteredTransformed_PSobjects.RData")





#### Filtering at each rank ####


# Run through each of the five taxonomic ranks, setting the dataset at the beginning of each rank, and saving at the end of all the ranks:
  
### Phylum
  
# set and subset PS object

PSobj <- PS_type_LocYr_phy   ##### set to whichever PS object for data type (16S or ITS), and location-year (LocYr) #####

PrevTax <- PrevType %>% filter(LocYr == "", TaxLevel == "Phylum")   ##### change Prev"Type" to either PrevBac or PrevFun and enter whichever LocYr #####

PrevTax <- PrevTax$Phylum

PSobj <- subset_taxa(PSobj, Phylum %in% PrevTax)



# transposes the OTU data frame and creates a column of the ASV row names

PSotu <- otu_table(PSobj) %>%
  data.frame %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ASV")


# creating a taxa df with only ASV name and taxa columns, and adding a rank suffix to the names

PStaxa <- tax_table(PSobj) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  select("ASV", "Phylum")%>% 
  mutate(Phylum = paste0(Phylum, "_phy"))


# Joining the otu and taxa dfs by ASV names, removing the ASV col, transposing, turning back into a data frame (which adds "x" in front of sample #s), and creating a SampleID col

phy_dataset <- PStaxa %>% 
  full_join(PSotu, by = "ASV") %>% 
  select(-"ASV") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("SampleID")


# Getting rid of the "Xs" in a new col called Sample

phy_dataset$Sample <- sub("X", "", phy_dataset$SampleID)


# Moving that Sample col to the front and deleting the SampleID col with the Xs

phy_dataset <- phy_dataset[ , c("Sample", names(phy_dataset)[names(phy_dataset) != "Sample"])] %>% 
  select(-"SampleID")


# Turning the taxa names that are in the top row into column names, which results in the taxa rank title now being the title for the Sample column

colnames(phy_dataset) <- phy_dataset[1, ]


# getting rid of the top row that was used for col names, and renaming the sample col as "Sample" instead of the taxa rank

phy_dataset <- phy_dataset[-1, ] %>% 
  rename(Sample = Phylum)




### Class

# set and subset PS object

PSobj <- PS_type_LocYr_cla   ##### set to whichever PS object for data type, LocYr, tax rank #####

PrevTax <- PrevFun %>% filter(LocYr == "", TaxLevel == "Class")   ##### set to either Bac or Fun and whichever LocYr #####

PrevTax <- PrevTax$Class

PSobj <- subset_taxa(PSobj, Class %in% PrevTax)



# transposes the OTU data frame and creates a column of the ASV row names

PSotu <- otu_table(PSobj) %>%
  data.frame %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ASV")


# creating a taxa df with only ASV name and taxa columns, and adding a rank suffix to the names

PStaxa <- tax_table(PSobj) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  select("ASV", "Class")%>% 
  mutate(Class = paste0(Class, "_cla"))


# Joining the otu and taxa dfs by ASV names, removing the ASV col, transposing, turning back into a data frame (which adds "x" in front of sample #s), and creating a SampleID col

cla_dataset <- PStaxa %>% 
  full_join(PSotu, by = "ASV") %>% 
  select(-"ASV") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("SampleID")


# Getting rid of the "Xs" in a new col called Sample

cla_dataset$Sample <- sub("X", "", cla_dataset$SampleID)


# Moving that Sample col to the front and deleting the SampleID col with the Xs

cla_dataset <- cla_dataset[ , c("Sample", names(cla_dataset)[names(cla_dataset) != "Sample"])] %>% 
  select(-"SampleID")


# Turning the taxa names that are in the top row into column names, which results in the taxa rank title now being the title for the Sample column

colnames(cla_dataset) <- cla_dataset[1, ]


# getting rid of the top row that was used for col names, and renaming the sample col as "Sample" instead of the taxa rank

cla_dataset <- cla_dataset[-1, ] %>% 
  rename(Sample = Class)



### Order

# set and subset PS object

PSobj <- PS_type_LocYr_ord   ##### set to whichever PS object for data type, LocYr, tax rank #####

PrevTax <- PrevFun %>% filter(LocYr == "", TaxLevel == "Order")   ##### set to either Bac or Fun and whichever LocYr #####

PrevTax <- PrevTax$Order

PSobj <- subset_taxa(PSobj, Order %in% PrevTax)



# transposes the OTU data frame and creates a column of the ASV row names

PSotu <- otu_table(PSobj) %>%
  data.frame %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ASV")


# creating a taxa df with only ASV name and taxa columns, and adding a rank suffix to the names

PStaxa <- tax_table(PSobj) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  select("ASV", "Order")%>% 
  mutate(Order = paste0(Order, "_ord"))


# Joining the otu and taxa dfs by ASV names, removing the ASV col, transposing, turning back into a data frame (which adds "x" in front of sample #s), and creating a SampleID col

ord_dataset <- PStaxa %>% 
  full_join(PSotu, by = "ASV") %>% 
  select(-"ASV") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("SampleID")


# Getting rid of the "Xs" in a new col called Sample

ord_dataset$Sample <- sub("X", "", ord_dataset$SampleID)


# Moving that Sample col to the front and deleting the SampleID col with the Xs

ord_dataset <- ord_dataset[ , c("Sample", names(ord_dataset)[names(ord_dataset) != "Sample"])] %>% 
  select(-"SampleID")


# Turning the taxa names that are in the top row into column names, which results in the taxa rank title now being the title for the Sample column

colnames(ord_dataset) <- ord_dataset[1, ]


# getting rid of the top row that was used for col names, and renaming the sample col as "Sample" instead of the taxa rank

ord_dataset <- ord_dataset[-1, ] %>% 
  rename(Sample = Order)



### Family

# set and subset PS object

PSobj <- PS_type_LocYr_fam   ##### set to whichever PS object for data type, LocYr, tax rank #####

PrevTax <- PrevFun %>% filter(LocYr == "", TaxLevel == "Family")   ##### set to either Bac or Fun and whichever LocYr #####

PrevTax <- PrevTax$Family

PSobj <- subset_taxa(PSobj, Family %in% PrevTax)



# transposes the OTU data frame and creates a column of the ASV row names

PSotu <- otu_table(PSobj) %>%
  data.frame %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ASV")


# creating a taxa df with only ASV name and taxa columns, and adding a rank suffix to the names

PStaxa <- tax_table(PSobj) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  select("ASV", "Family")%>% 
  mutate(Family = paste0(Family, "_fam"))


# Joining the otu and taxa dfs by ASV names, removing the ASV col, transposing, turning back into a data frame (which adds "x" in front of sample #s), and creating a SampleID col

fam_dataset <- PStaxa %>% 
  full_join(PSotu, by = "ASV") %>% 
  select(-"ASV") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("SampleID")


# Getting rid of the "Xs" in a new col called Sample

fam_dataset$Sample <- sub("X", "", fam_dataset$SampleID)


# Moving that Sample col to the front and deleting the SampleID col with the Xs

fam_dataset <- fam_dataset[ , c("Sample", names(fam_dataset)[names(fam_dataset) != "Sample"])] %>% 
  select(-"SampleID")


# Turning the taxa names that are in the top row into column names, which results in the taxa rank title now being the title for the Sample column

colnames(fam_dataset) <- fam_dataset[1, ]


# for Hawaii 2021 Bacteria there were two unknown families that needed their columns renamed:

# colnames(fam_dataset)[28] <- "Unknown1_fam"
# colnames(fam_dataset)[57] <- "Unknown2_fam"


# getting rid of the top row that was used for col names, and renaming the sample col as "Sample" instead of the taxa rank

fam_dataset <- fam_dataset[-1, ] %>% 
  rename(Sample = Family)



### Genus

# set and subset PS object

PSobj <- PS_type_LocYr_gen   ##### set to whichever PS object for data type, LocYr, tax rank #####

PrevTax <- PrevFun %>% filter(LocYr == "", TaxLevel == "Genus")   ##### set to either Bac or Fun and whichever LocYr #####

PrevTax <- PrevTax$Genus

PSobj <- subset_taxa(PSobj, Genus %in% PrevTax)



# transposes the OTU data frame and creates a column of the ASV row names

PSotu <- otu_table(PSobj) %>%
  data.frame %>% 
  t %>% 
  data.frame %>% 
  rownames_to_column("ASV")


# creating a taxa df with only ASV name and taxa columns, and adding a rank suffix to the names

PStaxa <- tax_table(PSobj) %>% 
  data.frame %>% 
  rownames_to_column("ASV") %>% 
  select("ASV", "Genus")%>% 
  mutate(Genus = paste0(Genus, "_gen"))


# Joining the otu and taxa dfs by ASV names, removing the ASV col, transposing, turning back into a data frame (which adds "x" in front of sample #s), and creating a SampleID col

gen_dataset <- PStaxa %>% 
  full_join(PSotu, by = "ASV") %>% 
  select(-"ASV") %>% 
  t %>% 
  as.data.frame %>% 
  rownames_to_column("SampleID")


# Getting rid of the "Xs" in a new col called Sample

gen_dataset$Sample <- sub("X", "", gen_dataset$SampleID)


# Moving that Sample col to the front and deleting the SampleID col with the Xs

gen_dataset <- gen_dataset[ , c("Sample", names(gen_dataset)[names(gen_dataset) != "Sample"])] %>% 
  select(-"SampleID")


# Turning the taxa names that are in the top row into column names, which results in the taxa rank title now being the title for the Sample column

colnames(gen_dataset) <- gen_dataset[1, ]


# getting rid of the top row that was used for col names, and renaming the sample col as "Sample" instead of the taxa rank

gen_dataset <- gen_dataset[-1, ] %>% 
  rename(Sample = Genus)



### Saving

# Save dataset dataframe

df_list <- list(phy_dataset, cla_dataset, ord_dataset, fam_dataset, gen_dataset)

PreppedTaxa_type_LocYr <- reduce(df_list, full_join, by = "Sample")   ##### set to whichever data type (16S/ITS) and LocYr #####



#### Saving after running through all the ranks for each dataset ####

# Save .RData of each dataset's dataframe of prepared (filtered, imputed, transformed, subsetted for prevalence) taxa values at each rank.

save(PreppedTaxa_16S_BZ21, PreppedTaxa_16S_BZ22, PreppedTaxa_16S_MO21, PreppedTaxa_16S_MO22, PreppedTaxa_16S_SD21, PreppedTaxa_16S_SD22, PreppedTaxa_16S_HI21,
     PreppedTaxa_ITS_BZ22, PreppedTaxa_ITS_MO21, PreppedTaxa_ITS_MO22, PreppedTaxa_ITS_SD21, PreppedTaxa_ITS_SD22, PreppedTaxa_ITS_HI21,
     file = "PreparedTaxaDFs.RData")











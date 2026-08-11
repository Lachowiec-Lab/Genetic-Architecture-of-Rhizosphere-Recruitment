
### The following script contains code for modeling prepared microbial traits (transformed taxa abundances), and agronomic traits
# to calculate heritability, block effect, and BLUPs where necessary,
# then preparing text files of phenotype data for input into GAPIT for GWAS.

# Here the swMT location is referred to as BZ and the cntrMT location is referred to as MO.


library(dplyr)
library(lmerTest)



##### Load and prepare data #####


# Load dataframes of prepared (filtered, imputed, transformed, subsetted for prevalence) taxa trait values generated in the "DataPrepForModeling.R" script.

load("PreparedTaxaDFs.RData")


# Make sure that all trait columns while numeric while leaving the first sample column as a character

PreppedTaxa_16S_BZ21[ , 2:ncol(PreppedTaxa_16S_BZ21)] <- lapply(PreppedTaxa_16S_BZ21[ , 2:ncol(PreppedTaxa_16S_BZ21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_BZ22[ , 2:ncol(PreppedTaxa_16S_BZ22)] <- lapply(PreppedTaxa_16S_BZ22[ , 2:ncol(PreppedTaxa_16S_BZ22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_MO21[ , 2:ncol(PreppedTaxa_16S_MO21)] <- lapply(PreppedTaxa_16S_MO21[ , 2:ncol(PreppedTaxa_16S_MO21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_MO22[ , 2:ncol(PreppedTaxa_16S_MO22)] <- lapply(PreppedTaxa_16S_MO22[ , 2:ncol(PreppedTaxa_16S_MO22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_SD21[ , 2:ncol(PreppedTaxa_16S_SD21)] <- lapply(PreppedTaxa_16S_SD21[ , 2:ncol(PreppedTaxa_16S_SD21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_SD22[ , 2:ncol(PreppedTaxa_16S_SD22)] <- lapply(PreppedTaxa_16S_SD22[ , 2:ncol(PreppedTaxa_16S_SD22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_16S_HI21[ , 2:ncol(PreppedTaxa_16S_HI21)] <- lapply(PreppedTaxa_16S_HI21[ , 2:ncol(PreppedTaxa_16S_HI21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_BZ22[ , 2:ncol(PreppedTaxa_ITS_BZ22)] <- lapply(PreppedTaxa_ITS_BZ22[ , 2:ncol(PreppedTaxa_ITS_BZ22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_MO21[ , 2:ncol(PreppedTaxa_ITS_MO21)] <- lapply(PreppedTaxa_ITS_MO21[ , 2:ncol(PreppedTaxa_ITS_MO21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_MO22[ , 2:ncol(PreppedTaxa_ITS_MO22)] <- lapply(PreppedTaxa_ITS_MO22[ , 2:ncol(PreppedTaxa_ITS_MO22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_SD21[ , 2:ncol(PreppedTaxa_ITS_SD21)] <- lapply(PreppedTaxa_ITS_SD21[ , 2:ncol(PreppedTaxa_ITS_SD21)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_SD22[ , 2:ncol(PreppedTaxa_ITS_SD22)] <- lapply(PreppedTaxa_ITS_SD22[ , 2:ncol(PreppedTaxa_ITS_SD22)],
                                                                function(x) as.numeric(as.character(x)))

PreppedTaxa_ITS_HI21[ , 2:ncol(PreppedTaxa_ITS_HI21)] <- lapply(PreppedTaxa_ITS_HI21[ , 2:ncol(PreppedTaxa_ITS_HI21)],
                                                                function(x) as.numeric(as.character(x)))





##### SET DATASET AND DO MODELING #####


#####
### Use this code to prepare dataframe for MICROBIAL trait datasets (skipping agronomic trait code below and jumping to "Add factor columns")

# Save whichever dataset's dataframe as Phen, run through all code through making the BLUPs dataframe, rerun for each dataset, then save .RData.

Phen <- PreppedTaxa_type_LocYr     ##### set to whichever data type (16S or ITS) and location-year (LocYr) #####

# Make sure that column names do not contain problematic characters (as sometimes happens with microbial taxa names)
names(Phen) <- make.names(names(Phen), unique = TRUE)


# Add Line and Block to the Phen data frame after setting for correct field trial.

### Pulling in plots, barley line names, and blocks from the field data file

field <- read.csv("RawFieldAgData.csv")

field <- field %>%
  filter(Year == "", Location == "") %>%      ##### set Year and Location #####
select(Plot, Line, Block)

field$Plot <- as.character(field$Plot)


Phen <- Phen %>% 
  mutate(Plot = sub("_.*", "", Sample)) %>% 
  select(-Sample) %>% 
  relocate(Plot, .before = 1)

Phen <- left_join(field, Phen, by = "Plot")
#####




#####
### Use this code to prepare dataframe for AGRONOMIC trait datasets (skipping code for microbial trait datasets above)

field <- read.csv("RawFieldAgData.csv")

Phen <- field %>% 
  select(-Year, -Location, -Row, -Range, -Plot, -Entry)     # remove unnecessary columns (keeping only LocYr, Block, Line, trait columns)

Phen <- Phen %>% 
  filter(LocYr == "MO22") %>%               ##### set LocYr #####
select(where(~ !all(is.na(.)))) %>%       # remove trait columns that are missing for the given LocYr
  rename_with(                              # add _LocYr suffixes to trait column names to help with file prep later
    ~ paste0(.x, "_LocYr"),                 ##### set LocYr #####
    .cols = -c(LocYr, Block, Line)
  )
#####



### Add factor columns

# Add columns and designate factors for modeling utilizing the randomized augmented design of the trial 
# in which check varieties are replicated at random positions in each block, 
# allowing for environmental effects captured by block to be accounted for in the un-replicated experimental barley lines.

### Adding a new column (CorE) to the data frame that designates lines as check (=0) or experimental (=1) and moving it toward the front of the data frame
Phen <- Phen %>%
  mutate(CorE = case_when(Line == "Hockett" | Line == "Lavina" | Line == "Merit57" | Line == "Odyssey"     # names of the replicated check varieties
                          ~ 0, TRUE ~ 1)) %>% 
  relocate("CorE", .after = "Line")


### Adding another column (CNames) to the data frame that gives the line name if it is a check, and "exp" if it is an experimental line and moving it toward the front of the data frame
Phen <- Phen %>% 
  mutate(CNames = case_when(CorE == 0 ~ Line, TRUE ~ "exp")) %>% 
  relocate("CNames", .after = "CorE")


### Designating explanatory variables for the models as factors
CNames <- as.factor(Phen$CNames)
Block <- as.factor(Phen$Block)
Line <- as.factor(Phen$Line)
CorE <- as.factor(Phen$CorE)



### BLUP Modeling

# Model each phenotype with a fixed effect for check varieties and random effects for block and genotype (represented by the interaction of line and the indicator for the line being a check or experimental)
# Save the BLUPs (genotype coefficients) from the phenotypes for which the model fit, and making note of which traits the model did not fit (due to lack of block effect and or lack of heritability).


# Extract a character vector of all of the trait column names
trait_cols <- names(Phen)[6:ncol(Phen)]     ##### double check the data frame to see which number col to start with #####

# Initializing a character vector in which to store the names of any traits that produce a message (probably the boundary singular fit message)
message_trait_cols <- character()

# Initializing a character vector in which to store the names of any traits that produce a warning
warning_trait_cols <- character()

# Fitting the model to each trait, storing a list of all these models, and filling the message and warning vectors
mods <- sapply(1:length(trait_cols), function(i){
  tryCatch(
    {
      lmer(Phen[[trait_cols[i]]] ~ CNames + (1|Block) + (1|Line:CorE), REML = TRUE, data = Phen)
    },
    message = function(msg){
      message_trait_cols <<- c(message_trait_cols, trait_cols[i])
      message(paste("Message in trait_col:", trait_cols[i]))
    },
    warning = function(wng){
      warning_trait_cols <<- c(warning_trait_cols, trait_cols[i])
      warning(paste("Warning in trait_col:", trait_cols[i]))
    }
  )
})



## The models that produced messages or warnings are stored as NULL or as characters (stating the warning) in the model list.
# These need to be removed in order to be able to pull the random effect coefficients out of the ones that did work all together.

# Identify by index number in the mods list, which models did not work (are NULL or a character)
invalid_indices <- which(sapply(mods, function(x) is.null(x) || is.character(x)))

# Identify which ones did work
valid_indices <- which(!sapply(mods, function(x) is.null(x) || is.character(x)))

# Subset the model list to only include the ones that worked
valid_mods <- mods[valid_indices]

# Make a character vector of the corresponding trait names for the models that worked
valid_traits <- trait_cols[valid_indices]



## Extract the random effect coefficients for the experimental lines and put them all into a data frame of BLUPs for the different traits.

# Make a vector of all the experimental line names to be able to pull out only the relevant BLUPs
exp_lines <- unique(Phen$Line[!Phen$Line %in% c("Hockett", "Lavina", "Merit57", "Odyssey")])

# Extract BLUPs for all valid models and combine into a list
blups_list <- lapply(valid_mods, function(model) {
  a <- as.data.frame(ranef(model))
  a %>%
    mutate(Line = gsub(":1", "", grp)) %>%
    filter(Line %in% exp_lines) %>%
    select(Line, condval)
})

# Join all BLUPs by Line into a data frame and give the columns their trait names
blups <- blups_list %>%
  reduce(full_join, by = "Line") %>%
  rename_with(~ paste0(valid_traits), starts_with("condval"))



### Heritability Calculations

# Calculating heritability (H2) of each trait by dividing the variance due to genotype by the total variance (each effect plus residuals divided by the number of replicates which is 1)

H2.Phen <- NULL

for(t in 1:length(trait_cols)){
  D1<-droplevels(Phen[!is.na(Phen[,colnames(Phen)==trait_cols[t]]),])
  mod<-lmer(eval(parse(text = paste(trait_cols[t],"~CNames + (1|Block) + (1|Line:CorE)",sep=""))),data = D1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  H2.Phen<-rbind(H2.Phen,data.frame(Trait=trait_cols[t],
                                    H2_LocYr=round(c(Var1[1]/sum(Var1[1], Var1[2], Var1[3]/1))     ##### Change H2 column suffix for LocYr #####
                                                  ,3)      # 3 is for the # of decimal points reported
  ))
}



#### Block Effect Calculations

# Calculating percent variance explained by the random effect of Block
# (like H2 but with Block as numerator- Var1[2] for Block as numerator (can run Var1[#] to check which # is for which factor)) 

PctBlock.Phen <- NULL

for(t in 1:length(trait_cols)){
  D1<-droplevels(Phen[!is.na(Phen[,colnames(Phen)==trait_cols[t]]),])
  mod<-lmer(eval(parse(text = paste(trait_cols[t],"~CNames + (1|Block) + (1|Line:CorE)",sep=""))),data = D1)
  Var1<-as.data.frame(VarCorr(mod))$vcov
  names(Var1)<-as.data.frame(VarCorr(mod))$grp
  PctBlock.Phen<-rbind(PctBlock.Phen,data.frame(Trait=trait_cols[t],
                                                PctBlock_LocYr=round(c(Var1[2]/sum(Var1[1], Var1[2], Var1[3]/1))      ##### Change PctBlock column suffix for LocYr #####
                                                                    ,3)      # the 3 is just for # of decimal points reported
  ))
}



### Making Effects dataframe

# Making a data frame with % var for random Block effect (PctBlock) and % var for random Line effect (H2) for each trait having used the full augmented design model for everything,
# as well as columns to indicate whether the trait produced a message or a warning when running the model for BLUPS.

# Make data frames from the message trait and warning trait lists
messages <- data.frame(
  Trait = message_trait_cols,
  Messages = "msg"
)

warnings <- data.frame(
  Trait = warning_trait_cols,
  Warnings = "wng"
)

##### If either or both of the above bits of code give errors and don't run because there aren't any traits that produced messages or warnings,
# run this code instead.

# messages <- data.frame(
#   Trait = trait_cols,
#   Warnings = "NA"
# )

# warnings <- data.frame(
#   Trait = trait_cols,
#   Warnings = "NA"
# )


# Join all the data frames
dfs <- list(PctBlock.Phen, H2.Phen, messages, warnings)

Effects <- reduce(dfs, left_join, by = "Trait")

Effects_type_LocYr <- Effects               ##### Change name to indicate data type (16S ot ITS or Ag) and LocYr #####



### Making BLUPs dataframe

# Remove any traits that had a 0 for either %var due to Block or %var due to genotype/Line from the data frame of "successful" BLUPs
# (a few of these might get through without a message or warning, but still should not have BLUPs calculated for them).
# If a trait had 0 block effect there is no reason to use BLUPs (which in this case are for the purpose of correcting for environmental variation captured by block) and the unBLUPed data should be used instead.

# Getting a list of all traits that do NOT have a 0 value for either %var due to block or genotype.
good_traits <- Effects$Trait[Effects$PctBlock != 0 & Effects$H2 != 0]


# Subsetting the dataframe of successful BLUPs to only include these traits.
BLUPs_type_LocYr <- blups %>%            ##### Change name to indicate data type (16S or ITS or Ag) and LocYr #####
select(Line, any_of(good_traits))



##### Now repeat the above modeling steps (starting at SET DATASET AND DO MODELING) for the different datasets #####



### Saving

save(Effects_16S_BZ21, Effects_16S_BZ22, Effects_16S_MO21, Effects_16S_MO22, Effects_16S_SD21, Effects_16S_SD22, Effects_16S_HI21,
     Effects_ITS_BZ22, Effects_ITS_MO21, Effects_ITS_MO22, Effects_ITS_SD21, Effects_ITS_SD22, Effects_ITS_HI21,
     Effects_Ag_BZ21, Effects_Ag_BZ22, Effects_Ag_MO21, Effects_Ag_MO22, Effects_Ag_SD21, Effects_Ag_SD22,
     file = "Effects.RData")

save(BLUPs_16S_BZ21, BLUPs_16S_BZ22, BLUPs_16S_MO21, BLUPs_16S_MO22, BLUPs_16S_SD21, BLUPs_16S_SD22, BLUPs_16S_HI21,
     BLUPs_ITS_BZ22, BLUPs_ITS_MO21, BLUPs_ITS_MO22, BLUPs_ITS_SD21, BLUPs_ITS_SD22, BLUPs_ITS_HI21,
     BLUPs_Ag_BZ21, BLUPs_Ag_BZ22, BLUPs_Ag_MO21, BLUPs_Ag_MO22, BLUPs_Ag_SD21, BLUPs_Ag_SD22,
     file = "BLUPs.RData")









####### Creating phenotype (Y) text files for input to GAPIT for GWAS #######

# Read in unBLUPed trait values and BLUPs, and field trial info

# Field data and raw agronomic trait values
field <- read.csv("RawFieldAgData.csv")

# For unBLUPed microbial trait values:
load("PreparedTaxaDFs.RData")

# For BLUPs:
load("BLUPs.RData")

# read-in genotyping data file for matching up barley lines ("taxa" in GAPIT lingo) to those that have genotyping data
geno <- read.table("BMB_GD.txt", head = TRUE)
taxa <- taxa %>% select("taxa")



##### SET TO A DATASET #####

# For trait BLUPs of microbial or agronomic datasets
BLUPs <- BLUPs_dataset     ##### change dataset name (eg. 16S_BZ21) ####


#####
# For unBLUPed trait values of MICROBIAL datasets and corresponding LocYr info (SKIP THIS FOR AGRONOMIC AND USE BELOW)
Raws <- PreppedTaxa_dataset     ##### change dataset name (eg. 16S_BZ21) ####

# Run names fn on Raws to get col names to match w/ BLUPs
names(Raws) <- make.names(names(Raws), unique = TRUE)

LocYr <- field %>%             ##### change name of LocYr #####
filter(LocYr == "") %>%      ##### enter name of LocYr #####
select(Plot, Line)
#####


#####
# For unBLUPed trait values of AGRONOMIC datasets (SKIP FOR MICROBIAL AND USE ABOVE)

Raws <- field %>% 
  select(-Year, -Location, -Row, -Range, -Plot, -Entry)     # remove unnecessary columns (keeping only LocYr, Block, Line, trait columns)

Raws <- Raws %>% 
  filter(LocYr == "MO22") %>%               ##### set LocYr #####
select(where(~ !all(is.na(.)))) %>%       # remove trait columns that are missing for the given LocYr
  rename_with(                              # add _LocYr suffixes to trait column names to help with file prep later
    ~ paste0(.x, "_LocYr"),                 ##### set LocYr #####
    .cols = -c(LocYr, Block, Line)
  )
#####




# Subset the raw traits dataframe to not include any that appear in the BLUPs dataframe (all others are retained)
# Thus mapping unBLUPed values for traits that did not have a Block effect
Raws <- Raws[ , !names(Raws) %in% names(BLUPs)]


# Setting up the taxa column (as GAPIT refers to Line) in the raw traits dataframe

Raws <- Raws %>% 
  mutate(Plot = sub("_.*", "", Sample)) %>% 
  select(-Sample) %>% 
  relocate(Plot, .before = 1)

LocYr$Plot <- as.character(LocYr$Plot)

Raws <- left_join(LocYr, Raws, by = "Plot")

Raws <- Raws %>% 
  select(-Plot) %>% 
  filter(!Line %in% c("Hockett", "Lavina", "Odyssey", "Merit57")) %>% 
  rename(taxa = Line)


# Setting up Taxa column for BLUPs
BLUPs <- BLUPs %>% rename(taxa = Line)


# Joining Raws and BLUPs
Traits <- left_join(Raws, BLUPs, by = "taxa")



# Save Y (GAPIT lingo for phenotype) text file.

write.table(Traits, "BMB_Y_dataset.txt", sep = "\t", row.names = FALSE)     #### change dataset name eg. 16S_BZ21) ####


### Repeat for other datasets




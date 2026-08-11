
### GWAS ###

library(GAPIT)


# In phenotype (BMB_Y) files, BZ = swMT and MO = cntrMT. (BMB is for barley microbiome project)


bmbY <- read.table("BMB_Y_.txt", head = TRUE)     # for phenotypes (Y), fill in .txt file name for either 16S, ITS, or Ag data and location-year
bmbGD <- read.table("BMB_GD.txt", head = TRUE)
bmbGM <- read.table("BMB_GM.txt" , head = TRUE)


bmbGAPIT <- GAPIT(
  Y = bmbY,
  GD = bmbGD,
  GM = bmbGM,
  PCA.total = 2,
  model = "FarmCPU",
  SNP.MAF = 0.05,
  cutOff = 0.5,                  # a cutoff of FDR-adjusted p < 0.2 or 0.05 was applied after mapping by sorting results in GAPIT output files
  Major.allele.zero = TRUE
)






### Wilcoxon's Tests of Ag traits at microbe-trait related SNPs from GWAS results ###

library(dplyr)


# open files

mic <- read.csv(".csv")     # a list of snps that had associations in the microbe GWAS after applying FDR-adjusted p < 0.2 threshold

geno <- read.table("BMB_GD.txt", head = TRUE)    # all the genotyping data

ag <- read.table("BMB_Y_Ag.txt", head = TRUE)    # the GAPIT phenotype input file of ag trait data


# create a dataframe (rest) that has genotype data for microbe snps and agronomic trait data

keep <- mic$SNP

geno <- geno %>% select(taxa, all_of(keep))

rest <- geno %>% left_join(ag, by = "taxa")


# define numbers for running tests

i <-  2   # first snp column in rest dataframe

j <-  770   # first ag trait column in rest dataframe



# create an empty matrix to store wilcox test

wilcoxMatrix <- matrix(data = NA, nrow = 768, ncol = 65)     # 768 is the number of microbe-associated snps, 65 is the number of ag traits

# loop through and run wilcox tests

for (i in 2:769){        # 2:769 is the range of snp columns
  
  print(colnames(rest)[i])
  
  for (j in 770:834){          # 770:834 is the range of agronomic trait columns
    
    temp <- as.data.frame(rest[,c(i,j)])
    
    temp[,1] <- as.factor(temp[,1])
    
    temp <- temp[temp[,1] %in% c(0,2),]
    
    temp2 <- wilcox.test(as.numeric(temp[,2]) ~ temp[,1])
    
    wilcoxMatrix[i-1,j-769] <- temp2$p.value            # 769 is the last snp column
    
  }
  
}



# add SNP and traits names to the matrix

rownames(wilcoxMatrix) <- colnames(rest)[2:769]

colnames(wilcoxMatrix) <- colnames(rest)[770:834]

summary(wilcoxMatrix)



# find the significant tests, based on Bonferonni correction

indexPadj <- which(wilcoxMatrix < 0.05/(768*65), arr.ind=TRUE)     # 768 is the number of snps tested, 65 is the number of traits tested (768*65 is the total number of tests)

d <- dim(indexPadj)  # the number of significant snp/trait associations by Wilcoxon's test over the 0.05 Bonferroni threshold



# pull the SNPs with significant p-values

sigSnp <- list()

for (k in 1:nrow(indexPadj)) {
  
  sigSnp[[k]] <- colnames(rest)[indexPadj[[k,1]]]
  
}



# pull the traits with significant p-values

sigTrait <- list()

for (k in 1:nrow(indexPadj)) {
  
  sigTrait[[k]] <- colnames(rest)[indexPadj[[k,2]] + 769]
  
}



# merge the SNPs and traits that are significant together

results <- matrix(data = NA, nrow = d, ncol = 2)

results[,1] <- unlist(sigSnp)

results[,2] <- unlist(sigTrait)

colnames(results) <- c("SNP", "agTrait")

results



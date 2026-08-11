setwd("~/Documents/PhD/Research Projects/Microbiome/BMB_R_Final")

rm(list = ls()) # Clears data from global environment


library(phyloseq)
library(microeco)
library(file2meco)
library(ggplot2)
library(ggpubr)
library(microViz)
library(RColorBrewer)
theme_set(theme_classic())
library(tidyverse)
library(readxl)
library(scales)


# Had to install some subpackages separately:

# library(BiocManager)

# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

# BiocManager::install("rhdf5")

# BiocManager::install("ggh4x")

# BiocManager::install("randomForest")

# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("david-barnett/microViz")



# saveRDS(my_plot, file = "my_plot.rds")

# save.image(file = "Fig1.RData")




# Based on analysis from  Riddley M, Hepp S, Hardeep F, Nayak A, Liu M, Xing X, Zhang H, Liao J (2025) 
# Differential roles of deterministic and stochastic processes in structuring soil bacterial ecotypes across 
# terrestrial ecosystems. Nat Commun 16: 2337. doi: 10.1038/s41467-025-57526-x

# Use updated phyloseq objects that I updated on 12.19.25 (Combined_Analysis_2025/Data Prep.R) to standardize names between 16S and ITS, add soil parameters, 
# and PCA line clusters. These phyloseq objects were from Erik's psField.RDS and ps.its.RDS


ps_16S <- readRDS("ps.filt.RDS")
ps_ITS <- readRDS("ps.q20.trunc200.RDS")

# Remove bulk samples
ps_16S <- subset_samples(ps_16S, Line != "Bulk")
ps_ITS <- subset_samples(ps_ITS, Line != "Bulk")

# Remove Bozeman 20221 fungal data because most samples failed
ps_ITS <- subset_samples(ps_ITS, !(Location == "BZ" & Year == "2021"))

m16S <- phyloseq2meco(ps_16S)
mITS <- phyloseq2meco(ps_ITS)

m16S$cal_abund()
mITS$cal_abund()

m16S$cal_betadiv()
mITS$cal_betadiv()

m16S$tidy_dataset()
mITS$tidy_dataset()

m16S$sample_table$Location <- factor(m16S$sample_table$Location, 
                                     levels = c("BZ","MO","SD","HI"),
                                     labels = c("swMT","cntrMT","SD","HI"))

mITS$sample_table$Location <- factor(mITS$sample_table$Location, 
                                     levels = c("BZ","MO","SD","HI"),
                                     labels = c("swMT","cntrMT","SD","HI"))
m16S$sample_table
saveRDS(m16S, "m16S.RDS")
saveRDS(mITS, "mITS.RDS")

m16S <- readRDS("m16S.RDS")
mITS <- readRDS("mITS.RDS")

# Bar plots
m16Sb <- trans_abund$new(dataset = m16S, taxrank = "Phylum", ntaxa = 8)
m16Sb$data_abund$Location <- factor(m16Sb$data_abund$Location, 
                                    levels = c("swMT","cntrMT","SD","HI"),
                                    labels = c("swMT","cntrMT","SD","HI"))

mITSb <- trans_abund$new(dataset = mITS, taxrank = "Phylum", ntaxa = 8)
mITSb$data_abund$Location <- factor(mITSb$data_abund$Location, 
                                    levels = c("swMT","cntrMT","SD","HI"),
                                    labels = c("swMT","cntrMT","SD","HI"))

bp1 <- m16Sb$plot_bar(others_color = "grey70", 
                      facet = c("Location","Year"), xtext_keep = FALSE, 
                      legend_text_italic = FALSE, barwidth = 1) + 
  theme(text = element_text(size=12), axis.title.y = element_text(size = 10),
        strip.text = element_text(size=12),
        panel.spacing = unit(0.1,"cm"))

bp2 <- mITSb$plot_bar(others_color = "grey70", 
                      facet = c("Location","Year"), xtext_keep = FALSE, 
                      legend_text_italic = FALSE, barwidth = 1) + 
  theme(text = element_text(size=12), axis.title.y = element_text(size = 10),
        strip.text = element_text(size=12),
        panel.spacing = unit(0.1,"cm"))
bp1
bp2

# Random Forest
m16Srf <- trans_diff$new(dataset = m16S, method = "rf", group = "Location", taxa_level = "Genus")
mITSrf <- trans_diff$new(dataset = mITS, method = "rf", group = "Location", taxa_level = "Genus")

# plot the MeanDecreaseGini bar
# group_order is designed to sort the groups
g1 <- m16Srf$plot_diff_bar(use_number = 1:20, group_order = c("swMT", "cntrMT", "SD", "HI")) +
  scale_fill_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.margin = margin(0,0.5,0,0, unit = "cm"))
# plot the abundance using same taxa in g1
g2 <- m16Srf$plot_diff_abund(group_order = c("swMT", "cntrMT", "SD", "HI"), select_taxa = m16Srf$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = TRUE, 
                             errorbar_addpoint = FALSE, errorbar_color_black = TRUE) + 
  scale_fill_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.margin = margin(0,1,0,0, unit = "cm")) 
# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g1 <- g1 + theme(legend.position = "none") + scale_x_discrete(labels = function(x) sub("^.*__", "", x))
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), axis.line.x  = element_line(color = "gray40", linewidth = 0.35),
                 axis.line.y  = element_blank(), legend.position = "none")
p1 <- ggarrange(g1, g2, ncol = 2, nrow = 1, widths = c(1, 1)) + theme(axis.title.y = element_text(size = 10))
p1

g3 <- mITSrf$plot_diff_bar(use_number = 1:20, group_order = c("swMT", "cntrMT", "SD", "HI"), ) +
  scale_fill_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.margin = margin(0,1,0,0, unit = "cm"))

# plot the abundance using same taxa in g1
g4 <- mITSrf$plot_diff_abund(group_order = c("swMT", "cntrMT", "SD", "HI"), select_taxa = mITSrf$plot_diff_bar_taxa, plot_type = "barerrorbar", add_sig = TRUE, 
                             errorbar_addpoint = FALSE, errorbar_color_black = TRUE, ytitle_size = 12) + 
  scale_fill_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title.x = element_text(size = 10), axis.text.y = element_text(size = 10), plot.margin = margin(0,1,0,0, unit = "cm"))

# now the y axis in g1 and g2 is same, so we can merge them
# remove g1 legend; remove g2 y axis text and ticks
g3 <- g3 + theme(legend.position = "none") + scale_x_discrete(labels = function(x) sub("^.*__", "", x)) 
g4 <- g4 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.border = element_blank(), axis.line.x  = element_line(color = "gray40", linewidth = 0.35),
                 axis.line.y  = element_blank(), legend.position = "none")
p2 <- ggarrange(g3, g4, ncol = 2, nrow = 1, widths = c(1, 1)) 
p2 

# MicroViz for PCA of CLR transformed data
library(microViz)

v16S <- tax_fix(meco2phyloseq(m16S))
vITS <- tax_fix(meco2phyloseq(mITS))

op1 <- v16S %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_transform(rank = "Genus", trans = "clr", zero_replace = "halfmin") %>%
  dist_calc(dist = "euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 2, color = "Location", auto_caption = NA) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title = element_text(size = 10)) +
  theme_classic(12) + coord_fixed(0.7) +  stat_ellipse(aes(color = Location))
op1

op2 <- vITS %>%
  tax_filter(min_prevalence = 2.5 / 100, verbose = FALSE) %>%
  tax_transform(rank = "Genus", trans = "clr", zero_replace = "halfmin") %>%
  dist_calc(dist = "euclidean") %>%
  ord_calc(method = "PCoA") %>%
  ord_plot(alpha = 0.6, size = 3, color = "Location", auto_caption = NA) +
  scale_color_manual(values = c(swMT ="skyblue", HI = "goldenrod1", cntrMT = "darkorange", SD = "seagreen3")) +
  theme(axis.title = element_text(size = 10)) +
  theme_classic(12) + coord_fixed(0.7) +  stat_ellipse(aes(color = Location))
op2




# List of taxa from RF
library(dplyr)
library(tidyr)
library(stringr)

rd_16S <- m16Srf$res_diff
ra_16S <- m16Srf$res_abund

ra_wide_16S <- ra_16S %>%
  pivot_wider(
    id_cols = Taxa,
    names_from = Location,
    values_from = c(Mean, SD, SE),
    names_sep = "_"
  )

merged_16S <- rd_16S %>%
  left_join(ra_wide_16S, by = "Taxa")

# Split taxa rankings
merged_16S <- merged_16S %>%
  separate(
    Taxa,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "\\|",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    across(
      Kingdom:Species,
      ~ str_remove(., "^[a-z]__")
    )
  ) %>%
  select(-Taxa)

rd_ITS <- mITSrf$res_diff
ra_ITS <- mITSrf$res_abund

ra_wide_ITS <- ra_ITS %>%
  pivot_wider(
    id_cols = Taxa,
    names_from = Location,
    values_from = c(Mean, SD, SE),
    names_sep = "_"
  )

merged_ITS <- rd_ITS %>%
  left_join(ra_wide_ITS, by = "Taxa")

# Split taxa rankings
merged_ITS <- merged_ITS %>%
  separate(
    Taxa,
    into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = "\\|",
    fill = "right",
    remove = FALSE
  ) %>%
  mutate(
    across(
      Kingdom:Species,
      ~ str_remove(., "^[a-z]__")
    )
  ) %>%
  select(-Taxa)

# Filter to remove taxa with MeanDecreaseGini < 0.05
merged_16S <- subset(merged_16S, MeanDecreaseGini > 0.5)
merged_ITS <- subset(merged_ITS, MeanDecreaseGini > 0.5)

write.csv(merged_16S, "rf_taxa_16S.csv")
write.csv(merged_ITS, "rf_taxa_ITS.csv")

imp_16S <- m16Srf$res_diff
imp_16S <- imp_16S[order(imp_16S$MeanDecreaseGini, decreasing = TRUE), ]

# Plot MDG and use a visual cutoff based on plot "elbow" as suggested in {Roguet, 2018 #2999}
imp_16S$Rank <- seq_len(nrow(imp_16S))  # create a rank column
ggplot(imp_16S, aes(x = Rank, y = MeanDecreaseGini)) +
  geom_line(color = "steelblue") + 
  geom_point(color = "steelblue") +
  labs(x = "Feature Rank (by importance)", y = "Mean Decrease Gini",
       title = "Scree Plot of Random Forest Feature Importances") +
  theme_minimal()


imp_ITS <- mITSrf$res_diff
imp_ITS <- imp_ITS[order(imp_ITS$MeanDecreaseGini, decreasing = TRUE), ]

imp_ITS$Rank <- seq_len(nrow(imp_ITS))  # create a rank column
ggplot(imp_ITS, aes(x = Rank, y = MeanDecreaseGini)) +
  geom_line(color = "steelblue") + 
  geom_point(color = "steelblue") +
  labs(x = "Feature Rank (by importance)", y = "Mean Decrease Gini",
       title = "Scree Plot of Random Forest Feature Importances") +
  theme_minimal()

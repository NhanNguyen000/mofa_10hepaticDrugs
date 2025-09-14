# set up enviroment
rm(list = ls())

library(MOFA2)
library(tidyverse)
library("viridis")

source("code/00_plotSettings.R")
# load data ====================================================================
#nameDat <- "normRNA_allEnsembleIds" # all measured Ensemble IDs
nameDat <- "normRNA_geneSymbols" # only measured Ensemble IDs converted into gene symbol

model <- load_model(paste0("processedDat/mofaModel_", nameDat, ".hdf5"))
# note: 8 factors for "normRNA_geneSymbols" and 10 factors for "normRNA_allEnsembleIds" 
# were found to explain no variance. They were automatically removed for downstream analysis.

# Overview of data and metadata -------------------------------------------------
plot_data_overview(model)
sum(model@dimensions$N) # Nsamples

# metadata
head(model@samples_metadata)

load("processedDat/inputDat.RData")
model@samples_metadata <- model@samples_metadata %>% 
  full_join(expDesignTab, by = c("sample" = "name"))
#View(model@samples_metadata)

rm(dat, expDesignTab)

# Visualisation of combinations of factors ------------------------------------
set.seed(42)
model@samples_metadata$condition3
model@samples_metadata <- model@samples_metadata %>%
  mutate(condition3 = ifelse(condition2 %in% c("CYC", "VPA", "ISO"), 
                             "CYC_VPA_ISO",
                             ifelse(condition2 %in% c("DIC", "RIF", "MTX", "APA"), 
                                    "DIC_RIF_MTX_APA",
                                    as.character(condition2)))) %>%
  mutate(condition3 = factor(condition3, 
                             levels = c("ConUNTR", "ConDMSO", "5FU", "AZA", "PHE", "CYC_VPA_ISO", "DIC_RIF_MTX_APA"))) %>%
  mutate(condition4 = ifelse(condition3 %in% c("ConUNTR", "ConDMSO"),
                             "Control (UNTR, DMSO)",
                             ifelse(condition3 %in% c("CYC_VPA_ISO", "DIC_RIF_MTX_APA"),
                                    as.character(condition3), 
                                    "other drugs (5FU, AZA, PHE)")))
# subset -------------
sample_5FU <- (model@samples_metadata %>% filter(drug == "5FU"))$sample
model_5FU <- subset_samples(model, sample_5FU)

plot_factor(model_5FU, factors=1, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=2, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=3, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=4, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=5, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=6, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=7, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))



plot_factor(model_5FU, factors=1, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=2, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=3, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=4, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=5, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=6, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_5FU, factors=7, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


# test set ------------- factor 5, seems to show the time
sample_subset <- (model@samples_metadata %>% 
                    filter(drug %in% c("DIC", "RIF", "MTX", "APA")))$sample
model_subset <- subset_samples(model, sample_subset)

plot_factor(model_subset, factors=1, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=2, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=3, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=4, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=5, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=6, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=7, group_by = "time") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))



plot_factor(model_subset, factors=1, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=2, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=3, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=4, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=5, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=6, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model_subset, factors=7, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


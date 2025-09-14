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

View(model@samples_metadata)

rm(dat, expDesignTab)

# Variance decomposition ------------------------------------------------
model@cache$variance_explained$r2_total[[1]] # Total variance explained per view and group
# RNAseq 69.72612 
model@cache$variance_explained$r2_per_factor[[1]] # Variance explained for every factor in per view and group

mofa_variance_factor <- plot_variance_explained(model, x="view", y="factor")
mofa_variance_factor

plot_variance_explained(model, x="group", y="factor", plot_total = T)[[1]]
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]
plot_variance_explained(model, x="view", y="group")

# save plots ------------------------------------------------------------
save_plot_png(mofa_variance_factor, 
              filename = paste0("output/03_", nameDat, "_MofaVarianceFactor.png"), 
              width = 60) 

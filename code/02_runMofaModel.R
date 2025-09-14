rm(list = ls())

library(tidyverse)
library(MOFA2)
# load data ====================================================================
load("processedDat/inputDat.RData")

#nameDat <- "normRNA_allEnsembleIds" # all measured Ensemble IDs
nameDat <- "normRNA_geneSymbols" # only measured Ensemble IDs converted into gene symbol

set.seed(42)
# prepare data for MOFA model -------------------------------------------------
inputDat <- list()
inputDat$RNAseq <- data.matrix(dat[[nameDat]])
lapply(inputDat, dim)

MOFAobject <- create_mofa(inputDat)
print(MOFAobject)

# prepare mofa model options ----------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

# run MOFA models --------------------------------------------------------------
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
print(MOFAobject)

#Train the MOFA model
outfile = file.path(getwd(), paste0("processedDat/mofaModel_", nameDat, ".hdf5"))
MOFAobject.trained <- run_mofa(MOFAobject, outfile, use_basilisk = TRUE)



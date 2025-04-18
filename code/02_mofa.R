rm(list = ls())

library(tidyverse)
library(MOFA2)
# load data =======================================================================
load("processedDat/inputDat.RData")
dim(norm_rna) #  60664   474
dim(norm_rna_v2) #  36244   474

# run mofa for all measure Ensemble ID gene -------------------------------------
## prepare data ------------------------------------------------------------
inputDat <- list()
inputDat$RNAseq <- data.matrix(norm_rna)
lapply(inputDat, dim)

MOFAobject <- create_mofa_from_matrix(inputDat)
print(MOFAobject)

## prepare moda options ----------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

## run MOFA --------------------------------------------------------------------
model_allRNAs <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

rm(inputDat, MOFAobject, data_opts, model_opts, train_opts)

# run MOFA for genes with converted symbols ------------------------------------
# prepare data --------------------------------------------------------------
inputDat <- list()
inputDat$RNAseq <- data.matrix(norm_rna_v2)
lapply(inputDat, dim)

MOFAobject <- create_mofa_from_matrix(inputDat)
print(MOFAobject)

## prepare moda options -----------------------------------------------------
data_opts <- get_default_data_options(MOFAobject)
head(data_opts)

model_opts <- get_default_model_options(MOFAobject)
head(model_opts)

train_opts <- get_default_training_options(MOFAobject)
head(train_opts)

## run MOFA ---------------------------------------------------------------
model_geneSymbols <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)


# save data -----------------------------------------------------
save(model_allRNAs, model_geneSymbols, file = "processedDat/mofa_models.RData")


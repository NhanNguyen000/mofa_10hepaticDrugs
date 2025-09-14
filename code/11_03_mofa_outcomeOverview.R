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

# Non-linear dimensionality reduction (UMAP, t-SNE) -----------------------
set.seed(42)
model <- run_tsne(model)

vars <- colnames(model@samples_metadata)[3:7]
plist_tsne <- list()
for(varName in vars) {
  plist_tsne[[varName]] <- plot_dimred(
    model, 
    method = "TSNE", 
    color_by = varName, 
    label = TRUE, 
    stroke=0.05, 
    dot_size = 5, 
    legend = TRUE)
}
plist_tsne$dose # separate control samples (DMSO and UNTR) vs. the rest
plist_tsne$drug # separate Con, PHE, 5FU, and group some drugs together (DIC, RIF, MTX, APA), AZA in between, (CYC, VPA, ISO)
plist_tsne$time  # not able to separate time points
plist_tsne$condition
plist_tsne$condition2 + plot_theme()

# UMAP plot
model <- run_umap(model)

plist_umap <- list()
for(varName in vars) {
  plist_umap[[varName]] <- plot_dimred(
    model, 
    method = "UMAP", 
    color_by = varName, 
    label = TRUE, 
    stroke=0.05, 
    dot_size = 2, 
    legend = TRUE) #+ scale_fill_manual(values=colors)
}
plist_umap$dose # separate control samples (DMSO and UNTR) vs. the rest - not good visualization as TSNE
plist_umap$drug # separate Con, PHE, 5FU, and group some drugs together (DIC, RIF, MTX, APA), AZA in between but more cloase to (CYC, VPA, ISO) - similar to TSNE
plist_umap$time # not able to separate time points
plist_umap$condition
plist_umap$condition2

# We can try to add some interpretatibility on the TSNE / UMAP by visualising the contribution of each Factor on the different groups of cells.

for (i in paste0("Factor",1:7)) {
  p <- plot_dimred(model, 
                   method = "TSNE", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 2
  )
  print(p)
}


for (i in paste0("Factor",1:7)) {
  p <- plot_dimred(model, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 2
  )
  print(p)
}

# Visualisation of combinations of factors ------------------------------------
plist_factors <- list()
for(varName in vars) {
  plist_factors[[varName]] <- plot_factors(model, factors = 1:7, color_by = varName)
}

plist_factors$dose # have some differences
plist_factors$drug # have some differences
plist_factors$time # not much differences (?)

# between 2 factors
plot_factors(model, 
             factors = 1:2,
             color_by = "drug")

plot_factors(model, 
             factors = 3:4,
             color_by = "dose")

p <- plot_factors(model, 
                  factors = c(2,5), 
                  color_by = "time",
                  shape_by = "dose",
                  dot_size = 2.5,
                  show_missing = T)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

# Correlation between factors ---------------------------------------------------------
plot_factor_cor(model)

# which Factors have a clear association with any of the covariates? 
corFactors_covariates <- correlate_factors_with_covariates(
  model, covariates = c("time","drug","dose"), 
  plot="log_pval")
corFactors_covariates

# Save plots ------------------------------------------------------------
png(file= paste0("output/03_", nameDat, "_mofa_variance_factors.png"), width = 320)
mofa_variance_factor + theme(text = element_text(size = 18))
dev.off()

for(varName in vars) {
  png(file= paste0("output/03_", nameDat, "_mofa_tsne_", varName, ".png"), width = 720)
  print(plist_tsne[[varName]] + theme(text = element_text(size = 18)))
  dev.off()
  
  png(file= paste0("output/03_", nameDat, "_mofa_combiFactors_", varName, ".png"), width = 576)
  print(plist_factors[[varName]])
  dev.off()
}

png(file= paste0("output/03_", nameDat, "_mofa_corFactors.png"))
plot_factor_cor(model)
dev.off()

png(file= paste0("output/03_", nameDat, "_mofa_corFactors_covariates.png"))
corFactors_covariates
dev.off()

# Extracting data for downstream analysis ---------------------------------------
# Extract "factors" - a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(model, factors = "all")
lapply(factors,dim)
View(factors$group1)

# Extract "weights" - a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(model, views = "all", factors = "all")
lapply(weights,dim)
View(weights$RNAseq)

# Extract "data" - a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
data <- get_data(model)
lapply(data, function(x) lapply(x, dim))[[1]]
View(data$RNAseq$group1)

# Extract the data in long data.frame format:
factors <- get_factors(model, as.data.frame = T)
head(factors, n=3)

# For convenience, we can extract the data in long data.frame format:
factors <- get_factors(model, as.data.frame = T)
head(factors, n=3)

weights <- get_weights(model, as.data.frame = T)
head(weights, n=3)

data <- get_data(model, as.data.frame = T)
head(data, n=3)
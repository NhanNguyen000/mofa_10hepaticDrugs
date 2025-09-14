# set up enviroment
rm(list = ls())

library(MOFA2)
library(tidyverse)
library("viridis")
library(PCAtools)

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

#View(model@samples_metadata)

#rm(dat, expDesignTab)

# PCA ----------------
factors <- get_factors(model, factors = "all")
inputDat <- factors$group1 %>% 
  t() %>% as.data.frame() %>% 
  dplyr::select(rownames(expDesignTab))
all(colnames(inputDat) == rownames(expDesignTab)) # TRUE, so can run the next code

pca_output <- pca(inputDat, metadata = expDesignTab, removeVar = 0.1)
pc_for75 <- which(cumsum(pca_output$variance) > 75)[1]
pc_for75
biplot_basic <- biplot(pca_output)
biplot_basic

biplot(pca_output, showLoadings = TRUE, lab = NULL)

condition2_palette <- c(
  "#3C5488", "#E64B35",  # 2 controls
  "#66A61E", "#00A087", "#F39B7F", "#A65628", "#F4A582", 
  "#80B1D3", "#DF65B0", "#7570B3", "#E6AB02",  "#4DBBD5" # 10 drugs
)

biplot_group <- biplot(
  pca_output,
  #colby = 'dose', colLegendTitle = 'Dose',
  colby = 'condition2', colLegendTitle = 'Condition', colkey = condition2_palette, 
  pointSize = 4, 
  hline = 0, vline = 0,
  # encircle = TRUE, encircleFill = FALSE, encircleAlpha = 0.7, encircleLineSize = 3,
  gridlines.major = FALSE, gridlines.minor = FALSE,
  legendPosition = 'right', legendLabSize = 7, legendIconSize = 6.0,
  #drawConnectors = FALSE,
  title = paste0('PCA bi-plot, ', nameDat ),
  subtitle = paste0("PC1 versus PC2, ", pc_for75, " PCs ≈ 75%")) #+ plot_theme()
biplot_group


# Non-linear dimensionality reduction (UMAP, t-SNE) -----------------------
## UMAP plot ---------------------------------------------------------------------
model <- run_umap(model)

plist_umap <- list()
for(varName in vars) {
  plist_umap[[varName]] <- plot_dimred(
    model, 
    method = "UMAP", 
    color_by = varName, 
    label = TRUE, 
    stroke=0.05, 
    dot_size = 5, 
    legend = TRUE) #+ scale_fill_manual(values=colors)
}
plist_umap$dose # separate control samples (DMSO and UNTR) vs. the rest - not good visualization as TSNE
plist_umap$drug # separate Con, PHE, 5FU, and group some drugs together (DIC, RIF, MTX, APA), AZA in between but more cloase to (CYC, VPA, ISO) - similar to TSNE
plist_umap$time # not able to separate time points
plist_umap$condition
plist_umap$condition2 + plot_theme()

## tSNE plot ---------------------------------------------------------------------
set.seed(42) # But t-SNE can give different results on each run, even with the same input data and parameters — unless you fix the random seed.
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


# # We can try to add some interpretatibility on the TSNE / UMAP by visualising the contribution of each Factor on the different groups of cells.
# 
# plist_factor_tSNE <- list()
# for (i in paste0("Factor",1:7)) {
#   plist_factor_tSNE[[i]] <- plot_dimred(model, 
#                    method = "TSNE", 
#                    color_by = i, 
#                    stroke = 0.05, 
#                    dot_size = 2)
# }
# plist_factor_tSNE$Factor1
# plist_factor_tSNE$Factor2
# plist_factor_tSNE$Factor3
# plist_factor_tSNE$Factor4
# plist_factor_tSNE$Factor5
# plist_factor_tSNE$Factor6
# 
# for (i in paste0("Factor",1:7)) {
#   p <- plot_dimred(model, 
#                    method = "UMAP", 
#                    color_by = i, 
#                    stroke = 0.05, 
#                    dot_size = 2
#   )
#   print(p)
# }

# t-SNE plot, with color and legend correponding to PCA plot ===============
set.seed(42)
model <- run_tsne(model)

condition2_palette <- c(
  "#3C5488", "#E64B35",  # 2 controls
  "#66A61E", "#00A087", "#F39B7F", "#A65628", "#F4A582", 
  "#80B1D3", "#DF65B0", "#7570B3", "#E6AB02",  "#4DBBD5" # 10 drugs
)

plot_tSNE_condition2 <- plot_dimred(
  model, method = "TSNE", 
  color_by = "condition2", 
  color_name = "Condition",
  label = TRUE, stroke=0.05,
  dot_size = 5, legend = TRUE) + 
  scale_fill_manual(values = condition2_palette) + plot_theme()
plot_tSNE_condition2

model@samples_metadata$condition3

model@samples_metadata <- model@samples_metadata %>%
  mutate(condition3 = ifelse(condition2 %in% c("CYC", "VPA", "ISO"), 
                             "CYC_VPA_ISO",
                             ifelse(condition2 %in% c("DIC", "RIF", "MTX", "APA"), 
                                    "DIC_RIF_MTX_APA",
                                    as.character(condition2))))

condition3_palette <- c(
  "ConDMSO"= "#3C5488", "ConUNTR" = "#E64B35",  # 2 controls
  "5FU" = "#66A61E", "AZA" = "#F39B7F",
  "CYC_VPA_ISO" = "#80B1D3", 
  "PHE" = "#7570B3", 
  "DIC_RIF_MTX_APA" = "#DF65B0" #"#66A61E" #
)

plot_tSNE_condition3 <- plot_dimred(
  model, method = "TSNE", 
  color_by = "condition3", 
  color_name = "Condition",
  label = TRUE, stroke=0.05,
  dot_size = 5, legend = TRUE) + 
  scale_fill_manual(values = condition3_palette) + plot_theme()
plot_tSNE_condition3


# save plots ------------------------------------------------------------
save_plot_png(biplot_group, 
              filename = paste0("output/03_", nameDat, "_PCAafterMOFA.png"), 
              width = 200, height = 130)

save_plot_png(plot_tSNE_condition2, 
              filename = paste0("output/03_", nameDat, "_tSNEcondition2.png"), 
              width = 200, height = 130)

save_plot_png(plot_tSNE_condition3, 
              filename = paste0("output/03_", nameDat, "_tSNEcondition3.png"), 
              width = 200, height = 130)


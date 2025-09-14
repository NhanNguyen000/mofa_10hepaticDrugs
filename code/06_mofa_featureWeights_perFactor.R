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

# metadata -------------------------------------------------
head(model@samples_metadata)

load("processedDat/inputDat.RData")
model@samples_metadata <- model@samples_metadata %>% 
  full_join(expDesignTab, by = c("sample" = "name"))
#View(model@samples_metadata)
head(model@samples_metadata)

rm(dat, expDesignTab)

# Visualisation of feature weights ------------------------------------
plot_weights(model, view = "RNAseq",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_top_weights(model, view = "RNAseq",
                 factor = 1,
                 nfeatures = 30
)


# more help on plot weight
plot_weights_fn <- function(mofa, factor=1, view=1, nfeatures=10) {
  p1 <- plot_weights(mofa, 
                     factors = factor, 
                     view = view,
                     nfeatures = nfeatures,
                     text_size = 4
  )
  
  p2 <- plot_top_weights(mofa, 
                         factors = factor, 
                         view = view,
                         nfeatures = nfeatures
  )
  
  p <- cowplot::plot_grid(plotlist=list(p1,p2), nrow=1)
  return(p)
}
plot_weights_fn(model, factor=1, view="RNAseq", nfeatures=20)

# Visualisation of patterns in the input data per factor -------------------------------

# Heatmap of observations: Top features are selected by its weight in the selected factor. 
plot_data_heatmap(model,
                  view = "RNAseq",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  #denoise = TRUE,
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE, #TRUE
                  #scale = "row"
)

plot_data_heatmap(model,
                  view = "RNAseq",
                  factor = 2, 
                  features = 50,
                  show_rownames = F, show_colnames = F, 
                  cluster_rows = T, cluster_cols = F,
                  annotation_samples = "drug",
                  denoise = TRUE
)

category.colors <- c(
  "The" = "#66C2A5", 
  "Tox" = "#8DA0CB",
  "Control" = "#E78AC3" #,
  #"?" = "#FC8D62"
)

plot_data_heatmap(model,
                  view = "RNAseq",
                  factor = 2, 
                  features = 20,
                  denoise = FALSE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  annotation_samples = "dose",  
                  annotation_colors = list("dose"=category.colors), 
                  annotation_legend = T,
                  scale = "row"
)

plot_data_heatmap(model,
                  view = "RNAseq",
                  factor = 4, 
                  features = 20,
                  denoise = FALSE,
                  cluster_rows = T, cluster_cols = F,
                  show_colnames = F, show_rownames = T,
                  annotation_samples = "dose",  
                  annotation_colors = list("dose"=category.colors), 
                  annotation_legend = T,
                  scale = "row"
)
# Scatter plots of observations vs factor values: add a linear regression estimate to visualise if the relationship between (top) features and factor values is linear.
plot_data_scatter(model,
                  view = "RNAseq",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "dose"
)
plot_data_scatter(model,
                  view = "RNAseq",
                  factor = 1,  
                  features = "ATP1A3",
                  sign = "positive",
                  color_by = "dose"
) + labs(y="RNA expression")

plot_data_scatter(model,
                  view = "RNAseq",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "dose"
) + labs(y="RNA expression")


plot_data_scatter(model,
                  view = "RNAseq",
                  factor = 1,  
                  features = 4,
                  sign =  "negative",
                  color_by = "dose"
) + labs(y="RNA expression")

# Save plots ------------------------------------------------------------
#png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_dose.png"), height = 320) # for Ensemble IDs data
png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_dose.png"), width = 576, height = 320) # for gene symbols data
factors_dose + theme(text = element_text(size = 18))
dev.off()

png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_drug.png"), width = 720, height = 320)
factors_drug + theme(text = element_text(size = 18))
dev.off()

png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_drug.png"), width = 720, height = 320)
factors_time + theme(text = element_text(size = 18))
dev.off()

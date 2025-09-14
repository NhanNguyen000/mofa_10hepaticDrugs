# PCA plot: https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
rm(list = ls())

library(tidyverse)
library(PCAtools)
library(ggplot2)
library(org.Hs.eg.db)

source("code/00_plotSettings.R")
# load data =======================================================================
load("processedDat/inputDat.RData")

#nameDat <- "normRNA_allEnsembleIds" # all measured Ensemble IDs
nameDat <- "normRNA_geneSymbols" # only measured Ensemble IDs converted into gene symbol

inputDat <- dat[[nameDat]]

# calculate PCA ------------------------------------------------------------
expDesignTab <- as.data.frame(expDesignTab)
rownames(expDesignTab) <- expDesignTab$name
all(colnames(inputDat) == rownames(expDesignTab)) # TRUE, so can run the next code

pca_output <- pca(inputDat, metadata = expDesignTab, removeVar = 0.1)
pc_for75 <-which(cumsum(pca_output$variance) > 75)[1]
pc_for75

# plotting --------------------------------------------------------------
## scree plot - determine the optimum number of PCs to retain -------------------------------
screeplot(pca_output, 
          axisLabSize = 18, 
          titleLabSize = 22) # scree plot - the accumulative proportion of explained variation 

# horn <- parallelPCA(inputDat) # Horn’s parallel analysis
# horn_value <- horn$n # horn$n = 16, for both case of the input data
horn_value <- 16

elbow_value <- findElbowPoint(pca_output$variance) # elbow method
elbow_value # 10, for both case of the input data

# In most cases, the identified values (horn vs. elbow) will disagree, 
# because finding the correct number of PCs is a difficult task and there is no correct answer.

# Taking these horn and elbow values, a new scree plot and mark these:
screeplot_output <-
  screeplot(pca_output,
            components = getComponents(pca_output, 1:20),
            vline = c(horn_value, elbow_value)) +
  geom_label(
    aes(x = horn_value + 1, y = 50,
        label = paste0('Horn\'s, ', 
                       round(cumsum(pca_output$variance)[horn_value], 2), "%"), 
        vjust = -1, size = 8)) +
  geom_label(
    aes(x = elbow_value + 1, y = 50,
        label = paste0('Elbow method, ',
                       round(cumsum(pca_output$variance)[elbow_value], 2), "%"), 
        vjust = -1, size = 8)) + 
  plot_theme()

screeplot_output

## biplot ---------------------------------------------------------------------
biplot_basic <- biplot(pca_output)
biplot_basic

biplot(pca_output, showLoadings = TRUE, lab = NULL)

biplot_group <- biplot(pca_output,
       #colby = 'dose',
       #colLegendTitle = 'Dose',
       colby = 'dose',
       colLegendTitle = 'Dose',
       hline = 0, vline = 0,
      # encircle = TRUE, encircleFill = FALSE, 
       encircleAlpha = 0.7, encircleLineSize = 3,
       gridlines.major = FALSE, gridlines.minor = FALSE,
       legendPosition = 'right', legendLabSize = 14, legendIconSize = 6.0,
      # drawConnectors = FALSE,
       title = 'PCA bi-plot',
       subtitle = paste0("PC1 versus PC2, ", pc_for75, " PCs ≈ 75%")) + 
  plot_theme()
biplot_group

## other plots ---------------------------------------------------------------------
pairsplot(pca_output)

plotloadings(pca_output, labSize = 2.5)

eigencorplot_output <-eigencorplot(
  pca_output,
  metavars = c("drug", "dose", "condition", "time"))

eigencorplot_output
# save plots ------------------------------------------------------------
png(file= paste0("output/01_", nameDat, "_screePlot.png"), width = 720)
screeplot_output
dev.off()

png(file= paste0("output/01_", nameDat, "_biplotBasic.png"), width = 720)
biplot_basic
dev.off()

png(file= paste0("output/01_", nameDat, "_biplotGroup.png"), width = 720)
biplot_group
dev.off()

png(file= paste0("output/01_", nameDat, "_eigencorplot.png"), width = 720)
eigencorplot_output
dev.off()

# how non-symboled genes contributed to the PCA loading ----------------------
symbols <- mapIds(org.Hs.eg.db, keys = rownames(dat$normRNA_allEnsembleIds),
                  column = c('SYMBOL'), keytype = 'ENSEMBL') %>% 
  as.data.frame() %>% rename("." = "symbol") %>% rownames_to_column("ensembl")

pca_geneSymbols_temp <- pca_output$loadings %>% 
  rownames_to_column("ensembl") %>% full_join(symbols)

pca_geneSymbols <- pca_geneSymbols_temp %>% filter(is.na(symbol)) %>% 
  as.data.frame() %>% dplyr::select(-symbol) %>% column_to_rownames("ensembl")

dim(dat$normRNA_allEnsembleIds) #[1] 60664   474
dim(pca_geneSymbols) #[1] 23956   474
sum(pca_geneSymbols, na.rm = TRUE)


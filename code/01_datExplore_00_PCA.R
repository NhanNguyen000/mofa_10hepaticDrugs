# PCA plot: https://www.bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html
rm(list = ls())

library(tidyverse)
library(PCAtools)
library(ggplot2)
# load data =======================================================================
load("processedDat/inputDat.RData")
dim(norm_rna) #  60664   474
dim(norm_rna_v2) #  36244   474

# PCA for all measured Ensemble IDs --------------------------------------------------------
rownames(expDesignTab) <- expDesignTab$name
all(colnames(norm_rna) == rownames(expDesignTab)) # TRUE, so can run the next code
#p <- pca(norm_rna, metadata = expDesignTab, removeVar = 0.1)

expDesignTab_v2 <- as.data.frame(expDesignTab)
rownames(expDesignTab_v2) <- expDesignTab_v2$name
p <- pca(norm_rna, metadata = expDesignTab_v2, removeVar = 0.1)

which(cumsum(p$variance) > 75)[1] # PC27
## A scree plot --------------------------------------------------------
screeplot(p, axisLabSize = 18, titleLabSize = 22) # A scree plot shows the accumulative proportion of explained variation, 

# # to determine the optimum number of PCs to retain:

# # horn <- parallelPCA(norm_rna) # Horn’s parallel analysis
# # horn$n # [1] 16
horn_value <- 16

# elbow <- findElbowPoint(p$variance) # elbow method
# elbow ## PC10
elbow_value <- 10
# In most cases, the identified values will disagree, because finding the correct number of PCs is a difficult task and there is no correct answer.


# Taking these values, we can produce a new scree plot and mark these:
screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn_value, elbow_value)) +
  geom_label(aes(x = horn_value + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow_value + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

## biplot ---------------------------------------------------------------------
biplot(p)
biplot(p, showLoadings = TRUE, lab = NULL)

biplot(p,
       colby = 'dose',
       colLegendTitle = 'Dose',
       hline = 0, vline = 0,
       encircle = TRUE, encircleFill = FALSE, 
       encircleAlpha = 0.7, encircleLineSize = 3,
       gridlines.major = FALSE, gridlines.minor = FALSE,
       legendPosition = 'right', legendLabSize = 14, legendIconSize = 6.0,
      # drawConnectors = FALSE,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2, 27 PCs ≈ 75%')

## Other ---------------------------------------------------------------------
pairsplot(p)

plotloadings(p, labSize = 3)

eigencorplot(p,
             metavars = c("drug", "dose", "condition", "time"))


# PCA for only measured Ensemble IDs converted into gene symbol ---------------------

## check how non-symboled genes contributed to the  PCA loading ----------------------
k <- p$loadings %>% rownames_to_column("ensembl") %>% full_join(symbols)
k2 <- k %>% filter(is.na(symbol))

## PCA ---------------------------------------------------------------------------
p <- pca(norm_rna_v2, metadata = expDesignTab_v2, removeVar = 0.1)

which(cumsum(p$variance) > 75)[1] # PC18
## A scree plot --------------------------------------------------------
screeplot(p, axisLabSize = 18, titleLabSize = 22) # A scree plot shows the accumulative proportion of explained variation, 

# to determine the optimum number of PCs to retain:

# horn <- parallelPCA(norm_rna) # Horn’s parallel analysis
# horn$n # [1] 16
horn_value <- 16

# elbow <- findElbowPoint(p$variance) # elbow method
# elbow ## PC10
elbow_value <- 10
# In most cases, the identified values will disagree, because finding the correct number of PCs is a difficult task and there is no correct answer.


# Taking these values, we can produce a new scree plot and mark these:
screeplot(p,
          components = getComponents(p, 1:20),
          vline = c(horn_value, elbow_value)) +
  geom_label(aes(x = horn_value + 1, y = 50,
                 label = 'Horn\'s', vjust = -1, size = 8)) +
  geom_label(aes(x = elbow_value + 1, y = 50,
                 label = 'Elbow method', vjust = -1, size = 8))

## biplot ---------------------------------------------------------------------
biplot(p)
biplot(p, showLoadings = TRUE, lab = NULL)

biplot(p,
       colby = 'dose',
       colLegendTitle = 'Dose',
       hline = 0, vline = 0,
       encircle = TRUE, encircleFill = FALSE, 
       encircleAlpha = 0.7, encircleLineSize = 3,
       gridlines.major = FALSE, gridlines.minor = FALSE,
       legendPosition = 'right', legendLabSize = 14, legendIconSize = 6.0,
       # drawConnectors = FALSE,
       title = 'PCA bi-plot',
       subtitle = 'PC1 versus PC2, 18 PCs ≈ 75%')

## Other ---------------------------------------------------------------------
pairsplot(p)

plotloadings(p, labSize = 3)

eigencorplot(p,
             metavars = c("drug", "dose", "condition", "time"))


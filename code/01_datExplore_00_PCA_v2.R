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
all(colnames(inputDat) == rownames(expDesignTab)) # TRUE, so can run the next code

pca_output <- pca(inputDat, metadata = expDesignTab, removeVar = 0.1)
pc_for75 <- which(cumsum(pca_output$variance) > 75)[1]
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
    aes(x = horn_value + 1, y = 70,
        label = paste0('Horn\'s, ', 
                       round(cumsum(pca_output$variance)[horn_value], 2), "%"), 
        vjust = -1, size = 7)) +
  geom_label(
    aes(x = elbow_value - 1, y = 45,
        label = paste0('Elbow method, ',
                       round(cumsum(pca_output$variance)[elbow_value], 2), "%"), 
        vjust = -1, size = 7)) + 
  ggtitle(paste0('PCA scree plot, ', nameDat )) +
  plot_theme()

screeplot_output

## biplot ---------------------------------------------------------------------
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

## other plots ---------------------------------------------------------------------
pairsplot(pca_output)

plotloadings(pca_output, labSize = 2.5)

eigencorplot_output <- eigencorplot(
  pca_output,
  metavars = c("drug", "dose", "condition", "time")) #+ 
 # labs( title = paste0('Correlate PCs to metadata, ', nameDat ))

# save plots ------------------------------------------------------------
save_plot_png(screeplot_output, 
              filename = paste0("output/01_", nameDat, "_screePlot.png"), 
              width = 190, height = 80) 

# save_plot_png(biplot_basic, 
#               filename = paste0("output/01_", nameDat, "_biplotBasic.png"), 
#               width = 150, height = 120) 

save_plot_png(biplot_group, 
              filename = paste0("output/01_", nameDat, "_biplotGroup.png"), 
              width = 200, height = 130)

save_plot_png(eigencorplot_output, 
              filename = paste0("output/01_", nameDat, "_eigencorplot.png"), 
              width = 200, height = 90)

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
sum(pca_geneSymbols, na.rm = TRUE) #[1] 72.23652


# make new heatmap (Cor PCs vs. metadta) ----------------------------
library(dplyr)
library(ggplot2)
library(tidyr)

# Select PCs and metadata variables
pcs <- 1:5
vars <- c("drug","dose","condition","time")

# Create empty list to store results
results_list <- list()

for (pc in pcs) {
  for (var in vars) {
    x <- pca_output$rotated[, pc]
    y <- metadata[[var]]
    
    if (is.numeric(y)) {
      # Numeric: Pearson correlation
      cor_val <- cor(x, y, method = "pearson")
      p_val <- cor.test(x, y)$p.value
      results_list <- append(results_list, list(data.frame(PC = paste0("PC", pc),
                                                           Variable = var,
                                                           Value = cor_val,
                                                           P = p_val)))
    } else {
      # Categorical: ANOVA F-statistic
      aov_res <- aov(x ~ y)
      f_val <- summary(aov_res)[[1]][["F value"]][1]
      p_val <- summary(aov_res)[[1]][["Pr(>F)"]][1]
      results_list <- append(results_list, list(data.frame(PC = paste0("PC", pc),
                                                           Variable = var,
                                                           Value = f_val,
                                                           P = p_val)))
    }
  }
}

# Combine into one data frame
results_df <- bind_rows(results_list)

# Optional: transform p-values for better display (e.g., -log10)
results_df <- results_df %>%
  mutate(P_label = signif(P, 2))

# Plot heatmap
ggplot(results_df, aes(x = Variable, y = PC, fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = P_label), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       name = "Correlation / F") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_blank(),
        plot.title = element_text(size = 16, face = "bold")) +
  ggtitle("PCs vs Metadata: Color = correlation/F, Number = p-value")


library(dplyr)

# Loop over PCs and categorical variables
results <- sapply(1:5, function(pc) {
  sapply(c("drug","dose","condition","time"), function(var) {
    aov_res <- aov(pca_output$rotated[,pc] ~ metadata[[var]])
    summary(aov_res)[[1]][["Pr(>F)"]][1]
  })
})
results

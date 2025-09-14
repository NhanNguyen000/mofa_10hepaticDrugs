# set up enviroments
rm(list = ls())

library(tidyverse)
library(pheatmap)

source("code/00_plotSettings.R") # plot settings

# load data =======================================================================
load("processedDat/inputDat.RData")

#nameDat <- "normRNA_allEnsembleIds" # all measured Ensemble IDs
nameDat <- "normRNA_geneSymbols" # only measured Ensemble IDs converted into gene symbol
inputDat <- dat[[nameDat]]

sample_metadata <- expDesignTab %>%
  mutate(treatment = ifelse(drug == "Con", "Control", "Treated"))
cor_matrix <- cor(inputDat, method = "pearson")

# Create annotation for samples
annotation_df <- data.frame(
  Treatment = sample_metadata$treatment,
  Condition = sample_metadata$condition2,
  Dose = sample_metadata$dose,
  Time = sample_metadata$time,
  row.names = sample_metadata$name
)

# Define colors for annotations
ann_colors = list(
  Treatment = c(Control = "#fff5f0", Treated = "#d6604d"),
  Condition = c("ConDMSO"= "#d1e5f0", "ConUNTR" = "#bcbd22", # 2 controls
           "5FU"= "#66A61E", "APA" = "#00A087", "AZA" = "#F39B7F", # 10 drugs
           "CYC" = "#A65628", "DIC" = "#F4A582", "ISO"  = "#80B1D3", 
           "MTX" = "#DF65B0", "PHE" = "#7570B3", "RIF" = "#E6AB02", 
           "VPA" =  "#4DBBD5"), 
  Dose = c("DMSO"= "#d1e5f0", "UNTR" = "#bcbd22", # 2 controls
           "The" =  "#7BB3DB", "Tox" = "#F8C471"), # The and Tox
   Time = c("000" = "#FFF5F0", "002" = "#F9DFD8", "008" ="#F3CAC1",
           "024" = "#EDB5AA", "072" = "#E79F92", "168" = "#E18A7B",
           "240" = "#DB7564", "336" = "#D6604D")
)

# colorRampPalette(c("#fff5f0", "#d6604d"))(8)
# [1] "#FFF5F0" "#F9DFD8" "#F3CAC1" "#EDB5AA" "#E79F92" "#E18A7B" "#DB7564" "#D6604D"

heatmap_plot <- pheatmap(
  cor_matrix,
  annotation_row = annotation_df,
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  color = colorRampPalette(journal_colors$reds)(100),
  #breaks = seq(-1, 1, length.out = 101),
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "ward.D2",
  fontsize = 8,
  fontsize_row = 7,
  fontsize_col = 7,
  angle_col = 45,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_legend = TRUE,
  legend = TRUE
)

heatmap_plot 

# save plots ------------------------------------------------------------
save_plot_png(heatmap_plot , 
              filename = paste0("output/01_heatmap.png"), 
              width = 250, height = 300, dpi = 600)
save_plot_png(heatmap_plot , 
              filename = paste0("output/01_heatmap_square.png"), 
              width = 210, height = 210, dpi = 300)

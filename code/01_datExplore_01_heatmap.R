# set up enviroments
rm(list = ls())

library(tidyverse)
library(pheatmap)

source("code/00_plotSettings.R") # plot settings

# load data, make heatmap =======================================================================
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

# subset of the dendrogram (hclust object) ---------------------------------------
library(dendextend)

row_dend <- as.dendrogram(heatmap_plot$tree_row)
col_dend <- as.dendrogram(heatmap_plot$tree_col)
all.equal(row_dend, col_dend) # TRUE

# make plot
plot(row_dend, main = "samples dendrogram")
abline(h = 3, lty = 2, lwd = 2, col = "red")

# cut tree (without order?)
clusters <- cutree(row_dend, h = 3)
table(clusters)
cluster_groups <- split(names(clusters), clusters)
cluster_groups

unique(substring(cluster_groups$`1`, 1, 7))
# [1] "5FU_Tox" "5FU_The" "ConDMSO" "ConUNTR"
unique(substring(cluster_groups$`2`, 1, 7))
# [1] "APA_Tox" "APA_The" "DIC_The" "DIC_Tox" "MTX_The" "MTX_Tox" "Rif_The" "Rif_Tox"
unique(substring(cluster_groups$`3`, 1, 7))
[1] "APA_The" "AZA_Tox" "AZA_The" "ConDMSO" "ConUNTR" "CYC_The" "CYC_Tox" "ISO_The"
[9] "ISO_Tox" "MTX_The" "MTX_Tox" "PHE_The" "PHE_Tox" "VPA_Tox" "VPA_The"


# reorder clusters by dendrogram order left-to-right -----------------------
clusters_ordered <- clusters[order.dendrogram(row_dend)]
cluster_newIds <- match(clusters_ordered, unique(clusters_ordered))
table(cluster_newIds)
# cluster_newIds
# 1   2   3 
# 154  66 254 

# which samples belong to which cluster
names(cluster_newIds) <- names(clusters)
cluster_groups_newIds <- split(names(cluster_newIds), cluster_newIds)
cluster_groups_newIds

unique(substring(cluster_groups_newIds$`1`, 1, 7))
# [1] "5FU_Tox" "5FU_The" "APA_Tox" "APA_The" "AZA_Tox" "AZA_The" "ConDMSO" "ConUNTR"

unique(substring(cluster_groups_newIds$`2`, 1, 7))
# [1] "ConUNTR" "CYC_The" "CYC_Tox" "DIC_The" "DIC_Tox"

unique(substring(cluster_groups_newIds$`3`, 1, 7))
# [1] "DIC_The" "DIC_Tox" "ISO_The" "ISO_Tox" "MTX_The" "MTX_Tox" "PHE_The" "PHE_Tox"
# [9] "Rif_The" "Rif_Tox" "VPA_Tox" "VPA_The"


#plot the dendrogram -  use WGCNA (work)----------------------
library(WGCNA)


# get the branch
cluster_samples <- cluster_groups$`1`
cluster_samples <- cluster_groups$`2`
cluster_samples <- cluster_groups$`3`
cluster_samples

sub_dend <- prune(row_dend, setdiff(labels(row_dend), cluster_samples))

# annotation
cluster_samples_treatment <- ifelse(grepl("Con", cluster_samples), "#fff5f0",  "#d6604d")
cluster_samples_condition <- ifelse(
  grepl("5FU", cluster_samples), "#66A61E",  
  ifelse(grepl("APA", cluster_samples), "#00A087", 
         ifelse(grepl("AZA", cluster_samples), "#F39B7F", 
                ifelse(grepl("CYC", cluster_samples), "#A65628", 
                       ifelse(grepl("DIC", cluster_samples), "#F4A582", 
                              ifelse(grepl("ISO", cluster_samples), "#80B1D3",
                                     ifelse(grepl("MTX", cluster_samples), "#DF65B0",
                                            ifelse(grepl("PHE", cluster_samples), "#7570B3", 
                                                   ifelse(grepl("RIF", cluster_samples), "#E6AB02",
                                                          ifelse(grepl("VPA", cluster_samples), "#4DBBD5",
                                                                 ifelse(grepl("ConDMSO", cluster_samples), "#d1e5f0", "#bcbd22")))))))))))

cluster_samples_dose <- ifelse(grepl("The", cluster_samples), "#7BB3DB", 
                               ifelse(grepl("Tox", cluster_samples), "#F8C471", 
                                      ifelse(grepl("ConDMSO", cluster_samples),"#d1e5f0", "#bcbd22")))
cluster_samples_time <- ifelse(grepl("000", cluster_samples), "#FFF5F0", 
                               ifelse(grepl("002", cluster_samples), "#F9DFD8", 
                                      ifelse(grepl("008", cluster_samples), "#F3CAC1",
                                             ifelse(grepl("024", cluster_samples), "#EDB5AA",
                                                    ifelse(grepl("072", cluster_samples), "#E79F92",
                                                           ifelse(grepl("168", cluster_samples), "#E18A7B",
                                                                  ifelse(grepl("240", cluster_samples), "#DB7564", "#D6604D")))))))



# plot
#sub_dend <- hang.dendrogram(sub_dend, hang_height = 0)  

png(file = "output/01_dedrogram_cluster1.png", width = 70, height = 90, unit = "mm", res = 600)
png(file = "output/01_dedrogram_cluster2.png", width = 160, height = 90, unit = "mm", res = 600)
png(file = "output/01_dedrogram_cluster3.png", width = 260, height = 90, unit = "mm", res = 600)

plotDendroAndColors(as.hclust(sub_dend) , 
                    cbind(cluster_samples_time,
                          cluster_samples_dose, 
                          cluster_samples_condition, 
                          cluster_samples_treatment),
                    c("Time", "Dose", "Condition", "Treatment"),
                    dendroLabels = FALSE, 
                    addGuide = TRUE
                    )

dev.off()

# save plots ------------------------------------------------------------
save_plot_png(heatmap_plot , 
              filename = paste0("output/01_heatmap.png"), 
              width = 250, height = 300, dpi = 600)
save_plot_png(heatmap_plot , 
              filename = paste0("output/01_heatmap_square.png"), 
              width = 210, height = 210, dpi = 300)

save_plot_png(heatmap_plot , 
              filename = paste0("output/01_heatmap_square_zoomIn.png"), 
              width = 500, height = 500, dpi = 2400)


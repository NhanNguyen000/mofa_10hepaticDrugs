# set up enviroment
rm(list = ls())

library(tidyverse)
library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(MOFA2)
library(ggrepel)


## correct the plot functions in MOFA2 source code ------------------------------
plot_enrichment <- function(enrichment.results, factor, alpha = 0.1, max.pathways = 25,
                            text_size = 1.0, dot_size = 5.0) {
  
  # Sanity checks
  stopifnot(is.numeric(alpha)) 
  stopifnot(length(factor)==1) 
  if (is.numeric(factor)) factor <- colnames(enrichment.results$pval.adj)[factor]
  if(!factor %in% colnames(enrichment.results$pval)) 
    stop(paste0("No gene set enrichment calculated for factor ", factor))
  
  # get p-values
  p.values <- enrichment.results$pval.adj
  
  # Get data  
  tmp <- data.frame(
    pvalues = p.values[,factor, drop=TRUE], 
    pathway = rownames(p.values)
  )
  
  # Filter out pathways
  tmp <- tmp[tmp$pvalue<=alpha,,drop=FALSE]
  if (nrow(tmp)==0) stop("No siginificant pathways at the specified alpha threshold")
  
  # If there are too many pathways enriched, just keep the 'max_pathways' more significant
  if (nrow(tmp)>max.pathways) tmp <- head(tmp[order(tmp$pvalue),],n=max.pathways)
  
  # Convert pvalues to log scale
  tmp$logp <- -log10(tmp$pvalue+1e-100)
  
  #order according to significance
  #tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$pathway <- factor(tmp$pathway <- tmp$pathway, levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$start <- 0
  
  p <- ggplot(tmp, aes_string(x="pathway", y="logp")) +
    geom_point(size=dot_size) +
    geom_hline(yintercept=-log10(alpha), linetype="longdash") +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes_string(xend="pathway", yend="start")) +
    ylab("-log pvalue") +
    coord_flip() +
    theme(
      axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

plot_enrichment_detailed <- function(enrichment.results, factor, 
                                     alpha = 0.1, max.genes = 5, max.pathways = 10, text_size = 3) {
  
  # Sanity checks
  stopifnot(is.list(enrichment.results))
  stopifnot(length(factor)==1) 
  if (!is.numeric(factor)) {
    if(!factor %in% colnames(enrichment.results$pval)) 
      stop(paste0("No feature set enrichment calculated for ", factor))
  }
  
  # Fetch and prepare data  
  
  # foo
  foo <- reshape2::melt(enrichment.results$feature.statistics[,factor], na.rm=TRUE, value.name="feature.statistic")
  #foo$feature <- rownames(foo)
  foo$feature <- rownames(enrichment.results$feature.statistics)
  
  # bar
  feature.sets <- enrichment.results$feature.sets %>% as.data.frame() %>% rownames_to_column("pathway")
  feature.sets[feature.sets==0] <- NA
  bar <- reshape2::melt(feature.sets, na.rm=TRUE)[,c(1,2)]
  colnames(bar) <- c("pathway","feature")
  bar$pathway <- as.character(bar$pathway)
  bar$feature <- as.character(bar$feature)
  
  # baz
  baz <- reshape2::melt(enrichment.results$pval.adj[,factor], value.name="pvalue", na.rm=TRUE)
  #baz$pathway <- rownames(baz)
  baz$pathway <- rownames(enrichment.results$pval.adj)
  
  # Filter out pathways by p-values
  baz <- baz[baz$pvalue<=alpha,,drop=FALSE]
  if(nrow(baz)==0) {
    stop("No siginificant pathways at the specified alpha threshold. \n
         For an overview use plot_enrichment_heatmap().")
  } else {
    if (nrow(baz)>max.pathways)
      baz <- head(baz[order(baz$pvalue),],n=max.pathways)
  }
  
  # order pathways according to significance
  baz$pathway <- factor(baz$pathway, levels = baz$pathway[order(baz$pvalue, decreasing = TRUE)])
  
  # Merge
  foobar <- merge(foo, bar, by="feature")
  tmp <- merge(foobar, baz, by="pathway")
  
  # Select the top N features with the largest feature.statistic (per pathway)
  tmp_filt <- top_n(group_by(tmp, pathway), n=max.genes, abs(feature.statistic))
  
  # Add number of features and p-value per pathway
  pathways <- unique(tmp_filt$pathway)
  
  # Add Ngenes and p-values to the pathway name
  feature.sets_v2 <- enrichment.results$feature.sets
  feature.sets_v2[feature.sets_v2==0] <- NA
  df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets_v2,na.rm=TRUE)[pathways])
  #df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets,na.rm=TRUE)[pathways])
  df <- merge(df, baz, by="pathway")
  df$pathway_long_name <- sprintf("%s\n (Ngenes = %d) \n (p-val = %0.2g)",df$pathway, df$nfeatures, df$pvalue)
  tmp <- merge(tmp, df[,c("pathway","pathway_long_name")], by="pathway")
  tmp_filt <- merge(tmp_filt, df[,c("pathway","pathway_long_name")], by="pathway")
  
  # sort pathways by p-value
  order_pathways <- df$pathway_long_name[order(df$pvalue,decreasing=TRUE) ]
  tmp$pathway_long_name <- factor(tmp$pathway_long_name, levels=order_pathways)
  tmp_filt$pathway_long_name <- factor(tmp_filt$pathway_long_name, levels=order_pathways)
  
  p <- ggplot(tmp, aes_string(x="pathway_long_name", y="feature.statistic")) +
    geom_text_repel(aes_string(x="pathway_long_name", y="feature.statistic", label="feature"), size=text_size, color="black", force=1, data=tmp_filt) +
    geom_point(size=0.5, color="lightgrey") +
    geom_point(aes_string(x="pathway_long_name", y="feature.statistic"), size=1, color="black", data=tmp_filt) +
    labs(x="", y="Weight (scaled)", title="") +
    coord_flip() +
    theme(
      axis.line = element_line(color="black"),
      axis.text.y = element_text(size=rel(0.75), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.0), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

plot_enrichment_v2 <- function(enrichment.results, factor, alpha = 0.1, selected.pathway = NA,
                            text_size = 1.0, dot_size = 5.0) {
  
  # Sanity checks
  stopifnot(is.numeric(alpha)) 
  stopifnot(length(factor)==1) 
  if (is.numeric(factor)) factor <- colnames(enrichment.results$pval.adj)[factor]
  if(!factor %in% colnames(enrichment.results$pval)) 
    stop(paste0("No gene set enrichment calculated for factor ", factor))
  
  # get p-values
  p.values <- enrichment.results$pval.adj
  
  # Get data  
  tmp <- data.frame(
    pvalues = p.values[,factor, drop=TRUE], 
    pathway = rownames(p.values)
  )
  
  # Filter out pathways
  tmp <- tmp[tmp$pvalue<=alpha,,drop=FALSE]
  if (nrow(tmp)==0) stop("No siginificant pathways at the specified alpha threshold")
  
  # select pathways
  tmp <- tmp[which(tmp$pathway %in% selected.pathway),]
  
  # Convert pvalues to log scale
  tmp$logp <- -log10(tmp$pvalue+1e-100)
  
  #order according to significance
  #tmp$pathway <- factor(tmp$pathway <- rownames(tmp), levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$pathway <- factor(tmp$pathway <- tmp$pathway, levels = tmp$pathway[order(tmp$pvalue, decreasing = TRUE)])
  tmp$start <- 0
  
  p <- ggplot(tmp, aes_string(x="pathway", y="logp")) +
    geom_point(size=dot_size) +
    geom_hline(yintercept=-log10(alpha), linetype="longdash") +
    scale_color_manual(values=c("black","red")) +
    geom_segment(aes_string(xend="pathway", yend="start")) +
    ylab("-log pvalue") +
    coord_flip() +
    theme(
      axis.text.y = element_text(size=rel(text_size), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.2), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}


plot_enrichment_detailed_v2 <- function(enrichment.results, factor, 
                                     alpha = 0.1, max.genes = 5, selected.pathway = NA, text_size = 3) {
  
  # Sanity checks
  stopifnot(is.list(enrichment.results))
  stopifnot(length(factor)==1) 
  if (!is.numeric(factor)) {
    if(!factor %in% colnames(enrichment.results$pval)) 
      stop(paste0("No feature set enrichment calculated for ", factor))
  }
  
  # Fetch and prepare data  
  
  # foo
  foo <- reshape2::melt(enrichment.results$feature.statistics[,factor], na.rm=TRUE, value.name="feature.statistic")
  #foo$feature <- rownames(foo)
  foo$feature <- rownames(enrichment.results$feature.statistics)
  
  # bar
  feature.sets <- enrichment.results$feature.sets %>% as.data.frame() %>% rownames_to_column("pathway")
  feature.sets[feature.sets==0] <- NA
  bar <- reshape2::melt(feature.sets, na.rm=TRUE)[,c(1,2)]
  colnames(bar) <- c("pathway","feature")
  bar$pathway <- as.character(bar$pathway)
  bar$feature <- as.character(bar$feature)
  
  # baz
  baz <- reshape2::melt(enrichment.results$pval.adj[,factor], value.name="pvalue", na.rm=TRUE)
  #baz$pathway <- rownames(baz)
  baz$pathway <- rownames(enrichment.results$pval.adj)
  
  # Filter out pathways by p-values
  baz <- baz[baz$pvalue<=alpha,,drop=FALSE]
  if(nrow(baz)==0) {
    stop("No siginificant pathways at the specified alpha threshold. \n
         For an overview use plot_enrichment_heatmap().")
  } else {
    # select pathways
    baz <- baz[which(baz$pathway %in% selected.pathway),]
  }
  
  # order pathways according to significance
  baz$pathway <- factor(baz$pathway, levels = baz$pathway[order(baz$pvalue, decreasing = TRUE)])
  
  # Merge
  foobar <- merge(foo, bar, by="feature")
  tmp <- merge(foobar, baz, by="pathway")
  
  # Select the top N features with the largest feature.statistic (per pathway)
  tmp_filt <- top_n(group_by(tmp, pathway), n=max.genes, abs(feature.statistic))
  
  # Add number of features and p-value per pathway
  pathways <- unique(tmp_filt$pathway)
  
  # Add Ngenes and p-values to the pathway name
  feature.sets_v2 <- enrichment.results$feature.sets
  feature.sets_v2[feature.sets_v2==0] <- NA
  df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets_v2,na.rm=TRUE)[pathways])
  #df <- data.frame(pathway=pathways, nfeatures=rowSums(feature.sets,na.rm=TRUE)[pathways])
  df <- merge(df, baz, by="pathway")
  df$pathway_long_name <- sprintf("%s\n (Ngenes = %d) \n (p-val = %0.2g)",df$pathway, df$nfeatures, df$pvalue)
  tmp <- merge(tmp, df[,c("pathway","pathway_long_name")], by="pathway")
  tmp_filt <- merge(tmp_filt, df[,c("pathway","pathway_long_name")], by="pathway")
  
  # sort pathways by p-value
  order_pathways <- df$pathway_long_name[order(df$pvalue,decreasing=TRUE) ]
  tmp$pathway_long_name <- factor(tmp$pathway_long_name, levels=order_pathways)
  tmp_filt$pathway_long_name <- factor(tmp_filt$pathway_long_name, levels=order_pathways)
  
  p <- ggplot(tmp, aes_string(x="pathway_long_name", y="feature.statistic")) +
    geom_text_repel(aes_string(x="pathway_long_name", y="feature.statistic", label="feature"), size=text_size, color="black", force=1, data=tmp_filt) +
    geom_point(size=0.5, color="lightgrey") +
    geom_point(aes_string(x="pathway_long_name", y="feature.statistic"), size=1, color="black", data=tmp_filt) +
    labs(x="", y="Weight (scaled)", title="") +
    coord_flip() +
    theme(
      axis.line = element_line(color="black"),
      axis.text.y = element_text(size=rel(0.75), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.0), vjust=0.5, color='black'),
      axis.title.y=element_blank(),
      legend.position='none',
      panel.background = element_blank()
    )
  
  return(p)
}

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

rm(dat, expDesignTab)

head(model@samples_metadata)
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


# GSEA analysis -------------------------------------------
# tutorial: https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/GSEA.html
# https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/CLL.html

# reactome
#data("reactomeGS") # reactome database provided by MOFAdata, but it is an old version 
load("processedDat/reactomeGS.RData") # up-to-date reactome database
head(rownames(reactomeGS_human), n=3)
head(colnames(reactomeGS_human), n=3)

# check gene name type
head(features_names(model)[["RNAseq"]])

# GSEA with default options and on positive/negative weights
res_list <- list()
# 
# res_list[["positive"]] <- run_enrichment(
#   model,
#   feature.sets = reactomeGS, 
#   view = "RNAseq",
#   sign = "positive", # GSEA on positive weights
#   alpha = 0.05
# )
# 
# res_list[["negative"]] <- run_enrichment(
#   model,
#   feature.sets = reactomeGS, 
#   view = "RNAseq",
#   sign = "negative" , # GSEA on positive weights
#   alpha = 0.05
# )

res_list[["all"]] <- run_enrichment(
  model,
  feature.sets = reactomeGS_human, 
  view = "RNAseq",
  alpha = 0.05
)

# check outcome
names(res_list)

# names(res_list$positive)
# res_list$positive$set.statistics[1:5,1]
# res_list$positive$pval.adj[1:5,1]
# 
# names(res_list$negative)
# res_list$negative$set.statistics[1:5,1]
# res_list$negative$pval.adj[1:5,1]

names(res_list$all)
res_list$all$set.statistics[1:5,1]
res_list$all$pval.adj[1:5,1]

# Plot results of GSEA analysis - heatmap ------------------------------------------
plot_enrichment_heatmap(res_list$all)
plot_enrichment_heatmap(res_list$all, alpha = 0.0001)

# plot_enrichment_heatmap(res_list$positive)
# plot_enrichment_heatmap(res_list$negative)

# modify the heatmap plot from MOFA
# add hierachy annotation
p <- plot_enrichment_heatmap(res_list$all)

# get p-values
a <- res_list$all$pval.adj
a <- a[,colMeans(is.na(a))<1]

# cap p-values 
a[a<1e-50] <- 1e-50

alpha = 0.01
# Apply Log transform
a <- -log10(a+1e-50)
alpha <- -log10(alpha)
col <- colorRampPalette(c("lightgrey","red"))(n=100)
pheatmap(a, color = col, cluster_cols = FALSE, show_rownames = FALSE)

# option
b <- a
b[b>0.05] <- NA
b <- b[rowSums(is.na(b)) < ncol(b),]

b <- -log10(b+1e-50)
b[is.na(b)] <- 0
alpha <- -log10(alpha)
col <- colorRampPalette(c("lightgrey","red"))(n=100)
pheatmap(b, color = col, cluster_cols = FALSE, show_rownames = FALSE)


# Define color palette: white for 0, then red gradient
col <- colorRampPalette(c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#de2d26", "#a50f15"))(100)

# Draw heatmap
pheatmap(b,
         color = col,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         legend_labels = c("–log10(p.adjvalue)"), # legend caption
         legend = TRUE)

# Make color palette for parent pathways
grouped_long <- grouped_matrix %>% unnest(TopGroup)

# Make wide format with TRUE/FALSE
grouped_matrix_v2 <- grouped_long %>%
  mutate(value = TRUE) %>%
  pivot_wider(names_from = TopGroup, values_from = value, values_fill = FALSE) %>%
  left_join(pathway_ids) %>% dplyr::select(-Species, -PathwayID) %>% 
  column_to_rownames("PathwayName")

# Use PathwayID as rownames
grouped_matrix_v2 <- as.data.frame(grouped_matrix_v2)


k<- grouped_matrix_v2[rownames(b), ]  %>% as.data.frame()
k[k=="TRUE"] <- "yes"
k[k=="FALSE"] <- "no"

k2 <- colnames(k) %>% as.data.frame() %>% rename("." = "PathwayID") %>% 
  left_join(pathway_ids) %>%
  mutate(PathwayName = ifelse(PathwayID == "R-HSA-937061",  
                              "TRIF (TICAM1)-mediated TLR4 signaling", PathwayName)) %>% 
  mutate(pathway_long_name = paste0(PathwayName, " (", PathwayID, ")"))

colnames(k3) <- k2$pathway_long_name

# Define colors for binary membership
ann_colors <- lapply(k, function(x) c("no" = "white", "yes" = "blue"))

pheatmap(b,
         color = col,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         annotation_row = k,
         annotation_colors = ann_colors,
         cutree_rows = 10)

ann_colors3 <- lapply(k3, function(x) c("no" = "white", "yes" = "blue"))

pheatmap(b,
         color = col,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         annotation_row = k3,
         annotation_colors = ann_colors3,
         cutree_rows = 10,
         annotation_legend = FALSE,
         annotation_width = 1.5,      # adjust width of annotation bar
         cellwidth = 18,
         fontsize_row = 5,
         angle_col = 315)


# Plot results of GSEA analysis per factors------------------------------------------
# plot top enriched pathways per factor
for (condition in names(res_list)) {
  res_list[[condition]]$pval.adj <- as.data.frame(res_list[[condition]]$pval.adj)
}

plot_enrichment(res_list$all, factor = 1, max.pathways = 30)
plot_enrichment(res_list$positive, factor = 1, max.pathways = 15)
plot_enrichment(res_list$negative, factor = 1, max.pathways = 15)

# plot selected enriched pathways per factor
selected_pathways <- res_list$all$sigPathways[[1]][1:30] # factor 1
selected_pathways <- k2$PathwayName # factor 1

plot_enrichment_v2(res_list$all, factor = 1, 
                   selected.pathway = selected_pathways)
plot_enrichment_v2(res_list$all, factor = 2, 
                   selected.pathway = selected_pathways)
plot_enrichment_v2(res_list$all, factor = 3, 
                   selected.pathway = selected_pathways)
plot_enrichment_v2(res_list$all, factor = 4, 
                   selected.pathway = selected_pathways)

plot_enrichment_v2(res_list$negative, factor = 1, 
                   selected.pathway = selected_pathways)

# plot top enriched pathways with top gene per factor
plot_enrichment_detailed(res_list$all, 
                         factor = 1, 
                         max.genes = 8, 
                         max.pathways = 15)

plot_enrichment_detailed(res_list$positive, 
                         factor = 1, 
                         max.genes = 8, 
                         max.pathways = 5)

plot_enrichment_detailed(res_list$negative, 
                         factor = 1, 
                         max.genes = 8, 
                         max.pathways = 5)

# plot selected enriched pathways with top gene per factor
plot_enrichment_detailed_v2(res_list$all, 
                         factor = 1, 
                         max.genes = 8, 
                         selected.pathway = selected_pathways)

plot_enrichment_detailed_v2(res_list$all, 
                            factor = 1, 
                            max.genes = 8, 
                            selected.pathway = selected_pathways)

# can run some GO pathway similarity analysis, and select the top GO pathways



# per gene
#genes <- list("ENSG00000134571","ENSF00000159173")
genes <- list("CYP1A1","APOA1")
genes <- list("VAMP2","TUBA1A")

genes %>% 
  map(~ plot_factors(
    model, factors = c(1,2), 
    color_by = "dose", scale = T, legend = T)) %>% 
  cowplot::plot_grid(plotlist=., nrow=1)

# Save plots ------------------------------------------------------------
png(file= paste0("output/05_", nameDat, "_mofa_pathwayEnrichment_heatmap_positive.png"), height = 720)
plot_enrichment_heatmap(res.positive)
dev.off()

png(file= paste0("output/05_", nameDat, "_mofa_pathwayEnrichment_heatmap_negative.png"), height = 720)
plot_enrichment_heatmap(res.negative)
dev.off()

png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_drug.png"), width = 720, height = 320)
factors_drug + theme(text = element_text(size = 18))
dev.off()

png(file= paste0("output/04_", nameDat, "_mofa_latentFactors_drug.png"), width = 720, height = 320)
factors_time + theme(text = element_text(size = 18))
dev.off()

# pathway database -----------------------------------------------------------------------
# # The Reactome database from MOFAdata package is an old version (v59)
# library(MOFAdata)
# data(reactomeGS)
# reactomeGS_MOFAdata <- reactomeGS
# # Reactome: https://reactome.org/download-data, current version v92, access date 22/07/2025.

library(GSEABase)
library(org.Hs.eg.db)

# Download the latest Reactome GMT file
temp_file <- tempfile(fileext = ".zip")
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", 
              temp_file)
unzip(temp_file, exdir = tempdir())
gmt_path <- list.files(tempdir(), pattern = "ReactomePathways\\.gmt$", full.names = TRUE)

# read file, and convert to a list (each element is a pathway, with a vector of gene symbols)
gene_sets <- getGmt(gmt_path)
gene_list <- lapply(gene_sets, geneIds)
names(gene_list) <- sapply(gene_sets, setName)

# Create binary membership matrix
gene_names <- sort(unique(unlist(gene_list)))
pathway_names <- names(gene_list)

binary_matrix <- matrix(0, nrow = length(gene_names), ncol = length(pathway_names),
                        dimnames = list(gene_names, pathway_names)) # Initialize empty binary matrix
for (p in pathway_names) {binary_matrix[gene_list[[p]], p] <- 1} # Fill in 1s where genes are in pathways

reactomeGS <- binary_matrix %>% t()
dim(reactomeGS) # [1]  2785 11886, more pathways, less gene (here are gene symbols)
dim(reactomeGS_MOFAdata) # [1]  1304 18818, here are gene Ensemble IDs

#save(reactomeGS, file = "processedDat/reactomeGS.RData")
# when we use this new reactomeGS feature, MOFA's functions still work well.

# map pathways hierachy -------------------------------
# pathway information
download.file("https://reactome.org/download/current/ReactomePathways.txt",
              file.path(tempdir(), "ReactomePathways.txt"))
download.file("https://reactome.org/download/current/ReactomePathwaysRelation.txt",
              file.path(tempdir(), "ReactomePathwaysRelation.txt"))

pathways_info <- read.delim(file.path(tempdir(), "ReactomePathways.txt"),
                            header = FALSE, stringsAsFactors = FALSE)
colnames(pathways_info) <- c("PathwayID", "PathwayName", "Species")
unique(pathways_info$Species)

pathway_info_human <- pathways_info %>% filter(Species ==  "Homo sapiens" )

pathway_ids_v1 <- pathway_info_human[pathway_info_human$PathwayName %in% rownames(reactomeGS), ]
pathway_ids_v2 <- pathway_info_human[pathway_info_human$PathwayName %in% rownames(reactomeGS_MOFAdata), ]

dim(reactomeGS) # [1]  2785 11886
dim(pathway_ids_v1) #[1] 2728    3
dim(pathway_ids_v2) #[1]1155    3, reactomeGS_MOFAdata has less pathways

pathway_ids <- pathway_ids_v1
# pathways hierachy
relations <- read.delim(file.path(tempdir(), "ReactomePathwaysRelation.txt"),
                        header = FALSE, stringsAsFactors = FALSE)
colnames(relations) <- c("ParentID", "ChildID")

# Map pathway names → IDs
# library(igraph)
valid_ids <- na.omit(unique(as.character(pathway_ids_v1$PathwayID))) # Remove NA pathway IDs

relations_subset <- relations[
  relations$ChildID %in% valid_ids | relations$ParentID %in% valid_ids, ]

g <- graph_from_data_frame(relations_subset, directed = TRUE) # (vertex names will be IDs)
get_top_parent <- function(node, g) { #get top-level parent(s)
  if (!(node %in% V(g)$name)) {
    return(NA_character_)  # skip if node not in graph
  }
  # traverse "upwards" (to parents)
  ancestors <- ego(g, order = vcount(g), nodes = node, mode = "in")[[1]]
  roots <- V(g)[degree(g, mode = "in") == 0]  # root pathways
  intersect(ancestors, roots)$name
}

top_parents <- sapply(valid_ids, get_top_parent, g = g)
head(top_parents)

# Add grouping info
grouped_matrix <- tibble(
  PathwayID = names(top_parents),
  TopGroup = top_parents
)

length(unique(unlist(top_parents))) # 30 pathways

reactomeGS_human <- reactomeGS[rownames(reactomeGS) %in% pathway_info_human$PathwayName,]
dim(reactomeGS) #[1]  2728 11886
dim(reactomeGS_human) #[1]  2728 11886
save(reactomeGS_human, grouped_matrix, pathway_ids, file = "processedDat/reactomeGS.RData")

# Compare different statistical tests -----------------------------
res_corAdj <- run_enrichment(
  model, view = "RNAseq", factors = 1:3,
  feature.sets = reactomeGS,
  # sign = "negative",
  statistical.test = "cor.adj.parametric")
# "cor.adj.parametric" is parametric test adjust for PCs, 
# but latent factors do not correlate with other, 
# so the defautl "parametric" test is MOFA is more suitable

# Compare the histogram of the p-values. 
# Clearly the correlation-adjusted test results in more conservative estimates:
dt <- rbind(res_list$all$pval[,1:3] %>% as.data.table %>% 
              .[,c("test","pathway"):=list("parametric",1:.N)],
            res_corAdj$pval[,1:3] %>% as.data.table %>% 
              .[,c("test","pathway"):=list("parametric.adj",1:.N)]) %>% 
  melt(id.vars=c("test","pathway"), variable.name="factor")

ggplot(dt, aes(x=value, fill=test)) +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  geom_histogram() +
  theme_bw() +
  theme(legend.position = "top",
        legend.title = element_blank())

# Yet, the resulting p-values are highlighy correlated between the two tests, as expected (???)
dt2 <- dt %>% dcast(factor+pathway~test)

ggplot(dt2, aes(x=parametric, y=parametric.adj)) +
  geom_point(size=0.5) +
  geom_abline(slope=1, intercept=0, color="orange") +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  labs(x="Parametric p-value", y="Adjusted parametric p-value") +
  theme_bw() +
  theme(
    legend.position = "top"
  )


# https://peder.quarto.pub/blog/posts/hm_clust/#comparison-with-the-pheatmap-package

rm(list = ls())

library(tidyverse) # For wrangling
library(ggdendro) # Required to make the dendrogram in ggplot
#library(survival) # Contains the pbc dataset
#library(mice) # Contains the function to calculate the Nelson-Aalen estimator
library(patchwork) # Required to assemble the plot composite
library(RColorBrewer) # Contains several color palettes
library(ggpubr) # An alternative to patchwork

font <- "Roboto"

# load data =======================================================================
load("processedDat/inputDat.RData")
dim(norm_rna) #  60664   474
dim(norm_rna_v2) #  36244   474

# clustering ------------------------------------------------------------------
# the data has been sclaed by DESeq2 vst(), so we countinue with clustering

# produce a dendrogram
dg <- hclust(dist(t(norm_rna_v2))) %>% as.dendrogram()

#dg <- hclust(dist(t(norm_rna_v2), method = "ward.D")) %>% # error, do not know why?
#  as.dendrogram() 

# It's not straight forward to extract data from a hclust object. To produce a ggplot dendrogram we use
# ggdendro::dendro_data() function
ddata_pts <- dendro_data(dg, type = "rectangle")

pts_dendrogram <- ggplot() +
  geom_segment(
    data = segment(ddata_pts),
    aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(
    data = label(ddata_pts),
    aes(x = x, y = -1, label = label),
    size = 2, family = font, color = "#444444", vjust = 0.5, angle = 90, hjust = 1) +
  #scale_x_continuous(limits = c(0, n + 1), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_dendro() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm")) +
  coord_cartesian(clip = "off")

ddata_analytes <- dendro_data(dg, type = "rectangle")
analytes_dendrogram <- ggplot() +
  geom_segment(data = segment(ddata_analytes),
    aes(x = x, y = y, xend = xend, yend = yend),
    position = position_nudge(x = -0.5)) +
  coord_flip(clip = "off") +
  scale_y_reverse() +
  #scale_x_continuous(limits = c(0, length(vars_to_clust)), expand = c(0, 0)) +
  theme_dendro() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "mm"))

# We can visualize this dendrogram
analytes_dendrogram +
  geom_text(
    data = label(ddata_analytes),
    aes(x = x, y = -1, label = label),
    size = 3.5, family = font, color = "#444444", vjust = 2, angle = 0, hjust = 0) +
  theme(plot.margin = unit(c(0, 15, 0, 0), "mm"))

# annotation bars ---------------------------
# Horizontal annotation bar: Histological stage
p_ann_stage <-
  label(ddata_pts) %>%
  rename(id = label) %>%
  left_join(expDesignTab, by = c("id"= "name")) %>%
  ggplot(aes(x = x, y = "drug", fill = drug)) +
  geom_tile(color = "white") +
 # scale_x_continuous(limits = c(0, n + 1), expand = c(0, 0)) +
  scale_y_discrete(position = "right") +
  theme(
    legend.direction = "horizontal",
    legend.title = element_text(family = font, size = 9, color = "#444444"),
    legend.text = element_text(family = font, size = 8, color = "#444444"),
    legend.key.height = unit(3, "mm")
  ) +
  scale_fill_manual(values = brewer.pal(12, "Set3")[c(2:12)]) +
  guides(fill = guide_legend(
    title = "drug",
    nrow = 1, title.position = "top"
  ))

# Get the legend
leg_stage <- ggpubr::get_legend(p_ann_stage)
# Remove legend from the original plot
p_ann_stage <-
  p_ann_stage +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    axis.text.y = element_text(hjust = 0, color = "#444444", family = font, size = 10)
  )

# Horizontal annotation bar: Nelson-Aalen estimator
p_ann_NelsonAalen <-
  label(ddata_pts) %>%
  rename(id = label) %>%
  left_join(expDesignTab, by = c("id"= "name")) %>%
  ggplot(aes(x = x, y = "time", fill = time)) +
  geom_tile(color = "white") +
 # scale_x_continuous(limits = c(0, n + 1), expand = c(0, 0)) +
  scale_y_discrete(position = "right") +
  #scale_fill_gradientn(colours = brewer.pal(9, "Greens"), na.value = "#F9F9F9") +
  theme_void() +
  theme(
    legend.direction = "horizontal",
    legend.title = element_text(family = font, size = 9, color = "#444444"),
    legend.text = element_text(family = font, size = 8, color = "#444444"),
    legend.key.height = unit(2, "mm")
  ) +
  scale_fill_manual(values = brewer.pal(9, "Greens")) +
  guides(fill = guide_legend(
    title = "time",
    title.position = "top"
  ))

# Get the legend
leg_NelsonAalen <- ggpubr::get_legend(p_ann_NelsonAalen)

# Remove legend from the original plot
p_ann_NelsonAalen <-
  p_ann_NelsonAalen +
  theme_void() +
  theme(
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    axis.text.y = element_text(hjust = 0, family = font, size = 10, color = "#444444")
  )

ggarrange(
  ggarrange(p_ann_stage, p_ann_NelsonAalen, nrow = 2, align = "v"),
  ggarrange(leg_stage, leg_NelsonAalen, nrow = 1, align = "h"),
  nrow = 2
)

legends <- ggarrange(leg_stage, leg_NelsonAalen, nrow = 1, align = "h")

# arranging the plot composite ----------------------------------------
# Define layout for the grid in the plot composite
layout <- "
AAAAAAA
BBBBBBB
CCCCCCC
DDDDDDD
"

# Now plot!
pts_dendrogram + p_ann_stage + p_ann_NelsonAalen +  legends +
  plot_layout(design = layout, heights = c(7, 1.6, 1.6, 5))

# old option
layout <- "
#AAAAAAA
#BBBBBBB
#CCCCCCC
DEEEEEEE
#FFFFFFF
"
pts_dendrogram + p_ann_stage + p_ann_NelsonAalen + analytes_dendrogram + legends +
  plot_layout(design = layout, heights = c(7, 1.6, 1.6, 20, 5))

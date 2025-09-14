# ============================================================================
# Figure Settings to suitable for journal
# ============================================================================

# Load required libraries
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(corrplot)
library(viridis)

# figure specifications (follow Nature/Cell journal)
plot_theme <- function() {
  theme_classic() +
    theme(
      # Text elements
      text = element_text(family = "Arial", color = "black"),
      axis.title = element_text(size = 8, color = "black"),
      axis.text = element_text(size = 7, color = "black"),
      legend.title = element_text(size = 8, color = "black"),
      legend.text = element_text(size = 7, color = "black"),
      plot.title = element_text(size = 10, hjust = 0.5, color = "black"),
      
      # Axes
      axis.line = element_line(color = "black", size = 0.5),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(1.5, "pt"),
      
      # Panel
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      
      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      legend.key.size = unit(3, "mm"),
      
      # Margins
      plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt")
    )
}

# PNG export function for journals
save_plot_png <- function(plot, filename, width = 85, height = 85, dpi = 300) {
  # Width/height in mm (Nature: single column = 85mm, double column = 170mm)
  png(filename, 
      width = width, 
      height = height, 
      units = "mm", 
      res = dpi)
  print(plot)
  dev.off()
}

# Color palettes suitable for journals (colorblind-friendly)
journal_colors <- list(
  # Sequential palettes
  blues = c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c"),
  reds = c("#fff5f0", "#fee0d2", "#fcbba1", "#fc9272", "#fb6a4a", "#ef3b2c", "#cb181d"),
  
  # Diverging palettes
  blue_red = c("#053061", "#2166ac", "#4393c3", "#92c5de", "#d1e5f0", 
               "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b", "#67001f"),
  
  # Categorical palettes (Nature/Cell preferred)
  categorical = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
                  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
)


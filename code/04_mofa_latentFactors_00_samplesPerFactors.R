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


# Visualisation of samples in factors 
# samples weighted across  factors side-by-size -------------------------------
condition2_palette <- c(
  "#3C5488", "#E64B35",  # 2 controls
  "#66A61E", "#00A087", "#F39B7F", "#A65628", "#F4A582", 
  "#80B1D3", "#DF65B0", "#7570B3", "#E6AB02",  "#4DBBD5" # 10 drugs
)

factors_condition2 <- plot_factor(
  model, factors = c(1:7), 
  color_by = "condition2", color_name = "Condition",
  add_violin = TRUE, violin_alpha = 0.1, 
  #add_boxplot = TRUE, boxplot_alpha = 0.1, 
  dodge = TRUE, #scale = TRUE
) +  scale_fill_manual(values = condition2_palette) #+ plot_theme()
factors_condition2

condition3_palette <- c(
  "ConDMSO"= "#3C5488", "ConUNTR" = "#E64B35",  # 2 controls
  "5FU" = "#66A61E", "AZA" = "#F39B7F",
  "CYC_VPA_ISO" = "#80B1D3", 
  "PHE" = "#7570B3", 
  "DIC_RIF_MTX_APA" = "#DF65B0" #"#66A61E" #
)

factors_condition3 <- plot_factor(
  model, factors = c(1:7), 
  color_by = "condition3", color_name = "Condition",
  add_violin = TRUE, violin_alpha = 0.1, 
  #add_boxplot = TRUE, boxplot_alpha = 0.1, 
  dodge = TRUE, #scale = TRUE
) +  scale_fill_manual(values = condition3_palette) #+ plot_theme()
factors_condition3

condition4_palette <- c(
  "Control (UNTR, DMSO)" = "#E64B35",  # 2 controls
  "other drugs (5FU, AZA, PHE)" = "#F39B7F",
  "CYC_VPA_ISO" = "#66A61E", 
  "DIC_RIF_MTX_APA" = "#3C5488" #"#66A61E" #
)
factors_condition4 <- plot_factor(
  model, factors = c(1:7), 
  color_by = "condition4", color_name = "Condition",
  add_violin = TRUE, violin_alpha = 0.1, 
  #add_boxplot = TRUE, boxplot_alpha = 0.1, 
  dodge = TRUE, #scale = TRUE
) +  scale_fill_manual(values = condition4_palette) #+ plot_theme()
factors_condition4

# samples weight per factors ---------------------------------------
plot_factor(model, 
            factors = 1, 
            color_by = "dose",
            add_violin = TRUE, violin_alpha = 0.1,
            dodge = TRUE)

plot_factor(model, 
            factor = 1:2,
            color_by = "time",
            shape_by = "dose")


plot_factor(model, 
            factors = 1, 
            color_by = "Factor1")

plot_factor(model, 
            factors = c(1:3), 
            color_by = "Factor1")

plot_factor(model, 
            factors = c(1:7), 
            color_by = "dose", 
            #scale = TRUE
)

plot_factor(model, 
            factors = c(1:7), 
            color_by = "drug", 
            #scale = TRUE
)

plot_factor(model, 
            factors = c(1:7), 
            color_by = "condition4", 
            #scale = TRUE
)

plot_factor(model, 
            factor = 2:3,
            color_by = "time",
            shape_by = "condition3")

# samples weight per group in one factor, side-by-side ---------------------------------------
plot_factor(model, factors = 1, color_by = "Factor1")
plot_factor(model, factors = 1, color_by = "drug")

plot_factor(model, factors=1, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=2, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=3, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=4, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=5, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=6, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

plot_factor(model, factors=7, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


plot_factor(model, factors=2, group_by = "dose") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))

# plot_factor(model, factors=1, group_by = "drug", color_by="ENSG00000160868") +
#   theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)) # for Ensemble ID data

plot_factor(model, factors=1, group_by = "drug", color_by="ATP1A3") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))# for gene symnbol data

# selected plot -------------------------

samples_factor1 <- plot_factor(model, factors=1, group_by = "drug") +
  theme(axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1))


#between 2 factors ------------------------------------------------------------
plot_factors(model, 
             factors = 1:2,
             color_by = "drug")

plot_factors(model, 
             factors = 1:2,
             color_by = "condition3")

factors_5FU_part1  <- plot_factors(model, 
             factors = 2:3,
             color_by = "condition3")
factors_5FU_part1 

factors_5FU_part2 <- plot_factors(model, 
             factors = c(2,4),
             color_by = "condition3")
factors_5FU_part2

plot_factors(model, 
             factors = c(3,5),
             color_by = "condition3")

plot_factors(model, 
             factors = c(2,5),
             color_by = "condition3")

plot_factors(model, 
             factors = c(2,7),
             color_by = "condition3")

plot_factors(model, 
             factors = 3:4,
             color_by = "condition3")

plot_factors(model, 
             factors = 3:4,
             color_by = "dose")

p <- plot_factors(model, 
                  factors = c(2,5), 
                  color_by = "time",
                  shape_by = "dose",
                  dot_size = 2.5,
                  show_missing = T)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

# Save plots ------------------------------------------------------------
save_plot_png(samples_factor1, 
              filename = paste0("output/04_", nameDat, "_MofaLatentFactor1.png"), 
              width = 95, height = 80)

save_plot_png(factors_5FU_part1, 
              filename = paste0("output/04_", nameDat, "_MofaLatentFactor_5FUpart1.png"), 
              width = 125, height = 80)

save_plot_png(factors_5FU_part2, 
              filename = paste0("output/04_", nameDat, "_MofaLatentFactor_5FUpart2.png"), 
              width = 125, height = 80)

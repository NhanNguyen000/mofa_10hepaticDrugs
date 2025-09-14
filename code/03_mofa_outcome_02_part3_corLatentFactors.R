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

# Overview of data and metadata -------------------------------------------------
plot_data_overview(model)
sum(model@dimensions$N) # Nsamples

# metadata
head(model@samples_metadata)

load("processedDat/inputDat.RData")
model@samples_metadata <- model@samples_metadata %>% 
  full_join(expDesignTab, by = c("sample" = "name"))
#View(model@samples_metadata)

rm(dat, expDesignTab)

# Visualisation of combinations of factors ------------------------------------
set.seed(42)
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

vars <- colnames(model@samples_metadata)[3:9]

plist_factors <- list()
for(varName in vars) {
  plist_factors[[varName]] <- plot_factors(model, factors = 1:7, color_by = varName)
}

plist_factors$dose # have some differences
plist_factors$drug # have some differences
plist_factors$time # not much differences (?)
plist_factors$condition2 # not much differences (?)
plist_factors$condition3
plist_factors$condition4 # combination of factor 1, with 4, 5, 6, 7

plist_factors$condition4
# save plots ------------------------------------------------------------
save_plot_png(plist_factors$condition4, 
              filename = paste0("output/03_", nameDat, "_00_corLatentFactors.png"), 
              width = 260, height = 180)

# Other plots ------------------------------------------------------------
# Correlation between factors ---------------------------------------------------------
plot_factor_cor(model)

# which Factors have a clear association with any of the covariates (but use pearson cor) -> need to modify calculation
corFactors_covariates <- correlate_factors_with_covariates(
  model, 
  covariates = c("time","drug","dose", "condition2", "condition3"), 
  plot="log_pval")
corFactors_covariates

# manual calculation

factors <- get_factors(model, factors = "all")
dat <- factors$group1 %>% as.data.frame() %>%
  rownames_to_column("sample") %>%
  full_join(model@samples_metadata)

levels(dat$drug)
levels(dat$dose)
levels(dat$time)
levels(dat$condition)
levels(dat$condition2)
levels(dat$condition3)

res.aov <- aov(Factor1 ~ condition3, data = dat)
res.aov
summary(res.aov)
TukeyHSD(res.aov)

a <- TukeyHSD(res.aov)
a$condition3
# all factors ~ condition 3
res.aov_list <- list()
res.aov_list[["Factor1"]] <- aov(Factor1 ~ condition3, data = dat)
res.aov_list[["Factor2"]] <- aov(Factor2 ~ condition3, data = dat)
res.aov_list[["Factor3"]] <- aov(Factor3 ~ condition3, data = dat)
res.aov_list[["Factor4"]] <- aov(Factor4 ~ condition3, data = dat)
res.aov_list[["Factor5"]] <- aov(Factor5 ~ condition3, data = dat)
res.aov_list[["Factor6"]] <- aov(Factor6 ~ condition3, data = dat)
res.aov_list[["Factor7"]] <- aov(Factor7 ~ condition3, data = dat)

summary(res.aov_list$Factor1)
summary(res.aov_list$Factor2)
summary(res.aov_list$Factor3)
summary(res.aov_list$Factor4)
summary(res.aov_list$Factor5)
summary(res.aov_list$Factor6)
summary(res.aov_list$Factor7)

TukeyHSD_list <- list()
for (factor_temp in names(res.aov_list)) {
  a <- (TukeyHSD(res.aov_list[[factor_temp]]))
  TukeyHSD_list[[factor_temp]] <- a$condition3
}

TukeyHSD_list$Factor1

TukeyHSD_dat <- TukeyHSD_list %>% 
  lapply(function(x) x %>% 
           as.data.frame() %>% rownames_to_column("comparision")) %>%
  bind_rows(.id = 'factor')
range(TukeyHSD_dat$`p adj`)


TukeyHSD_dat %>% ggplot(aes(x = factor, y = comparision))
TukeyHSD_dat %>%
  ggplot(aes(x = factor, y = comparision, size = log10(`p adj`))) +
  geom_point() +
  scale_size_continuous( trans = "reverse", # Smaller p-values → bigger dots
                         # range = c(0, 10), 
                         name = "Adjusted log10(p-value)") +  # Reverse size so small p-values = larger dots
  theme_minimal() +
  labs(x = "Factor", y = "Comparison", title = "Dot plot of p-values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + plot_theme()

unique(TukeyHSD_dat$comparision)
[1] "ConDMSO-ConUNTR"             "5FU-ConUNTR"                 "AZA-ConUNTR"                
[4] "PHE-ConUNTR"                 "CYC_VPA_ISO-ConUNTR"         "DIC_RIF_MTX_APA-ConUNTR"    
[7] "5FU-ConDMSO"                 "AZA-ConDMSO"                 "PHE-ConDMSO"                
[10] "CYC_VPA_ISO-ConDMSO"         "DIC_RIF_MTX_APA-ConDMSO"     "AZA-5FU"                    
[13] "PHE-5FU"                     "CYC_VPA_ISO-5FU"             "DIC_RIF_MTX_APA-5FU"        
[16] "PHE-AZA"                     "CYC_VPA_ISO-AZA"             "DIC_RIF_MTX_APA-AZA"        
[19] "CYC_VPA_ISO-PHE"             "DIC_RIF_MTX_APA-PHE"         "DIC_RIF_MTX_APA-CYC_VPA_ISO"

interested_group <- c("ConDMSO-ConUNTR", 
                      "5FU-ConUNTR", "AZA-ConUNTR", "PHE-ConUNTR",
                      "5FU-ConDMSO", "AZA-ConDMSO", "PHE-ConDMSO",
                      "CYC_VPA_ISO-ConUNTR", "DIC_RIF_MTX_APA-ConUNTR",
                      "CYC_VPA_ISO-ConDMSO", "DIC_RIF_MTX_APA-ConDMSO",
                      "DIC_RIF_MTX_APA-CYC_VPA_ISO")
b <- TukeyHSD_dat %>% 
  slice(which(comparision %in% interested_group)) %>%
  mutate(comparision = factor(comparision, levels = interested_group))

b %>%
  ggplot(aes(x = factor, y = comparision, size = log10(`p adj`))) +
  #ggplot(aes(x = factor, y = comparision, size = log10_padj)) +
  geom_point() +
  scale_size_continuous( trans = "reverse", # Smaller p-values → bigger dots
                         # range = c(0, 10), 
                         name = "Adjusted log10(p-value)") +  # Reverse size so small p-values = larger dots
  theme_minimal() +
  labs(x = "Factor", y = "Comparison", title = "Dot plot of p-values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + plot_theme()


a<- TukeyHSD_dat %>%
  mutate(log10_padj = log10(`p adj`)) # Take -log10 so small p-values become big positive values
range(a$log10_padj) # [1] -1.081474e+01 -3.943840e-06

## # all factors ~ condition 4 ------------------
res.aov_list <- list()
res.aov_list[["Factor1"]] <- aov(Factor1 ~ condition4, data = dat)
res.aov_list[["Factor2"]] <- aov(Factor2 ~ condition4, data = dat)
res.aov_list[["Factor3"]] <- aov(Factor3 ~ condition4, data = dat)
res.aov_list[["Factor4"]] <- aov(Factor4 ~ condition4, data = dat)
res.aov_list[["Factor5"]] <- aov(Factor5 ~ condition4, data = dat)
res.aov_list[["Factor6"]] <- aov(Factor6 ~ condition4, data = dat)
res.aov_list[["Factor7"]] <- aov(Factor7 ~ condition4, data = dat)

# summary(res.aov_list$Factor1)
# summary(res.aov_list$Factor2)
# summary(res.aov_list$Factor3)
# summary(res.aov_list$Factor4)
# summary(res.aov_list$Factor5)
# summary(res.aov_list$Factor6)
# summary(res.aov_list$Factor7)

TukeyHSD_list <- list()
for (factor_temp in names(res.aov_list)) {
  a <- (TukeyHSD(res.aov_list[[factor_temp]]))
  TukeyHSD_list[[factor_temp]] <- a$condition4
}

TukeyHSD_list$Factor1

TukeyHSD_dat <- TukeyHSD_list %>% 
  lapply(function(x) x %>% 
           as.data.frame() %>% rownames_to_column("comparision")) %>%
  bind_rows(.id = 'factor')
range(TukeyHSD_dat$`p adj`)


TukeyHSD_dat %>% ggplot(aes(x = factor, y = comparision))
TukeyHSD_dat %>%
  ggplot(aes(x = factor, y = comparision, size = log10(`p adj`))) +
  geom_point() +
  scale_size_continuous( trans = "reverse", # Smaller p-values → bigger dots
                         # range = c(0, 10), 
                         name = "Adjusted log10(p-value)") +  # Reverse size so small p-values = larger dots
  theme_minimal() +
  labs(x = "Factor", y = "Comparison", title = "Dot plot of p-values") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + plot_theme()

unique(TukeyHSD_dat$comparision)
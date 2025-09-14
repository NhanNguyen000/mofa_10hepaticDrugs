# set up enviroment
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(DESeq2)
library(org.Hs.eg.db)
# load files -----------------------------------------------------
rawDat_temp <- list()

liver1 <- list.files("./data/HeCaToS_liver1/")
for (fileName in liver1) {
  rawDat_temp[[fileName]] <- read_tsv(paste0("./data/HeCaToS_liver1/", fileName))
}

liver2 <- list.files("./data/HeCaToS_liver_2/")
for (fileName in liver2) {
  rawDat_temp[[fileName]] <- read_tsv(paste0("./data/HeCaToS_liver_2/", fileName))
}

rawDat <- rawDat_temp %>% purrr::reduce(full_join) %>% 
  column_to_rownames("...1") 
dim(rawDat) # [1] 60664   474

# VST normalization using DESeq2-----------------------------------------------------
# use of the concept of variance stabilizing transformations (VST) in DESeq2
# produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors.

expDesignTab <- colnames(rawDat) %>% 
  as_data_frame() %>%
  rename("value" = "name") %>% 
  mutate(drug = substring(name, 1, 3),
         dose = substring(name, 5, 7),
         condition = substring(name, 1, 7),
         time = substring(name, 9, 11)) %>%
  mutate(dose = gsub("NTR", "UNTR", dose)) %>%
  mutate(dose = gsub("MSO", "DMSO", dose)) %>%
  mutate(dose = factor(dose, level = c("UNTR", "DMSO", "The", "Tox" ))) %>%
  mutate(time = factor(time, levels = c("000", "002", "008", "024", "072", "168", "240", "336")))

dds <- DESeqDataSetFromMatrix(countData = round(rawDat),
                              colData = expDesignTab,
                              design= ~ condition + time)

vsd <- vst(dds, blind=FALSE)
norm_rna <- assay(vsd)

# convert ENSEMBL id to gene name ----------------------------
symbols <- mapIds(org.Hs.eg.db, keys = rownames(norm_rna),
                  column = c('SYMBOL'), keytype = 'ENSEMBL') %>% 
  as.data.frame() %>% rename("." = "symbol") %>% rownames_to_column("ensembl")

norm_rna_v2_temp <- norm_rna %>% as.data.frame() %>% 
  rownames_to_column("ensembl") %>% 
  full_join(symbols) %>% 
  dplyr::select(-ensembl) %>%
  group_by(symbol) %>%
  summarise(across(1:474, mean, na.rm = TRUE))

which(is.na(norm_rna_v2_temp$symbol)) # 36614

norm_rna_v2 <- norm_rna_v2_temp %>% 
  dplyr::slice(-which(is.na(norm_rna_v2_temp$symbol))) %>% 
  column_to_rownames("symbol")

dim(norm_rna_v2) # [1] 36613   474, lost haft of the gene

# save data -----------------------------------------------------
dim(norm_rna) #  60664   474
dim(norm_rna_v2) #  36613   474

# normalized data
dat <- list()
dat[["normRNA_allEnsembleIds"]] <- norm_rna # all measured Ensemble IDs
dat[["normRNA_geneSymbols"]] <- norm_rna_v2 # only measured Ensemble IDs converted into gene symbol

# metadata
expDesignTab <- as.data.frame(expDesignTab) %>% 
  mutate(drug = gsub("Rif", "RIF", drug),
         condition = gsub("Rif", "RIF", condition)) %>% 
  mutate(#drug = ifelse(drug == "Con", "Control", drug),
    condition2 = ifelse(drug == "Con", condition, drug)) %>%
  mutate(
    drug = factor(drug, levels = c("Con", 
                                   "5FU", "APA", "AZA", "CYC", "DIC", 
                                   "ISO", "MTX", "PHE", "RIF", "VPA")),
    condition2 = factor(condition2, 
                        levels = c("ConDMSO", "ConUNTR", 
                                   "5FU", "APA", "AZA", "CYC", "DIC",
                                   "ISO", "MTX", "PHE", "RIF", "VPA" )))
rownames(expDesignTab) <- expDesignTab$name

# save
save(dat, expDesignTab, file = "processedDat/inputDat.RData")

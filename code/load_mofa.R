# set up enviroment
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("MOFA2")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("rhdf5")

library(MOFA2)
#library(rhdf5)
library(tidyverse)
library("viridis")

# load files -----------------------------------------------------
getwd()

# h5ls("./Mofa models/01.10.24/mofa_model0110.hdf5")
# 
# mydata <- h5read("./Mofa models/01.10.24/mofa_model0110.hdf5", "/data")
# 
# mydata2 <- h5read("./Mofa models/01.10.24/mofa_model0110.hdf5", "/groups")
# 
# mydata3 <- h5read("./Mofa models/01.10.24/mofa_model0110.hdf5", "/variance_explained")

# use MOFA

#model <- load_model("./Mofa models/04.12.24/MOFA_model.hdf5")
model <- load_model("/Users/nhannguyen/Documents/02_work_VNU/project_withTung/mofa_10hepaticDrugs/data/04.12.24/MOFA_model.hdf5")

plot_data_overview(model)
head(model@samples_metadata)
model@samples_metadata <- model@samples_metadata %>% 
  mutate(name = gsub("X5FU", "5FU", sample)) %>%
  mutate(drug = substring(name, 1, 3),
         dose = substring(name, 5, 7),
         time = substring(name, 9, 11)) %>%
  mutate(dose = gsub("NTR", "UNTR", dose)) %>%
  mutate(dose = gsub("MSO", "DMSO", dose)) %>%
  mutate(dose = factor(dose, level = c("UNTR", "DMSO", "The", "Tox" )))

head(model@samples_metadata)
head(model@cache$variance_explained$r2_total[[1]])
head(model@cache$variance_explained$r2_per_factor[[1]])

plot_variance_explained(model, x="view", y="factor")
plot_variance_explained(model, x="group", y="factor", plot_total = T)[[2]]

plot_factor(model, 
            factor = 1:2,
            color_by = "time",
            shape_by = "dose"
)

p <- plot_factor(model, 
                 factors = c(1:3),
                 color_by = "dose",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
color_values <- viridis(4)
p <- p + 
  scale_color_manual(values= color_values) +
  scale_fill_manual(values=color_values)  +
  theme(legend.position = "bottom")

print(p)

plot_factors(model, 
             factors = 1:2,
             color_by = "dose"
)

plot_factors(model, 
             factors = 1:2,
             color_by = "time"
)

plot_weights(model,
           #view = "Gene Expression",
             view = "View1",
             factor = 1,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)


plot_weights(model,
             #view = "Gene Expression",
             view = "View1",
             factor =2,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_data_heatmap(model,
                  view = "Gene Expression",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE
)

plot_data_scatter(model,
                  view = "Gene Expression",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "dose"
)

plot_factor_cor(model)
correlate_factors_with_covariates(model, 
                                  covariates = c("time","drug","dose"), 
                                  plot="log_pval"
)

# group factor together
p <- plot_factors(model, 
                  factors = c(5,7), 
                  color_by = "dose",
                  shape_by = "time",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

plot_factors(model, 
             factors = c(1,4),
             color_by = "drug")

plot_factors(mofa, factors=1:7, color_by = "drug")

# other plot
plot_factor(model, factors=2, group_by = "drug") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )

plot_factor(model, factors=1, group_by = "drug", color_by="ENSG00000160868") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )


plot_factor(model, factors=1, group_by = "drug", color_by="ENSG00000114416") +
  theme(
    axis.text.x = element_text(color="black", angle=40, vjust=1, hjust=1)
  )

# other plot 
factors <- c(1:7)
mofa <- run_umap(model, 
                 factors = factors, 
                 n_neighbors = 15,  
                 min_dist = 0.30
)

plot_dimred(mofa, 
            method = "UMAP", 
            color_by = "drug", 
            label = TRUE, 
            stroke=0.05, 
            dot_size = 1, 
            legend = FALSE) #+ scale_fill_manual(values=colors)

# We can try to add some interpretatibility on the UMAP by visualising the contribution of each Factor on the different groups of cells.

for (i in paste0("Factor",1:3)) {
  p <- plot_dimred(mofa, 
                   method = "UMAP", 
                   color_by = i, 
                   stroke = 0.05, 
                   dot_size = 1
  )
  print(p)
}

# check data ---------
View(model@data$View1$group1)
View(model@intercepts$View1$group1)
View(model@samples_metadata)

#     expected values of the different variables of the model. A list of matrices, one per variable. The most relevant are "W" for weights and "Z" for factors.

View(model@expectations$W$View1)
View(model@expectations$Z$View1)
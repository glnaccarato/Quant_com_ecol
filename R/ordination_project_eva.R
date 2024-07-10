..........................................
# TASK 4: ORDINATION (PART 2)
..........................................


# Packages ----
library(tidyverse)
library(ggplot2)
library(vegan)


# Datasets ----
taxa_dat <- read_csv("../data/JenaExp_arthropod_taxa.csv")
treatments <- read_csv("../data/JenaExp_treatments.csv")


# Data Exploration ----
# checking for nans
anyNA(taxa_dat) # no none
anyNA(treatments) # no none

# checking for duplicates
sum(duplicated(taxa_dat)) # none
sum(duplicated(treatments)) #none

# min and max values taxa_dat
for (i in colnames(taxa_dat)[-1]) { # iterate on all but the plotcode column
  max_value <- max(taxa_dat[[i]])
  min_value <- min(taxa_dat[[i]])
  message_min <- sprintf("%s min: %s", i, min_value)
  print(message_min)
  message_max <- sprintf("%s max: %s", i, max_value)
  print(message_max)
}
# min is 0 (next is 2), max is 1917


# Transformation-based ordination ----

# as the values of species richness are up to 3 orders of magnitude apart, data transformation is most likely preferable
# the taxon with the highest species richness is to no surprise Coleoptera. If we do not scale, they will of course be weighted disproportionately
# however pca with raw data will be used as a first approach and as a reference point

# NOTE: Carabidae are a family of Coleoptera and that is probably why we have Coleoptera_rest as a group. We should probably add them in a new column and use that for our ordination


# creation of dataframe with merged coleoptera_rest + carabidae column
taxa_dat_new <- data.frame(taxa_dat)
taxa_dat_new$Coleoptera <- taxa_dat_new$Coleoptera_rest + taxa_dat_new$Carabidae
taxa_dat_new <- taxa_dat_new %>%
  select(-Coleoptera_rest, -Carabidae)

for (i in colnames(taxa_dat_new)[-1]) { # iterate on all but the plotcode column
  max_value <- max(taxa_dat_new[[i]])
  min_value <- min(taxa_dat_new[[i]])
  message_min <- sprintf("%s min: %s", i, min_value)
  print(message_min)
  message_max <- sprintf("%s max: %s", i, max_value)
  print(message_max)
}
# min 0 (next 2), max 2335


## 1. PCA on raw data ----
# from lucas' script
pca_sp_raw <- rda(taxa_dat[,-1], 
              scale = FALSE)

pca_sp_raw

screeplot(pca_sp_raw, bstick = TRUE, type = "l",
          main = NULL) # only PCA1 lay above the line!

# ..and 
summary(eigenvals(pca_sp_raw)) # our eigenvalues don't make to much sense

# they look fucking weird as well
ordiplot (pca_sp_raw, scaling = "symmetric") 
biplot(pca_sp_raw, scaling = "symmetric")


## 2. PCA on transformed data ----
taxa_hell <- decostand(taxa_dat[-1], method = "hellinger")

pca_sp_trans <- rda(taxa_hell[,-1])

pca_sp_trans

screeplot(pca_sp_trans, bstick = TRUE, type = "l",
          main = NULL) # PC1-3 above the line

summary(eigenvals(pca_sp_trans)) # Eigenvalues: PC1 = 0.020736, PC2 = 0.017359, Cumulative Proportion of PC1 + PC2 = 0.61758 (Coleoptera_rest & Carabidae)

# !!! Coleptera and Carabidae are PC1 and PC2 !!!
ordiplot (pca_sp_trans, scaling = "symmetric") 
biplot(pca_sp_trans, scaling = "symmetric")


## 2. PCA on transformed data REVISED ----
# using dataframe with coleoptera all in one column
taxa_new_hell <- decostand(taxa_dat_new[-1], method = "hellinger")

pca_sp_new <- rda(taxa_new_hell[,-1])

pca_sp_new

screeplot(pca_sp_new, bstick = TRUE, type = "l",
          main = NULL) # PC1-3 above the line

summary(eigenvals(pca_sp_new)) # Eigenvalues: PC1 = 0.019178, PC2 = 0.010799, Cumulative Proportion of PC1 + PC2 = 0.6790 (Coleoptera & Isopoda?)

# Seems better I think
ordiplot (pca_sp_new, scaling = "symmetric") 
biplot(pca_sp_new, scaling = "symmetric")



































































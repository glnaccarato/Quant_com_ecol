# Task 2 - Data analysis --------

## Preparation -----------------------------------------------------------------

# libraries needed
library(tidyverse)
library(vegan)
library(ggplot2)

# set working directory (only for my personal usage)
setwd("/Users/gianlucanaccarato/Documents/Dokumente/GitHub/Quant_com_ecol")

# import species data
spec_dat <- read.csv("data/JenaExp_arthropod_taxa.csv")

# import environmental data
env_dat <- read.csv("data/JenaExp_treatments.csv")

## Data exploration ------------------------------------------------------------

# check both data sets for NA‘s
anyNA(spec_dat) # no NA's
anyNA(env_dat) # no NA's

# checking for duplicates
sum(duplicated(spec_dat)) # none
sum(duplicated(env_dat)) #none

# graphical data exploration (abundances of different arthropod taxa)
par(mar = c(3, 7, 2, 2)) # adjusting margins
boxplot(spec_dat[, 2:12],
        horizontal = TRUE,
        las = 2,
        main = "Abundance") # large differences in abundances

# apply Hellinger-transformation
spec_dat_hell <- decostand(spec_dat[, -1], method = "hellinger")

# double check if down-weighting of dominant taxa worked
boxplot(spec_dat_hell,
        horizontal = TRUE,
        las = 2,
        main = "Abundance") # better

# double check gradient length of first DCA axis (only optional)
decorana(spec_dat_hell) # 0.875 < 3 SD = linear methods are applicable

# PCA on Hellinger-transformed data
pca_spec <- rda(spec_dat_hell, 
                scale = FALSE) # species data is already on the same scale

# Check relevance of PC axis
screeplot(pca_spec, bstick = TRUE, type = "l",
          main = NULL) # PC1:3 above the line 

# Check eigenvalues of each PC
summary(eigenvals(pca_spec)) # cumulative proportion of PC1 + PC2 = 58.44%

# PCA ordination plot
ordiplot (pca_spec, 
          scaling = "symmetric") # species and site scores are scaled symmetrically by square root of eigenvalues

# PCA biplot
biplot(pca_spec, scaling = "symmetric")

# I would go with this approach from here on ......

# extract species and sites scores for further plotting


# little outlook:


# fit environmental vectors 

set.seed(2)
fit <- envfit(pca_spec ~ ., data = env_dat[, -1], # environmental data
               choices = 1:2, # axes to plotted
               scaling = "symmetric",
               permutations = 1000)

fit


plot(pca_spec, display = "sites", type="none",
     scaling = "symmetric")
points(pca_spec, display = "sites",
       scaling = "symmetric")
plot(fit)














######## QUICK CONMPARIONS OF ALL APPROACHES TRIED #########################

# I also computed a PCA on raw data and on Hellinger transformed data 
# with subsequent scaling of species data in the PCA function

pca_raw <- rda(spec_dat[, -1], scale = FALSE) # raw data without scaling
pca_scaled <- rda(spec_dat_hell) # hell-trans data with scaling

# comparison of screeplots

# raw data
screeplot(pca_raw, bstick = TRUE, type = "l",
          main = NULL) # only PC1 above the line (makes no sense)

# scaled data
screeplot(pca_scaled, bstick = TRUE, type = "l",
          main = NULL) # PC1:3 above the line 

# our data
screeplot(pca_spec, bstick = TRUE, type = "l",
          main = NULL) # PC1:3 above the line 

# comparison of ordiplot

# raw data
ordiplot (pca_raw, scaling = "symmetric") 
# scaled data
ordiplot (pca_scaled, scaling = "symmetric") 
# our data
ordiplot (pca_spec, scaling = "symmetric") 

# comparion of biplots

# raw data
biplot(pca_raw, scaling = "symmetric")
# scaled data
biplot(pca_scaled, scaling = "symmetric")
# our data
biplot(pca_spec, scaling = "symmetric")

########################## END ################################







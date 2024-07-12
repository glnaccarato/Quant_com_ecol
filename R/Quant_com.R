# Task 2 - Data analysis --------

## Preparation -----------------------------------------------------------------

# libraries needed
library(tidyverse)
library(vegan)
library(ggplot2)
library("viridis")  # only for color palette

# set working directory (only for my personal usage)
setwd("/Users/gianlucanaccarato/Documents/Dokumente/GitHub/Quant_com_ecol")

# import species data
spec_dat <- read.csv("data/JenaExp_arthropod_taxa.csv")

# import environmental data
env_dat <- read.csv("data/JenaExp_treatments.csv")

## Data Analysis ---------------------------------------------------------------

### Exploration of data --------------------------------------------------------

# check both data sets for NA‘s
anyNA(spec_dat) # no NA's
anyNA(env_dat) # no NA's

# checking for duplicates
sum(duplicated(spec_dat)) # none
sum(duplicated(env_dat)) #none

# overview over the study design

# how many samples per plant-species richness treatment
summary <- tibble(
  treatment = unique(env_dat$Plant_SR),
  counts = sapply(treatment, 
                  function(value) sum(env_dat$Plant_SR == value)))

summary[order(summary$treatment), ] # unbalanced sampling design

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

### Computation of Principal Component Analaysis -------------------------------

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

### Plotting using the ggplot --------------------------------------------------

# extract species scores
species.scrs <- scores(pca_spec, display=c("sp")) %>% 
  as_tibble(rownames = "species_name") 
species.scrs

# extract site scores
site.scrs <- scores(pca_spec, display = c("sites")) %>% 
  as_tibble(rownames = "sites")
site.scrs

# add plot_code to site.scrs table

# add site column from environmental data to PCA site scores
site.scrs <- cbind(site.scrs, Plant_SR = env_dat$Plant_SR)
site.scrs

# Plot 1: PCA ordination plot of sites/ plots
# with plots colored by plant species richness
ggplot(site.scrs)+
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3, alpha = 0.9)+
  scale_color_viridis(option = "A", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="black", lty =1) +
  geom_vline(xintercept = 0, color="black", lty =1) +
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species")) # manual legend


# Plot 2: PCA ordination plot of species
# with plots colored by plant species richness
ggplot(species.scrs)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1, PC2,label = species_name), 
                  color = "#CD0000", size = 3)

# Plot 3: PCA ordination plot of species and plots
# with plots colored by plant species richness
ggplot(species.scrs)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_point(data = site.scrs, aes(PC1, PC2, color = factor(Plant_SR)), 
             size = 3, alpha = 0.5)+
  scale_color_viridis(option = "A", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1, PC2,label = species_name), 
                  color = "#CD0000", size = 3)

### Post-hoc analysis of environmental variables -------------------------------

# set seed to obtain equal results each time of analysis
set.seed(2) 

# fit environemntal data onto pca ordination
fit <- envfit(pca_spec ~ ., data = env_dat[, -1], 
               choices = 1:2, # axes to plotted
               scaling = "symmetric",
               permutations = 1000)
fit

# extract only significant p-values
env.scrs_sig <- as_tibble(scores(fit, display = "vectors")) %>% # extract scores
  mutate(p.value = fit$vectors$pvals, # extract p-values
         Env.variables = colnames(env_dat[2:6])) %>% # get variable names
  subset(p.value<=0.05) # exclude non-significant variables 

env.scrs_sig


### Plotting of significant environmental variables ----------------------------

# Plot 1: Only vectors of significant variables 
ggplot(env.scrs_sig)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1, PC2, label = Env.variables), 
                  color = "#008B45", size = 3)

# Plot 2: Samples and siginificant variables
ggplot(site.scrs)+
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3)+
  scale_color_viridis(option = "A", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_segment(data = env.scrs_sig, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+
  geom_text_repel(data = env.scrs_sig, aes(PC1, PC2, label = Env.variables), 
                  color = "#008B45", size = 3)
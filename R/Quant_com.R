# Task 2 -----------------------------------------------------------------------

## Preparation -----------------------------------------------------------------

# libraries needed
library(tidyverse)
library(vegan)
library(ggplot2)
library(viridis)  # only for color palette
library(ggrepel)
library(psych) # for (visual) inspection of multicollinearity

# set working directory (only for my personal usage)
setwd("/Users/gianlucanaccarato/Documents/Dokumente/GitHub/Quant_com_ecol")

# import species data
spec_dat <- read.csv("data/JenaExp_arthropod_taxa.csv")

# import environmental data
env_dat <- read.csv("data/JenaExp_treatments.csv")

### Exploration of data --------------------------------------------------------

# check both data sets for NA‘s
anyNA(spec_dat) # no NA's
anyNA(env_dat) # no NA's

# checking for duplicates
sum(duplicated(spec_dat)) # none
sum(duplicated(env_dat)) #none

# Overview over the study design

# how many samples per plant-species richness treatment
summary <- tibble(
  treatment = unique(env_dat$Plant_SR),
  counts = sapply(treatment, 
                  function(value) sum(env_dat$Plant_SR == value)))

summary[order(summary$treatment), ] # unbalanced sampling design

# graphical data exploration (abundances of different arthropod taxa)
par(mar = c(3, 7, 2, 2))

boxplot(spec_dat[, 2:12],
        horizontal = TRUE,
        col = "#CD0000",
        las = 2,
        main = "Absolute Abundance") # large differences in abundances

# apply Hellinger-transformation
spec_dat_hell <- decostand(spec_dat[, -1], method = "hellinger")

# double check if down-weighting of dominant taxa worked
boxplot(spec_dat_hell,
        horizontal = TRUE,
        col = "#CD0000",
        las = 2,
        main = " Hellinger-Transformed Abundance") # better

par(mfrow = c(1,1))

### Computation of Principal Component Analysis -------------------------------

# double check gradient length of first DCA axis (only optional)
decorana(spec_dat_hell) # 0.875 < 3 SD = linear methods are applicable

# PCA on "psych"# PCA on Hellinger-transformed data
pca_spec <- rda(spec_dat_hell, 
                scale = FALSE) # species data is already on the same scale

# Check relevance of PC axis via broken stick plot
screeplot(pca_spec, bstick = TRUE, type = "l",
          main = NULL) # PC1:3 above the line

# Double check axis relevance using the Kaiser-Guttman criterion

# extract eigenvalues
spec_ev <- pca_spec$CA$eig

# display only axis greater than the mean
spec_ev > mean(spec_ev) # PC1:3 are relevant

# Check eigenvalues of each PC
summary(eigenvals(pca_spec)) # cumulative proportion of PC1 + PC2 = 58.44%

# extract species scores with correct scaling (Bocard et al. 2011)
species.scrs <- scores(pca_spec, display=c("sp"),
                      scaling = "symmetric") %>% 
  as_tibble(rownames = "species_name") %>% 
  mutate(acronyms = c("Ara", "Car", "Chi", "Cic", "Col", "Dip", 
                    "Het", "Iso","Opi", "Ort", "Sta")) # acronyms for plotting 
species.scrs

# extract site scores with correct scaling (Bocard et al. 2011)
site.scrs <- scores(pca_spec, display = c("sites"),
                    scaling = "symmetric") %>% 
  as_tibble(rownames = "sites")
site.scrs

# add site column from environmental data to PCA site scores
site.scrs <- cbind(site.scrs, Plant_SR = env_dat$Plant_SR)
site.scrs

### Post-hoc analysis of environmental variables -------------------------------

# set seed to obtain equal results each time of analysis
set.seed(2) 

# fit environemntal data onto pca ordination
fit <- envfit(pca_spec ~ ., data = env_dat[, -1], 
              choices = 1:2, 
              scaling = "symmetric",
              permutations = 1000)
fit

# extract only significant p-values
env.scrs_sig <- as_tibble(scores(fit, display = "vectors")) %>% # extract scores
  mutate(p.value = fit$vectors$pvals, # extract p-values
         Env.variables = colnames(env_dat[2:6]), # get variable names
         Env_renamed = c("PSR","Grasses","Short Herbs", # rename variables
                        "Tall Herbs","Legumes")) %>% 
  subset(p.value<=0.05) # exclude non-significant variables 

head(env_dat)

### Plotting PCA results using the ggplot --------------------------------------

# Figure 1: PCA ordination plot of plots and significant environmental 
# variables with plots colored by plant species richness 
ggplot(site.scrs)+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3)+
  geom_segment(data = env.scrs_sig, aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  scale_colour_brewer(type="qual", palette="RdBu")+
  theme_bw()+
  labs(x = "PC1 (33.69 %)", y = "PC2 (24.74 %) ")+
  guides(color = guide_legend(title = "No. of Plant Species",
                              position = "bottom",
                              direction = "horizontal", nrow = 1))+
  theme(legend.background = element_rect(color = "black", 
                                   fill = "white", 
                                   size = 0.1, 
                                   linetype = "solid"))+
  geom_text_repel(data = env.scrs_sig, aes(PC1, PC2*1.2, label = Env_renamed), 
                  color = "#008B45", size = 3)

ggsave("PCA_plots.png", width = 15, height = 15, units = "cm")

# Figure 2: PCA ordination plot of species
ggplot(species.scrs)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "PC1", y = "PC2") +
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1, PC2,label = acronyms), 
                  cex = 3.5, direction = "both", color = "#CD0000")

ggsave("PCA_species.png", width = 15, height = 12, units = "cm")


### Data visualization --------------------------------------------------------- 

spec_dat_rel <- as_tibble(
  spec_dat[, -1] / rowSums(spec_dat[, -1])) %>% # calculate rel. abundances
  mutate(Plant_SR = as.factor(env_dat$Plant_SR)) %>% # add factorized Plat_SR
  pivot_longer(cols = 1:11, # bring taxa columnns into wider format
               names_to = "Taxa",
               values_to = "Abundances") %>% 
  group_by(Taxa, Plant_SR) %>% 
  summarise(tot_rel_abundance = mean(Abundances))

# Figure 3: Percent stacked bar plot
ggplot(spec_dat_rel, aes(x=Plant_SR, y=tot_rel_abundance, fill = Taxa))+
  geom_bar(stat = "identity", position = "fill", # create percent stacked barchart
           colour="white") + 
  scale_fill_brewer(type="qual", palette="RdBu")+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "right")+
  labs(x = "Plant Species Richness", y = "Mean Relative Abundance")

ggsave("bp_stacked.png", width = 15, height = 15, units = "cm")


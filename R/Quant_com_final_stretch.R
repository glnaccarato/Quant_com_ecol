# Task 2 -----------------------------------------------------------------------

## Preparation -----------------------------------------------------------------

# libraries needed
library(tidyverse)
library(vegan)
library(ggplot2)
library(viridis)  # only for color palette
library(ggrepel)
library(psych) # for (visual) inspection of multicollinearity
library(ggcorrplot)

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # to automatically set working directory to source file location

# import species data
spec_dat <- read.csv("../data/JenaExp_arthropod_taxa.csv") # adding "../" at the start of the path makes it possible to access the data file although it's in another folder

# import environmental data
env_dat <- read.csv("../data/JenaExp_treatments.csv")


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
        main = "Transformed Abundance") # better

par(mfrow = c(1,1))

### Computation of Principal Component Analysis -------------------------------

# double check gradient length of first DCA axis (only optional)
decorana(spec_dat_hell) # 0.875 < 3 SD = linear methods are applicable

# check for multicollinearity among predictor variables 
pairs.panels(env_dat[, -1], method = "spearman")

# PCA on Hellinger-transformed data
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
  as_tibble(rownames = "species_name") 
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
         Env.variables = colnames(env_dat[2:6])) %>% # get variable names
  subset(p.value<=0.05) # exclude non-significant variables 

env.scrs_sig

### Plotting PCA results using the ggplot --------------------------------------

# Figure 1: PCA ordination plot of sites/plots and significant environmental 
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
  geom_text_repel(data = env.scrs_sig, aes(PC1, PC2*1.2, label = Env.variables), 
                  color = "#008B45", size = 3)

ggsave("PCA_plots.png", width = 15, height = 15, units = "cm")

# Figure 2: PCA ordination plot of species
ggplot(species.scrs)+
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "PC1",
       y = "PC2") +
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1*1.3, PC2*1.3,label = species_name), 
                  cex = 3.5, direction = "both",
                  color = "#CD0000")

ggsave("PCA_species.png", width = 15, height = 12, units = "cm")


### Data visualization --------------------------------------------------------- 

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



# Post-hoc - PERMANOVA ---------------------------------------------------------

# calculation of different distance measures using transformed data, in order for comparisons to be "fairer"
dis_euc_taxa<-vegdist(dis_bc_taxa,"euclidean")
dis_bc_taxa<-vegdist(spec_dat_hell,"bray")
dis_jacc_taxa<-vegdist(dis_bc_taxa,"jaccard", binary=T)
dis_chi_taxa<-vegdist(dis_bc_taxa,"chisq")

# Distance vector of plot species richness differences
dis_plant_sr<-dist(env_dat$Plant_SR, 'euclidean')

# Relationship between distance in species space and environmental distance
par(mfrow=c(2,2))
plot(dis_euc_taxa,dis_plant_sr,main="Euclidean vs. Plant_SR")
plot(dis_bc_taxa,dis_plant_sr, main="Bray-Curtis vs. Plant_SR")
plot(dis_jacc_taxa,dis_plant_sr, main="Jaccard vs. Plant_SR")
plot(dis_chi_taxa,dis_plant_sr,main="Chi-Square vs. Plant_SR")

## Pearson correlation between distances
cor(dis_euc_taxa,dis_plant_sr)
# -0.02200944 --> very weak negative correlation, taxa composition not linearly related to changes in plant species richness

cor(dis_bc_taxa,dis_plant_sr)
# -0.02638225 --> very weak negative correlation, taxa composition not linearly related to changes in plant species richness

cor(dis_jacc_taxa,dis_plant_sr)
# -0.1120482 --> weak negative correlation, inversely proportionate, areas with higher plant species richness have less dissimilarity in plant species presence-absence

cor(dis_chi_taxa,dis_plant_sr)
# -0.03364025 --> very weak negative correlation, taxa composition not linearly related to changes in plant species richness

# NOTE: altogether, the weak correlations indicate little to no linear relationship between dissimilarity of taxa and variation in Plant_SR
# this in itself indicates the need for attempting to use non-linear methods to assess this relationship


# Perform PERMANOVA to test how significant the effect of each treatment is
# use of bray-curtis simply because it is generally well suited for abundance data, use of jaccard proved fruitless
set.seed(1)
permanova_result <- adonis2(dis_bc_taxa ~ Plant_SR + Grasses + Short.Herbs + Tall.Herbs + Legumes, data = env_dat, permutations = 999)

# Print the results
print(permanova_result)
#             Df SumOfSqs      R2       F Pr(>F)    
# Plant_SR     1  0.06867 0.05340  5.4448  0.001 ***  # 5.34 % of the variation explained by Plant_SR, however the relationship is statistically significant (p < 0.01)
# Grasses      1  0.12997 0.10107 10.3059  0.001 ***  # 10 % of the variation attributed to the presence or absence of grasses (highest), relationship significant
# Short.Herbs  1  0.05938 0.04618  4.7086  0.001 ***  # 4.6 % of variation attributed to short herbs, relationship significant
# Tall.Herbs   1  0.03213 0.02498  2.5473  0.022 *    # 2.5 % of variation, relationship significant
# Legumes      1  0.06260 0.04868  4.9642  0.002 ***  # 4.9 % ofvariation, relationship significant
# Residual    74  0.93323 0.72570                     # 73 % of the data variation is derived from factors that are not our vegetation types.
# Total       79  1.28598 1.00000                     #

# PUT TABLE IN REPORT
# ENVFIT FOR THE VECTORS AND POST-HOC WITH PERMANOVA

set.seed(1)
perm_stat_taxa <- permustats(permanova_result)
summary(perm_stat_taxa)
#             statistic     SES    mean lower  median   upper Pr(perm)    
# Plant_SR       5.4448  7.5222  0.9949        0.8819  2.1034    0.001 ***
# Grasses       10.3059 15.7386  0.9870        0.8825  2.1406    0.001 ***
# Short.Herbs    4.7086  5.9830  1.0168        0.9144  2.2320    0.001 ***
# Tall.Herbs     2.5473  2.6166  0.9911        0.8792  2.0707    0.022 *  
# Legumes        4.9642  5.7796  1.0331        0.8868  2.3755    0.002 ** 

densityplot(perm_stat_taxa)

# MUST ADD INTERPRETATION OF DENSITY PLOTS FOR ALL!!! ---------------------

# Plotting correlation matrix among vegetation types and plant_sr

# Zitat from https://www.quora.com/Can-Spearmans-be-used-for-both-binary-and-categorical-data):
# "For binary data, Spearman's rank correlation coefficient may not be the most appropriate measure
#  because binary data are categorical and do not have a natural ordering. In such cases, other measures
#  like the Point-Biserial correlation coefficient or the Phi coefficient may be more suitable for assessing
#  the relationship between two binary variables.
#  For categorical data with more than two categories, the use of Spearman's rank correlation coefficient 
#  is not recommended because it does not take into account the categorical nature of the data. In such cases, 
#  other measures like Kendall's tau or Cramer's V may be more appropriate for assessing the association between
#  categorical variables."
corr <- cor(env_dat[-1], method = "kendall")
corr
#              Plant_SR     Grasses Short.Herbs  Tall.Herbs     Legumes
# Plant_SR    1.0000000  0.34828653  0.31933356  0.31660786  0.28745902
# Grasses     0.3482865  1.00000000  0.22215277  0.16881931 -0.05593966
# Short.Herbs 0.3193336  0.22215277  1.00000000 -0.05534631  0.22215277
# Tall.Herbs  0.3166079  0.16881931 -0.05534631  1.00000000  0.16881931
# Legumes     0.2874590 -0.05593966  0.22215277  0.16881931  1.00000000

ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
# highest correlations observed between Plant_SR and Grasses (0.35), Tall.Herbs (0.32), Short.Herbs (0.32) and Legumes (0.29)
# none of them are high

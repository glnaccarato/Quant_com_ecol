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
# library(DescTools) # for Cochran-Armitage Trend test (trends between ordinal and binary variables)
# library(MASS) # for ordinal logistic regression (to model relationship between an ordinal dependent variable and independent variables)
# library(vcdExtra) # Goodman and Kruskal's Gamma (assessment of strength of association between ordinal and binary variables)



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

# Testing homogeneity of multivariate dispersion to decide whether raw or transformed data is preferable
#  with raw data
distance_matrix_raw <- vegdist(spec_dat[-1], method = "bray")
dispersion_raw <- betadisper(distance_matrix_raw, as.factor(env_dat$Plant_SR))
anova(dispersion_raw)
# p < 0.05, indicating significant heterogeneity

#  with Hellinger-transformed data
distance_matrix_hell <- vegdist(spec_dat_hell, method = "bray")
dispersion_hell <- betadisper(distance_matrix_hell, as.factor(env_dat$Plant_SR))
anova(dispersion_hell)
# p > 0.08 (marginally), indicating some heterogeneity in dispersion
# moving on with transformed data

# Calculation of Bray-Curtis dissimilarity
dis_bc_taxa<-vegdist(spec_dat_hell,"bray")

# Distance vector of plot species richness differences
dis_plant_sr<-dist(env_dat$Plant_SR, 'euclidean')

# Relationship between distance in species space and environmental distance
plot(dis_bc_taxa,dis_plant_sr, main="Bray-Curtis vs. Plant_SR")

## Pearson correlation between distances
cor(dis_bc_taxa,dis_plant_sr)
# -0.02638225 --> very weak negative correlation, taxa composition not linearly related to changes in plant species richness

# Altogether, the weak correlations indicate little to no linear relationship between dissimilarity of taxa and variation in Plant_SR
# This pionts to the need for attempting to use non-linear methods to assess this relationship

# Perform PERMANOVA to test how significant the effect of each treatment is
set.seed(1)
permanova_result <- adonis2(dis_bc_taxa ~ Plant_SR + Grasses + Short.Herbs + Tall.Herbs + Legumes, data = env_dat, permutations = 999)

# Print the results
print(permanova_result)
# All environmental factors significant (p < 0.05), with grasses accounting for the highest percentage of variance explained (10,3 %)
# more than 73 % of variance not explained by our model, indicating unaccounted variables

set.seed(1)
perm_stat_taxa <- permustats(permanova_result)
summary(perm_stat_taxa)
#             statistic     SES    mean median   upper Pr(perm)    
# Plant_SR       5.4448  7.5222  0.9949 0.8819  2.1034    0.001 ***
# Grasses       10.3059 15.7386  0.9870 0.8825  2.1406    0.001 ***
# --> extremely low chance that the variance attributed to these variables is due to chance

densityplot(perm_stat_taxa)

 main

# # Plotting correlation matrix among vegetation types and plant_sr
# 
# # Zitat from https://www.quora.com/Can-Spearmans-be-used-for-both-binary-and-categorical-data:
# # "For binary data, Spearman's rank correlation coefficient may not be the most appropriate measure
# #  because binary data are categorical and do not have a natural ordering. In such cases, other measures
# #  like the Point-Biserial correlation coefficient or the Phi coefficient may be more suitable for assessing
# #  the relationship between two binary variables.
# #  For categorical data with more than two categories, the use of Spearman's rank correlation coefficient 
# #  is not recommended because it does not take into account the categorical nature of the data. In such cases, 
# #  other measures like Kendall's tau or Cramer's V may be more appropriate for assessing the association between
# #  categorical variables."
# corr <- cor(env_dat[-1], method = "kendall")
# corr
# 
# ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
# # highest correlations observed between Plant_SR and Grasses (0.35), Tall.Herbs (0.32), Short.Herbs (0.32) and Legumes (0.29)
# # none of them are high


# # Cochran-Armitage Trend Test
# for (var in c("Grasses", "Short.Herbs", "Tall.Herbs", "Legumes")) {
#   trend_test <- CochranArmitageTest(table(env_dat$Plant_SR, env_dat[[var]]))
#   print(paste("Cochran-Armitage Trend Test for", var))
#   print(trend_test)
# }
# 
# 
# # Ordinal Logistic Regression
# model <- polr(as.factor(Plant_SR) ~ + Grasses + Short.Herbs + Tall.Herbs + Legumes, data = env_dat, Hess = TRUE)
# summary(model)
# 
# 
# # Goodman and Kruskal's Gamma
# for (var in c("Grasses", "Short.Herbs", "Tall.Herbs", "Legumes")) {
#   gamma <- GKgamma(table(env_dat$Plant_SR, env_dat[[var]]))
#   print(paste("Goodman and Kruskal's Gamma for", var))
#   print(gamma)
# }



# Chi-Square Test with Monte Carlo Simulation (make approximation more robust)
for (var in c("Grasses", "Short.Herbs", "Tall.Herbs", "Legumes")) {
  table <- table(env_dat$Plant_SR, env_dat[[var]])
  chisq_test <- chisq.test(table, simulate.p.value = TRUE, B = 2000)  # B is the number of simulations
  print(paste("Chi-Square Test with Monte Carlo Simulation for", var))
  print(chisq_test)
}

# Plant_SR and Grasses:
# 
# Chi-Square Test: p-value = 0.01399 (significant)

# Plant_SR and Short.Herbs:
# Chi-Square Test: p-value = 0.02249 (significant)


# Plant_SR and Tall.Herbs:
# Chi-Square Test: p-value = 0.04498 (significant)

# Plant_SR and Legumes:
# Chi-Square Test: p-value = 0.07746 (not significant)

# Plant_SR has significant associations with Grasses, Short.Herbs, and 
# Tall.Herbs in all tests, indicating moderate positive associations. Legumes has a significant
# trend but not in the Chi-Square Test, suggesting a weaker relationship.


# Binary vs Binary
binary_vars <- c("Grasses", "Short.Herbs", "Tall.Herbs", "Legumes")
for (i in 1:(length(binary_vars)-1)) {
  for (j in (i+1):length(binary_vars)) {
    var1 <- binary_vars[i]
    var2 <- binary_vars[j]
    
    table <- table(env_dat[[var1]], env_dat[[var2]])
    chisq_test <- chisq.test(table)
    phi <- assocstats(table)$phi
    
    print(paste("Chi-Square Test and Phi Coefficient for", var1, "and", var2))
    print(chisq_test)
    print(paste("Phi Coefficient:", phi))
  }
}

# Grasses and Short.Herbs:
#   
# Chi-Square Test: X-squared = 3.1063, p-value = 0.07799
# Phi Coefficient: 0.222
# p-value > 0.05, no statistically significant association. Phi Coefficient suggests a weak positive association

# Grasses and Tall.Herbs:
#   
# Chi-Square Test: X-squared = 1.6502, p-value = 0.1989
# Phi Coefficient: 0.169
# p-value > 0.05, no statistically significant association. Phi Coefficient suggests a weak positive association

# Grasses and Legumes:
#   
# Chi-Square Test: X-squared = 0.075883, p-value = 0.783
# Phi Coefficient: 0.056
# p-value > 0.05, no significant association. Phi Coefficient suggests a very weak association.

# Short.Herbs and Tall.Herbs:
#   
# Chi-Square Test: X-squared = 0.07291, p-value = 0.7871
# Phi Coefficient: 0.055
# p-value > 0.1, no significant association. Phi Coefficient suggests a very weak association.

# Short.Herbs and Legumes:
#   
# Chi-Square Test: X-squared = 3.1063, p-value = 0.07799
# Phi Coefficient: 0.222
# p-value > 0.05, no statistically significant association. Phi Coefficient suggests a weak positive association.

# Tall.Herbs and Legumes:
#   
# Chi-Square Test: X-squared = 1.6502, p-value = 0.1989
# Phi Coefficient: 0.169
# p-value > 0.05, no statistically significant association. Phi Coefficient suggests a weak positive association.


# Binary vs. Binary: Most pairs of binary variables do not show significant associations (p-value > 0.05).
# The Phi Coefficients suggest weak to very weak associations.

# Overall: The results suggest that while there are moderate positive associations between Plant_SR and some of 
# the binary variables, the binary variables themselves are mostly weakly associated with each other.





 main

# Task 2 - Data analysis --------

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

par(mfrow = c(1,1))


### Computation of Principal Component Analaysis -------------------------------

# double check gradient length of first DCA axis (only optional)
decorana(spec_dat_hell) # 0.875 < 3 SD = linear methods are applicable

# check for multicollinearity among predictor variables 
pairs.panels(env_dat[, -1])

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
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3, alpha = 0.9)+
  # stat_ellipse(aes(x=PC1,y=PC2,fill=factor(Plant_SR)),
  #              geom="polygon", level=0.95, alpha=0.2)+
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  labs(
    x = "PC1 (33.69 %)",
    y = "PC2 (24.74 %) ",
    title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
    subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species")) + # manual legend
  guides(color=guide_legend("Plant_SR"), fill = "none") # remove factorized legend 

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
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "PC1 (33.69 %)",
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
  labs(x = "PC1 (33.69 %)",
       y = "PC2 (24.74 %) ",
       title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(PC1, PC2, label = Env.variables), 
                  color = "#008B45", size = 3)

# Plot 2: Samples and significant variables
ggplot(site.scrs)+
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3)+
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_segment(data = env.scrs_sig, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  theme_bw()+
  labs(x = "PC1 (33.69 %)",
       y = "PC2 (24.74 %) ",
       title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+
  geom_text_repel(data = env.scrs_sig, aes(PC1, PC2, label = Env.variables), 
                  color = "#008B45", size = 3)


### Data visualization --------------------------------------------------------- 


# I was also thinking of visualizing the arthropod composition along
# the plant species gradient via a stacked bar plot or a bubble plot

# we need relative abundances as sampling design was unbalanced
# absolute abundances are not comparable e.g. plots with 60 species
# only got sampled 4 times -> overall lower absolute abundances

# add column of plant species richness and bring into long format,
# turn absolute into relative abundances, and calculate mean values
spec_dat_rel <- as_tibble(
  spec_dat[, -1] / rowSums(spec_dat[, -1])) %>% # calculate rel. abundances
  mutate(Plant_SR = as.factor(env_dat$Plant_SR)) %>% # add factorized Plat_SR
  pivot_longer(cols = 1:11, # bring taxa columnns into wider format
               names_to = "Taxa",
               values_to = "Abundances") %>% 
  group_by(Taxa, Plant_SR) %>% 
  summarise(tot_rel_abundance = mean(Abundances))

view(spec_dat_rel)

# percent stacked barplot
ggplot(spec_dat_rel, aes(x=Plant_SR, y=tot_rel_abundance, fill = Taxa))+
  geom_bar(stat = "identity", position = "fill", # create percent stacked barchart
           colour="white") + 
  scale_fill_viridis(option = "C", discrete = T)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        legend.position = "right")+
  labs(x = "Plant Species Richness", y = "Mean Relative Abundance")

# bubble plot
ggplot(spec_dat_rel, aes(x=Plant_SR, y=Taxa))+
  geom_point(aes(color=Taxa, size = tot_rel_abundance))+
  scale_size_continuous(range = c(1, 20))+ # Adjust the range of points size
  scale_color_viridis(option = "C", discrete = T, alpha = 0.7) +
  labs(x= "Plant Species Richness",
       y = "Taxonomic Arthtropod Groups",
       size = "Relative Mean Abundance (%)", 
       fill = "")+
  theme(panel.background = element_blank(), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5), 
        legend.position = "right")+
  guides(color = F, size = F) # removes irreleavnt legend


### Different approach using a NMDS instead ------------------------------------

# Rank correlations between dissimilarity indices and gradient separation
rankindex(env_dat[, -1], spec_dat[, -1]) 

# bray-curtis does not fit the data best, but is generally recommanded when 
# dealing with quantitative data

set.seed(2) # set seed to obtain similar results each time
nmds1 <- metaMDS(spec_dat[, -1], distance = "bray", k = 2,
                 autotransform = F, trymax = 20)

# check stress
nmds1$stress # 0.1649247 (< 0.1 – excellent; < 0.2 – good --> could be better)

# checking the influence of dimensionality 
set.seed(2)
nmds2 <- metaMDS(spec_dat[, -1], distance = "bray", k = 3,
                 autotransform = F, trymax = 20)
nmds3 <- metaMDS(spec_dat[, -1], distance = "bray", k = 4,
                 autotransform = F, trymax = 20)
nmds4 <- metaMDS(spec_dat[, -1], distance = "bray", k = 5,
                 autotransform = F, trymax = 20)
nmds5 <- metaMDS(spec_dat[, -1], distance = "bray", k = 6,
                 autotransform = F, trymax = 20)

# Saving stress values 
stress_df <- tibble(k = 2:6,
                    stress = c(nmds1$stress, nmds2$stress, nmds3$stress, 
                               nmds4$stress, nmds5$stress))

# plot stress against k
ggplot(stress_df, aes(k, stress))+ 
  geom_point()+
  geom_line(lty = 1)+
  theme_bw()

# 3 dimensions seem to have a workable level of stress

nmds2$stress # 0.1012874 

# Asses appropriateness of NMDS ordination
stressplot(nmds2, main = "Shepard plot") # high linear fit (0.94)

# which sample units contribute to the overall variation?
gof = goodness(object = nmds2)
plot(nmds2, display = "sites", type = "none")
points(nmds2, display = "sites", cex = 2*gof/mean(gof))
# sample units shown with larger circles account for more of the 
# variation in the overall fit of the data 


# Plot results
# extract species scores
spec.scrs.nmds <- scores(nmds2, display=c("sp")) %>% 
  as_tibble(rownames = "species_name") 
spec.scrs.nmds

# extract site scores
site.scrs.nmds <- scores(nmds2, display = c("sites")) %>% 
  as_tibble(rownames = "sites")
site.scrs.nmds

# add plot_code to site.scrs table

# add site column from environmental data to PCA site scores
site.scrs.nmds <- cbind(site.scrs.nmds, Plant_SR = env_dat$Plant_SR)
site.scrs.nmds

# Plot 1: NMDS ordination plot of sites/ plots
# with plots colored by plant species richness
ggplot(site.scrs.nmds)+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_point(aes(NMDS1, NMDS2, color = factor(Plant_SR)), size = 3, alpha = 0.9)+
#  stat_ellipse(aes(x=NMDS1,y=NMDS2,fill=factor(Plant_SR)),
#              geom="polygon", level=0.6, alpha=0.2)+
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       title = "NMDS on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on Bray-Curtis Dissimilarity Matrix; Stress = 0.1012874")+
  guides(color = guide_legend(title = "No. of Plant Species")) + # manual legend
  guides(color=guide_legend("Plant_SR"), fill = "none") # remove factorized legend 

# Plot 2: NMDS ordination plot of species
# with plots colored by plant species richness
ggplot(spec.scrs.nmds)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       title = "NMDS on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on Bray-Curtis Dissimilarity Matrix; Stress = 0.1012874")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(NMDS1, NMDS2,label = species_name), 
                  color = "#CD0000", size = 3)

# Plot 3: NMDS ordination plot of species and plots
# with plots colored by plant species richness
ggplot(spec.scrs.nmds)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#CD0000")+
  geom_point(data = site.scrs.nmds, aes(NMDS1, NMDS2, color = factor(Plant_SR)), 
             size = 3, alpha = 0.5)+
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       title = "NMDS on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on Bray-Curtis Dissimilarity Matrix; Stress = 0.1012874")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(NMDS1, NMDS2,label = species_name), 
                  color = "#CD0000", size = 3)



# results are similar to the PCA

# fit environmental vectors

set.seed(4) # random number generator to obtain equal results each time

fit.nmds <- envfit(nmds2 ~ ., data = env_dat[, -1], 
                   choices = 1:2, # axes to plotted
                   scaling = "symmetric",
                   permutations = 1000)
fit.nmds

# extract only significant p-values
env.scrs.nmds_sig <- as_tibble(scores(fit.nmds, display = "vectors")) %>% # extract scores
  mutate(p.value = fit.nmds$vectors$pvals, # extract p-values
         Env.variables = colnames(env_dat[2:6])) %>% # get variable names
  subset(p.value<=0.05) # exclude non-significant variables 

env.scrs.nmds_sig


### Plotting of significant environmental variables ----------------------------

# Plot 1: Only vectors of significant variables 
ggplot(env.scrs.nmds_sig)+
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), # arrows for environmental variables starting at 0
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  theme_bw()+
  labs(x = "NMDS1",
       y = "NMDS2",
       title = "NMDS on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on Bray-Curtis Dissimilarity Matrix; Stress = 0.1012874")+
  guides(color = guide_legend(title = "No. of Plant Species"))+ # manual legend
  geom_text_repel(aes(NMDS1, NMDS2, label = Env.variables), 
                  color = "#008B45", size = 3)

# Plot 2: Samples and siginificant variables
ggplot(site.scrs)+
  geom_point(aes(PC1, PC2, color = factor(Plant_SR)), size = 3)+
  scale_color_viridis(option = "C", # color pallet with different options
                      discrete = TRUE)+ 
  scale_fill_viridis(discrete = TRUE) +
  geom_hline(yintercept = 0, color="grey", lty =1) +
  geom_vline(xintercept = 0, color="grey", lty =1) +
  geom_segment(data = env.scrs_sig, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), color = "#008B45")+
  theme_bw()+
  labs(x = "PC1 (33.69 %)",
       y = "PC2 (24.74 %) ",
       title = "PCA on Arthropod Composition along a Plant Species Richness Gradient",
       subtitle = "Performed on the Hellinger-transformed species matrix")+
  guides(color = guide_legend(title = "No. of Plant Species"))+
  geom_text_repel(data = env.scrs_sig, aes(PC1, PC2, label = Env.variables), 
                  color = "#008B45", size = 3)


# 3. Approach - Rudeundanca Analysis (RDA) -------------------------------------

# Check for multicollinearity

# Perform RDA on Hellinger-transformed data

# apply Hellinger-transformation
spec_dat_hell <- decostand(spec_dat[, -1], method = "hellinger")

# RDA
rda1 <- rda(spec_dat_hell ~ ., env_dat[, -1]) # all variables except Plot_ID

# RDA results
summary(rda1)

# check eigenvalues of RDA axis
summary(eigenvals(rda1))

ordiplot(rda1)

# Model evaluation 
gof = goodness(object = rda1)
plot(rda1, display = "sites", type = "none")
points(rda1, display = "sites", cex = 2*gof/mean(gof))
# larger circles indicate more contribution to variation

# choose most simplified model (principlle of parsimony)

# check for multicolinearity
vif.cca(rda1) # VIF 

RsquareAdj(rda1) # adjusted R-squared: 0.23

# test models for their significance
perm_stat <- permustats(anova(rda1))
summary(perm_stat) # the model is significant 






?step

test <- step(rda1, )








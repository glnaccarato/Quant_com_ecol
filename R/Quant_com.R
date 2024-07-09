# Overview over Task4 ------------------------------------------------------------------------

# Part1 – Biodiversity of ant communities along an elevational gradient --------

# Dataset: dat4_ants_elevation.xlsx
# Description of the data set:
# This data set includes samples of ant communities along an elevational 
# gradient. The sampling used a hierarchically nested design. The basic sampling 
# unit are 1m2 plots for which abundance counts for all species are available.
# # Sites ranged in elevation from 379-1828 m and were in relatively undisturbed 
# deciduous forest sites, away from trails. At each site, there is one randomly 
# placed 50 × 50 m plot, from which 16 1-m2 quadrats were arranged in a 
# nested design: 10 x 10 m subplots were placed in the corners of 
# each 50 x 50 m plot, and 1-m2 quadrats were placed in the corners of 
# each 10 x 10 m subplot, for a total of 16 1-m2 quadrats per site. Ants were 
# sampled by collecting all the leaflitter within each quadrat and sifting 
# through it with a coarse mesh screen (1-cm grid) to remove the largest 
# fragments and concentrate the fine litter. Concentrated litter from each 
# quadrat was then put in its own mini-Winkler sacks for 2 days in the lab. 
# After 2 days, all worker ants were extracted and enumerated.

# The main question is now, how ant diversity at different scales 
# (alpha, gamma, and beta) changes along the elevational gradient

# Tasks:
#   • How many samples (1m2 plots) were taken at each elevation?
#   • Provide appropriate measures of alpha (1m2) and gamma-diversity
#     (16 x 1 m2 plot) for each elevation.
#   • Apply standardized biodiversity comparisons for species richness 
#     along the elevational gradient using:
#        o Equalsamplingeffort(i.e.equalno.ofsamples)
#        o Equalno.ofindividuals(i.e.individual-basedrarefaction) 
#        o Equalcoverage(i.e.coverage-basedrarefaction)
#   • Calculate beta-diversity for each elevation and its partitioning into 
#     turnover and nestedness components.
#        o Usemultisitemeasureofbeta-diversityNOTpair-wisemeasures
# • Visualise the relationship of all biodiversity indices with elevation

# Part 2 – Taxonomic composition of arthropods along plant diversity gradient in the Jena Experiment ----

# Datasets:
#   o JenaExp_arthropod_taxa.csv o JenaExp_treatments.csv

# Description of the data set:
#   The Jena Experiment (https://the-jena-experiment.de) is the largest and the 
#   longest grassland biodiversity experiment in Europe. It is located in city 
#   Jena in Germany. On 80 experimental plots (each of 20 x 20 m in size) the 
#   number of plant species and the presence of four different plant functional 
#   groups (legumes, grasses, tall herbs and short herbs) were experimentally 
#   manipulated. In 2010 the arthropod communities were sampled on each plot.
#   The dataset called JenaExp_arthropod_taxa.csv includes community composition
#   for each plot (called “plotecode”) of arthropods with the arthropod taxa and
#   their abundances. The dataset called JenaExp_treatments.csv includes the 
#   treatment variables of the Jena Experiment:

#    “Plant_SR” – species richness (number of species) of plants sawn on 
#                 the plot: 1, 2, 4, 8, 16, 60 species
#    “Grasses” – presence of grasses (0-absent, 1-present); 
#               “Legumes” – presence of legumes (0-absent, 1-present); 
#               “Short.Herbs” – presence of short herbs (0-absent, 1-present); 
#               “Tall.Herbs” – presence of tall herbs (0-absent, 1-present);

# Main questions:
#   • Are there patterns in the taxonomic composition of arthropods along 
#     the plant diversity gradient?
#   • What are the effects of plant species richness and of the presence of each
#     plant functional groups on the taxonomic composition of arthropod communities?

# Tasks:
#   • Explore if there are patterns in the taxa composition of arthropods across
#     the plant diversity gradient. Create an ordination diagram of taxa and 
#     plots for the first two ordination axes. Colour the plots based on the 
#     level of plant species richness of each plot. Investigate any patterns 
#     in the taxonomic composition of arthropods.
#   • Post-hoc, analyse the effects of the plant-diversity treatments 
#     (i.e., plant species richness, presence of grasses, legumes, small herbs
#     and tall herbs) on the community composition of arthropods. Add vectors 
#     representing these effects to the ordination diagram.
#   • Describe the results of the ordination and of the post-hoc analysis:
#     o Identify the plant-diversity treatments with the strongest effect. 
#     o Determine which drivers have similar and different effects on community composition.
#     o Identify the arthropod taxa that occur most frequently and least 
#       frequently at high and low plant diversity levels.
#     o Identify the arthropod taxa that associate most with grasses.


# Task 2 - Data analysis --------

# libraries needed
library(vegan)
library(ggplot2)


## Preparation -----

# Import data

# set working directory (only for my personal usage)
setwd("/Users/gianlucanaccarato/Documents/Dokumente/GitHub/Quant_com_ecol")

# species data
spec_dat <- read.csv("data/JenaExp_arthropod_taxa.csv")

# environmental data
env_dat <- read.csv("data/JenaExp_treatments.csv")

# check both data sets for NA‘s
anyNA(spec_dat) # no NA's
anyNA(env_dat) # no NA's

# no zero values in species matrix indicates a homogeneous dataset

# Heterogenity of species data can be measured by DCA (the length of the 
# first DCA axis) = measure of heterogeneity;

# 1. Using the first Detrended Component Analysis (DCA) axis to detect,
#    if species response is rather linear or unimodal along the 
#    environmental gradient 

# perform DCA to check for first axis length
decorana(spec_dat[,-1]) # 1.5048 = short gradient 

# indicates a homogeneous dataset
#   - linear (< 3 SD) as well as distance-based methods are suitable, but 
#     linear methods are more powerful

# And here the confusion starts because:

# if I understood it correctly we do not need to use a hellinger transformation
# as our species data already response linear to the environmental gradient
# even though its abundance data (compare to oksanas final task script)

# but if we continue like this, that happens:

# 1. Approach:

# perform PCA on untransformed abundance data
pca_sp <- rda(spec_dat[,-1], 
              scale = FALSE) # species data does not need scaling

pca_sp

screeplot(pca_sp, bstick = TRUE, type = "l",
          main = NULL) # only PCA1 lay above the line!

# ..and 
summary(eigenvals(pca_sp)) # our eigenvalues don't make to much sense

# 2. Approach 

# would be to use a hellinger-transformation, but if we do, there is 
# no point in performing a DCA at the beginning, because we are using a 
# linear method on transformed data 

# Quote from Oksanas source:
# https://www.davidzeleny.net/anadat-r/doku.php/en:confusions

# Use transformation-based ordination method. You do not need to make a decision 
# between linear or unimodal methods, simply you opt for transformation-based 
# method. Transform the species composition data by Hellinger transformation 
# (there are also other transformations, but this is the most often used), 
# and analyse data either by unconstrained (tb-PCA) or constrained (tb-RDA) 
# approach, depending on your purpose. No DCA is done before, it is not 
# necessary, since you are not making decision between linear or unimodal 
# ordination method (you use linear method applied on Hellinger transformed data).

# so we would than justify our transformation only because we're dealing with
# abundance data. The output with transformed data makes also more sense

# hellinger transformation
spec_dat_hell <- decostand(spec_dat[, -1], method = "hellinger")

# PCA
pca_sp_hell <- rda(spec_dat_hell, 
              scale = FALSE) # no scaling needed, species data

pca_sp_hell

# check relevance of ordination axis
screeplot(pca_sp_hell, bstick = TRUE, type = "l",
          main = NULL) # now PCA1, PCA2, and PCA3 lay above the line 

# eigenvalues of axis now make more sense 
summary(eigenvals(pca_sp_hell))

# Explained proportions of the axis (should be included in the final plot)
# PCA1: 0.33690
# PCA2: 0.24746

# So should we go with Approach 2?

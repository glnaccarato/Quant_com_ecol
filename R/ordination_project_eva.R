#..........................................
# TASK 4: ORDINATION (PART 2)
#..........................................



# Packages ----
library(tidyverse)
library(ggplot2)
library(vegan)
library(car)
library(ggcorrplot)
library(ggvegan)
library(ggrepel)



# Datasets ----
taxa_dat <- read_csv("data/JenaExp_arthropod_taxa.csv")
treatments <- read_csv("data/JenaExp_treatments.csv")



#..........................................
# Data Exploration ----
#..........................................

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
rm(i, max_value, min_value, message_min, message_max)

# DCA to determine length of gradient in our data
decorana(taxa_dat[-1])
# < 3 --> linear ordination



#........................................
# Data Transformation (Hellinger) ----
#........................................

# as the values of species richness are up to 3 orders of magnitude apart, data transformation is most likely preferable
# the taxon with the highest species richness is to no surprise Coleoptera. If we do not scale, they will of course be weighted disproportionately
# however pca with raw data will be used as a first approach and as a reference point

# Transformation
community_matrix <- taxa_dat[,-1] %>%
  as.matrix() %>%
  decostand(method = "hellinger")



#............................
# PCA - Unconstrained ---- 
#............................

## 1. on raw data ----
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


## 2. on transformed data ----
set.seed(1)

pca_sp_trans <- rda(community_matrix[,-1])

pca_sp_trans

screeplot(pca_sp_trans, bstick = TRUE, type = "l",
          main = NULL) # PC1-3 above the line

summary(eigenvals(pca_sp_trans))
# Eigenvalues: PC1 = 0.020736, PC2 = 0.017359, Cumulative Proportion of PC1 + PC2 = 0.61758 (Coleoptera_rest & Carabidae)

## PCA scores/Plots of PC1 & PC2 ----
scores(pca_sp_trans, choices = 1:2, display = c("sites"))
#                PC1         PC2
# sit1  -0.186705378  0.06418547
# sit2   0.018085154 -0.13481872
# sit3   0.105135188  0.10917365
# sit4  -0.308882906  0.09115304
# sit5  -0.149175784 -0.30953832
# sit6   0.035866437 -0.12237521
# sit7   0.131679208 -0.27252811
# sit8   0.115704917 -0.19252997
# sit9   0.126334546 -0.17428549
# sit10 -0.040591980 -0.16133299
# sit11  0.042576790 -0.15282056
# sit12  0.061540544 -0.04658713
# sit13  0.017097470 -0.28706882
# sit14 -0.229384690 -0.06537081
# sit15  0.069072404 -0.21003616
# sit16  0.119236586 -0.29420839
# sit17  0.122382470  0.17196991
# sit18  0.006256004  0.03538089
# sit19  0.003435056 -0.10808886
# sit20 -0.013646202  0.03247214
# sit21  0.028892280 -0.01185861
# sit22  0.121673189  0.05435465
# sit23  0.044409649  0.14070324
# sit24 -0.051674864 -0.13143408
# sit25  0.193450592 -0.17530373
# sit26 -0.523077235  0.05274663
# sit27 -0.074228141 -0.10790953
# sit28 -0.270398054 -0.12410684
# sit29  0.069311622 -0.01534822
# sit30 -0.114420753 -0.07282007
# sit31 -0.426769096 -0.03492671
# sit32 -0.040211957  0.06861553
# sit33  0.167999572 -0.12216154
# sit34 -0.166862774 -0.04877407
# sit35 -0.205702569 -0.02317180
# sit36 -0.098419121  0.05454056
# sit37 -0.284167738  0.01252650
# sit38 -0.594338165 -0.08864262
# sit39 -0.104166100  0.06853572
# sit40  0.005049315  0.21472816
# sit41  0.048699270 -0.23006682
# sit42 -0.114508139 -0.22171473
# sit43 -0.008573942 -0.03978803
# sit44  0.153402175  0.17117782
# sit45  0.105831954  0.03620708
# sit46  0.199169147 -0.01930446
# sit47  0.111191250 -0.02170382
# sit48  0.197189306 -0.03903269
# sit49  0.152698821  0.03405989
# sit50  0.118671516  0.06998222
# sit51  0.145537220 -0.16932331
# sit52  0.171221291  0.17047749
# sit53  0.134432328  0.30573551
# sit54  0.204681371  0.08937474
# sit55  0.190439760 -0.24665039
# sit56  0.228034647 -0.18215850
# sit57 -0.011536688  0.02418312
# sit58  0.014604309 -0.28098541
# sit59  0.093661308 -0.02052391
# sit60  0.039840421 -0.08173128
# sit61  0.306192150  0.05462706
# sit62  0.098615912  0.22151507
# sit63 -0.059019369  0.23065833
# sit64 -0.124997816  0.31922625
# sit65  0.017983693  0.21561411
# sit66  0.067171269  0.22685385
# sit67 -0.104611827  0.08752550
# sit68  0.032907145  0.14846271
# sit69 -0.075932135  0.25602283
# sit70  0.170480792  0.38532807
# sit71  0.185446451  0.01131288
# sit72 -0.028472122 -0.24404289
# sit73 -0.310655841  0.07726064
# sit74 -0.054348898 -0.01199550
# sit75  0.127302321  0.18339174
# sit76 -0.041364997 -0.15733219
# sit77  0.055667203  0.30804911
# sit78 -0.132497222  0.31573200
# sit79  0.006660241  0.17086034
# sit80 -0.033579765  0.16967691

scores(pca_sp_trans, choices = 1:2, display = c("species"))
#                         PC1         PC2
# Carabidae        0.20870280 -0.51198667
# Chilopoda        0.03892654  0.01397785
# Cicadina         0.20202933  0.33059858
# Coleoptera_rest -0.80018082 -0.02044252
# Diplopoda        0.02554667 -0.02819436
# Heteroptera     -0.01144545  0.09112284
# Isopoda          0.05536583  0.48770884
# Opiliones        0.06675238 -0.02279022
# Orthoptera       0.05952568 -0.03417802
# Staphylinidae    0.06365365  0.02320699





#...............................
## PCA - Ordination Plots ----
#...............................

# Coleoptera_rest and Carabidae the main contributors as their vectors run the longest across the axes of PC1 and PC2 respectively
ordiplot (pca_sp_trans, scaling = "symmetric") 
biplot(pca_sp_trans, scaling = "symmetric")

set.seed(1)
ordipointlabel(pca_sp_trans,
               display = "sites",
               scaling = "symmetric")


set.seed(1)
plot(pca_sp_trans, display = "sites",
     scaling = "symmetric", type = "n")
points(pca_sp_trans, display = "sites",
        scaling = "symmetric", pch = 19,
        col = "orange")
ordipointlabel(pca_sp_trans,
                display = "sites",
                scaling = "symmetric",
                add = TRUE)


# Including categorical data
length(unique(treatments$Plant_SR))
# 6 as expected

### Colors based on plant_SR value ----
col_vec <- c("green", "pink", "blue", "magenta", "orange", "turquoise")

plot(pca_sp_trans, type = "n", scaling = "symmetric",
     display = "sites")
cols <- with(treatments, col_vec[as.factor(Plant_SR)]) # as.factor to not consider them numeric

points(pca_sp_trans, display = "sites", scaling = "symmetric",
       pch = 19, col = cols, cex = 2)

lvl <- levels(as.factor(treatments$Plant_SR))

legend("topright", legend = lvl,
       bty = "n", col = col_vec, pch = 19)


### Convex hulls around groups of data --> TOO MUCH GOING ON HERE? ----
plot(pca_sp_trans, type = "n", scaling = "symmetric",
     display = "sites")

ordihull(pca_sp_trans, groups = as.factor(treatments$Plant_SR),
         col = col_vec,
         scaling = "symmetric", lwd = 2)

ordispider(pca_sp_trans, groups = as.factor(treatments$Plant_SR),
           col = col_vec,
           scaling ="symmetric", label = TRUE)

points(pca_sp_trans, display = "sites", scaling = "symmetric",
       pch = 21, col = "black", bg ="gray")


### Elipsoid --> AGAIN TOO MUCH? ----

# first plot coordinates
plot(pca_sp_trans, type = "n", scaling = "symmetric",
     display = "sites")

## ellipsoid hull
ordiellipse(pca_sp_trans, groups = as.factor(treatments$Plant_SR),
            kind = "ehull", col = col_vec,
            scaling = "symmetric", lwd = 2)

## standard error of centroid ellipse (statard error calculated for the center)
ordiellipse(pca_sp_trans, groups = as.factor(treatments$Plant_SR),
            draw = "polygon", col = col_vec,
            scaling = "symmetric", lwd = 2)

ordispider(pca_sp_trans, groups = as.factor(treatments$Plant_SR),
           col = col_vec,
           scaling = "symmetric", label = TRUE)

points(pca_sp_trans, display = "sites", scaling = "symmetric",
       pch = 21, col = "black", bg ="gray")



## Data across environmental/diversity gradients ----

# Calculate Shannon diversity index
shannon_diversity <- diversity(taxa_dat[, -1], index = "shannon")

# Calculate Simpson diversity index
simpson_diversity <- diversity(taxa_dat[, -1], index = "simpson")

# Calculate Gamma Diversity (overall diversity)
# gamma diversity for the entire dataset can be calculated by combining all plots
gamma_diversity_shannon <- diversity(colSums(taxa_dat[, -1]), index = "shannon")
gamma_diversity_simpson <- diversity(colSums(taxa_dat[, -1]), index = "simpson")

# Calculate Beta Diversity
# using Bray-Curtis dissimilarity as a measure of beta diversity
beta_diversity <- vegdist(taxa_dat[, -1], method = "bray")

# Combine the site identifier and calculated diversity indices
site_diversity <- data.frame(plotcode = taxa_dat$plotcode, shannon_diversity = shannon_diversity, simpson_diversity = simpson_diversity)

# Print the results
print(site_diversity)

gamma_diversity_shannon
# 1.800359

gamma_diversity_simpson
# 0.797958

beta_diversity
# Calculate summary statistics for beta diversity
mean(beta_diversity) # average
# 0.329235778643104

min(beta_diversity) # minimum
# 0.0553459119496855

max(beta_diversity) # maximum
# 0.790218394340203


# # Visualize beta diversity with a heatmap
# heatmap(as.matrix(beta_diversity), main = "Heatmap of Bray-Curtis Dissimilarity", symm = TRUE)
# 
# # Alternatively, visualize with a dendrogram
# dendrogram <- hclust(beta_diversity, method = "average")
# plot(dendrogram, main = "Dendrogram of Bray-Curtis Dissimilarity", xlab = "", sub = "")



#.......................................................
# Visualizing data patterns across gradients ----
#.......................................................

# Ordisurf  with shannon gradient
par(mfrow=c(2,2))
fit_shannon <- ordisurf(pca_sp_trans ~ shannon_diversity,
                   data = site_diversity,
                   main = NULL)
title("Shannon")
# points(pca_sp_trans, display = "sites")
# text(pca_sp_trans, display = "sites", labels = site_diversity$plotcode, cex = 0.7, pos = 3)
summary(fit_shannon)


# rdisurf  with simpson gradient
fit_simpson <- ordisurf(pca_sp_trans ~ simpson_diversity,
                        data = site_diversity,
                        main = NULL)
title("Simpson")
# points(pca_sp_trans, display = "sites")
# text(pca_sp_trans, display = "sites", labels = site_diversity$plotcode, cex = 0.7, pos = 3)
summary(fit_simpson)

# plant_sr gradient
fit_sr <- ordisurf(pca_sp_trans ~ Plant_SR,
                 data = treatments,
                 main = NULL)
title("Plant_SR")
# points(pca_sp_trans, display = "sites")
# text(pca_sp_trans, display = "sites", labels = site_diversity$plotcode, cex = 0.7, pos = 3)
summary(fit_sr)


#  If these lines are aligned with the PCA axes, it suggests a strong relationship between species richness.


# Fit environmental variables (plant-diversity treatments) to the PCA
envfit_result <- envfit(pca_sp_trans ~ Plant_SR + Grasses + Legumes + Short.Herbs + Tall.Herbs, data = treatments)

# Add the environmental fit vectors to the PCA biplot
plot(pca_sp_trans, scaling = 2, main = "PCA with Environmental Vectors")
points(pca_sp_trans, display = "sites", col = treatments$Plant_SR, pch = 19)
plot(envfit_result, p.max = 0.05, col = "red")



# Plotting ordination diagram using ggplot  --> TOTALLY TOO MUCH BUT MIGHT BE NICE TO HAVE ALL CONCENTRATED ----
# Extract PCA site scores
site_scores <- scores(pca_sp_trans, display = "sites", scaling = 2)
site_scores_df <- as.data.frame(site_scores)
site_scores_df$Plant_SR <- as.factor(treatments$Plant_SR)
site_scores_df$plotcode <- treatments$plotcode

# Extract PCA species scores
species_scores <- scores(pca_sp_trans, display = "species", scaling = 2)
species_scores_df <- as.data.frame(species_scores)
species_scores_df$taxa <- rownames(species_scores_df)

# Fit environmental variables (plant-diversity treatments) to the PCA
envfit_result <- envfit(pca_sp_trans ~ Plant_SR + Grasses + Legumes + Short.Herbs + Tall.Herbs, data = treatments)
envfit_vectors <- as.data.frame(scores(envfit_result, display = "vectors"))
envfit_vectors$variable <- rownames(envfit_vectors)

# Define custom color vector
col_vec <- c("green", "pink", "blue", "magenta", "orange", "turquoise")

# Create a combined biplot with sites, taxa, and environmental vectors
p <- ggplot() +
  # Plot site scores
  geom_point(data = site_scores_df, aes(x = PC1, y = PC2, color = Plant_SR), size = 3) +
  geom_text_repel(data = site_scores_df, aes(x = PC1, y = PC2, label = plotcode), size = 3, max.overlaps = 20) +
  # Plot species scores with ggrepel to avoid overlap and larger font size
  geom_text_repel(data = species_scores_df, aes(x = PC1, y = PC2, label = taxa), color = "red", size = 3, max.overlaps = 20) +
  # Plot environmental vectors with larger font size
  geom_segment(data = envfit_vectors, aes(x = 0, xend = PC1, y = 0, yend = PC2), arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text_repel(data = envfit_vectors, aes(x = PC1, y = PC2, label = variable), color = "blue", size = 5, max.overlaps = 20) +
  # Additional plot styling
  scale_color_manual(values = col_vec) +
  labs(title = "PCA of Hellinger-transformed Arthropod Data with Environmental Vectors",
       x = "PC1", y = "PC2", color = "Plant Species Richness") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        plot.title = element_text(size = 16, face = "bold"),
        panel.background = element_blank()) +
  coord_fixed(ratio = 1, xlim = c(min(site_scores_df$PC1, species_scores_df$PC1) - 0.1, max(site_scores_df$PC1, species_scores_df$PC1) + 0.1),
              ylim = c(min(site_scores_df$PC2, species_scores_df$PC2) - 0.1, max(site_scores_df$PC2, species_scores_df$PC2) + 0.1)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed")

print(p)



#...................................
# Distance measures/PERMANOVA ----
#...................................
# calculation of different distance measures using transformed data, in order for comparisons to be "fairer"
dis_euc_taxa<-vegdist(community_matrix,"euclidean")
dis_bc_taxa<-vegdist(community_matrix,"bray")
dis_jacc_taxa<-vegdist(community_matrix,"jaccard", binary=T)
dis_chi_taxa<-vegdist(community_matrix,"chisq")

# Distance vector of plot species richness differences
dis_plant_sr<-dist(treatments$Plant_SR, 'euclidean')

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


# Perform PERMANOVA to test how significant the effect of each treatment is, (SIMPLE ANOVA DOESN'T WORK HERE, BUT PERMANOVA IS A ROBUST TOOL FOR MULTIVARIATE AND NON-PARAMETRIC ANALYSIS, MEANING I DON'T EVEN NEED TO CHECK NORMALITY OF RESIDUALS)
# use of bray-curtis simply because it is generally well suited for abundance data, use of jaccard proved fruitless
set.seed(1)
permanova_result <- adonis2(dis_bc_taxa ~ Plant_SR + Grasses + Short.Herbs + Tall.Herbs + Legumes, data = treatments, permutations = 999)

# Print the results
print(permanova_result)
#             Df SumOfSqs      R2       F Pr(>F)    
# Plant_SR     1  0.06867 0.05340  5.4448  0.001 ***  # only 5.34 % of the variation explained by Plant_SR, however the relationship is statistically significant (p < 0.01)
# Grasses      1  0.12997 0.10107 10.3059  0.001 ***  # 10% of the valiation attributed to the presence or absence of grasses (highest )
# Short.Herbs  1  0.05938 0.04618  4.7086  0.001 ***
# Tall.Herbs   1  0.03213 0.02498  2.5473  0.022 *  
# Legumes      1  0.06260 0.04868  4.9642  0.002 ***
# Residual    74  0.93323 0.72570                     # 73% of the data variation is derived from factors that are not our vegetation types          
# Total       79  1.28598 1.00000 


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
corr <- cor(treatments[-1], method = "spearman") # spearman can be used with mixed data (continuous and boolean)
corr
#              Plant_SR     Grasses Short.Herbs  Tall.Herbs     Legumes
# Plant_SR    1.0000000  0.39227902  0.35966898  0.35659899  0.32376832
# Grasses     0.3922790  1.00000000  0.22215277  0.16881931 -0.05593966
# Short.Herbs 0.3596690  0.22215277  1.00000000 -0.05534631  0.22215277
# Tall.Herbs  0.3565990  0.16881931 -0.05534631  1.00000000  0.16881931
# Legumes     0.3237683 -0.05593966  0.22215277  0.16881931  1.00000000

ggcorrplot(corr, hc.order = TRUE, lab = TRUE)
# highest correlations observed between Plant_SR and Grasses (0.39), Tall.Herbs (0.36), Short.Herbs (0.36) and Legumes (0.32)
# none of them are high


# To plot a density plot only with Plant_SR
set.seed(1)
permanova_result_dp <- adonis2(dis_bc_taxa ~ Plant_SR, data = treatments, permutations = 999)

perm_stat_taxa_dp <- permustats(permanova_result_dp)
summary(perm_stat_taxa_dp)
#          statistic    SES   mean lower median  upper Pr(perm)    
# Plant_SR    4.3998 5.7820 0.9939       0.8776 2.1022    0.001 ***

densityplot(perm_stat_taxa_dp)
# What to take away:
# 1. The observed F-statistic is 4.4, which is significantly higher (p < 0.01) than the mean 0.99. Only an outlier is found where the observed F-statistic is.
# 2. p-value 0.001 means that there is a very low probability of obtaining such an extreme F-statistic by chance. Additionally, SES > 0 and high (5.78)
# 3. According to the above, even though the correlation function used above gave a weak relationship between PLant_SR and arthropod taxa abundance and composition, that is not the case
# 4. The density plot provides a visual confirmation that the observed F-statistic is an outlier, thus supporting the opinion that it is an outlier.





#........................
# RDA - Constrained? ----
#........................
set.seed(1)
rda_taxa <- rda(community_matrix ~ Plant_SR + Grasses,
                data = treatments)

rda_taxa

summary(eigenvals(rda_taxa))
#                          RDA1     RDA2
# Eigenvalue            0.01104 0.005351
# Proportion Explained  0.15734 0.076281
# Cumulative Proportion 0.15734 0.233618

# permutations default
set.seed(1)

perm_stat_rda <- permustats(anova(rda_taxa))

summary(perm_stat_rda)

densityplot(perm_stat_rda)






















































































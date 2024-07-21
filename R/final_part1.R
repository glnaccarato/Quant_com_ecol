# 0 - Setup ----

# Load needed libraries and set default theme for ggplot
library(betapart)
library(iNEXT)
library(readxl)
library(tidyverse)
library(vegan)

theme_set(theme_bw())

# Read in data
data_ants <- read_excel("data/dat4_ants_elevation.xlsx")


# 1 - Calculate number of samples for each elevation ----
samples <- data_ants %>% 
  count(elevation_m, name = "n_samples")
view(samples)
# 2 sites at same elevation (= 941m)
# therefore 16 samples for each elevation, except 32 samples at elevation = 941m


# 2 - Alpha and gamma diversity for each elevation ----

# Create table with species abundance data only
spec_table <- data_ants %>% 
  select(-(1:4))

  ## 2.1 - Alpha diversity ----

S <- specnumber(spec_table)
hill_shannon <- as.numeric(renyi(spec_table, scales = 1, hill = T))
hill_simpson <- as.numeric(renyi(spec_table, scales = 2, hill = T))
evenness <- log(hill_shannon)/log(S) # log(hill-shannon) = shannon

data_alpha <- (tibble(elevation_m = data_ants$elevation_m,
                      site = data_ants$site,
                      S = S,
                      hill_shannon = hill_shannon,
                      hill_simpson = hill_simpson,
                      evenness = evenness))
rm(S, hill_shannon, hill_simpson, evenness)

# Now we can summarise the measures taking the mean for each,
# so that we get one value for each elevation
data_alpha <- data_alpha %>% 
  group_by(elevation_m) %>% 
  summarise(across(S:evenness, ~ mean(., na.rm = T)))

# Create a long table with only one column for the index-values
data_alpha_l <- data_alpha %>% 
  pivot_longer(S:evenness,
               names_to = "index",
               values_to = "value") %>% 
  # we can treat index as a factor, so we get facets in the order we choose
  mutate(index = factor(index, levels = 
                          c("S", "hill_shannon", "hill_simpson", "evenness")))

# The long table can now be used for a single plot with facet wrap
ggplot(data_alpha_l, aes(elevation_m, value, 
                         color = index, fill = index)) +
  geom_point(shape = 21, color = "black") +
  geom_smooth(se = T, method = "lm") +
  facet_wrap(~index, scales = "free") +
  labs(title = "Alpha diversity",
       x = "Elevation (m)",
       y = "Value") +
  guides(color = "none", fill = "none")


  ## 2.2 - Gamma diversity ----

# This time we pool together abundances before assessing indices
spec_pooled <- data_ants %>% 
  # we keep elevation_m for grouping until we don't need it anymore
  select(-(2:4)) %>% 
  group_by(elevation_m) %>%
  summarise(across(everything(), sum)) %>%
  # now we can get rid of elevation_m
  column_to_rownames(var = "elevation_m")

# Like this we get a table similar to the full spec_table,
# but with just one row for each elevation

S <- specnumber(spec_pooled)
hill_shannon <- as.numeric(renyi(spec_pooled, scales = 1, hill = T))
hill_simpson <- as.numeric(renyi(spec_pooled, scales = 2, hill = T))
evenness <- log(hill_shannon)/log(S) # log(hill-shannon) = shannon

data_gamma <- (tibble(elevation_m = samples$elevation_m,
                      n_samples = samples$n_samples,
                      S = S,
                      hill_shannon = hill_shannon,
                      hill_simpson = hill_simpson,
                      evenness = evenness))
rm(S, hill_shannon, hill_simpson, evenness)

# Create a long table with only one column for the index-values
data_gamma_l <- data_gamma %>% 
  pivot_longer(cols = S:evenness,
               names_to = "index",
               values_to = "value") %>% 
  mutate(index = factor(index, 
                        levels = c("S", "hill_shannon", "hill_simpson", "evenness")))

# The long table can now be used for a single plot with facet wrap
ggplot(data_gamma_l, aes(elevation_m, value, 
                         color = index, fill = index)) +
  geom_point(shape = 21, color = "black") +
  geom_smooth(se = T, method = "lm") +
  facet_wrap(~index, scales = "free") +
  labs(title = "Gamma diversity",
       x = "Elevation (m)",
       y = "Value") +
  guides(color = "none", fill = "none")


# 3 - Standardized comparisons ----

  ## 3.1 - Same no. of samples ----
# Sampling was standardized, but one elevation had two sites and
# therefore double the number of samples

elevations <- unique(data_ants$elevation_m)
dat_inc_freq <- list()
for (lvl in elevations) {
  subset <- spec_table[data_ants$elevation_m == lvl,]
  subset[subset > 0] <- 1
  n_samples <- nrow(subset)
  vector <- c(n_samples, sort(colSums(subset), decreasing = T))
  vector <- vector[vector > 0] 
  dat_inc_freq[[as.name(lvl)]] <- vector
  
  rm(lvl, n_samples, vector)
}

sbr_dat <- estimateD(dat_inc_freq, q = 0, datatype = "incidence_freq", 
                     level = 16)

sbr_dat

data_div <- data_gamma %>% 
  select(elevation_m)
data_div$S_t16 <- sbr_dat$qD


  ## 3.2 - Individual-based rarefaction ----

inext_table <- t(spec_pooled)
colnames(inext_table) <- elevations

inext <- iNEXT(inext_table, se = F)
ggiNEXT(inext) +
  scale_color_discrete(breaks = elevations) +
  scale_shape_manual(breaks = elevations, values = rep(23, times = 28)) +
  labs(title = "Rarefaction/extrapolation of species richness by sample size")
# Too many levels for plotting rarefaction curves
rm(inext)

ibr_dat <- estimateD(inext_table, q = 0, base = "size") 
ibr_dat

data_div$S_n12 <- ibr_dat$qD


  ## 3.3 - Coverage-based rarefaction ----
 
ggiNEXT(inext, type = 3)
cbr_dat <- estimateD(inext_table, q = 0, base = "coverage")
cbr_dat

# If base="coverage" and level=NULL, then this function computes 
# the diversity estimates for the minimum among all the coverage values 
# for samples extrapolated to double the reference sample sizes.

data_div$S_cov95 <- cbr_dat$qD

data_div_l <- data_div %>% 
  pivot_longer(S_t16:S_cov95, names_to = "index", values_to = "value") %>% 
  mutate(index = factor(index, levels = c("S_t16", "S_n12", "S_cov95")))

ibm_pal <- c("#648FFF", "#DC267F", "#FFB000")
ggplot(data_div_l, aes(elevation_m, value, color = index, fill = index)) +
  geom_point(shape = 21, color = "black") +
  geom_smooth(se = T, method = "lm") +
  facet_wrap(~index, scales = "free") +
  guides(color = "none", fill = "none") +
  scale_color_manual(values = ibm_pal) +
  scale_fill_manual(values = ibm_pal) +
  labs(x = "Elevation (m)",
       y = "Value")

# 4 - Beta diversity ----

# Change the species abundance table to presence-absence data
data_pa <- spec_table
data_pa[data_pa > 0] <- 1

# Initialize a tibble to store information about beta diversity 
data_beta <- tibble(elevation_m = numeric(),
                    beta.SIM = numeric(),
                    beta.SNE = numeric(),
                    beta.SOR = numeric())

# Iterate over elevation levels.
# For each elevation:
for (lvl in elevations) {
  # We extract the corresponding subset from the new PA-tibble.
  subset <- data_pa[data_ants$elevation_m == lvl,]
  
  # We assess the index for the row of data_beta to add current entries
  i <- match(lvl, elevations)
  
  # We add the value for elevation_m
  data_beta[i,]$elevation_m <- lvl

  # For cases where we have more than 16 samples per elevation (see h941m)
  if (nrow(subset) > 16) {
    # We use beta-sample to rarefy to 16 samples
    sampled <- beta.sample(subset, sites = 16, samples = 100)
    partitioning <- sampled$mean.values
    # We add the mean values for SIM, SNE and SOR to our data_beta tibble
    data_beta[i,]$beta.SIM <- partitioning[1]
    data_beta[i,]$beta.SNE <- partitioning[2]
    data_beta[i,]$beta.SOR <- partitioning[3]
    rm(sampled)
  }
  # In all other cases
  else {
    # We have 16 samples, so we can use beta.multi
    partitioning <- beta.multi(subset)
    # and directly add the output to our data_beta tibble
    data_beta[i,]$beta.SIM <- partitioning$beta.SIM
    data_beta[i,]$beta.SNE <- partitioning$beta.SNE
    data_beta[i,]$beta.SOR <- partitioning$beta.SOR
  }
  rm(i, lvl, partitioning, subset)
}

# Again, we change tibble to the long format:
data_beta_l <- data_beta %>%
  # We want to plot nestedness and turnover as proportions of SOR
  mutate(prop_nest = beta.SNE/beta.SOR,
         prop_turn = beta.SIM/beta.SOR) %>% 
  select(-c(beta.SIM, beta.SNE)) %>% 
  pivot_longer(beta.SOR:prop_turn,
               names_to = "index",
               values_to = "value")

# Now we can use facet_wrap:
ggplot(data_beta_l, aes(elevation_m, value, color = index, fill = index)) +
  geom_point(color = "black", shape = 21) +
  geom_smooth(se = T, method = "lm") +
  facet_wrap(~index, scales = "free") +
  labs(title = "Beta diversity along the elevational gradient",
       x = "Elevation (m)",
       y = "Value") +
  guides(color = "none", fill = "none")

# We can see that overall beta diversity (Sørensen dissimilarity) shows a
# slight decrease at higher elevations, even though variance gets much higher.
# What is interesting here is the transition from turnover to nestedness
# at higher elevations, with complete nestedness at elevations >= 1651m!
# (except for 1828 m <- outlier?).

# Plot without potential outlier
ggplot(data_beta_l[1:81,], aes(elevation_m, value, color = index, fill = index)) +
  geom_point(color = "black", shape = 21) +
  geom_smooth(se = T, method = "lm") +
  facet_wrap(~index, scales = "free") +
  labs(title = "Beta diversity along the elevational gradient",
       x = "Elevation (m)",
       y = "Value") +
  guides(color = "none", fill = "none")
# Doesn't really change the outcome, even though the decrease in overall beta
# is more clear - keeping it for the sake of completeness.

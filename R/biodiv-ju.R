# Setup ----

# Load needed libraries and set default theme for ggplot
library(betapart)
library(readxl)
library(tidyverse)
library(vegan)


theme_set(theme_bw())


# Read in data
data_ants <- read_excel("data/dat4_ants_elevation.xlsx")

# Calculate number of samples for each elevation ----
samples <- data_ants %>% 
  count(elevation_m, name = "n_samples")
view(samples)
# 2 sites at same elevation (= 941m)
# therefore 16 samples for each elevation, except 32 samples at elevation = 941m


# alpha and gamma diversity for each elevation ----

# Create table with species abundance data only
spec_table <- data_ants %>% 
  select(-(1:4))

  ## alpha diversity ----
S <- specnumber(spec_table)
N <- rowSums(spec_table)
shannon <- diversity(spec_table, index = "shannon")
simpson <- diversity(spec_table, index = "simpson")
evenness <- shannon / log(S)
  
data_alpha <- (tibble(elevation_m = data_ants$elevation_m,
                      site = data_ants$site,
                      S = S,
                      N = N,
                      shannon = shannon,
                      simpson = simpson,
                      evenness = evenness))
rm(S, N, shannon, simpson, evenness, site)
# In many cases only one species was found in a sample, 
# leading to shannon = 0 and evenness = NA

# Create a long table with only one column for the index-values
data_alpha_l <- data_alpha %>% 
  pivot_longer(S:evenness,
               names_to = "index",
               values_to = "value")

# The long table can now be used for a single plot with facet wrap
ggplot(data_alpha_l, aes(elevation_m, value, fill = elevation_m)) +
  geom_point(shape = 21) +
  facet_wrap(~index, scales = "free") +
  labs(title = "Alpha diversity along elevation gradient") +
  guides(fill = "none")
rm(data_alpha_l)

  ## gamma diversity ----

# List of samples Sites
elevations <- unique(data_ants$elevation_m) 

# Initialize empty variable to store values for S per elevation
S <- NULL

# Iterate through sample sites and calculate corresponding values for S
for (lvl in elevations) {
  # subset for current site
  subset <- spec_table[data_ants$elevation_m == lvl,]
  S <- append(S, specnumber(colSums(subset)))
  rm(lvl, subset)
}

data_gamma <- data_alpha %>% 
  group_by(elevation_m) %>%
  summarise(N = sum(N),
            shannon = mean(shannon),
            simpson = mean(simpson),
            evenness = mean(evenness, na.rm = T)) %>% 
  mutate(S = S,
         n_samples = samples$n_samples)
rm(S)
# Warning: many NAs removed for assessing the mean evenness

# Create a long table with only one column for the index-values
data_gamma_l <- data_gamma %>% 
  pivot_longer(cols = shannon:S,
               names_to = "index",
               values_to = "value")

# The long table can now be used for a single plot with facet wrap
ggplot(data_gamma_l, aes(elevation_m, value, fill = elevation_m)) +
  geom_point(shape = 21) +
  facet_wrap(~index, scales = "free") +
  labs(title = "Gamma diversity along elevation gradient") +
  guides(fill = "none")
rm(data_gamma_l)


# Rarefaction ----

data_inext <- data_ants %>% 
  select(-(2:4)) %>% 
  group_by(elevation_m) %>% 
  summarise(across(everything(), sum)) %>%
  column_to_rownames(var = "elevation_m") %>% 
  t() 

inext1 <- iNEXT(data_inext)
ggiNEXT(inext1)

# groups <- NULL
# for (lvl in elevations) {
#   name <- paste0("h",lvl,"m")
#   groups <- groups %>% append(name)
#   rm(lvl, name)
# }
# 
# i = 1
# for (lvl in elevations) {
#   S_sorted <- sort(data_alpha[data_alpha$elevation_m == lvl,]$S, decreasing = T)
#   assign(groups[i], S_sorted)
#   i <- i + 1
# }
# # 
# # data_inext <- list


# Beta-Diversity ----

# Whittaker's beta function (from Felix)
beta_w <- function(community)
{
  require(vegan)
  gamma <- specnumber(colSums(community))
  alpha <- mean(specnumber(community))
  return(gamma/alpha)
}

# Change the species abundance table to presence-absence 
data_pa <- spec_table
data_pa[data_pa > 0] <- 1

# Initialize a tibble to store information about beta diversity 
data_beta <- tibble(elevation_m = numeric(),
                    whittaker = numeric(),
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
  
  # We add values for elevation_m and for Whittaker's beta diversity
  data_beta[i,]$elevation_m <- lvl
  data_beta[i,]$whittaker <- beta_w(subset)
  
  # For cases where we have more than 16 samples per elevation (see h941m)
  if (nrow(subset) > 16) {
    # We use beta-sample to rarefy to 16 samples
    sampled <- beta.sample(subset, index.family = "sor", sites = 16,
                           samples = 100)
    partitioning <- sampled$mean.values
    # We add the mean values for SIM, SNE and SOR to our data_beta tibble
    data_beta[i,]$beta.SIM <- round(partitioning[1], digits = 2)
    data_beta[i,]$beta.SNE <- round(partitioning[2], digits = 2)
    data_beta[i,]$beta.SOR <- round(partitioning[3], digits = 2)
  }
  # In all other cases
  else {
    # We have 16 samples, so we can use beta.multi
    partitioning <- beta.multi(subset)
    # and directly add the output to our data_beta tibble
    data_beta[i,]$beta.SIM <- round(partitioning$beta.SIM, digits = 2)
    data_beta[i,]$beta.SNE <- round(partitioning$beta.SNE, digits = 2)
    data_beta[i,]$beta.SOR <- round(partitioning$beta.SOR, digits = 2)
  }
}

# First we can plot Whittaker's beta diversity along the elevation gradient:
ggplot(data_beta, aes(elevation_m, whittaker, fill = elevation_m)) +
  geom_point(shape = 21) +
  guides(fill = "none")

# We can see that we have an outlier at max elevation.
# In this group only 2 species were observed in total.
# Thus, the mean alpha div is extremely low (0.125) and the resulting
# beta diversity extremely high.
# We could plot it without it, to better observe the general trend..
ggplot(filter(data_beta, elevation_m < 1828), 
       aes(elevation_m, whittaker, fill = elevation_m)) +
  geom_point(shape = 21) +
  guides(fill = "none")
# Trend seems slightly decreasing, but not clearly 

# In order to plot partitioning together, we need to make adjustments
# to the structure of the tibble (same as we did for alpha diversity):
data_beta_l <- data_beta %>%
  # We want to plot nestedness and turnover as proportions of SOR
  mutate(prop_nest = beta.SNE/beta.SOR,
         prop_turn = beta.SIM/beta.SOR) %>% 
  select(-c(beta.SIM, beta.SNE)) %>% 
  pivot_longer(beta.SOR:prop_turn,
               names_to = "index",
               values_to = "value")

# Now we can use facet_wrap:
ggplot(data_beta_l, aes(elevation_m, value, fill = elevation_m)) +
  geom_point(shape = 21) +
  facet_wrap(~index, scales = "free") +
  labs(title = "Partitioning of beta diversity into nestedness and spatial turnover along the elevation gradient") +
  guides(fill = "none")

# We can see how the variance of Sørensen dissimilarity increases
# at higher elevations!
# What is interesting here is the transition from turnover to nestedness
# at higher elevations, with complete nestedness at elevations >= 1651m!
# (except for the outlier h1828m).


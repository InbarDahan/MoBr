

   # Calculating aggregations using the Negative binomial dispersion
     # parameter k for each species in each treatment (depth_)

# libraries:
library(dplyr)
library(tidyr)
library(tidyverse)
library(MASS) # For negative binomial model
# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________

# read data:
setwd(wd_processed_data)
fish <- read.csv("wide_red.csv")
sponges <- read.csv("sponges_data.csv")
stony <- read.csv("stony_wide_data.csv")
soft <- read.csv("soft_data.csv")

fish = fish %>% dplyr::select(-X) 
stony  = stony %>% dplyr::select(-X) 
soft = soft %>% dplyr::select(-X) 
# _______________________________________________________________

      # data preparations for this analysis:

# renaming columns:
fish <- fish %>% rename(sample = OpCode)
sponges <- sponges %>% rename(sample = quad)
stony <- stony %>% rename(sample = plot, depth_group = depth)
soft <- soft %>% rename(depth_group = depth)

# ------------

# rename
fish <- fish %>% rename(lon = lon_x)
fish <- fish %>% rename(lat = lat_y)
# ------------

# adding a column for the taxonomic groups id:
fish    <- fish    %>% mutate(taxon = 'fish',           .before = 1)
sponges <- sponges %>% mutate(taxon = 'sponges',        .before = 1)
stony   <- stony   %>% mutate(taxon = 'stony corals',   .before = 1)
soft    <- soft    %>% mutate(taxon = 'soft corals',    .before = 1)
# ------------

# binning the fish data:
# 8.1 - 149 m to - 15 m bins (10 breaks)
fish <- fish %>% mutate(depth_group = cut(fish$depth,                                         
         breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
          labels=c('7.5', '22.5', '37.5', '52.5', '67.5','82.5', '97.5', '112.5', '127.5', '142.5'))) %>% 
           relocate(depth_group,.after = depth) 
# ------------

fish    <- fish %>% dplyr::select(-sea, -temperature, - depth)
sponges <- sponges %>% dplyr::select(-depth, -Habitat, -Site)
soft    <- soft %>% dplyr::select(-site)
stony   <- stony %>% dplyr::select(-site)
# ------------

taxon_list <- list(fish = fish, sponges = sponges, stony = stony, soft = soft)

# _______________________________________________________________

# General steps:
# for each taxon, and each species (as columns) and each depth_group - 
# 1 - fit negative binomial distribution
# 2 - calc k (dispersion parameter) * do it only for species with sufficient data
# 3 - plot k against depth
# 4 - run senlm model to catch the trend
# add it to the main measures vs depth plot

# _______________________________________________________________
 
# a loop to fit NB:

# Create a list to store model results
nb_model_results <- list()

# Loop over each taxon (each dataset in the list)
for (taxon_name in names(taxon_list)) {
  
  df <- taxon_list[[taxon_name]]  # Get the taxon's data
  depth_groups <- unique(df$depth_group)  # Get unique depth groups
  
  # Create a list to store results per taxon
  taxon_results <- list()
  
  # Loop over depth groups
  for (depth in depth_groups) {
    
    df_depth <- df[df$depth_group == depth, ]  # Subset data for current depth
    
    # Create a list to store results per species
    species_results <- list()
    
    # Loop over species (columns after the first 5 explanatory columns)
    for (species in names(df)[6:ncol(df)]) {
      
      # Get abundance data for the species
      abundance <- df_depth[[species]]
      
      # Check if species has enough data to fit a model (avoid models with too many zeros)
      if (sum(abundance > 0) >= 3) {  # Ensure at least some nonzero values
        
        # Fit Negative Binomial model (abundance ~ 1 means no predictor except intercept)
        model <- tryCatch(
          glm.nb(abundance ~ 1, data = df_depth),
          error = function(e) return(NULL)  # Handle errors
        )
        
        # Store model if successful
        if (!is.null(model)) {
          species_results[[species]] <- model
        }
      }
    }
    
    # Store species results in taxon-depth structure
    taxon_results[[as.character(depth)]] <- species_results
  }
  
  # Store taxon results
  nb_model_results[[taxon_name]] <- taxon_results
}

# _______________________________________________________________

# - interpertation: 
# after I filtered for species with more then 3 indv per depth group
# assuming less then 3 is not enough to determind grouping 
# and runnig the binomial models on thoes species left, I see that the variations
# in the number of species per depth group that worked is great
# it might mean that some depths had a lot of rare species, that where not included

# q - is it ok? or will it biass the next analyses? 
############ Ask Yoni
# _______________________________________________________________

# A loop to extract k:

# Create a list to store extracted k values
k_values <- list()

for (taxon_name in names(nb_model_results)) {
  
  taxon_data <- nb_model_results[[taxon_name]]
  taxon_k <- list()
  
  for (depth in names(taxon_data)) {
    
    depth_data <- taxon_data[[depth]]
    depth_k <- list()
    
    for (species in names(depth_data)) {
      
      model <- depth_data[[species]]
      
      if (!is.null(model)) {
        k_value <- 1 / model$theta  # Extract dispersion parameter k
        depth_k[[species]] <- k_value
      }
    }
    
    taxon_k[[depth]] <- depth_k
  }
  
  k_values[[taxon_name]] <- taxon_k
}

# _______________________________________________________________

# convert results to data frame and plot

# Convert the list into a data frame
k_df <- do.call(rbind, lapply(names(k_values), function(taxon) {
  do.call(rbind, lapply(names(k_values[[taxon]]), function(depth) {
    data.frame(
      taxon = taxon,
      depth_group = as.numeric(depth),
      species = names(k_values[[taxon]][[depth]]),
      k_value = unlist(k_values[[taxon]][[depth]])
    )
  }))
}))

# Plot k against depth
library(ggplot2)

ggplot(k_df, aes(x = depth_group, y = k_value, color = taxon)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  labs(x = "Depth Group", y = "Dispersion Parameter (k)",
       title = "Intraspecific Aggregation (k) Along Depth Gradient") +
  theme_minimal()

# smaller k -> higher aggregations


























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
wd_plots <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/output/plots/agg_plot" 
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

          # method 1: k dispersion\clumping parameter:

# General steps:
# for each taxon, and each species (as columns) and each depth_group - 
# 1 - fit negative binomial distribution
# 2 - calc k (dispersion parameter) 
# * do it only for species with sufficient data - more than 5 ind
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
      if (sum(abundance) >= 5) {  # Ensure at least some nonzero values
        
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

# - interpretation: 
# after I filtered for species with more then 5 indv per depth group
# assuming less then 5 is not enough to determine grouping 
# and running the binomial models on the species left, I see that the variations
# in the number of species per depth group that worked is great
# it might mean that some depths had a lot of rare species, that where not included

# q - is it ok? or will it bias the next analyses? 
# Yoni:
# answer - check other methods and see if the trend is consistent + increase to 5 
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
        k_value <- model$theta # Extract dispersion parameter k
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

# _______________________________________________________________


# Plot k against depth, faceted by taxon
library(ggplot2)

agg_plot <- ggplot(k_df, aes(x = depth_group, y = k_value, color = taxon)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~ taxon) +
  labs(
    x = "Depth Group",
    y = "Dispersion Parameter (k)",
    title = "Intraspecific Aggregation (k) Along Depth Gradient (n >= 5)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")  # remove redundant legend if each facet is a taxon

# smaller k -> higher aggregations

# _______________________________________________________________

# save plot
setwd(wd_plots)
ggsave("plot_k_dispersion_5_indv.png", plot = agg_plot, width = 8, height = 6, dpi = 300, bg = "white")

# _______________________________________________________________
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# _______________________________________________________________


# _____________________________________
# _____________________________________

        # method 2: j index of ives (j = 1/k)
   # b - keeping species with more than 5 indv:

# Create a list to store J index results
j_results <- list()

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
      
      # Get abundance vector for the species at this depth group
      abundance <- df_depth[[species]]
      
      # Filtering step: ensure sufficient data to avoid instability
      if (sum(abundance > 0) >= 5 && mean(abundance) > 1) {
        
        xj <- abundance
        N <- length(xj)
        X <- mean(xj)
        
        # Compute J index (Ives 1988, 1991)
        J <- ((sum(xj * (xj - 1)) / (X * N)) - X) * (1 / X)
        
        # Store the J index and supporting info
        species_results[[species]] <- list(
          J = J,
          mean_abundance = X,
          nonzero_sites = sum(xj > 0)
        )
      }
    }
    
    # Store species-level J results in depth group
    taxon_results[[as.character(depth)]] <- species_results
  }
  
  # Store taxon-level results
  j_results[[taxon_name]] <- taxon_results
}

# _____________________________________

# Create an empty data frame to collect all results
j_summary <- data.frame()

for (taxon in names(j_results)) {
  for (depth in names(j_results[[taxon]])) {
    species_data <- j_results[[taxon]][[depth]]
    for (species in names(species_data)) {
      entry <- species_data[[species]]
      j_summary <- rbind(j_summary, data.frame(
        taxon = taxon,
        depth_group = depth,
        species = species,
        J = entry$J,
        mean_abundance = entry$mean_abundance,
        nonzero_sites = entry$nonzero_sites
      ))
    }
  }
}

# _____________________________________

ggplot(j_summary, aes(x = as.numeric(depth_group), y = J, color = taxon)) +
  geom_point(aes(shape = J < 0), size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ taxon, scales = "free_y") +
  labs(
    x = "Depth Group",
    y = "J (Intraspecific Aggregation)",
    shape = "Negative J",
    title = "Intraspecific Aggregation (J) Across Depths by Taxon"
  ) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16)) +
  theme_minimal()

# _______________________________________________________________
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# _______________________________________________________________

        # check alignment between the two measures

# Ensure species names are character, not factors (if needed)
j_summary$species <- as.character(j_summary$species)
k_df$species <- as.character(k_df$species)

# Ensure depth is numeric
j_summary$depth_group <- as.numeric(as.character(j_summary$depth_group))
k_df$depth_group <- as.numeric(as.character(k_df$depth_group))

# Merge by taxon, depth group, and species
comparison_data <- merge(j_summary, k_df, 
                         by = c("taxon", "depth_group", "species"))

# Remove zeros or NAs to avoid issues with log()
comparison_data <- comparison_data[comparison_data$k_value > 0 & comparison_data$J > 0, ]

# Add log-transformed variables
comparison_data$log_1_k <- log(1 / comparison_data$k_value)
comparison_data$log_J <- log(comparison_data$J)

# Now you can plot it
ggplot(comparison_data, aes(x = log_1_k, y = log_J)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~ taxon)

#####################################################################
#####################################################################
#####################################################################
 
# I calculated aggregations using 2 matrices: J and K. 
# They showed some consistencies and some inconsistencies, although they are 
# supposed to be an inverse measure of each other - 
#  --> understand why

# something is wrong,  i a not convinced

# also, I stopped in the middle of chat gpt explainig about the next measure

#####################################################################
#####################################################################
#####################################################################

# _______________________________________________________________
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# _______________________________________________________________

       # method 3: standardized version of the Morisita index:








# method 2: j index of ives (j = 1/k)
# a - keeping species with more than 1 indv:

# Create a list to store J index results
j_index_results <- list()

# Loop over each taxon (each dataset in the list)
for (taxon_name in names(taxon_list)) {
  
  df <- taxon_list[[taxon_name]]  # Get the taxon's data
  depth_groups <- unique(df$depth_group)  # Get unique depth groups
  
  # Create a list to store results per taxon
  taxon_results <- list()
  
  for (depth in depth_groups) {
    
    df_depth <- df[df$depth_group == depth, ]  # Subset data for current depth
    species_results <- list()
    
    for (species in names(df)[6:ncol(df)]) {
      
      abundance <- df_depth[[species]]
      
      if (sum(abundance > 1)) {  # Only consider non-singleton species
        
        X <- mean(abundance)
        N <- length(abundance)
        
        if (X > 0) {
          numerator <- sum(abundance * (abundance - 1))
          J <- ((numerator / (X * N)) - X) / X
          species_results[[species]] <- J
        }
      }
    }
    
    taxon_results[[as.character(depth)]] <- species_results
  }
  
  j_index_results[[taxon_name]] <- taxon_results
}

# _____________________________________

# Create an empty data frame to collect all results
j_summary_2 <- data.frame()

for (taxon in names(j_index_results)) {
  for (depth in names(j_index_results[[taxon]])) {
    species_data <- j_index_results[[taxon]][[depth]]
    for (species in names(species_data)) {
      j_value <- species_data[[species]]
      j_summary_2 <- rbind(j_summary_2, data.frame(
        taxon = taxon,
        depth_group = depth,
        species = species,
        J = j_value
      ))
    }
  }
}

# _____________________________________

ggplot(j_summary_2, aes(x = as.numeric(depth_group), y = J, color = taxon)) +
  geom_point(aes(shape = J < 0), size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ taxon, scales = "free_y") +
  labs(
    x = "Depth Group",
    y = "J (Intraspecific Aggregation)",
    shape = "Negative J",
    title = "Intraspecific Aggregation (J) Across Depths by Taxon"
  ) +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = 16)) +
  theme_minimal()










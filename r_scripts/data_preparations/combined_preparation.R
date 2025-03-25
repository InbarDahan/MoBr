      

          # creating a combined data frame for the plots

# libraries:
library(dplyr)
library(tidyr)
library(poilog) # calculate evenness metric (sigma)
# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________

# read data:
setwd(wd_processed_data)
fish <- read.csv("wide_red.csv")
sponges <- read.csv("sponges_data.csv")
stony <- read.csv("stony_wide_data.csv")
soft_corals <- read.csv("soft_data.csv")
#stony_s <- read.csv("stony_wide_species_data.csv") # unused for now - just species and not other taxonomic levels

fish = fish %>% dplyr::select(-X) 
stony  = stony %>% dplyr::select(-X) 
soft_corals = soft_corals %>% dplyr::select(-X) 
#stony_s  = stony_s %>% dplyr::select(-X) 
# _______________________________________________________________

         # calculate sample abundance and richness 
              # for each taxonomic groups:

## fish:
first_species_f <- 7

fish_sum <- fish %>% 
  mutate(abundance = rowSums(fish[,first_species_f:length(fish)])) %>%
  mutate(richness = rowSums(fish[, first_species_f:ncol(fish)] > 0)) %>%
  dplyr::select(OpCode, depth, lon_x, lat_y, richness, abundance)

## sponges:
first_species_sp <- 8

sponges_sum <- sponges %>% 
  mutate(abundance = rowSums(sponges[,first_species_sp:length(sponges)])) %>%
  mutate(richness = rowSums(sponges[, first_species_sp:ncol(sponges)] > 0)) %>%
  dplyr::select(quad, depth_group, depth, Site, lon, lat, richness, abundance)

## stony corals:
first_species_st <- 6

stony_sum <- stony %>% 
  mutate(abundance = rowSums(stony[,first_species_st:length(stony)])) %>%
  mutate(richness = rowSums(stony[, first_species_st:ncol(stony)] > 0)) %>%
  dplyr::select( site,depth, plot, lon, lat, richness, abundance)

## soft corals:
first_species_sf <- 6

soft_sum <- soft_corals %>% 
  mutate(abundance = rowSums(soft_corals[,first_species_sf:length(soft_corals)])) %>%
  mutate(richness = rowSums(soft_corals[, first_species_sf:ncol(soft_corals)] > 0)) %>%
  dplyr::select(site,depth, sample, lon, lat, richness, abundance)
# _______________________________________________________________
# _______________________________________________________________

       # evenness using Poisson distribution (sigma)
                      # per sample:

## fish:

species_matrix_f <- fish[,first_species_f:length(fish)] 

# Create an empty vector to store the sigma values
sigma_values_f <- numeric(nrow(species_matrix_f))

# Loop over each sample (row) in the matrix
for (i in 1:nrow(species_matrix_f)) {
  # Extract the abundance vector for the i-th sample
  sample_abundance <- species_matrix_f[i, ]
  
  # Remove species with zero abundance (optional, based on your needs)
  sample_abundance <- sample_abundance[sample_abundance > 0]
  
  # Fit the Poisson Lognormal model to the sample
  fit <- poilogMLE(sample_abundance)
  
  # Extract the sigma (evenness) value and store it
  sigma_values_f[i] <- fit$par[2]
}

# Add the sigma values as a new column in your original data frame
fish_sum$sigma <- sigma_values_f

# - - - 

## sponges:

species_matrix_sp <- sponges[,first_species_sp:length(sponges)] 

# Create an empty vector to store the sigma values
sigma_values_sp <- numeric(nrow(species_matrix_sp))

# Loop over each sample (row) in the matrix
for (i in 1:nrow(species_matrix_sp)) {
  # Extract the abundance vector for the i-th sample
  sample_abundance <- species_matrix_sp[i, ]
  
  # Remove species with zero abundance
  sample_abundance <- sample_abundance[sample_abundance > 0]
  
  # Check if the sample has any species left (after removing zeros)
  if (length(sample_abundance) > 0) {
    # Fit the Poisson Lognormal model to the sample
    fit <- poilogMLE(sample_abundance)
    
    # Extract the sigma (evenness) value and store it
    sigma_values_sp[i] <- fit$par[2]
  } else {
    # If the row has no species (all zeros), assign NA or another placeholder
    sigma_values_sp[i] <- NA
  }
}
# Add the sigma values as a new column in your original data frame
sponges_sum$sigma <- sigma_values_sp

# I have many samples with abundance = 0 so there is no distribution
# to evaluate and hence the sigma is defined NA.

# - - - 

## stony corals:

species_matrix_st <- stony[,first_species_st:length(stony)] 

# Create an empty vector to store the sigma values
sigma_values_st <- numeric(nrow(species_matrix_st))

# Loop over each sample (row) in the matrix
for (i in 1:nrow(species_matrix_st)) {
  # Extract the abundance vector for the i-th sample
  sample_abundance <- species_matrix_st[i, ]
  
  # Remove species with zero abundance (optional, based on your needs)
  sample_abundance <- sample_abundance[sample_abundance > 0]
  
  # Fit the Poisson Lognormal model to the sample
  fit <- poilogMLE(sample_abundance)
  
  # Extract the sigma (evenness) value and store it
  sigma_values_st[i] <- fit$par[2]
}

# Add the sigma values as a new column in your original data frame
stony_sum$sigma <- sigma_values_st

# - - - 

## soft corals:

species_matrix_sf <- soft_corals[,first_species_sf:length(soft_corals)] 

# Create an empty vector to store the sigma values
sigma_values_sf <- numeric(nrow(species_matrix_sf))

# Loop over each sample (row) in the matrix
for (i in 1:nrow(species_matrix_sf)) {
  # Extract the abundance vector for the i-th sample
  sample_abundance <- species_matrix_sf[i, ]
  
  # Remove species with zero abundance (optional, based on your needs)
  sample_abundance <- sample_abundance[sample_abundance > 0]
  
  # Fit the Poisson Lognormal model to the sample
  fit <- poilogMLE(sample_abundance)
  
  # Extract the sigma (evenness) value and store it
  sigma_values_sf[i] <- fit$par[2]
}

# Add the sigma values as a new column in your original data frame
soft_sum$sigma <- sigma_values_sf

# _______________________________________________________________
# _______________________________________________________________

# renaming columns before combining the data frames:
fish_sum <- fish_sum %>% rename(sample = OpCode)
sponges_sum <- sponges_sum %>% rename(sample = quad)
stony_sum <- stony_sum %>% rename(sample = plot)

fish_sum <- fish_sum %>% rename(lon = lon_x)
fish_sum <- fish_sum %>% rename(lat = lat_y)

# adding a column for the taxonomic groups id:
fish_sum$taxon <- 'fish'
sponges_sum$taxon <- 'sponges'
stony_sum$taxon <- 'corals stony'
soft_sum$taxon <- 'corals soft'
# _______________________________________________________________

# Select the relevant columns from each data frame
fish_selected <- fish_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance", "sigma")]
sponges_selected <- sponges_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance", "sigma")]
stony_selected <- stony_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance", "sigma")]
soft_selected <- soft_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance", "sigma")]
# _______________________________________________________________

# Combine the selected columns from each data frame into a new data frame
combined_data <- rbind(fish_selected, sponges_selected, stony_selected, soft_selected)
combined_data$lon[combined_data$lon == 0] <- NA
combined_data$lat[combined_data$lat == 0] <- NA
 # _______________________________________________________________

## save:
setwd(wd_processed_data)
write.csv(combined_data, file = "combined_data.csv")

# ________________________________________________________________

# ________________________________________________________________



























# creating another df format to allow faceted graphs:

combined_long <- combined_data %>% pivot_longer(cols = c(abundance, richness, sigma), 
               names_to = "diversity_measure", 
               values_to = "value")

setwd(wd_processed_data)
write.csv(combined_long, file = "combined_long.csv")

# ________________________________________________________________

# ________________________________________________________________

             # filter data - by location (NR SOUTH)

library(sf)
library(raster)       
library(mapview)

# - - -
# df coor
df_with_coords <- combined_long %>% filter(!is.na(lon) & !is.na(lat))  # Keep rows with valid coordinates
# df no coo - will be joined later
df_without_coords <- combined_long %>% filter(is.na(lon) | is.na(lat))  # Keep rows with missing coordinates

# - - -

# df to sf object
taxon_data_sf <- st_as_sf(df_with_coords, coords = c('lon', 'lat'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# - - -

# creat polygon for the NR SOUTH
polygon_coords <- matrix(c(
  34.923950, 29.511717,  # Point 1: Coral Beach (replace with actual coordinates)
  34.900030, 29.489906,  # Point 2: South (adjust latitude to southern border)
  34.899837, 29.467587,  # Point 3: Adjust west longitude
  34.937480, 29.503077,  # Point 4: West of Coral Beach
  34.923950, 29.511717   # Closing the polygon back to Point 1
), byrow = TRUE, ncol = 2)

# Convert coordinates to an sf polygon object
polygon <- st_polygon(list(polygon_coords))
polygon_sf <- st_sfc(polygon, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")  # Set CRS to WGS84 (latitude/longitude)

# map the polygon and my points to check I included everything:
mapview(polygon_sf) + mapView(taxon_data_sf, zcol = "taxon", legend = TRUE, layer.name = 'Taxon')

# Filter out the points that are outside the polygon
df_filtered_sf <- taxon_data_sf[st_intersects(taxon_data_sf, polygon_sf, sparse = FALSE), ]

# map the polygon and my point to see that I included everything:
mapview(polygon_sf) + mapView(df_filtered_sf, zcol = "taxon", legend = TRUE, layer.name = 'Taxon')

df_filtered <- st_drop_geometry(df_filtered_sf)


           ### might be unneeded when Tom coordinates are added:
# Combine the filtered spatial data with the rows that had missing coordinates
df_final <- bind_rows(as.data.frame(df_filtered), df_without_coords)

setwd(wd_processed_data)
write.csv(df_final, file = "NR_SOUTH_data.csv")
# ________________________________________________________________

# ________________________________________________________________

           # filter data - by features 

## what do we want: from sponges, corals and fish? 

## how do I separate juveniles from adults by size? Does Shahar has these data?









# identifying outliers in the data:

# stony corals:
# Identify the row index with sigma value 18 (replace `row_with_sigma_18` if known)
row_with_sigma_18 <- 20  # Replace with the actual index if you have it

# Extract the abundance data for this sample
abundance_data <- as.numeric(species_matrix_st[row_with_sigma_18, ])

# Plot a histogram of species abundances
hist(abundance_data, 
     breaks = 20, 
     main = "Histogram of Species Abundances (Sigma = 18)", 
     xlab = "Species Abundance", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black")



# Identify the row index with sigma value 18 (replace `row_with_sigma_18` if known)
row_with_sigma_6 <- 323  # Replace with the actual index if you have it

# Extract the abundance data for this sample
abundance_data <- as.numeric(species_matrix_st[row_with_sigma_6, ])

# Plot a histogram of species abundances
hist(abundance_data, 
     breaks = 20, 
     main = "Histogram of Species Abundances (Sigma = 18)", 
     xlab = "Species Abundance", 
     ylab = "Frequency", 
     col = "lightblue", 
     border = "black")











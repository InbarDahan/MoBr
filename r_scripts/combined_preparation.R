      

          # creating a combined data frame for the plots

# libraries:
library(dplyr)
library(tidyr)
# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________

# read data:
setwd(wd_processed_data)
fish <- read.csv("wide_red.csv")
sponges <- read.csv("sponges_data.csv")
stony <- read.csv("stony_wide_data.csv")

fish = fish %>% dplyr::select(-X) 
stony  = stony %>% dplyr::select(-X) 
# _______________________________________________________________

# calculate abundance and richness per sample for each taxonomic groups:
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
stony_sum$taxon <- 'stony corals'
# _______________________________________________________________

# Select the relevant columns from each data frame
fish_selected <- fish_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance")]
sponges_selected <- sponges_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance")]
stony_selected <- stony_sum[, c("taxon", "sample", "depth", "lon", "lat", "richness", "abundance")]
# _______________________________________________________________

# Combine the selected columns from each data frame into a new data frame
combined_data <- rbind(fish_selected, sponges_selected, stony_selected)
# _______________________________________________________________

## save:
setwd(wd_processed_data)
write.csv(combined_data, file = "combined_data.csv")













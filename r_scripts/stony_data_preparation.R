
                      # stony corals data preparations

library(dplyr)
library(tidyr)
library(ggplot2)
library(mapview)

wd_raw_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw" 
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 

# read data
setwd(wd_raw_data)
stony <- read.csv("stony_corals_iui.csv")

# remove soft corals
stony_rm_octa <- stony %>% filter(!grepl("Octocorallia", Comment, ignore.case = TRUE))

# make a column for the species and genus combined
stony_rm_octa <- stony_rm_octa %>% mutate(species_id = paste(Genus, Species, sep = " "))

# long to wide format
stony_wide <- stony_rm_octa %>%
  group_by(Site, Depth, Plot, species_id) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = list(count = 0))

# change the cloumn name from Depth to depth: 
stony_wide <- stony_wide %>% rename(depth = Depth)
stony_wide <- stony_wide %>% rename(site = Site)
stony_wide <- stony_wide %>% rename(plot = Plot)

# change the rows with NA to 0:
stony_wide[is.na(stony_wide)] <- 0
 
# reconstruct coordinates:

# add two columns:lon and lat - filled with NA:
stony_wide$lon <- NA
stony_wide$lat <- NA

# move this columns after the plot column:
stony_wide <- stony_wide %>% relocate(lon, lat, .after = plot)

######################################
# reconstruct the coordinates values

# save:
setwd(wd_processed_data)
write.csv(stony_wide, file = "stony_wide_data.csv")

# for Tom:
species <- as.data.frame(unique(stony_rm_octa$species_id))
species <- species %>% rename (species_list = `unique(stony_rm_octa$species_id)`)
setwd(wd_processed_data)
write.csv(species, file = "stony_species.csv")

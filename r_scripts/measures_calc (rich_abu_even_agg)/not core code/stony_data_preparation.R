
                     
                 # stony corals data preparations

library(dplyr)
library(tidyr)
library(ggplot2)
library(mapview)

wd_raw_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw" 
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_coor <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw/stony_coor"

# read data
setwd(wd_raw_data)
stony <- read.csv("stony_corals_iui.csv")

# change the column name from Depth to depth: 
stony <- stony %>% rename(depth = Depth)
stony <- stony %>% rename(site = Site)
stony <- stony %>% rename(plot = Plot)

# reconstruct coordinates:
# add two columns:lon and lat - filled with NA:
stony$lon <- NA
stony$lat <- NA

# move this columns after the plot column:
stony <- stony %>% relocate(lon, lat, .after = plot)

###################################################
setwd(wd_coor)
# read coor:
coor <- read.csv("stony_coor_f.csv")
##################################################

# remove soft corals
stony_rm_octa <- stony %>% filter(!grepl("Octocorallia", Comment, ignore.case = TRUE))

# filter data to leave just species level data (removing family and genus data):
stony_rm_fam_gen <- stony_rm_octa %>% filter(!is.na(Species) & Species != "" & !Species == "sp.")

# make a column for the species and genus combined
# species:
stony_rm_fam_gen <- stony_rm_fam_gen %>% mutate(species_id = paste(Genus, Species, sep = " "))
# all data:
stony_rm_octa <- stony_rm_octa %>% mutate(species_id = paste(Genus, Species, sep = " "))

# clean species_id column from mistakes (spaces that create 2 species from one species)
# just species:
stony_rm_fam_gen$species_id[stony_rm_fam_gen$species_id == "Cynarina  lacrymalis"] <- "Cynarina lacrymalis"
stony_rm_fam_gen$species_id[stony_rm_fam_gen$species_id == "lobophyllia hemprichii"] <- "Lobophyllia hemprichii"
# all data:
stony_rm_octa$species_id[stony_rm_octa$species_id == "Cynarina  lacrymalis"] <- "Cynarina lacrymalis"
stony_rm_octa$species_id[stony_rm_octa$species_id == "lobophyllia hemprichii"] <- "Lobophyllia hemprichii"

# long to wide format
# just species:
stony_wide_s <- stony_rm_fam_gen %>%
  group_by(site, depth, plot, species_id, lon, lat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = list(count = 0))
# all data:
stony_wide <- stony_rm_octa %>%
  group_by(site, depth, plot, species_id, lon, lat) %>%
  summarise(count = n(), .groups = 'drop') %>%
  pivot_wider(names_from = species_id, values_from = count, values_fill = list(count = 0))

# change the rows with NA to 0:
# just species:
stony_wide_s[is.na(stony_wide_s)] <- 0
# all data:
stony_wide[is.na(stony_wide)] <- 0
 
# save:
# species
setwd(wd_processed_data)
write.csv(stony_wide_s, file = "stony_wide_species_data.csv")
# all data:
setwd(wd_processed_data)
write.csv(stony_wide, file = "stony_wide_data.csv")

# for Tom:
# species <- as.data.frame(unique(stony_rm_octa$species_id))
# species <- species %>% rename (species_list = `unique(stony_rm_octa$species_id)`)
# setwd(wd_processed_data)
# write.csv(species, file = "stony_species.csv")

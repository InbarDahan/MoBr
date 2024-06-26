

            # combining Shahar's data sets

# packages:
library(dplyr)

# _______________________________________

# wd's:
wd_raw_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw"
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed"

# _______________________________________

# set wd:
setwd(wd_raw_data)

# _______________________________________

# read:
fish_data_no_lon_lat <- read.csv("fish_data_no_lon_lat.csv")
lon_lat_red <- read.csv("shahar_lon_lat_red.csv") # just for the red sea, but includes habitat
lon_lat_both <- read.csv("shahar_lon_lat_both.csv") # for both seas but without habitat

# _______________________________________

# delete irrelevant columns:
fish_data_no_lon_lat <- fish_data_no_lon_lat %>% dplyr::select(OpCode, species, MaxN, depth, temperature, sea)
lon_lat_red <- lon_lat_red %>% dplyr::select(OpCode, lat_y, lon_x, substrate)

# _______________________________________

# left_join the two data sets by the OpCode column:

# only red lon\lat, with habitat:
fish_data_lon_lat_red <- merge(fish_data_no_lon_lat, lon_lat_red, by = "OpCode", all.x =TRUE)

# all lon\lat, without habitat:
fish_data <- merge(fish_data_no_lon_lat, lon_lat_both, by = c("OpCode", "sea"), all.x =TRUE)

#removing samples from Sdot Yam from wide med:
fish_data <- fish_data[-which(fish_data$OpCode == "isrsdot201119A" | fish_data$OpCode == "isrsdot201119B"),]
# _______________________________________

# number of species at each sea:

med <- filter(fish_data, sea == "med")
n_distinct(med$species) # --- 70

red <- filter(fish_data, sea == "red")
n_distinct(red$species) # --- 173

length(intersect(med$species, red$species)) # --- 13 intersecting species

# _______________________________________

setwd(wd_processed_data)
write.csv(red, file = "long_red.csv")
write.csv(med, file = "long_med.csv")
write.csv(fish_data, file = "long_opCode.csv")

# _______________________________________

# transform into a wide format

# just red sea species
wide_red <- spread(data = red, key= species, value= MaxN, fill=0)
# just med sea species
wide_med <- spread(data = med, key= species, value= MaxN, fill=0)
# combined species
wide_opCode <- spread(data = fish_data, key= species, value= MaxN, fill=0)

# _______________________________________

# save_csv - wide formate:
setwd(wd_processed_data)
write.csv(wide_red, file = "wide_red.csv")
write.csv(wide_med, file = "wide_med.csv")
write.csv(wide_opCode, file = "wide_opCode.csv")
# _______________________________________

           # df coordinates red sea only, with habitat:

# transform into a wide format
wide_opCode_coo_red <- spread(data = fish_data_lon_lat_red, key= species, value= MaxN, fill=0)

# save_csv - without med coordinates:
setwd(wd_processed_data)
write.csv(wide_opCode_coo_red, file = "wide_opCode_coo_red.csv")

# - - - - - - - - - 



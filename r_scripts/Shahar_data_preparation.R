

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
  
# Shahar's data:
fish_data_no_lon_lat <- read.csv("fish_data_no_lon_lat.csv")
lon_lat_red <- read.csv("shahar_lon_lat_red.csv") # just for the red sea, but includes habitat
lon_lat_both <- read.csv("shahar_lon_lat_both.csv") # for both seas but without habitat

# _______________________________________

# delete irrelevant columns:
fish_data_no_lon_lat <- fish_data_no_lon_lat %>% select(OpCode, species, MaxN, depth, temperature, sea)
lon_lat_red <- lon_lat_red %>% select(OpCode, lat_y, lon_x, substrate)

# _______________________________________

# left_join the two data sets by the OpCode column:

# only red lon\lat, with habitat:
fish_data_lon_lat_red <- merge(fish_data_no_lon_lat, lon_lat_red, by = "OpCode", all.x =TRUE)

# all lon\lat, without habitat:
fish_data <- merge(fish_data_no_lon_lat, lon_lat_both, by = c("OpCode", "sea"), all.x =TRUE)

# _______________________________________


# df coordinates red sea only, with habitat:

# transform into a wide format
wide_opCode_coo_red = spread(data = fish_data_lon_lat_red, key= species, value= MaxN, fill=0)

# save_csv - without med coordinates:
setwd(wd_processed_data)

write.csv(wide_opCode_coo_red, file = "wide_opCode_coo_red.csv")

# - - - - - - - - - 

# df coordinates red sea only, with habitat:

# transform into a wide format
wide_opCode = spread(data = fish_data, key= species, value= MaxN, fill=0)

# save_csv - without med coordinates:
setwd(wd_processed_data)

write.csv(wide_opCode, file = "wide_opCode.csv")






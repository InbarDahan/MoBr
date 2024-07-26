
                 # sponges data preparations

library(dplyr)        
library(ggplot2)
library(mapview)

wd_raw_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw" 
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 

# read data
setwd(wd_raw_data)
sponges <- read.csv("sponges_speceis_data.csv")

# add coordination for the IUI and oil Getty:
sponges$lon[sponges$Site == "IUI"] <- "34.9174723"
sponges$lat[sponges$Site == "IUI"] <- "29.5018611"
sponges$lon[sponges$Site == "OJ"] <- "34.9323889"
sponges$lat[sponges$Site == "OJ"] <- "29.5226944"

# add a grouping var ("site") for the analysis (might change after the exploration of the data part:
sponges$depth_group <- sponges$Depth

# move the env columns to the beginning:
sponges <- sponges %>% dplyr::select(depth_group, lon, lat, everything())

# change the cloumn name from Depth to depth: 
sponges <- sponges %>% rename(depth = Depth)

# change the rows with NA to 0:
sponges[is.na(sponges)] <- 0

# save:
setwd(wd_processed_data)
write_csv(sponges, file = "sponges_data.csv")


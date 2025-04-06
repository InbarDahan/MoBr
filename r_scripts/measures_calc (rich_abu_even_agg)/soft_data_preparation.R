

                        # soft data preparations

library(dplyr)
library(ggplot2)
library(mapview)

# set wd:
wd_raw_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw" 
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 

# read data:
setwd(wd_raw_data)
soft_corals <- read.csv("ErezShoham_shallow_transects.csv")  

# add sample (id), lon and lat columns with NA:
soft_corals$lon <- NA
soft_corals$lat <- NA
soft_corals$sample <- 1:nrow(soft_corals)

# move this columns after the Site column:
soft_corals <- soft_corals %>% relocate(lon, lat, sample, .after = Site)

# rename Depth and Site to depth and site to fit common format:
soft_corals <- soft_corals %>% rename(depth = Depth)
soft_corals <- soft_corals %>% rename(site = Site)

#save processed data:
setwd(wd_processed_data)
write.csv(soft_corals, file = "soft_data.csv")

















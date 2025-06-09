

             # Plotting the sampling area

# -----

library(tidyverse)    # organizing data
library(raster)       # for raster object
library(sf)           # creating sf objects
library(ggplot2)      # creating plots
library(RColorBrewer) # for color palate

# _____________________________________

# wd:
wd_bay <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw/aqaba_bay_plot"   # map of Aqaba bay - for plotting
wd_occ <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_soft_coor <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/raw/soft_coor"

# _____________________________________

           # read and prepare data
# -----

# ensemble model for boundary
setwd(wd_bay)
pred_ens <- raster("final_sdm_percent.tif") 

# -----

# Aqaba bay map for the simple plots

aqaba_bay <- st_read(layer = "iho", dsn = wd_bay)
box <- st_bbox(pred_ens) + c(0.025, 0.01, 0, 0)                                  # creating the study area based on the final model and a buffer zone - numbers are how much to ass or crop in degrees
eilat <- st_crop(aqaba_bay, box)                                                 # cropping the Aqaba bay map to zoom into our area - just the sea
box_poly <- st_polygon(list(cbind(c(box[1], box[3], box[3], box[1], box[1]),     # creating  square polygon around the sea
                                  c(box[2], box[2], box[4], box[4], box[2]))))
box_poly <- st_sfc(box_poly, crs = crs(eilat))  

# -----

setwd(wd_occ)
occ <- read.csv("combined_long.csv")

occ <- occ %>% dplyr::select(-X)

occ_clean <- occ[!is.na(occ$lon) & !is.na(occ$lat), ]                            # removing NA's from lon lat cloumns

occ_clean$taxon[occ_clean$taxon == 'fish'] <- 'Fish'
occ_clean$taxon[occ_clean$taxon == 'sponges'] <- 'Porifera'
occ_clean$taxon[occ_clean$taxon == 'stony'] <- 'Scleractinia'  
occ_clean$taxon[occ_clean$taxon == 'soft'] <- 'Octocorallia'

occ_clean$taxon <- factor(occ_clean$taxon, levels = c("Fish", "Porifera", "Scleractinia", "Octocorallia"))

occ_sf <- st_as_sf(occ_clean, coords = c('lon', 'lat'), crs = crs(box_poly)) 

# _____________________________________

# plots:

# -----

# all occ together:
ggplot(eilat) + 
  geom_sf(fill = "#8ccfeb") + 
  geom_sf(data = occ_sf, size = 2)

# -----

# adding 4 general location to the soft corals data:

# seperating to the 4 taxonomic groups
#p_methods <- 
  
# plots:
# -----
# VERSION 1: Legend on top with fewer coordinate labels (grid lines remain)
ggplot(eilat) + 
  geom_sf(fill = "#8ccfeb") + 
  geom_sf(data = occ_sf, aes(color = taxon), size = 2) +
  scale_color_manual(values = c('Fish' = "#4B0082", 
                                'Porifera' = "#FF4500", 
                                'Scleractinia' = "#FFA500", 
                                'Octocorallia' = "#FF69B4")) +
  scale_x_continuous(breaks = c(34.92, 34.96, 35.00)) +
  scale_y_continuous(breaks = c(29.48, 29.52)) +
  theme_minimal() +
  theme(
       legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),  
       plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) +
  labs(color = "taxon")

# VERSION 2: More similar to original - just fewer coordinate labels, no other changes
ggplot(eilat) + 
  geom_sf(fill = "#8ccfeb") + 
  geom_sf(data = occ_sf, aes(color = taxon), size = 2) +
  scale_color_manual(values = c("Fish" = "#4B0082", 
                                "Porifera" = "#FF4500", 
                                "Scleractinia" = "#FFA500", 
                                "Octocorallia" = "#FF69B4")) +
  scale_x_continuous(breaks = c(34.92, 34.96, 35.00)) +
  scale_y_continuous(breaks = c(29.48, 29.52))

# ggplot(eilat) + 
#   geom_sf(fill = "#8ccfeb") + 
#   geom_sf(data = occ_sf,  aes(color = taxon), size = 2.5) +
#   scale_color_manual(values = c("#FEE0B6", "#FB6A4A", "#9970AB", "#E7A6C7")) +
#   theme_bw() 

# -----
pal <-  mapviewPalette("mapviewSpectralColors") #from mapview doc. example
fun_color_range <- colorRampPalette(c('lightblue', '#e31a1c' ,'#f768a1', '#feb24c')) 
mapView(occ_sf, zcol = "taxon", legend = TRUE, layer.name = 'Taxon', cex = 5, alpha.regions = 0.2, col.regions = fun_color_range)

   




# ensemble model for boundary
setwd(wd_bay)
pred_ens <- raster("final_sdm_percent.tif") 
# -----
# Aqaba bay map for the simple plots
aqaba_bay <- st_read(layer = "iho", dsn = wd_bay)
box <- st_bbox(pred_ens) + c(0.025, 0.01, 0, 0)                                  # creating the study area based on the final model and a buffer zone - numbers are how much to ass or crop in degrees
eilat <- st_crop(aqaba_bay, box)                                                 # cropping the Aqaba bay map to zoom into our area - just the sea
box_poly <- st_polygon(list(cbind(c(box[1], box[3], box[3], box[1], box[1]),     # creating  square polygon around the sea
                                  c(box[2], box[2], box[4], box[4], box[2]))))
box_poly <- st_sfc(box_poly, crs = crs(eilat))  
# -----
setwd(wd_occ)
occ <- read.csv("combined_long.csv")
occ <- occ %>% dplyr::select(-X)
occ_clean <- occ[!is.na(occ$lon) & !is.na(occ$lat), ]                            # removing NA's from lon lat cloumns
occ_clean$taxon <- factor(occ_clean$taxon, levels = c("Fish", "Porifera", "Scleractinia", "Octocorallia"))
occ_sf <- st_as_sf(occ_clean, coords = c('lon', 'lat'), crs = crs(box_poly)) 
# _____________________________________
# plots:
# -----
# VERSION 1: ggplot with compact layout
ggplot(eilat) + 
  geom_sf(fill = "#8ccfeb") + 
  geom_sf(data = occ_sf, aes(color = taxon), size = 2) +
  scale_color_manual(values = c("Fish" = "#4B0082", 
                                "Porifera" = "#FF4500", 
                                "Scleractinia" = "#FFA500", 
                                "Octocorallia" = "#FF69B4")) +
  scale_x_continuous(breaks = c(34.92, 34.96, 35.00)) +
  scale_y_continuous(breaks = c(29.48, 29.52)) +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 12),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
  labs(color = "taxon")



# VERSION 2: More compact ggplot version
ggplot(eilat) + 
  geom_sf(fill = "#8ccfeb") + 
  geom_sf(data = occ_sf, aes(color = taxon), size = 2) +
  scale_color_manual(values = c("Fish" = "#4B0082", 
                                "Porifera" = "#FF4500", 
                                "Scleractinia" = "#FFA500", 
                                "Octocorallia" = "#FF69B4")) +
  scale_x_continuous(breaks = c(34.92, 34.96, 35.00)) +
  scale_y_continuous(breaks = c(29.48, 29.52)) +
  theme(axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.margin = margin(0.2, 0.2, 0.2, 0.2, "cm")) +
  labs(color = "taxon")




# -----

# This code was moves after changes to the "soft data preperation" r script

# setwd(wd_soft_coor)
# soft_coor <- read.csv("soft_coor.csv")
# 
# # Step 1: Filter the relevant subset
# soft_rows <- occ %>% filter(taxon == "soft")
# 
# # Step 2: Repeat the 4 coordinates to match the number of soft_rows (80)
# set.seed(123)  # For reproducibility
# coords_to_assign <- soft_coor[sample(1:4, size = nrow(soft_rows), replace = TRUE), ]
# 
# # If you want **equal quarters**, not random sampling with replacement:
# repeated_indices <- rep(1:4, length.out = nrow(soft_rows))
# coords_to_assign <- soft_coor[sample(repeated_indices), ]
# 
# # Step 3: Add the new coordinates to the soft_rows
# soft_rows$lon <- coords_to_assign$lon
# soft_rows$lat <- coords_to_assign$lat
# 
# # Step 4: Replace the old soft rows in occ with the new ones
# occ_updated <- occ
# occ_updated[occ$taxon == "soft", ] <- soft_rows
# 
# occ_up_clean <- occ_updated[!is.na(occ_updated$lon) & !is.na(occ_updated$lat), ]                            # removing NA's from lon lat cloumns
# 
# 
# occ_updated_sf <- st_as_sf(occ_up_clean, coords = c('lon', 'lat'), crs = crs(box_poly)) 




























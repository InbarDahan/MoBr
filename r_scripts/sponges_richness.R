

                  ###  sponges Richness and abundance   ###

# Does the sponges richness changes across depths layers ?

library(vegan)        # running rarefaction analysis
library(tidyverse)    # organizing the data 
library(tidyr)        # organizing the data
library(dplyr)        # organizing the data
library(plotrix)      #
library(rareNMtests)  # running rarefaction analysis
library(mobr)         # running rarefaction analysis - for gradients
library(mapview)      # visualization
library(sf)           # visualization
library(raster)       # for raster object
library(RColorBrewer)    # color palettes
require(ggmisc)

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 

# read data:
setwd(wd_processed_data)
sponges <- read.csv("sponges_data.csv")

# _______________________________________________________________

# preliminary exploration #
# _______________________________________________________________


# define the first column of species
first_species <- 8

# create species metric
species_metric <- sponges[,first_species:length(sponges)]

# _______________________________________________________________

# depth-ranges
range(sponges$depth, na.rm=TRUE) # 5 - 70 m: 5,10,20,50,60,70 m

# depths histogram:
sponges %>% 
  ggplot()+
  aes(x = depth)+
  geom_histogram(bins = 20) + ggtitle("histogram of samples at each depth")  # most of the depths have 3-5 observations\opcodes

# _______________________________________________________________

# map view
# make sf objects:
sponges_sf <- st_as_sf(sponges, coords = c('lon', 'lat'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# visualize:
mapView(sponges_sf, zcol = "depth", legend = TRUE, layer.name = 'Depth')

# _______________________________________________________________

# sampling effort

# red: add 2 additional columns: abundance and richness per sample:
sponges_sum <- sponges %>% 
  mutate(abundance = rowSums(sponges[,first_species:length(sponges)])) %>%
  mutate(richness = rowSums(sponges[, first_species:ncol(sponges)] > 0)) %>%
  dplyr::select(quad, depth_group, depth, Site, lon, lat, richness, abundance)

# check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(sponges_sum$abundance) # 1
max(sponges_sum$abundance) # 252
# 252 folds variation in the alpha scale (sample scale)

# sum richness:
min(sponges_sum$richness) # 0
max(sponges_sum$richness) # 12

# _______________________________________________________________

# plot abundance per sample ~ depth (wide_data)
# Fit the linear model
model_a <- lm(abundance ~ depth, data = sponges_sum)
# Extract the equation and R-squared
eq_a <- paste0("y = ", round(coef(model_a)[1], 2), " ", round(coef(model_a)[2], 2), "x")
r_squared_a <- paste0("R² = ", round(summary(model_a)$r.squared, 2))

# Create the plot
ggplot(sponges_sum, aes(x = depth, y = abundance)) +
  labs(y = "Abundance", x = "Depth") +
  ggtitle("                                   Abundance Across Depths (Sponges)") +
  geom_point(size = 3, color = "orange") + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 17.5, y = Inf, label = eq_a, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 11, y = 37, label = r_squared_a, hjust = 1.1, vjust = 3, size = 4.5, color = "black")

# -----------

# plot richness per sample ~ depth (wide_data)

# Fit the linear model
model_r_s <- lm(richness ~ depth, data = sponges_sum)
# Extract the equation and R-squared
eq_r_s <- paste0("y = ", round(coef(model_r_s)[1], 2), " ", round(coef(model_r_s)[2], 2), "x")
r_squared_r_s <- paste0("R² = ", round(summary(model_r_s)$r.squared, 2))

# Create the plot
ggplot(sponges_sum, aes(x = depth, y = richness)) +
  labs(y = "Richness", x = "Depth") +
  ggtitle("                                   Richness Across Depths (Sponges)") +
  geom_point(size = 3, color = "orange") + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 18, y = Inf, label = eq_r_s, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 12.5, y = 11.5, label = r_squared_r_s, hjust = 1.1, vjust = 3, size = 4.5, color = "black")

# _______________________________________________________________
# alternative plots:
        # Box plot - abundance over depth layers

      # shows mean and quarters, all samples distribution

# red:                                                  # we have outliers (step 1 from zuur 2010)
boxplot(abundance~depth,                               # we do not have homogeneity of variance (step 2 from zuur 2010)
        data=sponges_sum,                             # we do not have normal distributions (step 3 from zuur 2010)
        main= "Abundance Across Depth Layers",       # we have zero inflated data (step 4)
        xlab="Depth ranges (m)",                             # if I will add more explanatory variable except depth - I will need to check for co-linearity (the relationship between a few x) (step 5)
        ylab="Abundance",                    # the relationship between the abu ~ depth is  (step 6 from zuur 2010)
        col="indianred2",
        border="black", add = F
)

# - - -
      # Bar plot  -  mean abundance over depth layers

fun_color_range_red <- brewer.pal(6, "Reds") # doensn't work for more than 9

ggplot(sponges_sum, aes(x = depth, y = abundance)) + 
  stat_summary(fun = "mean", geom = "bar", fill = fun_color_range_red)+
  ggtitle("avr abu ~ depth")

# _______________________________________________________________

# Box plot - richness 

# shows mean and quarters, all samples distribution
boxplot(richness~depth,
        data=sponges_sum,
        main="Richness Across Depth Layers",
        xlab="Depth layer",
        ylab="Richness per sample",
        col="indianred2",
        border="black"
)


# -----------

# Bar plot - mean richness

ggplot(sponges_sum, aes(x = depth, y = richness)) + 
  stat_summary(fun = "mean", geom = "bar", fill = fun_color_range_red)+
  ggtitle("avr rich ~ depth")
# _______________________________________________________________



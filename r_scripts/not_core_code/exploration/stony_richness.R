

           ###  Stony Corals Richness and abundance - NOT UNUPDATED ###

# Does the stony richness changes across depths layers ?

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
stony <- read.csv("stony_wide_data.csv")
stony  = stony %>% dplyr::select(-X) 

# _______________________________________________________________

                  # preliminary exploration #
# _______________________________________________________________


# define the first column of species
first_species <- 6 

# create species metric
species_metric <- stony[,first_species:length(stony)]

# _______________________________________________________________

# depth-ranges
range(stony$depth, na.rm=TRUE) # 2 - 55 m: 2, 5, 15, 30, 45, 55 m

# depths histogram:
stony %>% 
  ggplot()+
  aes(x = depth)+
  geom_histogram(bins = 20) + ggtitle("IUI - histogram of samples at each depth")  # most of the depths have 3-5 observations\opcodes

# _______________________________________________________________
# 
# # map view
# # make sf objects:
# stony_sf <- st_as_sf(stony, coords = c('lon', 'lat'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# 
# # visualize:
# mapView(stony_sf, zcol = "depth", legend = TRUE, layer.name = 'Depth')

# _______________________________________________________________

# general sampling effort

# red: add 2 additional columns: abundance and richness per sample:
stony_sum <- stony %>% 
  mutate(abundance = rowSums(stony[,first_species:length(stony)])) %>%
  mutate(richness = rowSums(stony[, first_species:ncol(stony)] > 0)) %>%
  dplyr::select(site ,depth, plot, richness, abundance)

# sum abundance:
min(stony_sum$abundance) # 74
max(stony_sum$abundance) # 325
# 252 folds variation in the alpha scale (sample scale)

# sum richness:
min(stony_sum$richness) # 29
max(stony_sum$richness) # 54

# _______________________________________________________________

# plot abundance per sample ~ depth (wide_data)
# Fit the linear model
model_a <- lm(abundance ~ depth, data = stony_sum)
# Extract the equation and R-squared
eq_a <- paste0("y = ", round(coef(model_a)[1], 2), " ", round(coef(model_a)[2], 2), "x")
r_squared_a <- paste0("R² = ", round(summary(model_a)$r.squared, 2))

# Create the plot
ggplot(stony_sum, aes(x = depth, y = abundance)) +
  labs(y = "Abundance", x = "Depth") +
  ggtitle("                                   Abundance Across Depths (stony)") +
  geom_point(size = 3, color = "orange") + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 17.5, y = Inf, label = eq_a, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 11, y = 37, label = r_squared_a, hjust = 1.1, vjust = 3, size = 4.5, color = "black")

# -----------

# plot richness per sample ~ depth (wide_data)

# Fit the linear model
model_r_s <- lm(richness ~ depth, data = stony_sum)
# Extract the equation and R-squared
eq_r_s <- paste0("y = ", round(coef(model_r_s)[1], 2), " ", round(coef(model_r_s)[2], 2), "x")
r_squared_r_s <- paste0("R² = ", round(summary(model_r_s)$r.squared, 2))

# Create the plot
ggplot(stony_sum, aes(x = depth, y = richness)) +
  labs(y = "Richness", x = "Depth") +
  ggtitle("                                   Richness Across Depths (stony)") +
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
        data=stony_sum,                             # we do not have normal distributions (step 3 from zuur 2010)
        main= "Stony Corals - Abundance Across Depth Layers",       # we have zero inflated data (step 4)
        xlab="Depth ranges (m)",                             # if I will add more explanatory variable except depth - I will need to check for co-linearity (the relationship between a few x) (step 5)
        ylab="Abundance",                    # the relationship between the abu ~ depth is  (step 6 from zuur 2010)
        col="yellow",
        border="black", add = F
)

# - - -
# Bar plot  -  mean abundance over depth layers

fun_color_range_red <- brewer.pal(5, "Reds") # doensn't work for more than 9

ggplot(stony_sum, aes(x = depth, y = abundance)) + 
  stat_summary(fun = "mean", geom = "bar", fill = fun_color_range_red)+
  ggtitle("avr abu ~ depth")

# _______________________________________________________________

# Box plot - richness 

# shows mean and quarters, all samples distribution
boxplot(richness~depth,
        data=stony_sum,
        main="Stony Corals - Richness Across Depth Layers",
        xlab="Depth layer",
        ylab="Richness per sample",
        col="yellow",
        border="black"
)


# -----------

# Bar plot - mean richness

ggplot(stony_sum, aes(x = depth, y = richness)) + 
  stat_summary(fun = "mean", geom = "bar", fill = fun_color_range_red)+
  ggtitle("avr rich ~ depth")
# _______________________________________________________________



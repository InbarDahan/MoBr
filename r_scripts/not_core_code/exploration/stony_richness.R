

           ###  Stony Corals Richness and abundance - NOT UNUPDATED ###

# Does the stony richness changes across depths layers ?

library(vegan)        # running rarefaction analysis
library(tidyverse)    # organizing the data 
library(tidyr)        # organizing the data
library(dplyr)        # organizing the data
library(plotrix)   
library(ggplot2)
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
first_species <- 7 

# create species metric
species_metric <- stony[,first_species:length(stony)]

# _______________________________________________________________

# depth-ranges
range(stony$depth, na.rm=TRUE) # 2 - 55 m: 2, 5, 15, 30, 45, 55 m

dev.off() 
# depths histogram:
ggplot(stony, aes(x = depth, fill = site)) +
  geom_histogram(bins = 20, position = "dodge", alpha = 0.5) +
  ggtitle("Histogram of Samples at Each Depth by Site") +
  theme_minimal()
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
min(stony_sum$abundance[stony_sum$site == "iui"]) #70
max(stony_sum$abundance[stony_sum$site == "iui"]) #301
min(stony_sum$abundance[stony_sum$site == "nr"])  #41
max(stony_sum$abundance[stony_sum$site == "nr"])  #256
# big variation among the alpha and jamma scales

# sum richness:
min(stony_sum$richness[stony_sum$site == "iui"]) #24
max(stony_sum$richness[stony_sum$site == "iui"]) #48
min(stony_sum$richness[stony_sum$site == "nr"])  #23
max(stony_sum$richness[stony_sum$site == "nr"])  #60
# _______________________________________________________________

# plot abundance per sample ~ depth (wide_data)

# Fit linear models for each site
model_iui <- lm(abundance ~ depth, data = filter(stony_sum, site == "iui"))
model_nr  <- lm(abundance ~ depth, data = filter(stony_sum, site == "nr"))
model_a <- lm(abundance ~ depth, data = stony_sum) # together

# Extract the equation and R-squared
#iui:
eq_iui <- paste0("iui: y = ", round(coef(model_iui)[1], 2), " + ", round(coef(model_iui)[2], 2), "x")
r2_iui <- paste0("R² = ", round(summary(model_iui)$r.squared, 2))
#nr
eq_nr <- paste0("nr: y = ", round(coef(model_nr)[1], 2), " + ", round(coef(model_nr)[2], 2), "x")
r2_nr <- paste0("R² = ", round(summary(model_nr)$r.squared, 2))
# together
eq_a <- paste0("y = ", round(coef(model_a)[1], 2), " ", round(coef(model_a)[2], 2), "x")

# Create the plot
ggplot(stony_sum, aes(x = depth, y = abundance, color = site)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("iui" = "orange", "nr" = "blue")) +
  labs(y = "Abundance", x = "Depth", title = "Abundance Across Depths (Stony Corals)") +
  annotate("text", x = Inf, y = Inf, label = eq_iui, hjust = 1.1, vjust = 2, size = 4.5, color = "orange") +
  annotate("text", x = Inf, y = Inf, label = r2_iui, hjust = 1.1, vjust = 4, size = 4.5, color = "orange") +
  annotate("text", x = -Inf, y = Inf, label = eq_nr, hjust = -0.1, vjust = 2, size = 4.5, color = "blue") +
  annotate("text", x = -Inf, y = Inf, label = r2_nr, hjust = -0.1, vjust = 4, size = 4.5, color = "blue") +
  theme_minimal()

ggplot(stony_sum, aes(x = depth, y = abundance)) +
  labs(y = "Abundance", x = "Depth") +
  ggtitle("                                   Abundance Across Depths (stony)") +
  geom_point(size = 3, color = "darkgreen") + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = 17.5, y = Inf, label = eq_a, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 11, y = 37, label = r_squared_a, hjust = 1.1, vjust = 3, size = 4.5, color = "black")
# -----------

# plot richness per sample ~ depth (wide_data)

# Fit the linear model
model_r_iui <- lm(richness ~ depth, data = filter(stony_sum, site == "iui"))
model_r_nr  <- lm(richness ~ depth, data = filter(stony_sum, site == "nr"))
model_r <- lm(richness ~ depth, data = stony_sum)
# Extract the equation and R-squared

eq_r_iui <- paste0("y = ", round(coef(model_r_iui)[1], 2), " ", round(coef(model_r_iui)[2], 2), "x")
r_squared_r_iui <- paste0("R² = ", round(summary(model_r_iui)$r.squared, 2))

eq_r_nr <- paste0("y = ", round(coef(model_r_nr)[1], 2), " ", round(coef(model_r_nr)[2], 2), "x")
r_squared_r_nr <- paste0("R² = ", round(summary(model_r_nr)$r.squared, 2))

eq_r_s <- paste0("y = ", round(coef(model_r)[1], 2), " ", round(coef(model_r)[2], 2), "x")
r_squared_r_s <- paste0("R² = ", round(summary(model_r)$r.squared, 2))

# Create the plot
ggplot(stony_sum, aes(x = depth, y = richness, color = site)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("iui" = "darkred", "nr" = "darkblue")) +
  labs(y = "Richness", x = "Depth", title = "Richness Across Depths (Stony Corals)") +
  annotate("text", x = Inf, y = Inf, label = eq_r_iui, hjust = 1.1, vjust = 2, size = 4.5, color = "orange") +
  annotate("text", x = Inf, y = Inf, label = r_squared_r_iui, hjust = 1.1, vjust = 4, size = 4.5, color = "orange") +
  annotate("text", x = -Inf, y = Inf, label = eq_r_nr, hjust = -0.1, vjust = 2, size = 4.5, color = "blue") +
  annotate("text", x = -Inf, y = Inf, label = r_squared_r_nr, hjust = -0.1, vjust = 4, size = 4.5, color = "blue") +
  theme_minimal()

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


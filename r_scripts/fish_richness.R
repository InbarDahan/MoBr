

               ###  Fish Richness and abundance   ###

     # Does the fish richness changes across depths layers ?

                    # Red and Med  

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
library(RColorBrewer) # color palettes
library(ggpmisc)      # trade line

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 

# read data:
setwd(wd_processed_data)
wide_med <- read.csv("wide_med.csv")
wide_red <- read.csv("wide_red.csv")
wide_opCode <- read.csv("wide_opCode.csv")

# removing the wired column x:
wide_med = wide_med %>% dplyr::select(-X)
wide_red = wide_red %>% dplyr::select(-X) 
wide_opCode = wide_opCode %>% dplyr::select(-X) 
# _______________________________________________________________

                # preliminary exploration #
# _______________________________________________________________

# define the first column of species
first_species <- 7

# create species metric
species_metric_med <- wide_med[,first_species:length(wide_med)]
species_metric_red <- wide_red[,first_species:length(wide_red)]

# _______________________________________________________________

                    # depth-ranges

# med
range(wide_med$depth, na.rm=TRUE) # 5.6 - 147 m
# depths histogram for med sea
wide_med %>% 
  ggplot()+
  aes(x = depth)+     
  geom_histogram() + ggtitle("med sea - samples count at each depth") 
 # some depths have up to 6 observations\opcode

# -----------

# red:
range(wide_red$depth, na.rm=TRUE) # 8.1 - 149 m
# depths histogram for red sea
wide_red %>% 
  ggplot()+
  aes(x = depth)+
  geom_histogram(bins = 6) + ggtitle("red sea - samples count at each depth")  # most of the depths have 3-5 observations\opcodes

# _______________________________________________________________

                      # map view
# make sf objects:
wide_med_sf <- st_as_sf(wide_med, coords = c('lon_x', 'lat_y'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) # WGS84 is in decimal degrees
wide_red_sf <- st_as_sf(wide_red, coords = c('lon_x', 'lat_y'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# visualize:
mapView(wide_med_sf, zcol = "depth", legend = TRUE,
        layer.name = 'Depth')
mapView(wide_red_sf, zcol = "depth", legend = TRUE,
        layer.name = 'Depth')

# _______________________________________________________________

          # sampling effort at each sea

wide_med_sum <- wide_med %>% 
  mutate(abundance = rowSums(wide_med[,first_species:length(wide_med)])) %>%
  mutate(richness = rowSums(wide_med[, first_species:ncol(wide_med)] > 0)) 

# med: add 2 additional columns: abundance and richness per sample:
wide_med_sum <- wide_med %>% 
  mutate(abundance = rowSums(wide_med[,first_species:length(wide_med)])) %>%
  mutate(richness = rowSums(wide_med[, first_species:ncol(wide_med)] > 0)) %>%
  dplyr::select(OpCode, depth, lon_x, lat_y, richness, abundance)


# check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(wide_med_sum$abundance) # 0
max(wide_med_sum$abundance) # 338
# 338 folds variation in the alpha scale (sample scale)

# sum richness:
min(wide_med_sum$richness) # 0
max(wide_med_sum$richness) # 15
# -----------

# red: add 2 additional columns: abundance and richness per sample:
wide_red_sum <- wide_red %>% 
  mutate(abundance = rowSums(wide_red[,first_species:length(wide_red)])) %>%
  mutate(richness = rowSums(wide_red[, first_species:ncol(wide_red)] > 0)) %>%
  dplyr::select(OpCode, depth, lon_x, lat_y, richness, abundance)

# sf:
sample_red_sf <- st_as_sf(wide_red_sum, coords = c('lon_x', 'lat_y'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# map of sample richness:
mapView(sample_red_sf, zcol = "richness", legend = TRUE,
        layer.name = 'Richness')

# check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(wide_red_sum$abundance) # 1
max(wide_red_sum$abundance) # 252
# 252 folds variation in the alpha scale (sample scale)

# sum richness:
min(wide_red_sum$richness) # 1
max(wide_red_sum$richness) # 42

# _______________________________________________________________

           # plot abundance per sample ~ depth (wide_data)

# med - points: 
ggplot(wide_med_sum,aes(x=depth,y=abundance))+
  geom_point(size = 1.5, col = "cyan") + ggtitle("med sea abundance ~ depth") 

# med - line: 
ggplot(wide_med_sum,aes(x=depth,y=abundance))+                               
  geom_line(size = 1.5, col = "cyan") + ggtitle("med sea abundance ~ depth") 

# -----------
# Fit the linear model
model <- lm(abundance ~ depth, data = wide_red_sum)
# Extract the equation and R-squared
eq <- paste0("y = ", round(coef(model)[1], 2), " ", round(coef(model)[2], 2), "x")
r_squared <- paste0("R² = ", round(summary(model)$r.squared, 2))

# Create the plot
ggplot(wide_red_sum, aes(x = depth, y = abundance)) +
  labs(y = "Abundance", x = "Depth") +
  ggtitle("                         Fish - Abundance Across Depths") +
  geom_point(size = 3, color = "indianred2") + 
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text", x = Inf, y = 250, label = eq, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 145, y = 230, label = r_squared, hjust = 1.1, vjust = 3, size = 4.5, color = "black")

# red - line: 
ggplot(wide_red_sum,aes(x=depth,y=abundance))+
  geom_line(size = 1.5, col = "indianred2") + ggtitle("Red Sea abundance ~ depth") 

# -----------

           # plot richness per sample ~ depth (wide_data)

# med - points: 
ggplot(wide_med_sum,aes(x=depth,y=richness))+
  geom_point(size = 1.5, col = "cyan") + ggtitle("med sea richness ~ depth") 

# med - line: 
ggplot(wide_med_sum,aes(x=depth,y=richness))+
  geom_line(size = 1.5, col = "cyan") + ggtitle("med sea richness ~ depth") 

# -----------

# Fit the linear model
model_r <- lm(richness ~ depth, data = wide_red_sum)
# Extract the equation and R-squared
eq_r <- paste0("y = ", round(coef(model_r)[1], 2), " ", round(coef(model_r)[2], 2), "x")
r_squared_r <- paste0("R² = ", round(summary(model)$r.squared, 2))

# red - points: 
ggplot(wide_red_sum,aes(x=depth,y=richness))+
  labs(y= "Richness", x = "Depth")+
  geom_point(size = 3, col = "indianred2") + ggtitle("                                       Richness Across Depths (Fish)")+
  geom_smooth(method = "lm", se = TRUE, color = "black") + 
  annotate("text", x = Inf, y = 50, label = eq_r, hjust = 1.1, vjust = 2, size = 4.5, color = "black") +
  annotate("text", x = 145, y = 46, label = r_squared_r, hjust = 1.1, vjust = 3, size = 4.5, color = "black")

# red - line: 
ggplot(wide_red_sum,aes(x=depth,y=richness))+
  geom_line(size = 1.5, col = "indianred2") + ggtitle("red sea richness ~ depth")

# _______________________________________________________________

# I noticed I have samples with very low abundances (0, 1, 2...) in the two seas
# since we are interested in exploring the effect of depth, I wanted to merge samples from the same depth strip 
# but what should be the width of the bin? 10m? 20m?

# _______________________________________________________________

              # Convert the data into 20 m depth bins
# _______________________________________________________________

# I chose 20m because I saw other papers group in bins of 20 m, but 10 is also an option

med_depth_bins <- wide_med
red_depth_bins <- wide_red

# 5.6 - 147 m to 8 layers
med_depth_bins$depth <- cut(med_depth_bins$depth,                                          
                       breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 147),
                       labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-147'))

# 15 meters bins:
red_depth_bins$depth <- cut(red_depth_bins$depth,                                         
                            breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
                            labels=c('0-15', '16-30', '31-45', '46-60', '61-75','76-90', '91-105', '106-120', '121-135', '136-149'))

# 20 meters bins - 8.1 - 149 m to 8 layers
# red_depth_bins$depth <- cut(red_depth_bins$depth,                                         
#                        breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 149),
#                        labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

# _______________________________________________________________

          # df of richness and abundance for each sample at different depths

# to get the richness in each depth bin, without recounting the same species twice
# I will move to long format and then will group by depth and species and count them:

#med:
med_bins_long <- med_depth_bins %>%                                                      # the df to convert
  dplyr::select(-c(2, 4, 5, 6)) %>%                                                      # columns to remove
  tidyr::pivot_longer(cols = !c(depth, OpCode), names_to = "species", values_to = "count") %>%  # transforming to wide
  filter(!count == 0) %>%                                                                #removing rows with 0 abundance so species will not be counted when they are not present
  group_by(depth, OpCode) %>%                                                            # grouping by sea and species
  summarise(sample_r = n_distinct(species), sample_abu = sum(count))                     # richness per depth bin - no repeated species
  
#red:
red_bins_long <- red_depth_bins %>%                                                      # the df to convert
  dplyr::select(-c(2, 4, 5, 6)) %>%                                                      # columns to remove
  tidyr::pivot_longer(cols = !c(depth, OpCode), names_to = "species", values_to = "count") %>%  # transforming to wide
  filter(!count == 0) %>%                                                                #removing rows with 0 abundance so species will not be counted when they are not present
  group_by(depth, OpCode) %>%                                                            # grouping by sea and species
  summarise(sample_r = n_distinct(species), sample_abu = sum(count))                     # richness per depth bin - no repeated species

# _______________________________________________________________

            # Box plot - abundance over depth layers

      # shows mean and quarters, all samples distribution
                                                                      # charectaristics:

# med:                                                 # we have outliers (step 1 from zuur 2010)
boxplot(sample_abu~depth,                              # we do not have homogeneity of variance (step 2 from zuur 2010)
        data=med_bins_long,                            # we do not have normal distributions (step 3 from zuur 2010)
        main= "Boxplot: sample abundance ~ depth",     # we have zero inflated data (step 4)
        xlab="depth layer",                            # if I will add more explanatory variable except depth - I will need to check for co-linearity (the relationship between a few x) (step 5) - for example richness is dependent on abundance
        ylab="abundance per sample",                   # the relationship between the abu ~ depth is  (step 6 from zuur 2010)
        col="cyan",                                    # 
        border="black"
)

par(mfrow=c(1,1))

# Remark - to add the n samples above the graph
# red:                                                  # we have outliers (step 1 from zuur 2010)
boxplot(sample_abu~depth,                               # we do not have homogeneity of variance (step 2 from zuur 2010)
        data=red_bins_long,                             # we do not have normal distributions (step 3 from zuur 2010)
        main= "Fish - Abundance Across Depth Layers",       # we have zero inflated data (step 4)
        xlab="Depth ranges (m)",                             # if I will add more explanatory variable except depth - I will need to check for co-linearity (the relationship between a few x) (step 5)
        ylab="Abundance",                    # the relationship between the abu ~ depth is  (step 6 from zuur 2010)
        col="indianred2",
        border="black", add = F
)

# -----------

            # Bar plot  -  mean abundance over depth layers

# blue color ramp:
fun_color_range_med <- brewer.pal(8, "Blues")
fun_color_range_red <- brewer.pal(8, "Reds") # doensn't work for more than 9

ggplot(med_bins_long, aes(x = depth, y = sample_abu)) + 
  stat_summary(fun = "mean", geom = "bar",fill = fun_color_range_med, col = "black")+
  ggtitle("Med: avr abu ~ depth")

ggplot(red_bins_long, aes(x = depth, y = sample_abu)) + 
  stat_summary(fun = "mean", geom = "bar", col = "black")+
  ggtitle("Red: avr abu ~ depth")

# _______________________________________________________________

                     # Box plot - richness 

      # shows mean and quarters, all samples distribution

# med:
boxplot(sample_r~depth,
        data=med_bins_long,
        main="Richness Across Depth Layers",
        xlab="Depth layer",
        ylab="Richness per sample",
        col="cyan",
        border="black"
)

# red:
boxplot(sample_r~depth,
        data=red_bins_long,
        main="Fish - Richness Across Depth Layers",
        xlab="Depth layer",
        ylab="Richness per sample",
        col="indianred2",
        border="black"
)


# -----------

                   # Bar plot - mean richness

ggplot(med_bins_long, aes(x = depth, y = sample_r)) + 
  stat_summary(fun = "mean", geom = "bar",fill = fun_color_range_med, col = "black")+
  ggtitle("Med: avr rich ~ depth")

ggplot(red_bins_long, aes(x = depth, y = sample_r)) + 
  stat_summary(fun = "mean", geom = "bar", col = "black")+
  ggtitle("Red: avr rich ~ depth")
# _______________________________________________________________

     # total richness over each depth range using species identities

#med:
long_med_bins <- med_depth_bins %>%     # the df to convert
  dplyr::select(-c(1, 2, 4, 5, 6)) %>%     # columns to remove
  pivot_longer(cols = !c(depth), names_to = "species", values_to = "count") %>% # transforming to wide
  filter(!count == 0) %>%               #removing rows with 0 abundance so species will not be counted when they are not present
  group_by(depth) %>%                   # grouping by depth layer
  summarise(depth_layer_r = n_distinct(species), layer_abu = sum(count)) # richness per depth bin - no repeated species

#red:
long_red_bins <- red_depth_bins %>%     # the df to convert
  dplyr::select(-c(2, 4, 5, 6)) %>%     # columns to remove
  pivot_longer(cols = !c(OpCode, depth), names_to = "species", values_to = "count") %>% # transforming to wide
  filter(!count == 0) %>%               #removing rows with 0 abundance so species will not be counted when they are not present
  group_by(depth) %>%                   # grouping by depth
  summarise(depth_layer_r = n_distinct(species), layer_abu = sum(count), sample = n_distinct(OpCode)) # richness per depth bin - no repeated species

# - - - - - - - - - - - - 

          # Bar plot  - total richness at each depth layers

ggplot(long_med_bins, aes(x = depth, y = depth_layer_r)) + 
  stat_summary(geom = "bar",fill = fun_color_range_med, col = "black")+
  ggtitle("Med: total rich ~ depth")

ggplot(long_red_bins, aes(x = depth, y = depth_layer_r)) + 
  stat_summary(geom = "bar",col = "black")+
  ggtitle("Red: total rich ~ depth")

# _______________________________________________________________

                # The problem with our plots:

# we still don't know if the diff' between the depth bins
    # richness and abundance arise from sampling effort biases or real ecological reasons
       # to account for it we will create rarefaction plots
# - Should be done for both binned depth and continuous depth

# _______________________________________________________________

               # individual based rarefaction
              
                     # OpCode scale 

# med:

ind_based_rare_all_med = list() # empty list

for (i in 1:nrow(med_depth_bins)) { # for each row...
  
  one_opCode = med_depth_bins[i,] # filter

  one_opCode_sp_matrix = one_opCode[,first_species:ncol(one_opCode)] # create sp_matrix

  if(rowSums(one_opCode_sp_matrix > 0)){

  rarefaction = rarefaction.individual(one_opCode_sp_matrix, method = "sample-size", q = 0) # apply rarefaction function

  rarefaction$depth = one_opCode$depth
  
  rarefaction$OpCode = one_opCode$OpCode # add the opcode_id to the data

  ind_based_rare_all_med[[i]] = rarefaction } # save in my list

}

ind_based_rare_all_med = bind_rows(ind_based_rare_all_med) # convert from list to data frame

colnames(ind_based_rare_all_med) <- c("Individuals","Richness","Depth_layer", "sample") # change columns names

# -----------

# plot

ggplot(ind_based_rare_all_med, aes(x=Individuals, y=Richness, group = sample ,color = Depth_layer))+
  geom_line()

# facet:
ggplot(ind_based_rare_all_med, aes(x=Individuals, y=Richness, group = sample, color = Depth_layer))+
  geom_line() +
  facet_wrap(~Depth_layer)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# red:
ind_based_rare_all_red = list() # empty list

for (i in 1:nrow(red_depth_bins)) { # for each row...
  
  one_opCode = red_depth_bins[i,] # filter
  
  one_opCode_sp_matrix = one_opCode[,first_species:ncol(one_opCode)] # create sp_matrix
  
  if(rowSums(one_opCode_sp_matrix > 0)){
    
    rarefaction = rarefaction.individual(one_opCode_sp_matrix, method = "sample-size", q = 0) # apply rarefaction function
    
    rarefaction$depth = one_opCode$depth
    
    rarefaction$OpCode = one_opCode$OpCode # add the opcode_id to the data
    
    ind_based_rare_all_red[[i]] = rarefaction } # save in my list
  
}

ind_based_rare_all_red = bind_rows(ind_based_rare_all_red) # convert from list to data frame

colnames(ind_based_rare_all_red) <- c("Individuals","Richness","Depth_layer", "sample") # change columns names

# -----------

# plot
ggplot(ind_based_rare_all_red, aes(x=Individuals, y=Richness, group = sample ,color = Depth_layer))+
  geom_line()

# facet:
ggplot(ind_based_rare_all_red, aes(x=Individuals, y=Richness, group = sample, color = Depth_layer))+
  geom_line() +
  facet_wrap(~Depth_layer)

# _______________________________________________________________

                 # individual based rarefaction

                       # depth bin scale 

# Now we will sum the total individuals of each species in each `depth bin`. 

# med:
group_data_med = med_depth_bins %>% 
  dplyr::select(depth, 7:ncol(med_depth_bins)) %>% 
  group_by(depth) %>% #this tells R to view the df in terms of groups of depth bins
  summarize(across(.fns = sum), .groups = "drop") # summarize all other values by summing all the rows in each group. (.groups = "drop" is to ungroup the data after we are done)

group_data_med$order <- 1:8

ind_based_rare_med = list()

for (i in unique(group_data_med$depth)) {
  
  one_depth = group_data_med %>% filter(depth == i)
  
  depth_sp_matrix = one_depth[,2:ncol(one_depth)]
  
  rarefaction = rarefaction.individual(depth_sp_matrix, method = "sample-size", q = 0)
  
  rarefaction$depth = i
  
  rarefaction$order = one_depth$order
  
  ind_based_rare_med[[i]] = rarefaction
  
}

ind_based_rare_med = bind_rows(ind_based_rare_med)

colnames(ind_based_rare_med) = c("Individuals","Richness","Depth_bin", "order")

# plot:
ggplot(ind_based_rare_med, aes(x=Individuals, y=Richness, color = reorder(Depth_bin, order)))+
  geom_line(size = 1.5)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# red:

group_data_red = red_depth_bins %>% 
  dplyr::select(depth, 7:ncol(red_depth_bins)) %>% 
  group_by(depth) %>% #this tells R to view the df in terms of groups of depth bins
  summarize(across(.fns = sum), .groups = "drop") # summarize all other values by summing all the rows in each group. (.groups = "drop" is to ungroup the data after we are done)

group_data_red$order <- 1:8

ind_based_rare_red = list()

for (i in unique(group_data_red$depth)) {
  
  one_depth = group_data_red %>% filter(depth == i)
  
  depth_sp_matrix = one_depth[,2:ncol(one_depth)]
  
  rarefaction = rarefaction.individual(depth_sp_matrix, method = "sample-size", q = 0)
  
  rarefaction$depth = i
  
  rarefaction$order = one_depth$order
  
  ind_based_rare_red[[i]] = rarefaction
  
}

ind_based_rare_red = bind_rows(ind_based_rare_red)

colnames(ind_based_rare_red) = c("Individuals","Richness","Depth_bin", "order")

# plot:
ggplot(ind_based_rare_red, aes(x=Individuals, y=Richness, color = reorder(Depth_bin, order)))+
  geom_line(size = 1.5)

# _______________________________________________________________

                   # Sample based rarefaction

# When we use Sample-based rarefaction we summing up the number of new species that added with each 
  # new sample and ignoring the number of individuals belonging to each species. 
    #Therefore we *neutralize the effect of large schools* (herds/flocks etc.) on the rarefaction curve shape.
# package: rareNMtests

# med:

sample_based_rare_med = list()

for (i in unique(med_depth_bins$depth)) {
  
  one_depth = med_depth_bins %>% filter(depth == i)
  
  depth_sp_matrix = one_depth[,first_species:ncol(one_depth)]
  
  rarefaction = rarefaction.sample(depth_sp_matrix, method = "sample-size", q = 0)
  
  rarefaction$depth = i
  
  sample_based_rare_med[[i]] =  rarefaction
  
}

sample_based_rare_med = bind_rows(sample_based_rare_med)

colnames(sample_based_rare_med) = c("Samples","Richness","Depth_bin")

# reordering a factor :
sample_based_rare_med$Depth_bin = factor(sample_based_rare_med$Depth_bin, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-147'))

ggplot(sample_based_rare_med, aes(x=samples, y=richness, color = Depth_bin))+geom_line(size = 1.2)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# red:

sample_based_rare_red = list()

for (i in unique(red_depth_bins$depth)) {
  
  one_depth = red_depth_bins %>% filter(depth == i)
  
  depth_sp_matrix = one_depth[,first_species:ncol(one_depth)]
  
  rarefaction = rarefaction.sample(depth_sp_matrix, method = "sample-size", q = 0)
  
  rarefaction$depth = i
  
  sample_based_rare_red[[i]] =  rarefaction
  
}

sample_based_rare_red = bind_rows(sample_based_rare_red)

colnames(sample_based_rare_red) = c("Samples","Richness","Depth_bin")

# reordering a factor :
sample_based_rare_red$Depth_bin = factor(sample_based_rare_red$Depth_bin, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

ggplot(sample_based_rare_red, aes(x=samples, y=richness, color = Depth_bin))+geom_line(size = 1.2)

# _______________________________________________________________

# model for the relationship between richness and depth, accounting for abundance


# _______________________________________________________________


# clean environment:
rm(list=setdiff(ls(), c("red_depth_bins","first_species")))

# _______________________________________________________________

################################################################


# - continue Tal's cod

# -  run the spatial sample based rarefaction (understand the 2 mobr approaches)

# - rum the coverage based rarefaction - understand it




                             ### Richness ###
     # Does the fish community richness changes across seas and depths?
                               # sea: med
      # depths: (for now it's values but maybe I shoud trasform to bins of 10 m)

library(vegan)        # running rarefaction analysis
library(tidyverse)    # organizing the data 
library(dplyr)        # organizing the data
library(plotrix)      #
library(rareNMtests)  # running rarefaction analysis
library(mobr)         # running rarefaction analysis - for gradients
library(mapview)      # visualization
library(sf)           # visualization
library(raster)       # for raster object

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed"

# read data:
setwd(wd_processed_data)
wide_opCode <- read.csv("wide_opCode.csv")
wide_opCode = wide_opCode %>% dplyr::select(-X) # rm weird column

# _______________________________________________________________

                # preliminary exploration #
# _______________________________________________________________

# keep med sea observation only
wide_opCode_med <- filter(wide_opCode, sea == "med")

# keep red sea observations only
wide_opCode_red <- filter(wide_opCode, sea == "red")

# define the first column of species
first_species <- 7

# create species metric
species_metric_med <- wide_opCode_med[,first_species:length(wide_opCode_med)]
species_metric_red <- wide_opCode_red[,first_species:length(wide_opCode_red)]

# _______________________________________________________________

                    # depth-ranges

# med
range(wide_opCode_med$depth, na.rm=TRUE) # 5.6 - 147 m
# depths histogram for med sea
wide_opCode_med %>% 
  ggplot()+
  aes(x = depth)+     
  geom_histogram() + ggtitle("med sea - samples count at each depth") 
 # some depths have up to 6 observations\opcode

# -----------

# red:
range(wide_opCode_red$depth, na.rm=TRUE) # 8.1 - 149 m
# depths histogram for red sea
wide_opCode_red %>% 
  ggplot()+
  aes(x = depth)+
  geom_histogram() + ggtitle("red sea - samples count at each depth")  # most of the depths have 3-5 observations\opcodes

# _______________________________________________________________

                      # map view
# make sf objects:
wide_opCode_med_sf <- st_as_sf(wide_opCode_med, coords = c('lon_x', 'lat_y'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")) # WGS84 is in decimal degrees
wide_opCode_red_sf <- st_as_sf(wide_opCode_red, coords = c('lon_x', 'lat_y'), crs = crs("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

# visualize:
mapView(wide_opCode_med_sf, zcol = "depth")
mapView(wide_opCode_red_sf, zcol = "depth")

# _______________________________________________________________

          # explore sampling effort at each sea

# add 2 additional columns: abundance and richness per sample:
wide_opCode_med_sum<- wide_opCode_med %>% 
  mutate(abundance = rowSums(wide_opCode_med[,first_species:length(wide_opCode_med)])) %>%
  mutate(richness = rowSums(wide_opCode_med[, first_species:ncol(wide_opCode_med)] > 0)) %>%
  dplyr::select(OpCode, depth, lon_x, lat_y, richness, abundance)


# check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(wide_opCode_med_sum$abundance) # 0
max(wide_opCode_med_sum$abundance) # 338
# 338 folds variation in the alpha scale (sample scale)

# -----------

# add 2 additional columns: abundance and richness per sample:
wide_opCode_red_sum <- wide_opCode_red %>% 
  mutate(abundance = rowSums(wide_opCode_red[,first_species:length(wide_opCode_red)])) %>%
  mutate(richness = rowSums(wide_opCode_red[, first_species:ncol(wide_opCode_red)] > 0)) %>%
  dplyr::select(OpCode, depth, lon_x, lat_y, richness, abundance)

# check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(wide_opCode_red_sum$abundance) # 1
max(wide_opCode_red_sum$abundance) # 252
# 252 folds variation in the alpha scale (sample scale)

# _______________________________________________________________

# I noticed I have samples with very low abundances (0, 1, 2...) in the two seas
# since we are interested in exploring the effect of depth, I wanted to merge samples from the same depth strip 
# but what should be the width of the bin? 10m? 20m?
# _______________________________________________________________

           # sum abundance and richness for the wide data_sets:

# med:
med_matrix_r <- wide_opCode_med %>% 
  mutate(abundance = rowSums(wide_opCode_med[first_species:length(wide_opCode_med)])) %>%
  mutate(richness = rowSums(wide_opCode_med[, first_species:ncol(wide_opCode_med)] > 0)) 

# red:
red_matrix_r <- wide_opCode_red %>% 
  mutate(abundance = rowSums(wide_opCode_red[first_species:length(wide_opCode_red)])) %>%
  mutate(richness = rowSums(wide_opCode_red[, first_species:ncol(wide_opCode_red)] > 0))

# _______________________________________________________________

                 # plot individuals ~ depth (wide_data)

# med - points: 
ggplot(med_matrix_r,aes(x=depth,y=abundance))+
  geom_point(size = 1.5, col = "cyan") + ggtitle("med sea abundance ~ depth") 

# med - line: 
ggplot(med_matrix_r,aes(x=depth,y=abundance))+
  geom_line(size = 1.5, col = "cyan") + ggtitle("med sea abundance ~ depth") 

# red - points: 
ggplot(red_matrix_r,aes(x=depth,y=abundance))+
  geom_point(size = 1.5, col = "indianred2") + ggtitle("red sea abundance ~ depth") 

# red - line: 
ggplot(red_matrix_r,aes(x=depth,y=abundance))+
  geom_line(size = 1.5, col = "indianred2") + ggtitle("red sea abundance ~ depth") 

# -----------

                  # plot richness ~ depth (wide_data)

# med - points: 
ggplot(med_matrix_r,aes(x=depth,y=richness))+
  geom_point(size = 1.5, col = "cyan") + ggtitle("med sea richness ~ depth") 

# med - line: 
ggplot(med_matrix_r,aes(x=depth,y=richness))+
  geom_line(size = 1.5, col = "cyan") + ggtitle("med sea richness ~ depth") 

# red - points: 
ggplot(red_matrix_r,aes(x=depth,y=richness))+
  geom_point(size = 1.5, col = "indianred2") + ggtitle("red sea richness ~ depth") 

# red - line: 
ggplot(red_matrix_r,aes(x=depth,y=richness))+
  geom_line(size = 1.5, col = "indianred2") + ggtitle("red sea richness ~ depth")

# _______________________________________________________________

              # Convert the data into 20 m depth bins
# _______________________________________________________________

# I chose 20m because I saw other papers group in bins of 20 m, but 10 is also an option

med_depth_bins <- med_matrix_r
red_depth_bins <- red_matrix_r

# 5.6 - 147 m to 8 layers
med_depth_bins$depth <- cut(med_depth_bins$depth,                                          
                       breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 147),
                       labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-147'))

# 8.1 - 149 m to 8 layers
red_depth_bins$depth <- cut(red_depth_bins$depth,                                         
                       breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 149),
                       labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

# _______________________________________________________________

                       # summary table 

# med:
med_table <- med_depth_bins %>% 
  group_by(depth) %>% 
    summarize(abundance = sum(abundance), richness = sum(richness), samples = n())
  
# red:
red_table <- red_depth_bins %>% 
  group_by(depth) %>% 
    summarize(abundance = sum(abundance), richness = sum(richness), samples = n())

# for the binned depths groups: check the difference between the smallest n sample and the largest n sample in each sea (the variance in abundance):
min(med_table$abundance) # 40 - deepest band
max(med_table$abundance) # 1064 - shallowest band

min(red_table$abundance) # 30 - deepest band
max(red_table$abundance) # 1169 - second to shallowest band

# _______________________________________________________________


             # plot individuals ~ depth (10 m bins)

# med - points: 
ggplot(med_table,aes(x=depth,y=n))+
  geom_point(size = 1.5, col = "blue") + ggtitle("med sea abundance ~ depth (binned)") 

# red - points: 
ggplot(red_table,aes(x=depth,y=n))+
  geom_point(size = 1.5, col = "red") + ggtitle("red sea abundance ~ depth (binned)") 

# -----------

                # plot richness ~ depth (10 m bins)

# med - points: 
ggplot(med_table,aes(x=depth,y=s))+
  geom_point(size = 1.5, col = "blue") + ggtitle("med sea richness ~ depth (binned)") 

# red - points: 
ggplot(red_table, aes(x=depth,y=s))+
  geom_point(size = 1.5, col = "red") + ggtitle("red sea richness ~ depth (binned)") 

# _______________________________________________________________



                     # MOBR package 



# _______________________________________________________________
 
               # check for linear relationship 
           # between the # samples and # individuals

# ( The MoB methods assume a linear relationship between the number of plots and the number of individuals. 
#  This function provides a means of verifying the validity of this assumption )

# med:
plot_N(species_metric_med, n_perm = 1000) # - V - indeed linear relationship
# red: 
plot_N(species_metric_red, n_perm = 1000) # - V - indeed linear relationship

# _______________________________________________________________

                    # "mob_in" object 

# `comm` = community matrix (only species)
# `plot_attr` =  the meta-data for each sample (only the groups we want to compare. A character string) and coordinate if we have and want to use them.
# `coord_names`  = column names of longitude and latitude 

# -----------

            # for discrete \ binned depth groups:

# create env dataset
env_data_med_b <- med_depth_bins[, c(3, 5, 6)]
env_data_red_b <- red_depth_bins[, c(3, 5, 6)]

# create the mob_in objects:
mob_in_med_b <- make_mob_in(species_metric_med, plot_attr = env_data_med_b, coord_names = c("lon_x", "lat_y"))
mob_in_red_b <- make_mob_in(species_metric_red, plot_attr = env_data_red_b, coord_names = c("lon_x", "lat_y"))

# -----------

                 # for continuous depths

# create env dataset
env_data_med <- med_matrix_r[, c(3, 5, 6)]
env_data_red <- med_matrix_r[, c(3, 5, 6)]

# create the mob_in objects:
mob_in_med <- make_mob_in(species_metric_med, plot_attr = env_data_med, coord_names = c("lon_x", "lat_y"))
mob_in_red <- make_mob_in(species_metric_red, plot_attr = env_data_red, coord_names = c("lon_x", "lat_y"))

# _______________________________________________________________
 
                  # ' get_delta_stats '

# Conduct the MoB tests on drivers of biodiversity across scales: shape of the SAD, treatment\group density - N and aggregations degree:

# -----------
 
         # for discrete \ binned depth groups:

# med:
delta_med_b = get_delta_stats(mob_in = mob_in_med_b, 
                             env_var = 'depth',
                             type='discrete', log_scale=FALSE, n_perm = 20)

plot(delta_med_b, 'b1', )

# red:
delta_red_b = get_delta_stats(mob_in = mob_in_red_b, 
                            env_var = 'depth',
                            type='discrete', log_scale=FALSE, n_perm = 20)

plot(delta_red_b, 'b1')

# -----------

              # for continuous depths
# med:
delta_med = get_delta_stats(mob_in = mob_in_med, 
                             env_var = 'depth',
                             type='continuous', log_scale=FALSE, n_perm = 20)

plot(delta_med_b, 'b1', )

# red:
delta_red = get_delta_stats(mob_in = mob_in_red, 
                            env_var = 'depth',
                            type='continuous', log_scale=FALSE, n_perm = 20)

plot(delta_red_b, 'b1')









# delta_med <- get_delta_stats(mob_in_med, mob_in_med$env$depth, ref_level = NULL, tests = c("SAD", "N", "agg"), spat_algo = 'kNN', type = "continuous",
#                  stats = NULL, inds = NULL, log_scale = FALSE, min_plots = NULL, density_stat = c("mean", "max", "min"), n_perm = 1000, overall_p = FALSE )
# 
# delta_red <- get_delta_stats( mob_in, env_var, group_var = NULL, ref_level = NULL, tests = c("SAD", "N", "agg"), spat_algo = NULL, type = c("continuous", "discrete"),
#                               stats = NULL, inds = NULL, log_scale = FALSE, min_plots = NULL, density_stat = c("mean", "max", "min"), n_perm = 1000, overall_p = FALSE )



     # plot rarefaction
# plot_rarefaction(mob_in_med, group_var = 'depth', ref_level = NULL, method = c('IBR', 'SBR', 'nsSBR', 'sSBR'), pooled = false, spat_algo = 'kNN', col = NULL, lwd = 3, log = "X", leg_loc = "topleft") 
# 
# plot_rarefaction(mob_in_med, group_var = 'depth', method = 'IBR', dens_ratio = 1, pooled = TRUE, col = NULL, lwd = 3, log = NULL, leg_loc = "topleft")
# 
# 
# plot_abu(mob_in_med, group_var = mob_in_med$env[1], ref_level = NULL, type = c("sad", "rad"), pooled = FALSE, col = NULL, lwd = 3, log = "", leg_loc = "topleft" )
# 





# next steps:

# after we created richness\depth profiles - now lets create richness rarefactions: per sea, per depth bins (for each sea separated), per sample 

# then ask myself what syntax the most complex function in mob r (the one that calculated the difference between the rarefactions) needs
# correct the rarefaction syntax I used to what is required
# created the other rarefactions

# at the same time - set a date and a staff to my committee

# and later later I think of adding analysis that are connected to species compositions - although this was done before










             # measuring diversity - rarefactions
# 
# 
# 
# # _______________________________________________________________
# 
#              # rarefied max richness - bar plot #
# 
# # a variable that indicate from which column the species list start:
# first_species = 7 #the column form which the species list begin in wide_day2 dataset. 4,if we use the regular wide_day
# 
# # -----------
# 
# # subset our data to create the *species matrix data frame*, which we will call sp_matrix:
# sp_matrix = wide_opCode[,first_species:length(wide_opCode)]
# 
# # -----------
# 
# # check the minimum number of records in a sample (opcode)
# raremax = sp_matrix %>% rowSums(.) %>% min() # 0 
# 
# # -----------
# 
#                 # samples abundance histogram
# 
# # when we rarefy by sample, we may see an extremely low individual count. Rarefying to a low number such as *0* isn’t really helpful.
# # lets observe how abundance varies in our data set. The x- axis is plotted on a log10 scale for better clarity:
# sp_matrix %>% 
#   mutate(abundance = rowSums(.)) %>% 
#   ggplot()+
#   aes(x = abundance)+
#   geom_histogram()+
#   scale_x_log10() # for clarity
# 
# # -----------
# 
#                # filter low abundance samples
# 
# # some samples have extremely low abundances (o). 
# # We will remove these samples and stay only with samples that have more than *9* individuals
# wide_opCode_clean = wide_opCode %>% 
#   mutate(abundance = rowSums(wide_opCode[first_species:length(wide_opCode)])) %>% 
#   filter(abundance > 4) %>%  # set a threshold (based on the histogram)
#   mutate(abundance = NULL) # remove this column so we will have a fresh start
# 
# # -----------
# 
#           # new sp_matrix based on the filtered df:
# sp_matrix_filter = wide_opCode_clean[,first_species:length(wide_opCode_clean)]
# 
# # check lowest abundance count now:
# raremax_n = sp_matrix_filter %>% rowSums(.) %>% min() # = 5
# 
# # -----------
# 
#             # "rarefy" function (vegan package)
# 
# # - *x* = the data frame we will use 
# # - *sample* = how many individuals to use in the rarefraction.
# # - We will take the *minimum* number of records in a single site (7)
# 
# # species richness in each sample when we sample randomly 5 individuals:
# rare_index_raremax = rarefy(x=sp_matrix_filter, sample=raremax_n)
# 
# # reorder the values in the sea column:
# wide_opCode_clean$sea <- factor(as.character(wide_opCode_clean$sea), levels=c("red", "med"))
# 
# # -----------
#                       # Bar plot #
# 
# # mutate the data by adding a rarefied_richness column, and plot a richness barplot
# rare_richness_plot <- wide_opCode_clean %>%
#   mutate(rarefied_richness = rare_index_raremax) %>%
#   ggplot() + aes(x = sea, y = rarefied_richness, fill = sea) +
#   stat_summary(geom = "bar",fun.data = mean_se)+
#   stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge",width = 0.3)
# 
# print(rare_richness_plot) # the red sea has higher richness per sample if we rarefy to 7 individuals
#  
# # _______________________________________________________________
# 
#      ######################## richness! ##############################
# 
# # exploring how richness changes under different treatments and env gradients scenarios
# 
# # _______________________________________________________________
# 
#            # individual based rarefaction curve - one sample #
#  
# # But maybe the cut off\sampling effort I chose is wrong? (samples > 6 individuals?)
# # Better to plot the individual rarefaction curve using `rarefaction.individual` function:
# 
# # Accumulating richness with accumulating individuals at a single sample. q=0 is richness:
# # The sample size for a given sample is the N - number of individuals:
# one_sample_rare = rarefaction.individual(sp_matrix[1,], method = "sample-size", q = 0)
# 
# head(one_sample_rare) # richness Hill numbers for each number of individuals
# 
# # plot:
# ggplot(data = one_sample_rare,aes(x = `sample-size`,y = `Hill (q=0)`)) +
#   geom_line() +
#   xlab("Individuals")+
#   ylab("Richness")
# 
# # _______________________________________________________________
#  
#                                   # 1 #
# 
#       # individual based rarefaction curve - all samples + both seas #
# 
# # A loop to apply the rarefaction curve to each row in the data
# # p.s. the "rarefaction.individual" function doesn't work for abundance = 0
# # inside the loop we remove the 3 samples from the med sea with 0 abundance
# # we don't have to do it when we group by a treatment. it's even better to leave them
# 
# ind_based_rare_all = list() # empty list
# 
# for (i in unique(wide_opCode$OpCode)) { # for each row...
#   
#   one_opCode = wide_opCode %>% filter(OpCode == i) # filter
#   
#   one_opCode_sp_matrix = one_opCode[,first_species:ncol(one_opCode)] # create sp_matrix
#   
#   if(rowSums(one_opCode_sp_matrix > 0)){
#   
#   rarefaction = rarefaction.individual(one_opCode_sp_matrix, method = "sample-size", q = 0) # apply rarefaction function
#   
#   rarefaction$sea = rep(unique(one_opCode$sea)) # add which sea
#   
#   rarefaction$OpCode = i # add the opcode_id to the data
# 
#   ind_based_rare_all[[i]] = rarefaction } # save in my list
#   
# }
# 
# ind_based_rare_all = bind_rows(ind_based_rare_all) # convert from list to data frame
# colnames(ind_based_rare_all) <- c("Individuals","Richness","sea", "sample") # change columns names
# 
# # -----------
# 
# # plot
# ggplot(ind_based_rare_all,aes(x=Individuals,y=Richness,group = sample ,color = sea))+
#   geom_line()
# 
# # _______________________________________________________________
# 
#                                 # 2 #
# 
#           # individual based rarefaction curve - sea scale # 
# 
# # # plot                     # shows all samples and both seas in an organise way
# # ggplot(ind_based_rare_all,aes(x=Individuals,y=Richness,group = sea ,color = sea))+
# #   geom_line()
# 
# # create the data set:
# group_data = wide_opCode %>% 
#   select(sea, 7:ncol(wide_opCode)) %>% 
#   group_by(sea) %>% #this tells R to view wide_day2 in terms of groups of habitat
#   summarize(across(.fns = sum), .groups = "drop") #summarize all other values by summing all the rows in each group. 
#                                                   # (.groups = "drop" is to ungroup the data after we are done)
# 
# # -----------
# 
# # The result is a much shorter dataset, where each row is the sum abundance of each species in each sea.  
# group_data
# 
# # -----------
# 
# # changing the order of the seas for visualization purposes:
# group_data$sea <- factor(as.character(group_data$sea), levels=c("red", "med"))
# 
# # -----------
# 
# # loop to create IBR for the 2 seas:
# ind_based_rare = list()
# 
# for (i in unique(group_data$sea)) {
#   
#   one_sea = group_data %>% filter(sea == i)
#   
#   sea_sp_matrix = one_sea[,first_species:ncol(one_sea)]
#   
#   rarefaction = rarefaction.individual(sea_sp_matrix, method = "sample-size", q = 0)
#   
#   rarefaction$sea = i
#   
#   ind_based_rare[[i]] = rarefaction
#   
# }
# 
# ind_based_rare = bind_rows(ind_based_rare)
# colnames(ind_based_rare) = c("Individuals","Richness","Survey")
# 
# # -----------
# 
# #plot: 
# ggplot(ind_based_rare,aes(x=Individuals,y=Richness,color = Survey))+
#   geom_line(size = 1.5) +
#   scale_color_manual(values=c("cyan", "indianred2"))
# 
# # _______________________________________________________________
# 
#                           #  mobr package  #
# 
#     # Creating richness rarefaction curves over treatments and env gradients #  
# 
# # -----------
#  
# # calc_biodiv - Calculating biodiversity statistics from sites by species table for each sea separately: red and med
# bio_stat <- calc_biodiv(abund_mat = wide_opCode[,first_species:ncol(wide_opCode)], 
#                         groups = wide_opCode$sea, 
#                         index = c("N", "S", "S_n", "S_asymp", "f_0", "pct_rare", "PIE", "S_PIE"), 
#                         effort = 10, ########################################### change it to something like the min abundance in a sea\ somethinf else that makes sense and is not arbitrary. 
#                         extrapolate = TRUE, 
#                         return_NA = FALSE, 
#                         rare_thres = "N/S")
# # I can use this table to extract any possible indices I am interested in, for each sea. The indices are calculated per sapmle.
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   
  
  
  




















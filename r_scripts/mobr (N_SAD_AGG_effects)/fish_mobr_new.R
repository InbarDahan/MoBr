

                   # Fish MOBR - generation 2

# Variables:
# OpCode - the sample id
# lat_y - the sample lat
# lon_x - the sample lon
# depth - num (continuous - 29.5 )
# temperature
# -- there is no repetition in the exact depth and so I need to group depth into depth layers in order to 
# -- have repetitions. 

                       # Mobr Analysis:

library(mobr)# running rarefaction analysis - for gradients
library(dplyr)        
library(ggplot2)
# -----------------------------------------------------------------------------

# read
# wd
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_plots <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/output/plots"
wd_results <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/output/files"

# read data
setwd(wd_processed_data)
wide_red <- read.csv("wide_red.csv")

# removing the wired column x
wide_red = wide_red %>% dplyr::select(-X)
# -----------------------------------------------------------------------------

# prepare the data

# define the first column of species
first_species <- 7

# create species metric
sp_matric <- wide_red[,first_species:length(wide_red)]

# binning the data:
# 8.1 - 149 m to - 15 m bins (10 breaks)
red_depth_bins <-  wide_red %>% mutate(depth_group = cut(wide_red$depth,                                         
                   breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
                   labels=c('7.5', '22.5', '37.5', '52.5', '67.5','82.5', '97.5', '112.5', '127.5', '142.5'))) %>% 
                   relocate(depth_group,.after = depth) 
  
# prepare the plot attributes data frame
plot_attr <- red_depth_bins %>% select(, c(1:7))

# Create the mob_in object
mob_in <- make_mob_in(sp_matric, plot_attr, coord_names = c('lon_x', 'lat_y'))
# -----------------------------------------------------------------------------

                   # Run the MOBr analysis

# Perform the MoB analysis using the depth groups
# 1 - env_var = depth_group (needs to be continuous)
delta_stats_log <- get_delta_stats(mob_in, 'depth_group', type = 'continuous', log_scale = TRUE, n_perm = 199)
delta_stats_not_log <- get_delta_stats(mob_in, 'depth_group', type = 'continuous', n_perm = 199)
# ------------

# ###
# # wrong alternatives - understand why with Daniel:
# # 2 - env_var = depth
# delta_red <- get_delta_stats(mob_in, env_var = 'depth',group_var = 'depth_group',type = 'continuous', n_perm = 199)
# ### env_var = depth_group or env_var = depth completely reverses the results.
# 
# # 3 - env_var = depth
# mob_in_rev <- mob_in
# # Reversing Depth Order to check if it fixes the problem:
# mob_in_rev$env$depth <- -mob_in_rev$env$depth  # Reverse depth sign1
# 
# delta_red_rev_not_log <- get_delta_stats(mob_in_rev, env_var = 'depth',group_var = 'depth_group', type = 'continuous', n_perm = 199)
# delta_red_rev_log <- get_delta_stats(mob_in_rev, env_var = 'depth',  group_var = 'depth_group', type = 'continuous', n_perm = 199, log_scale = TRUE)
# ###
# ------------

# results
dev.off() 
plot(delta_stats_log, 'b1')          # good option - continuous depth_group log
plot(delta_stats_not_log, 'b1')      #good option - continuous depth_group not log
# plot(delta_red, 'b1')                # wrong - reversed - probably wrong
# plot(delta_red_rev_log,'b1')         # wrong - also reversed - although looks otherwise
# plot(delta_red_rev_not_log, 'b1')    # wrong - also reversed - although looks otherwise
# ------------

# in summary:
# - rev do not fixes anything - its just the result of the delta_red in reverse
# - log scale doesn't change much
# - env_var = depth group make the result more distinct and also changes the order in the depth rarefaction curves

# alternative # 2 - env_var = depth
# delta_red <- get_delta_stats(mob_in, env_var = 'depth',  # define the env gradient
#                              group_var = 'depth_group', stats = c('betas', 'r'), # define the site\grouping var
#                              type = 'continuous', spat_algo = 'kNCN', # define the type of var # define the spatial rarefaction method
#                              n_perm = 199, overall_p = FALSE) #, density_stat = c("mean", "max", "min")) # 

# -----------------------------------------------------------------------------

   # N vs depth for 3 sampling efforts: 10, 20, 30 ind
                  # delta_red_rev 

delta_red_rev_3_eff <- get_delta_stats(mob_in, 'depth_group', type = 'continuous', n_perm = 199, inds = c(10,20,30))

dev.off() 
plot(delta_red_rev_3_eff, 'b1') 

plot(delta_red_rev_3_eff, stat = 'b1',  # stat can be changed based on the wanted effect size method 
     eff_sub_effort = T,          # subset\all samples   #   
     eff_disp_pts = F,                              # T\F - show the raw effect points
     eff_disp_smooth = T)

# -----------------------------------------------------------------------------

### Here I am trying to flip the color of the gradient but I am also wondering if
### it will change the results of the analysis - The right results should be that
### evenness reduce with depth, N reduces with depth and agg increases with depth

### Doesn't work because I can not - a factor object

### env_var = depth_group or env_var = depth completely reverses the results.
mob_in_rev <- mob_in
# Reversing Depth Order to check if it fixes the problem:
mob_in_rev$env$depth_group <- -mob_in_rev$env$depth_group  # Reverse depth sign
#### cant do it because it a factor, not an numeric.
# check the levels of the factors

delta_stats_rev <- get_delta_stats(mob_in_rev, 'depth_group', 
                          type = 'continuous', log_scale = TRUE, n_perm = 199)

dev.off() 
plot(delta_stats, 'b1') 
plot(delta_stats_rev, 'b1') 

# -----------------------------------------------------------------------------
 
       # creating 2 analysis - assuming no linear trend

   # The delta effect is assumed to be linear but is not necessarily so, 
# hence we need to run two analysis, before and after the change in slopes

red_depth_bins$depth_group <- as.numeric(as.character(red_depth_bins$depth_group))
# ------------

   # try once before and after 37.5 m (create 2 data sets), and once for 52.5 m

# - Define breakpoints
breakpoints <- c(37.5, 52.5)

# - create a general function to split the data 
split_data <- function(df, cutoff) {
  list(
    under = df %>% filter(as.numeric(as.character(depth_group)) <= cutoff),
    above = df %>% filter(as.numeric(as.character(depth_group)) > cutoff)
  )
}

# - Apply function to all breakpoints (x takes each break-point in the list)
data_subsets <- lapply(breakpoints, function(x) split_data(red_depth_bins, x))

# - Flatten list (reduce from 2 lists to a united one) and name the elements
data_subsets <- unlist(data_subsets, recursive = FALSE)
names(data_subsets) <- c("under_37", "above_37", "under_52", "above_52")
# ------------

       # define env_attr + mob_in for the subsets

delta_stats <- list()

for(i in 1:length(data_subsets))
{
  # pull one data set:
  sub_gr <- data_subsets[[i]]
  
  # prepare the plot attributes data frame
  env_attr <- sub_gr %>% select(, c(1:7))
  
  # prepare species matrix
  first_species <- 8
  sp_matric <- sub_gr[,first_species:length(sub_gr)]
  
  # Create the mob_in object
  mob_in <- make_mob_in(sp_matric, env_attr, coord_names = c('lon_x', 'lat_y'))
  
  # calculate delta
  one_delta_stats <- get_delta_stats(mob_in, 'depth_group', type = 'continuous', n_perm = 199)
  
  # add the stats results to the list of stats:
  delta_stats[[i]] <- one_delta_stats
}

# add the name of the groups to the results:
names(delta_stats) <- names(data_subsets)

plot(delta_stats$under_37, 'b1') 
plot(delta_stats$above_37, 'b1') 
plot(delta_stats$under_52, 'b1') 
plot(delta_stats$above_52, 'b1') 

# seems like the cut of at 37 splits the data better
# into 2 linear trends:

# - under_37: between 8- 37 m
### agg - as we deepen the agg effect on richness become more negative 
# meaning that there are more agg at 30 m then 8 m. 
### N - as we deepen the N is reduced and hence we get a reduction in richness.
# meaning that there are more individuals at 8 m versus 30 m.
### SAD - as we deepen the SAD goes up slightly
# meaning that the community is more even at 30 m then shallowr.

# - above_37: between 38 - 149
### agg - as we deepen the agg effect on richness is reduced
# meaning that there are less aggregations in deep waters (100) than 40 m. 
### N - as we deepen the N is reduced and hence we get a reduction in richness.
# meaning that there are more individuals at 40 m versus 125 m.
### SAD - as we deepen the SAD is reduced
# meaning that the community is more even at 30 m then deeper.


























                          
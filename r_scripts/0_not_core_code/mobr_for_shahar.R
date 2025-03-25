

                ###  MOBr - fish (here marked as f) ###

  # what drivers shape the change in fish richness over the depth gradient ?

# -----------------------------------------------------------------------------

# codes from the papers:

# discrete explanatory variables:
# "Measurement of Biodiversity (MoB): A method to separate the scale-dependent effects of species abundance distribution, density, and aggregation on diversity change"

# continuous explanatory variables:
# "A multiscale framework for disentangling the roles of evenness, density, and aggregation on diversity gradients"

#' @author Dan McGlinn and Xiao Xiao

# -----------------------------------------------------------------------------
install.packages(c('devtools', 'mobr', 'leaflet', 'mapview', 'tidyr',
                   'vegan', 'dplyr', 'ggplot2', 'egg'))

library(mobr)         # running rarefaction analysis - for gradients
library(vegan)
library(dplyr)        
library(ggplot2)
library(egg)
devtools::install_version("broom", version = "0.5.4", repos = "http://cran.us.r-project.org")
library(broom)
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

# calculate samples abundances:


# -----------------------------------------------------------------------------

# prepare the data

# define the first column of species
first_species <- 7

# create species metric
sp_matric <- wide_red[,first_species:length(wide_red)]

#if binned data required:
red_depth_bins <- wide_red

# 8.1 - 149 m to - 15 m bins (8 layers)
red_depth_bins <-  wide_red %>% mutate(depth_group = cut(wide_red$depth,                                         
                                                         breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
                                                         labels=c('0-15', '16-30', '31-45', '46-60', '61-75','76-90', '91-105', '106-120', '121-135', '136-149')))

# -----------------------------------------------------------------------------

# make the "mob in" object
# create env var
env_red <- red_depth_bins[, c(1, 3:6, 180)]

# make mob in
red_bin_mob_in <- make_mob_in(sp_matric, env_red, # define the sp matrix and the variables   # no wornings apear - v
                              coord_names = c('lon_x', 'lat_y'))   # define coor

# -----------------------------------------------------------------------------

# compute mob stats: 

# binned:
# here the α and γ scales are different. One sample vs depth bin
stats_f <- get_mob_stats(red_bin_mob_in, group_var = "depth_group", extrapolate = FALSE) 

plot(stats_f)

#______________________________________________________________________________

#                  Multi-scale analysis -  gradient
#                       N, SAD and aggregations
#______________________________________________________________________________

# ' get_delta_stats '

# Conduct the MoB tests on drivers of biodiversity across scales: 
# shape of the SAD, treatment\group density (N) and aggregations degree:

# red:
delta_red = get_delta_stats(mob_in = red_bin_mob_in, env_var = 'depth',  # define the env gradient
                            group_var = 'depth_group', stats = c('betas', 'r'), # define the site\grouping var
                            type = 'continuous', spat_algo = 'kNCN', # define the type of var # define the spatial rarefaction method
                            n_perm = 199, overall_p = FALSE, density_stat = c("mean", "max", "min")) # 

plot(delta_red, stat = 'b1', scale_by = 'indiv',  # stat can be changed based on the wanted effect size method 
     eff_sub_effort = F, eff_log_base = 1,          # subset\all samples   #   
     eff_disp_pts = T,                              # T\P - show the raw effect points
     eff_disp_smooth = F) # T\F - show the linear effect of the expl var on the effect sizes

# in the upper plotting scheme it delets 492 rowa

plot(delta_red)

# in this plot it delets only 42 rows




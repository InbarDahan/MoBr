

                 "MOB MULTIMETRIC ANALAYSIS"
       "2 scales - calculating each effect separetly "

library(mobr)
library(dplyr)                 
                 
# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________
                 
# read data:
setwd(wd_processed_data)
fish <- read.csv("wide_red.csv")
# removing the wired column x
fish = fish %>% dplyr::select(-X)
# _______________________________________________________________

# define the first column of species
first_species <- 7

# create species metric
sp_matric <- fish[,first_species:length(fish)]


# if binned data required:
# 8.1 - 149 m to - 15 m bins (10 breaks)
red_depth_bins <-  fish %>% mutate(depth_group = cut(fish$depth,                                         
                   breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
                   labels=c('7.5', '22.5', '37.5', '52.5', '67.5','82.5', '97.5', '112.5', '127.5', '142.5'))) %>% 
                   relocate(depth_group,.after = depth) 

# prepare the plot attributes data frame
plot_attr <- red_depth_bins %>% dplyr::select(, c(1:7))

# _______________________________________________________________
                 
# make mob_in:
mob_in <- make_mob_in(sp_matric, plot_attr, coord_names = c('lon_x','lat_y'))
                 
# _______________________________________________________________                 
                 
# count samples abundances:
species_abundances_a <- rowSums(red_depth_bins[ ,8:length(red_depth_bins)])                 

# add the abundance column to the data frame
red_depth_bins$sample_abu <- species_abundances_a
red_depth_bins <- red_depth_bins %>% relocate(sample_abu,.after = lon_x)

# check the number of samples with more then 5 individuals (sampling effort)
test_a <- red_depth_bins %>% group_by(depth_group) %>%
summarise(high_samples = sum(sample_abu >= 5))

length(which(red_depth_bins$sample_abu >= 5)) # 110\123 # this might be the problam in the graph

# _______________________________________________________________
                 
              # 2 scale analysis - continuous
                    # alpha + delta

# calculating the abundances (N), richness (S), rarefied richness (S_n) 
# and the effective number of species for the probability of interspecific 
# encounter index (S_PIE), for two scales (alpha, gamma) and thier ratio, 
# which is the beta diversity, estimates of species turnover.

indices <- c('N', 'S', 'S_n', 'S_PIE')
inv_div <- tibble(sp_matric) %>%
  group_by(group = plot_attr$depth_group) %>%
  group_modify(~ calc_comm_div(.x, index = indices, effort = c(10,20,30),
                               extrapolate = TRUE,
                               scales = c("alpha", "gamma", "beta")))

# head:
head(inv_div)

# plot richness:
plot_comm_div(inv_div, 'S') # α - boxplot sp per sample| γ - total per group                
plot_comm_div(inv_div, 'N')  # α - boxplot indv per sample| γ - total per group                    
plot_comm_div(inv_div, 'S_n') # rarefied S - ??????????
plot_comm_div(inv_div, 'S_PIE')

plot_comm_div(inv_div, index = c('S', 'S_n', 'S_PIE'), multi_panel = TRUE)
# _______________________________________________________________

                    # beta diversity

betas <- calc_beta_div(sp_matric, c('S', 'S_n', 'S_PIE'), effort = c(10, 20,30))
betas
                 
# _______________________________________________________________
                                                                                                                                
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
                 
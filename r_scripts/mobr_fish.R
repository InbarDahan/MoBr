
    
                           ###  MOBr - fish   ###

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

# create env var
env_red_b <- red_depth_bins[, c(1, 3:6, 180)]

# -----------------------------------------------------------------------------

                     # make the "mob in" object

# binned
red_bin_mob_in <- make_mob_in(sp_matric, env_red_b, # define the sp matrix and the variables
                          coord_names = c('lon_x', 'lat_y'))   # define coor
# no warnings appear - we answer the assumptions required!

# -----------------------------------------------------------------------------

                         # compute mob stats: 

# binned:
 # here the α and γ scales are different. One sample vs depth bin
stats_b <- get_mob_stats(red_bin_mob_in, group_var = "depth_group", extrapolate = FALSE) 

plot(stats_b)

# Multiple graphs the binned depths var. It seems like there is a decrease in abundance and richness with depth.
# but not in rarefied richness. Scale dependent.

#______________________________________________________________________________

#                  Multiscale analysis - continuous 
#                       N, SAD and aggregations
#______________________________________________________________________________

                        # ' get_delta_stats '

# Conduct the MoB tests on drivers of biodiversity across scales: 
# shape of the SAD, treatment\group density (N) and aggregations degree:

# red:
delta_red = get_delta_stats(mob_in = red_bin_mob_in, env_var = 'depth',  # define the env gradient
                            group_var = 'depth_group', stats = c('betas', 'r'), # define the site\grouping var
                            type = 'continuous', spat_algo = 'kNCN', # define the type of var # define the spatial rarefaction method
                            n_perm = 100, overall_p = FALSE) # 

plot(delta_red, stat = 'b1', scale_by = 'indiv',
     eff_sub_effort = F, eff_log_base = 2,
     eff_disp_pts = F,
     eff_disp_smooth = T)

plot(delta_red, stat = 'b1') 

setwd(wd_results)
save(delta_red, file = 'deltas_red.Rdata')

# (1) richness over depth - there is no clear gradient. Some very deep depths have high richness, Yet shallow is the highest.
# (2.a) agg over depth - agg has strong negative effect in intermediate to deep depths (60-90 m), meaning there are more agg that reduce the richness
# (2.b)  N over depth  - N counter-reflects the agg effect - it has a positive effect (60-90 m) meaning that more individuals in the deep support higher richness
# (2.c) SAD over depth - has positive effect all throughout, with a bit higher effect in the same depths (60-90 m). The evenness of the community support richness in the same way across depths
# (3.a) agg over scale - as scale increase, the negative effect of agg increase
# (3.b)   N over scale - as scale increase, the positive effect of N increase
# (3.c) SAD over scale - as scale increase, SAD starts to gain negative effect, but just from 30 ind

# summary - there are changes in the behavior of different underlying components on richness over the depth gradient
#           the components are also influenced by scale.

#______________________________________________________________________________

                    # 2 scale analysis - continuous
                             # alpha + delta
#______________________________________________________________________________

# binned:
alphas_b <- stats_b$samples_stats
# alphas_b$OpCode <- red_depth_bins$OpCode[match(row.names(alphas_b), row.names(red_depth_bins))]
gammas_b <- stats_b$groups_stats  

# -----------------------------------------------------------------------------

                           # fit linear models 
# binned
lm_alpha_b = alphas_b %>% group_by(index, effort) %>%
  do(mod_b = lm(value ~ group, data = .))

lm_gamma_b = gammas_b %>% group_by(index, effort) %>%
  do(mod_b = lm(value ~ group, data = .))

# -----------------------------------------------------------------------------

                            # get model coefs
# binned
mod_coef_alpha_b = broom::tidy(lm_alpha_b, mod_b)
mod_coef_gamma_b = broom::tidy(lm_gamma_b, mod_b)

# -----------------------------------------------------------------------------

                           # get model summary
# binned
mod_sum_alpha_b = broom::glance(lm_alpha_b, mod)
mod_sum_gamma_b = broom::glance(lm_gamma_b, mod)

# -----------------------------------------------------------------------------

                      # binding the α and γ data stats
# binned
mob_met_b = rbind(data.frame(scale = 'alpha', alphas_b),
                  data.frame(scale = 'gamma', gammas_b))

# -----------------------------------------------------------------------------

                 # reordering the levels of the "index" column
# binned
mob_met_b$index = factor(mob_met_b$index, 
                         levels = levels(mob_met_b$index)[c(2:1, 3:7)])

# -----------------------------------------------------------------------------

           # plotting the effect of depth on total abundance
# binned
colnames(mob_met_b)[colnames(mob_met_b) == "group"] <- "depth"

  N1_b = mob_met_b %>%
    subset(index %in% 'N') %>%
    ggplot(aes(x = depth, y = value)) + 
    geom_point(aes(color = scale)) +
    geom_smooth(aes(color = scale), method = 'lm', se = F) +
    labs(x = "binned depth (m)", y = "Total abundance (N)")

  ggsave("binned_grad_vs_N.pdf", plot = N1_b, path = wd_plots, width = 15, height = 12, 
         units = "cm")

# -----------------------------------------------------------------------------
  
     # plotting the effect of depth on α scale abundance & richness

# binned  
  SN1_b = mob_met_b %>% 
    subset(index %in% c('S', 'N')) %>%
    subset(scale == 'alpha') %>%
    ggplot(aes(x = depth, y = value )) + 
    geom_point(col = "indianred2") +
    geom_smooth(method = 'lm', se = T) +
    facet_wrap(~ index, scales = "free",  
               labeller = as_labeller(c(S = "α Species richness (α S)",
                                        N = "α Total abundance (α N)"))) +
                                        labs(x = "binned_depth (m)")
  
 ggsave("binned_grad_vs_S&N.pdf", plot = SN1_b, path = wd_plots, width = 15, height = 12, 
        units = "cm")

# -----------------------------------------------------------------------------

      # plotting the effect of depth on 'S', 'S_n', 'S_PIE'

 # binned 
 p1_b = mob_met_b %>% 
   subset(abs(value) < 1000) %>%
   subset(index %in% c('S', 'S_n', 'S_PIE')) %>% 
   ggplot(aes(x = depth, y = value, col = scale)) + 
   geom_point() +
   geom_smooth(method = 'lm', se = F) +
   facet_wrap(. ~ index, scales = "free")

# -----------------------------------------------------------------------------

     # plotting the effect of depth on 'S', 'S_n', 'S_PIE'
 
 # binned 
 p2_b = mob_met_b %>% 
   subset(abs(value) < 1000) %>%
   subset(index %in% c('beta_S', 'beta_S_n', 'beta_S_PIE')) %>% 
   ggplot(aes(x = depth, y = value)) + 
   geom_point() +
   geom_smooth(method = 'lm', se = F) +
   facet_wrap(. ~ index, scales = "free")
 
 g_b = ggarrange(p1_b, p2_b)
 
ggsave("ENS.pdf", plot = g_b, path = wd_plots, width = 20, height = 15, units = "cm")
  
###############################################################################

                # ?????????????????????????? # 

 # the graphs work but something is very very weird with them. Try to find the 
 # problem
  


  
  
  
  
  
  
  
  
  

# calculate biodiversity statistics: 

biodiv <- calc_biodiv(sp_matric, groups = "depth", index = c('N', 'S', 'S_n', 'S_PIE'), effort = 5, extrapolate = FALSE, return_NA =  FALSE) 

# ----------------------------------------------------------------------------- 
  
  
  
  
  
  
  
##################################################

# my code - comes from the mobR package instructions:

 #################################################
# _______________________________________________________________

# check for linear relationship 
# between the # samples and # individuals

# (The MoB methods assume a linear relationship between the number of plots (samples\opcode) and the number of individuals. 
# This function provides a mean of verifying the validity of this assumption )
# red: 
plot_N(species_metric_red, n_perm = 1000) # - V - indeed linear relationship

# ***************************************************************

# red:

for (i in unique(red_depth_bins$depth)) {
  
  one_depth = red_depth_bins %>% filter(depth == i)
  
  depth_sp_matrix = one_depth[,first_species:ncol(one_depth)]
  
  mobr::plot_N(depth_sp_matrix, n_perm = 100)
  
}

# _______________________________________________________________

# mob_in object 

# `comm` = community matrix (only species)
# `plot_attr` =  the meta-data for each sample (only the groups we want to compare. A character string) and coordinate if we have and want to use them.
# `coord_names`  = column names of longitude and latitude  

# _______________________________________________________________








#                           #  mobr package  #
# 
#     # Creating richness rarefaction curves over treatments and env gradients #  
# 
# -----------


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


# # _______________________________________________________________

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
# # y_lim <- max(ind_based_rare_med$Richness) + max(ind_based_rare_med$Richness) * 0.05
# x_lim <- max(ind_based_rare_med$Individuals) + max(ind_based_rare_med$Individuals) * 0.05
# 
# par(mar = c(5,5,1,10))
# 
# plot(1, type = 'n', xlab = '', ylab = '',
#      ylim = c(0, y_lim),
#      xlim = c(0, x_lim))
# 
# lines(ind_based_rare_med$Individuals[], ind_based_rare_med$Richness)







# next steps:



# then ask myself what syntax the most complex function in mob r (the one that calculated the difference between the rarefactions) needs
# correct the rarefaction syntax I used to what is required
# created the other rarefactions

# and later later I think of adding analysis that are connected to species compositions - although this was done before



# if needed:



#if binned data required:
red_depth_bins <- wide_red

# 8.1 - 149 m to 8 layers
red_depth_bins$depth <- cut(red_depth_bins$depth,                                         
                            breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 149),
                            labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))



















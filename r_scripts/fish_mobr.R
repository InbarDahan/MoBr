
    
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
                   labels=c('0-15', '16-30', '31-45', '46-60', '61-75','76-90', '91-105', '106-120', '121-135', '136-149'))) %>% 
  relocate(depth_group,.after = lon_x) 

# count samples abundances:
species_abundances <- rowSums(red_depth_bins[, 8:ncol(red_depth_bins)])

# add the abundance column to the data frame
red_depth_bins$sample_abu <- species_abundances


red_depth_bins<-red_depth_bins %>% relocate(sample_abu,.after = depth_group)


test<-red_depth_bins %>% group_by(depth_group) %>% summarise(high_samples = sum(sample_abu >= 5))

length(which(red_depth_bins$sample_abu >= 5)) # 110\123

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

# the output is multiple graphs for the binned depths var. It seems like there is a decrease in abundance and richness with depth.
# but not in rarefied richness. Scale dependent.

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

plot(delta_red)

setwd(wd_results)
save(delta_red, file = 'deltas_red.Rdata')

#  Insights on the graphs:
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
                                # α and γ
#______________________________________________________________________________

                       # extract α and γ data
 
alphas_f <- stats_f$samples_stats
alphas_f$depth <- red_depth_bins$depth[match(alphas_f$group, red_depth_bins$depth_group)]

gammas_f <- stats_f$groups_stats  
gammas_f$depth <- red_depth_bins$depth[match(gammas_f$group, red_depth_bins$depth_group)]

# -----------------------------------------------------------------------------

                        # fit linear models 

lm_alpha_f = alphas_f %>% group_by(index, effort) %>%
  do(mod_f = lm(value ~ depth, data = .))

lm_gamma_f = gammas_f %>% group_by(index, effort) %>%
  do(mod_f = lm(value ~ depth, data = .))

# -----------------------------------------------------------------------------

                        # get model coefs

mod_coef_alpha_f= broom::tidy(lm_alpha_f, mod_f)
mod_coef_gamma_f= broom::tidy(lm_gamma_f, mod_f)

# -----------------------------------------------------------------------------

                      # get model summary

mod_sum_alpha_f= broom::glance(lm_alpha_f, mod_f)
mod_sum_gamma_f= broom::glance(lm_gamma_f, mod_f)

# -----------------------------------------------------------------------------

                  # binding the α and γ data stats

mob_met_f = rbind(data.frame(scale = 'alpha', alphas_f),
                  data.frame(scale = 'gamma', gammas_f))

# -----------------------------------------------------------------------------

             # reordering the levels of the "index" column
# binned
mob_met_f$index = factor(mob_met_f$index, 
                         levels = levels(mob_met_f$index)[c(2:1, 3:7)])

# -----------------------------------------------------------------------------

colnames(mob_met_f)[colnames(mob_met_f) == "group"] <- "depth_group"

           # plotting the effect of depth on total abundance
                            # α and γ

N_f = mob_met_f %>%
    subset(index %in% 'N') %>%
    ggplot(aes(x = depth, y = value)) + 
    geom_point(aes(color = scale)) +
    geom_smooth(aes(color = scale), method = 'lm', se = F) +
    labs(x = "depth (m)", y = "Total abundance (N)")

  ggsave("fish_grad_vs_N.pdf", plot = N_f, path = wd_plots, width = 15, height = 12, 
         units = "cm")

# -----------------------------------------------------------------------------
  
     # plotting the effect of depth on abundance & richness 
                           # α scale

SN1_f = mob_met_f%>% 
    subset(index %in% c('S', 'N')) %>%
    subset(scale == 'alpha') %>%
    ggplot(aes(x = depth, y = value )) + 
    geom_point(col = "indianred2") +
    geom_smooth(method = 'lm', se = T) +
    facet_wrap(~ index, scales = "free",  
               labeller = as_labeller(c(S = "α Species richness (α S)",
                                        N = "α Total abundance (α N)"))) +
                                        labs(x = "depth (m)")
  
 ggsave("fish_grad_alpha_vs_S&N.pdf", plot = SN1_f, path = wd_plots, width = 15, height = 12, 
        units = "cm")
 
 # ------
                          # γ scale
 
 SN_f = mob_met_f%>% 
   subset(index %in% c('S', 'N')) %>%
   subset(scale == 'gamma') %>%
   ggplot(aes(x = depth, y = value )) + 
   geom_point(col = "indianred2") +
   geom_smooth(method = 'lm', se = T) +
   facet_wrap(~ index, scales = "free",  
              labeller = as_labeller(c(S = "γ Species richness (γ S)",
                                       N = "γ Total abundance (γ N)"))) +
   labs(x = "depth (m)")
 
 ggsave("fish_grad_gamma_vs_S&N.pdf", plot = SN_f, path = wd_plots, width = 15, height = 12, 
        units = "cm")

# -----------------------------------------------------------------------------

      # plotting the effect of depth on 'S', 'S_n', 'S_PIE'

 p_f = mob_met_f %>% 
   subset(abs(value) < 1000) %>%
   subset(index %in% c('S', 'S_n', 'S_PIE')) %>% 
   ggplot(aes(x = depth, y = value, col = scale)) + 
   geom_point() +
   geom_smooth(method = 'lm', se = F) +
   facet_wrap(. ~ index, scales = "free")

 ggsave("fish_grad_vs_S&Sn&SPIE.pdf", plot = p_f, path = wd_plots, width = 15, height = 12, 
        units = "cm")
 
  # -----------------------------------------------------------------------------

     # plotting the effect of depth on 'S', 'S_n', 'S_PIE'
 
 # binned 
 p2_f= mob_met_f%>% 
   subset(abs(value) < 1000) %>%
   subset(index %in% c('beta_S', 'beta_S_n', 'beta_S_PIE')) %>% 
   ggplot(aes(x = depth, y = value)) + 
   geom_point() +
   geom_smooth(method = 'lm', se = F) +
   facet_wrap(. ~ index, scales = "free")
 
 g_f= ggarrange(p1_f, p2_f)
 
ggsave("ENS.pdf", plot = g_f, path = wd_plots, width = 20, height = 15, units = "cm")
  


  
  
  
  
  
  
  
  
  

# calculate biodiversity statistics: 

# biodiv <- calc_biodiv(sp_matric, groups = "depth", index = c('N', 'S', 'S_n', 'S_PIE'), effort = 5, extrapolate = FALSE, return_NA =  FALSE) 

# ----------------------------------------------------------------------------- 
# my code - comes from the mobR package instructions:
# ___________________________________________________

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






# and later later I think of adding analysis that are connected to species compositions - although this was done before

#if binned data required:
red_depth_bins <- wide_red

# 8.1 - 149 m to 8 layers
red_depth_bins$depth <- cut(red_depth_bins$depth,                                         
                            breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 149),
                            labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))



# Trying to chnage the colors
#values = c('white','yellow','orange','#f768a1','red','brown','purple','#3f007d', 'blue','black'
# fun_color_range <- colorRampPalette(c('#3f007d', '#e31a1c' ,'#f768a1', '#feb24c', 'yellow'))    
# my_colors <- fun_color_range(100) 















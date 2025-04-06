

            ###  MOBr - sponges (here marked as s) ###

   # what drivers shape the change in sponges richness over the depth gradient ?

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
sponges <- read.csv("sponges_data.csv")

# -----------------------------------------------------------------------------

                   # prepare the data

# define the first column of species
first_species <- 8

# create species metric
sp_matric <- sponges[,first_species:length(sponges)]

# -----------------------------------------------------------------------------

                       # make the "mob in" object
# create env var
env_red <- sponges[, c(1:7)]

# make mob in
sponges_mob_in <- make_mob_in(sp_matric, env_red, # define the sp matrix and the variables   # no wornings apear - v
                              coord_names = c('lon', 'lat'))   # define coor

# -----------------------------------------------------------------------------

# compute mob stats: 

# binned:
# here the α and γ scales are different. One sample vs depth bin
stats_f <- get_mob_stats(sponges_mob_in, group_var = "depth_group", extrapolate = FALSE) 

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
delta_red = get_delta_stats(mob_in = sponges_mob_in, env_var = 'depth',  # define the env gradient
                            group_var = 'depth_group', stats = c('betas', 'r'), # define the site\grouping var
                            type = 'continuous', spat_algo = 'kNCN', # define the type of var # define the spatial rarefaction method
                            n_perm = 100, overall_p = FALSE) # 

setwd(wd_plots)
pdf('sponges_mobr.pdf')
plot(delta_red, stat = 'b1', scale_by = 'indiv',    # stat can be changed based on the wanted effect size method 
     eff_sub_effort = F, eff_log_base = 2,          # subset\all samples   #   
     eff_disp_pts = T,                              # T\P - show the raw effect points
     eff_disp_smooth = F)                           # T\F - show the linear effect of the expl var on the effect sizes
dev.off()

setwd(wd_results)
save(delta_red, file = 'deltas_sponges.Rdata')


#  Insights on the graphs:
# (1) richness over depth - 
# (2.a) agg over depth - 
# (2.b)  N over depth  - 
# (2.c) SAD over depth - 
# (3.a) agg over scale - 
# (3.b)   N over scale - 
# (3.c) SAD over scale - 

# summary - there are changes in the behavior of different underlying components on richness over the depth gradient
#           the components are also influenced by scale.

#______________________________________________________________________________

# 2 scale analysis - continuous
# α and γ
#______________________________________________________________________________

# extract α and γ data

alphas_f <- stats_f$samples_stats
alphas_f$depth <- sponges$depth[match(alphas_f$group, sponges$depth_group)]

gammas_f <- stats_f$groups_stats  
gammas_f$depth <- sponges$depth[match(gammas_f$group, sponges$depth_group)]

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

ggsave("sponges_grad_vs_N.pdf", plot = N_f, path = wd_plots, width = 15, height = 12, 
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

ggsave("sponges_grad_alpha_vs_S&N.pdf", plot = SN1_f, path = wd_plots, width = 15, height = 12, 
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

ggsave("sponges_grad_gamma_vs_S&N.pdf", plot = SN_f, path = wd_plots, width = 15, height = 12, 
       units = "cm")

# -----------------------------------------------------------------------------

# plotting the effect of depth on 'S', 'S_n', 'S_PIE'

p1_f = mob_met_f %>% 
  subset(abs(value) < 1000) %>%
  subset(index %in% c('S', 'S_n', 'S_PIE')) %>% 
  ggplot(aes(x = depth, y = value, col = scale)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(. ~ index, scales = "free")

ggsave("sponges_grad_vs_S&Sn&SPIE.pdf", plot = p1_f, path = wd_plots, width = 15, height = 12, 
       units = "cm")

# -----------------------------------------------------------------------------

# plotting the effect of depth on 'S', 'S_n', 'S_PIE'

p2_f= mob_met_f%>% 
  subset(abs(value) < 1000) %>%
  subset(index %in% c('beta_S', 'beta_S_n', 'beta_S_PIE')) %>% 
  ggplot(aes(x = depth, y = value)) + 
  geom_point() +
  geom_smooth(method = 'lm', se = F) +
  facet_wrap(. ~ index, scales = "free")

g_f= ggarrange(p1_f, p2_f)

ggsave("sponges_ENS.pdf", plot = g_f, path = wd_plots, width = 20, height = 15, units = "cm")






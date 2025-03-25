

# General information: the samples were distributed in a nested way. 
# at each elevation - 50 m quadrature - in its corners - 10 m quadrature - in its corners 1 m quadrate
# total of 16 quadrates, with coordinated for each of them - to know which is closer to which

install.packages(c('devtools', 'mobr', 'leaflet', 'mapview', 'tidyr',
                   'vegan', 'dplyr', 'ggplot2', 'egg'))

library(mobr)
library(vegan)
library(dplyr)
library(ggplot2)
library(egg)
devtools::install_version("broom", version = "0.5.4", repos = "http://cran.us.r-project.org")
library(broom)

wd <-  "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr"

setwd(wd)

# read in data file
ant_dat <- read.csv('./suppl/gradients/smokies_all.csv')
names(ant_dat)    # gives the names of the columns
ant_mob_in <- make_mob_in(ant_dat[ , 5:42], ant_dat[ c(1:4, 43:50)], # define the sp matrix and the variables
                          coord_names = c('sample_x', 'sample_y'))   # define coor - 

# drop site with low number of individuals
ant_mob_in = subset(ant_mob_in, site != 'NODI', drop_levels = T) ##### what is the limit? 5 individuals? ####

# count samples abundances:
species_abundances_a <- rowSums(ant_dat[ , 5:42])

# add the abundance column to the data frame
ant_dat$sample_abu <- species_abundances_a
ant_dat<-ant_dat %>% relocate(sample_abu,.after = sample)

# check the number of samples with more then 5 individuals (sampling effort)
test_a <- ant_dat %>% group_by(elevation_m, site) %>% summarise(high_samples = sum(sample_abu >= 5))
length(which(ant_dat$sample_abu >= 5)) # 110\123 # this might be the problam in the graph

#______________________________________________________________________________

                   # 2 scale analysis - continuous
                            # alpha + delta
#______________________________________________________________________________


## compute mob stats-------------------

 ### They chanrged the function from get_mob_stat to that and now it doesn't work for me ###

stats_a <- calc_comm_div(ant_dat[ , 5:42],
                      index = c('N', 'S' ,'S_n', 's_C'),
                       effort = c(15, 80))

plot(stats_a)

alphas_a <- stats_a$samples_stats
alphas_a$elevation <- ant_dat$elevation_m[match(alphas_a$group, ant_dat$site)]

gammas_a <- stats_a$groups_stats  
gammas_a$elevation <- ant_dat$elevation_m[match(gammas_a$group, ant_dat$site)]

# fit linear models
lm_alpha_a = alphas_a %>% group_by(index, effort) %>%
    do(mod_a = lm(value ~ elevation, data = .))

lm_gamma_a = gammas_a %>% group_by(index, effort) %>%
    do(mod_a = lm(value ~ elevation, data = .))

# get model coefs

mod_coef_alpha_a = broom::tidy(lm_alpha_a, mod_a)
mod_coef_gamma_a = broom::tidy(lm_gamma_a, mod_a)

# get model summary
mod_sum_alpha_a = broom::glance(lm_alpha_a, mod_a)
mod_sum_gamma_a = broom::glance(lm_gamma_a, mod_a)

mob_met_a = rbind(data.frame(scale = 'alpha', alphas_a),
                data.frame(scale = 'gamma', gammas_a))
mob_met_a$index = factor(mob_met_a$index, 
                       levels = levels(mob_met_a$index)[c(2:1, 3:7)])

N1_a = mob_met_a %>%
     subset(index %in% 'N') %>%
     ggplot(aes(x = elevation, y = value)) + 
         geom_point(aes(color = scale)) +
         geom_smooth(aes(color = scale), method = 'lm', se = F) +
         labs(x = "Elevation (m)", y = "Total abundance (N)")
#ggsave("grad_vs_N.pdf", plot = N1, path = "./figs", width = 15, height = 12, 
#       units = "cm")

       
       
SN1_a = mob_met_a %>% 
     subset(index %in% c('S', 'N')) %>%
     subset(scale == 'alpha') %>%
     ggplot(aes(x = elevation, y = value)) + 
         geom_point() +
         geom_smooth(method = 'lm', se = T) +
         facet_wrap(~ index, scales = "free",  
                    labeller = as_labeller(c(S = "Species richness (S)",
                                             N = "Total abundance (N)"))) +
         labs(x = "Elevation (m)")

#ggsave('./figs/grad_vs_S&N.pdf', SN1)


p1_a = mob_met_a %>% 
  subset(abs(value) < 1000) %>%
  subset(index %in% c('S', 'S_n', 'S_PIE')) %>% 
  ggplot(aes(x = elevation, y = value, col = scale)) + 
    geom_point() +
    geom_smooth(method = 'lm', se = F) +
    facet_wrap(. ~ index, scales = "free")


p2_a = mob_met_a %>% 
    subset(abs(value) < 1000) %>%
    subset(index %in% c('beta_S', 'beta_S_n', 'beta_S_PIE')) %>% 
    ggplot(aes(x = elevation, y = value)) + 
    geom_point() +
    geom_smooth(method = 'lm', se = F) +
    facet_wrap(. ~ index, scales = "free")

g_a = ggarrange(p1_a, p2_a)

#ggsave("ENS.pdf", plot = g, path = "./figs", width = 20, height = 15, units = "cm")

#______________________________________________________________________________

                # Multiscale analysis - continuous 
                #     N, SAD and aggregations

#______________________________________________________________________________

deltas = get_delta_stats(ant_mob_in, env_var = 'elevation_m',
                         group_var = 'site', stats = c('betas', 'r'),
                         type = 'continuous', spat_algo = 'kNCN',
                         n_perm = 199, overall_p = FALSE)

#save(deltas, file = './results/deltas.Rdata')
#load('./results/deltas.Rdata')

                                                # The figure from the paper!!!
#pdf('./figs/deltas_b1.pdf')
plot(deltas, stat = 'b1', scale_by = 'indiv',
     eff_sub_effort = F, eff_log_base = 2,
     eff_disp_pts = F,
     eff_disp_smooth = T)
dev.off()

          # changing the measure of the effect size changes the graph completly

#pdf('./figs/deltas_r.pdf')
plot(deltas, stat = 'r', scale_by = 'indiv',
     eff_sub_effort = T, eff_log_base = 2,
     eff_disp_pts = T,
     eff_disp_smooth = F)
dev.off()

# exploration of centering effect size
deltas$S_df <- deltas$S_df %>%
    group_by(test, effort) %>%
    mutate(effect = effect - mean(effect)) %>%
    ungroup()

deltas$S_df[deltas$S_df$test == 'SAD', ] <-
    deltas$S_df %>%
    subset(test == 'SAD') %>%
    group_by(effort) %>%
    mutate(effect = effect - mean(effect))

#______________________________________________________________________________
 
      # compositional change over the elevation gradient 

#______________________________________________________________________________

## non-published supplemental analysis of compositional change ------------

# indirect gradient analysis -----------
empty_sites = rowSums(ant_mob_in$comm) == 0
ant_comm_no0 = ant_mob_in$comm[!empty_sites, ]

ant_mds = metaMDS(ant_comm_no0)
stressplot(ant_mds)

plot(ant_mds, display='sp', type='n')
orditorp(ant_mds, display='sp', priority = colSums(ant_comm_no0))
fit = envfit(ant_mds, ant_mob_in$env[!empty_sites, c('elevation_m', 'utm_n', 'utm_e')])
plot(fit)

# at the site rather than trap scale
ant_comm_agg = aggregate(ant_mob_in$comm, by = list(ant_mob_in$env$site), sum)[ , -1]
ant_dat_agg = aggregate(ant_mob_in$env, by = list(ant_mob_in$env$site),
                              function(x) x[1])[ , -1]

ant_mds = metaMDS(ant_comm_agg)

pdf('./figs/ant_mds.pdf')
plot(ant_mds, display='sp', type='n')
orditorp(ant_mds, display='sp', priority = colSums(ant_comm_no0))
fit = envfit(ant_mds, ant_dat_agg[, c('elevation_m', 'utm_n', 'utm_e')])
plot(fit)
stressplot(ant_mds)
dev.off()


ant_rda = rda(sqrt(ant_comm_agg) ~ elevation_m, data = ant_dat_agg)
ant_rda
RsquareAdj(ant_rda)
anova(ant_rda, by = 'terms', permutations = 2e3)
ant_mso = mso(ant_rda, ant_dat_agg[ , c('utm_e', 'utm_n')], grain = 5000, 
              permutations = 1000)

pdf('./figs/rda_results.pdf')
plot(ant_rda, display = c('sp', 'bp'), type='n')
orditorp(ant_rda, display = 'sp')
text(ant_rda, display = 'bp', col='red')
msoplot(ant_mso)
dev.off()

#elev = ant_dat_agg$elevation_m[rowSums(ant_comm) > 0]
#xy = ant_dat[rowSums(ant_comm) > 0, c('UTM_E', 'UTM_N')]
#ants = ant_comm[rowSums(ant_comm) > 0, ]
ant_cca = cca(sqrt(ant_comm_agg) ~ elevation_m, data = ant_dat_agg)
RsquareAdj(ant_cca)
anova(ant_cca, by = 'terms', permutations = 2e3)
ant_mso_cca = mso(ant_cca, ant_dat_agg[ , c('utm_e', 'utm_n')], grain = 5000, permutations = 1000)

pdf('./figs/cca_results.pdf')
plot(ant_cca, display = c('sp', 'bp'), type='n')
orditorp(ant_cca, display = 'sp')
text(ant_cca, display = 'bp', col='red')
msoplot(ant_mso_cca)
dev.off()

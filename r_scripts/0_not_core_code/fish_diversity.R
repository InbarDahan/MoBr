

                  ###  Diversity  ###

   # How does beta diversity changes within and among depth layers?

# libraries:
library(vegan)
library(tidyverse)
library(plotrix)
library(betapart)

# _______________________________________________________________

wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed"

# read data:
setwd(wd_processed_data)
wide_red <- read.csv("wide_red.csv")

# removing the wired column x:
wide_red = wide_red %>% dplyr::select(-X) 

# _______________________________________________________________

                 # data preperations

# define the first column of species
first_species <- 7

# create species metric
sp_matrix <- wide_red[,first_species:length(wide_red)]

# transform depth to bins:
red_depth_bins <- wide_red
red_depth_bins$depth <- cut(wide_red$depth,                                         
                            breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 149),
                            labels=c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

# _______________________________________________________________

          # Diversity Profile - Hill Numbers

# calculating *Hill number* values for each *depth layer*. We will use the `renyi` function.
#   x = data  
#   scales = desired q values  
#   hill = whether to calculate Hill numbers or not (T/F)  
# When the value is *0* you will see the *species richness*, when it is *1* you will get the *Shannon diversity index*.

# _______________________________________________________________

# a table of the Hill number for each sample (opcode) by the corresponding q number (scale).
renyi_profile = renyi(sp_matrix,  
                      scales = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2, 4, 8, 16, 32, 64, Inf),
                      hill = T)


# keep the ID and depth columns:
meta_data = red_depth_bins[,c(1,3)]


# combine the two data sets: Hill numbers for each q value + information on each row\sample:
renyi_df = bind_cols(meta_data,renyi_profile) 


# summarize the Q values and Hill number columns:
renyi_df = gather(renyi_df, "Q", "Value",3:length(renyi_df))  


# changing the type of variable to numeric
renyi_df$Q = as.numeric(renyi_df$Q) 

# _______________________________________________________________

             # Hill numbers over q

          # colored by depth per sample:  

ggplot(data = renyi_df ) +
  aes(x = Q, y = Value, group = OpCode , color = depth)+
  geom_line() +
  geom_point()+
  scale_x_continuous(limits = c(0,64))

# _______________________________________________________________

       # preparing hill numbers per depth layer

group_data = red_depth_bins %>% 
  dplyr::select(depth, first_species:length(red_depth_bins)) %>% 
  group_by(depth) %>% # use this unit for analyses
  summarise(across(.fns = sum),.groups = "keep") #summarize all other values by summing all the rows in each group

renyi_profile_group = renyi(group_data[,-c(1)],
                            scales = c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2, 4, 8, 16, 32, 64, Inf), 
                            hill = T)

renyi_df_group = bind_cols(group_data[,1],renyi_profile_group)


renyi_df_group = gather(renyi_df_group,"Q","Value",2:length(renyi_df_group))


renyi_df_group = renyi_df_group %>% arrange(depth)


renyi_df_group$Q = as.numeric(renyi_df_group$Q)

# _______________________________________________________________

             # Hill numbers over q

           # colored by depth layers:  

ggplot(renyi_df_group,aes(x=Q, y=Value, group = depth, color= depth))+
  geom_line(size = 1.5)+
  scale_x_continuous(limits = c(0,3)) + ggtitle ("Hill numbers for each depth layer - rare to common")

# _______________________________________________________________

                   # Pairwise Beta Diversity
              
# This analyses is good to evaluate turnover among samples of the same depth layer
# We will use both abundance and incidence data
# Bray-Curtis Dissimilarity Index  - uses abundance data
# Sorenson (and Simpson), Jaccard - use incidence data

bray_scores = list() # create an empty list to store my newly created data

for (i in unique(red_depth_bins$depth)){
  
  depth_data = red_depth_bins %>% filter(depth == i) # keep only the observation of sample i
  
  depth_sp_matrix = depth_data[,first_species:length(depth_data)] # create species matrix
  
  depth_bray_part = bray.part(depth_sp_matrix) # apply the function that calculate bray Curtis distances 
  
  depth_bray_results = depth_bray_part[[3]] # keep only the bray Curtis results
  
  depth_bray_results = as.numeric(depth_bray_results) # convert to numeric object
  
  mean_bray = mean(depth_bray_results) # calculate the mean bray curtis distance
  
  se_bray = std.error(depth_bray_results)# calculate SE
  
  depth = i # argument with the the name of the site
  
  bray_data = data.frame(depth,mean_bray,se_bray) # create data frame that save those variables 
  
  bray_scores[[i]] = bray_data # save it in my list
  
}

bray_scores = bind_rows(bray_scores) # convert from list to data frame

bray_scores$depth = factor(bray_scores$depth, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))


# plot:

ggplot(bray_scores,aes(x = depth,
                       y = mean_bray,
                       color = depth)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= mean_bray - se_bray,
                    ymax= mean_bray + se_bray),size =1.2,width = 0.2)
# We learn that the shallow samples are more similar to each other than the deep samples.

# _______________________________________________________________

              # Data transformations - log

        # gives less weight to the common species            

depth_sp_matrix_log_2 = decostand(sp_matrix, method = 'log',logbase = 2)
depth_sp_matrix_log_10 = decostand(sp_matrix, method = 'log',logbase = 10)

# _______________________________________________________________

depth_bray_scores_log_10 = list() # create an empty list to store my newly created data

for (i in unique(red_depth_bins$depth)){
  
  depth_data = red_depth_bins %>% filter(depth == i) # keep only the observation of sample i
  
  depth_sp_matrix = depth_data[,first_species:length(depth_data)] # create species matrix
  
  depth_sp_matrix_log = decostand(depth_sp_matrix,method = 'log',logbase = 10)
  
  depth_bray_part_log = bray.part(depth_sp_matrix_log) # apply the function that calculate bray Curtis distances 
  
  depth_bray_results_log = depth_bray_part_log[[3]] # keep only the bray Curtis results
  
  depth_bray_results_log = as.numeric(depth_bray_results_log) # convert to numeric object
  
  depth_mean_bray = mean(depth_bray_results_log) # calculate the mean bray curtis distance
  
  depth_se_bray = std.error(depth_bray_results_log)# calculate SE
  
  depth_Sample = i # argument with the the name of the site
  
  depth_bray_data_log = data.frame(depth_Sample,depth_mean_bray,depth_se_bray) # create data frame that save those variables 
  
  depth_bray_scores_log_10[[i]] = depth_bray_data_log # save it in my list
  
}

depth_bray_scores_log_10 = bind_rows(depth_bray_scores_log_10) # convert from list to data frame


#lets plot with transformation - reducing the effect of common species:

depth_bray_scores_log_10$depth_Sample = factor(depth_bray_scores_log_10$depth_Sample, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

ggplot(depth_bray_scores_log_10,aes(x = depth_Sample,
                                      y = depth_mean_bray,
                                      color = depth_Sample)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= depth_mean_bray - depth_se_bray,
                    ymax= depth_mean_bray + depth_se_bray),size =1.2,width = 0.2)+
  ggtitle("depth_log 10 - data transformation")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# and compare to the results without transformation:

bray_scores$depth = factor(bray_scores$depth, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

ggplot(bray_scores,aes(x = depth,
                       y = mean_bray,
                       color = depth)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= mean_bray - se_bray,
                    ymax= mean_bray + se_bray),size =1.2,width = 0.2)+
  ggtitle("depth_No data transformation")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# not too big of a diff between log and un_logged - the common species don't affect that much
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# cleaning the environment:
rm(list=setdiff(ls(), c("red_depth_bins","first_species","sp_matrix")))

# _______________________________________________________________

         # Different Components of Beta Diversity
  # Nestedness and turnover - analyses on incidence data

# We get 3 distance matrixes in one list using 'beta.pair':
# 1. The distance which is due to species *turnover* between sites - b_sim
# 2. The distance which is due to species *nestedness* between sites - b_nes
# 3. The total distance (*nestedness + turnonver*) - b_sor 

beta_pair_data = list()

for (i in unique(red_depth_bins$depth)){
  
  one_sample = red_depth_bins %>% filter(depth == i)
  
  one_sample_sp_matrix = one_sample[,7:length(one_sample)] # keep only species data  
  
  one_sample_incidence = one_sample_sp_matrix %>% replace(one_sample_sp_matrix > 0, 1) # convert to presence-absence data
  
  one_sample_beta_pairs = beta.pair(one_sample_incidence) # calculate the beta values
  
  one_sample_beta_pairs = bind_cols(one_sample_beta_pairs) # tie the list element toghter
  
  one_sample_beta_pairs$depth <- rep(i) # add site
  
  one_sample_beta_pairs = as.data.frame(one_sample_beta_pairs) # convert to data frame 
  
  one_sample_beta_pairs$beta.sim = as.numeric(one_sample_beta_pairs$beta.sim) # change from distance object to numeric
  one_sample_beta_pairs$beta.sne = as.numeric(one_sample_beta_pairs$beta.sne)
  one_sample_beta_pairs$beta.sor = as.numeric(one_sample_beta_pairs$beta.sor)
  
  beta_pair_data[[i]] = one_sample_beta_pairs # save to list
  
}

beta_pair_data = bind_rows(beta_pair_data) # convert list to data frame

beta_pair_data = gather(beta_pair_data,"component","beta",1:3) # long format

beta_pair_data = beta_pair_data %>% filter(component != "beta.sor")  # filter beta.sor (sum of sim and sne) 

# plot:

beta_pair_data$depth = factor(beta_pair_data$depth, levels = c('0-20', '21-40', '41-60', '61-80', '81-100', '101-120', '121-140', '141-149'))

ggplot(beta_pair_data)+
  aes(x=depth,y=beta,fill=component)+
  stat_summary(geom = "bar",fun.data = mean_se)+
  stat_summary(geom = "errorbar",fun.data = mean_se,width = 0.3)+
  ylab(expression(paste(beta,"-diversity")))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# I think - We can see what causes the differences among samples in each depth layer
# red is the turnover - change in species among samples
# blue is the nestedness 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 




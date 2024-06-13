

                  ###  Diversity  ###

   # How does beta diversity changes over depth layers?

        # Is ther a nestedness pattern in the data?

# libraries:
library(vegan)
library(tidyverse)
library(plotrix)
library(betapart)

# _______________________________________________________________

# wd:
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
              
# Good for analyses of turnover. We will use both abundance and incidence data.
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

sp_matrix_log_2 = decostand(sp_matrix, method = 'log',logbase = 2)
sp_matrix_log_10 = decostand(sp_matrix, method = 'log',logbase = 10)

# _______________________________________________________________

habitat_bray_scores_log_10 = list() # create an empty list to store my newly created data

for (i in unique(wide_day2$habitat)){
  
  habitat_data = wide_day2 %>% filter(habitat == i) # keep only the observation of sample i
  
  habitat_sp_matrix = habitat_data[,first_species:length(habitat_data)] # create species matrix
  
  habitat_sp_matrix_log = decostand(habitat_sp_matrix,method = 'log',logbase = 10)
  
  habitat_bray_part_log = bray.part(habitat_sp_matrix_log) # apply the function that calculate bray Curtis distances 
  
  habitat_bray_results_log = habitat_bray_part_log[[3]] # keep only the bray Curtis results
  
  habitat_bray_results_log = as.numeric(habitat_bray_results_log) # convert to numeric object
  
  habitat_mean_bray = mean(habitat_bray_results_log) # calculate the mean bray curtis distance
  habitat_se_bray = std.error(habitat_bray_results_log)# calculate SE
  habitat_Sample = i # argument with the the name of the site
  
  habitat_bray_data_log = data.frame(habitat_Sample,habitat_mean_bray,habitat_se_bray) # create data frame that save those variables 
  
  habitat_bray_scores_log_10[[i]] = habitat_bray_data_log # save it in my list
  
}

habitat_bray_scores_log_10 = bind_rows(habitat_bray_scores_log_10) # convert from list to data frame

# lets plot:

ggplot(habitat_bray_scores_log_10,aes(x = habitat_Sample,
                                      y = habitat_mean_bray,
                                      color = habitat_Sample)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= habitat_mean_bray - habitat_se_bray,
                    ymax= habitat_mean_bray + habitat_se_bray),size =1.2,width = 0.2)+
  ggtitle("habitat_log 10 - data transformation")

# and compare to the results without transformation:

ggplot(habitat_bray_scores,aes(x = habitat_Sample,
                               y = habitat_mean_bray,
                               color = habitat_Sample)) +
  geom_point(size = 4)+
  geom_errorbar(aes(ymin= habitat_mean_bray - habitat_se_bray,
                    ymax= habitat_mean_bray + habitat_se_bray),size =1.2,width = 0.2)+
  ggtitle("habitat_No data transformation")


























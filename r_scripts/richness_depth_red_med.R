

                           ### Richness ###
       # Does the fish community richness changes across seas and depths?
                            # sea: red
    # depths: (for now it's values but maybe I shoud trasform to bins of 10 m)

library(vegan)
library(tidyverse)
library(plotrix)
library(rareNMtests)
library(mobr)

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed"

setwd(wd_processed_data)
wide_opCode <- read.csv("wide_opCode.csv")
wide_opCode = wide_opCode %>% select(-X) # rm weird column

# _______________________________________________________________


# *Individual Based Rarefaction*

# a variable that indicate from which column the species list start:
first_species = 7 #the column form which the species list begin in wide_day2 dataset. 4,if we use the regular wide_day

# _______________________________________________________________

# subset our data to create the *species matrix data frame*, which we will call sp_matrix:
sp_matrix = wide_opCode[,first_species:length(wide_opCode)]

# _______________________________________________________________

# checking the minimum number of records in a sample (opcode)
raremax = sp_matrix %>% rowSums(.) %>% min() # 0 

# _______________________________________________________________

# when we rarefy by sample, we may see an extremely low individual count. Rarefying to a low number such as *2* isn’t really helpful.
# lets observe how abundance varies in our data set:
sp_matrix %>% 
  mutate(abundance = rowSums(.)) %>% 
  ggplot()+
  aes(x = abundance)+
  geom_histogram()+
  scale_x_log10() # for clarity

# _______________________________________________________________

# We plotted the x-axis on a log10 scale for better clarity. You can see that some knolls have extremely low abundances. 
# We will remove these samples for this demonstration and stay only with samples that have more than *10* individuals
wide_opCode_clean = wide_opCode %>% 
  mutate(abundance = 
           rowSums(wide_opCode[first_species:length(wide_opCode)])) %>% 
  filter(abundance > 10) %>%  # set a threshold
  mutate(abundance = NULL) # remove this column so we will have a fresh start

sp_matrix = wide_opCode_clean[,first_species:length(wide_opCode_clean)]

# - - - - - - 

# filter low abundance samples

# some samples have extremely low abundances (o). 
# We will remove these samples and stay only with samples that have more than *9* individuals
wide_opCode_clean = wide_opCode %>% 
  mutate(abundance = rowSums(wide_opCode[first_species:length(wide_opCode)])) %>% 
  filter(abundance > 4) %>%  # set a threshold (based on the histogram)
  mutate(abundance = NULL) # remove this column so we will have a fresh start

# -----------

# create a new sp_matrix based on the filtered df:
sp_matrix = wide_opCode_clean[,first_species:length(wide_opCode_clean)]

# check lowest abundance count now:
raremax_n = sp_matrix %>% rowSums(.) %>% min() # = 5
raremax_n

# -----------

# after filtering #
# new depth ranges and distributions #

# original:
length(which(wide_opCode$sea == "red")) # 123
length(which(wide_opCode$sea == "med")) # 126

# after filtering samples with less than 7 samples:
length(which(wide_opCode_clean$sea == "red")) # 110 - we lost 13 samples - 10 %
length(which(wide_opCode_clean$sea == "med")) # 89 - we lost 37 samples - 29 %

# - - - - -

# separate each sea observations:
wide_opCode_med <-  filter(wide_opCode_clean, sea == "med")
wide_opCode_red <-  filter(wide_opCode_clean, sea == "red")

# new depth ranges;
print(range(wide_opCode_med$depth, na.rm=TRUE)) # 5.6 - 141 m (was 5.6 - 141 m)
print(range(wide_opCode_red$depth, na.rm=TRUE)) # 8.1 - 144 m (was 8.1 - 144 m)

# new depth histograms
# med:
wide_opCode_med %>% 
  ggplot()+
  aes(x = depth)+     
  geom_histogram() # some depths have up to 6 observations\opcode

# red:
wide_opCode_red %>% 
  ggplot()+
  aes(x = depth)+
  geom_histogram() # most of the depths have 3-5 observations\opcodes



                    # united plots for all taxonomic groups

# libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(mapview)
# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________

# read data:
setwd(wd_processed_data)
taxon_data <- read.csv('combined_data.csv')
# _______________________________________________________________

             # Create abundance\depth plot:

ggplot(taxon_data, aes(x = depth, y = abundance, color = taxon)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # Add a trend line for each taxon
  labs(
    title = "Abundance vs. Depth",
    x = "Depth (m)",
    y = "Abundance per sample",
    color = "Taxonomic Group"
  )  +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )
# _______________________________________________________________

            # Create richness\depth plot:

ggplot(taxon_data, aes(x = depth, y = richness, color = taxon)) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE) +  # Add a trend line for each taxon
  labs(
    title = "Richness vs. Depth",
    x = "Depth (m)",
    y = "Richness per sample",
    color = "Taxonomic Group"
  )  +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )













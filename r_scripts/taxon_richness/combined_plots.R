

                    # united plots for all taxonomic groups

# libraries:
library(dplyr)
library(tidyr)
library(ggplot2)
library(mapview)
library(ggh4x)
# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _______________________________________________________________

                          # Full data 

# read data:
setwd(wd_processed_data)
taxon_data <- read.csv('combined_data.csv')
taxon_data <- taxon_data %>% dplyr::select(-X)

taxon_long <- read.csv('combined_long.csv')
taxon_long <- taxon_long %>% dplyr::select(-X)
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

# Create abundance\depth plot:

ggplot(taxon_data, aes(x = depth, y = abundance, color = taxon)) +
  geom_point() +
  #geom_smooth(method = "loess", se = FALSE) +  # Add a trend line for each taxon
  labs(
    title = "Abundance vs. Depth",
    x = "Depth (m)",
    y = "Abundance per sample",
    color = "Taxonomic Group"
  )  + 
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5)  # Center the title
  )+ facet_wrap(~ taxon)




            # Create richness\depth plot:

ggplot(taxon_data, aes(x = depth, y = richness, color = taxon)) +
  geom_point() +
  #geom_smooth(method = "loess", se = FALSE) +  # Add a trend line for each taxon
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
# _______________________________________________________________

            # facet plot measures\depth for each taxon:

# without colors
ggplot(taxon_long, aes(x = depth, y = value)) +  
  geom_point() +  
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(taxon ~ diversity_measure, scales = "free") +  # Facet by diversity measure and taxonomic group
  theme_minimal() +
  labs(x = "Depth", y = "Value", title = "Diversity Measures by Taxonomic Groups")

# with colors
ggplot(taxon_long, aes(x = depth, y = value, color = taxon)) +  # Color by taxonomic group
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggh4x::facet_grid2(taxon ~ diversity_measure  ,  scales = "free_y", independent = "y") +  # Facet and move categories to the left
  theme_minimal() +
  theme(
    strip.placement = "outside",  # Puts the facet labels outside
    strip.text.y.left = element_text(angle = 0),  # Keeps the facet labels horizontal
    strip.text.y = element_blank(),
    plot.title = element_text(hjust = 0.5),  # Centers the plot title
    legend.title = element_blank(),
    panel.spacing = unit(2, "lines")
  ) +
  labs(x = "depth", y = "Value", title = "Diversity Measures by Taxonomic Groups") +
  scale_color_manual(values = c("purple", "red", "orange", "yellow3"))  # Customize colors for taxonomic groups

# ________________________________________________________________

           # filter data - by location (NR SOUTH)

# read data:
setwd(wd_processed_data)
NR_SOUTH_data <- read.csv('NR_SOUTH_data.csv')
NR_SOUTH_data <- NR_SOUTH_data %>% dplyr::select(-X)
# _______________________________________________________________

         # facet plot measures\depth for each taxon:

ggplot(NR_SOUTH_data, aes(x = depth, y = value, color = taxon)) +  
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggh4x::facet_grid2(taxon ~ diversity_measure ,scales = "free_y", independent = "y") +   
  theme_minimal() +
  theme(
    strip.placement = "outside",  
    strip.text.y.left = element_text(angle = 0),  
    plot.title = element_text(hjust = 0.5),  
    legend.title = element_blank()
  ) +
  labs(x = "depth", y = "Value", title = "Diversity Measures by Taxonomic Groups (NR SOUTH)") +
  scale_color_manual(values = c("purple","darkred", "darkorange", "yellow"))  


# ________________________________________________________________

# ________________________________________________________________

        # run just for the species level and not higher levels - stony









# other design options:

# all data facet:
ggplot(taxon_long, aes(x = depth, y = value, color = taxon)) +
  geom_point(size = 1.5, alpha = 0.8) +  # Moderate point size and transparency
  geom_smooth(method = "lm", se = FALSE, size = 0.8, linetype = "dashed") +  # Dashed lines for linear trends
  ggh4x::facet_grid2(taxon ~ diversity_measure, scales = "free_y", independent = "y") +
  theme_bw() +  # Classic, clean theme suitable for scientific articles
  theme(
    strip.placement = "outside",  # Keep facet labels outside
    strip.text.y = element_blank(),  # Remove facet labels on the y-axis
    strip.background = element_blank(),  # Remove background color from facet labels
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),  # Bold, centered title
    plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 10)),  # Subtle subtitle with spacing
    axis.title = element_text(size = 11, face = "bold"),  # Bold axis titles for emphasis
    axis.text = element_text(size = 9),  # Slightly smaller axis text for clarity
    legend.position = "right",  # Position legend on the right for a cleaner layout
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 9),  # Moderate legend text size
    panel.grid.major = element_line(color = "grey80", linetype = "dotted"),  # Subtle grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.border = element_rect(color = "black", fill = NA)  # Black border for panels
  ) +
  labs(
    x = "Depth (m)",  # Include units
    y = "Value",
    title = "Diversity Measures Across Depth Gradients",
    subtitle = "Comparing Taxonomic Groups by Measure"
  ) +
  scale_color_manual(values = c("purple","darkred", "steelblue", "darkgreen"))  # Softer, professional colors


# Facet plot of diversity measures by depth for each taxon with updated colors
ggplot(NR_SOUTH_data, aes(x = depth, y = value, color = taxon)) +  
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  ggh4x::facet_grid2(taxon ~ diversity_measure, scales = "free_y", independent = "y") +   
  theme_minimal() +
  theme(
    strip.placement = "outside",  
    strip.text.y.left = element_text(angle = 0),  
    plot.title = element_text(hjust = 0.5),  
    legend.title = element_blank()
  ) +
  labs(x = "Depth (m)", y = "Value", title = "Diversity Measures by Taxonomic Groups (NR SOUTH)") +
  scale_color_manual(values = c("purple","#1B4F72", "#28B463", "#F1C40F"))  # Updated colors

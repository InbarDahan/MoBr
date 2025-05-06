

                      # The final agg code

# Calculating aggregations using 3 methods for each species in each treatment (depth):

     # 1 - the Negative binomial dispersion parameter (k)
                   # 2 - Ives J index
                # 3 - Morisita Index (Iᵈ)

# Required packages
library(dplyr)
library(tidyr)
library(purrr)
library(fitdistrplus)  # for Negative Binomial k
library(vegan)         # for Morisita index
library(ggplot2)
library(patchwork)     # to show plots side by side
library(MASS)          # for negative binomial 

# _______________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_plots <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/output/plots/agg_plot" 
# _______________________________________________________________

# read data:
setwd(wd_processed_data)
fish <- read.csv("wide_red.csv")
sponges <- read.csv("sponges_data.csv")
stony <- read.csv("stony_wide_data.csv")
soft <- read.csv("soft_data.csv")

fish = fish %>% dplyr::select(-X) 
stony  = stony %>% dplyr::select(-X) 
soft = soft %>% dplyr::select(-X) 
# _______________________________________________________________

# data preparations for this analysis:

# renaming columns:
fish <- fish %>% rename(sample = OpCode)
sponges <- sponges %>% rename(sample = quad)
stony <- stony %>% rename(sample = plot, depth_group = depth)
soft <- soft %>% rename(depth_group = depth)

# ------------

# rename
fish <- fish %>% rename(lon = lon_x)
fish <- fish %>% rename(lat = lat_y)
# ------------

# adding a column for the taxonomic groups id:
fish    <- fish    %>% mutate(taxon = 'fish',           .before = 1)
sponges <- sponges %>% mutate(taxon = 'sponges',        .before = 1)
stony   <- stony   %>% mutate(taxon = 'stony',   .before = 1)
soft    <- soft    %>% mutate(taxon = 'soft',    .before = 1)
# ------------

# binning the fish data:
# 8.1 - 149 m to - 15 m bins (10 breaks)
fish <- fish %>% mutate(depth_group = cut(fish$depth,                                         
                                          breaks=c(0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 149),
                                          labels=c('7.5', '22.5', '37.5', '52.5', '67.5','82.5', '97.5', '112.5', '127.5', '142.5'))) %>% 
  relocate(depth_group,.after = depth) 
# ------------

fish    <- fish %>% dplyr::select(-sea, -temperature)
sponges <- sponges %>% dplyr::select(-Habitat, -Site)
soft    <- soft %>% dplyr::select(-site)
stony   <- stony %>% dplyr::select(-site)
# ------------

taxon_list <- list(fish = fish, sponges = sponges, stony = stony, soft = soft)

# _______________________________________________________________
 
        # Step 1: Define aggregation functions:

# J index (Ives)
calc_J <- function(ab_vec) {
  if (length(ab_vec) < 2 || sum(ab_vec) < 5 || all(ab_vec == 0)) return(NA)
  xj <- ab_vec
  X <- mean(xj)
  N <- length(xj)
  J <- (sum(xj * (xj - 1)) / (X * N) - X) / X
  return(J)
}

# - - - 

# Negative binomial k using glm.nb
calc_k <- function(ab_vec) {
  if (length(ab_vec) < 2 || sum(ab_vec) < 5 || var(ab_vec) == 0) return(NA)
  tryCatch({
    fit <- MASS::fitdistr(ab_vec, densfun = "Negative Binomial")
    return(fit$estimate["size"])  # 'size' is k
  }, error = function(e) return(NA))
}

# calc_k <- function(ab_vec) {
#   if (length(ab_vec) < 2 || sum(ab_vec) < 5 || var(ab_vec) == 0) return(NA)
#   tryCatch({
#     df <- data.frame(count = ab_vec, x = seq_along(ab_vec))  # Dummy predictor
#     model <- MASS::glm.nb(count ~ 1, data = df)  # Intercept-only model
#     return(model$theta)  # theta is the dispersion parameter (k)
#   }, error = function(e) return(NA))
# }

# - - - 

# Morisita index (standardized)
calc_morisita <- function(ab_vec) {
  if (length(ab_vec) < 2 || sum(ab_vec) < 5 || all(ab_vec == 0)) return(NA)
  mat <- matrix(ab_vec, ncol = 1)
  mori <- vegan::dispindmorisita(mat)[1, 1]
  return(mori)
}

# _______________________________________________________________

     # Step 2: Apply Functions to Each Taxon by Depth Group

agg_results <- map_dfr(names(taxon_list), function(taxon_name) {
  df <- taxon_list[[taxon_name]]
  
  df_long <- df %>%
    pivot_longer(cols = 6:ncol(.), names_to = "species", values_to = "abundance") %>%
    filter(abundance > 0)  # exclude species not observed
  
  df_grouped <- df_long %>%
    mutate(depth_group = as.character(depth_group)) %>%
    group_by(depth_group, species) %>%
    summarise(
      abund_vec = list(abundance),
      .groups = "drop"
    ) %>%
    mutate(
      J = map_dbl(abund_vec, calc_J),
      k = map_dbl(abund_vec, calc_k),
      morisita = map_dbl(abund_vec, calc_morisita),
      taxon = taxon_name
    )
  
  return(df_grouped)
})

# _______________________________________________________________


        # Step 3: Plotting the Results Side by Side

# Convert depth_group to numeric if possible (for better spacing on x-axis)
 agg_results <- agg_results %>%
   mutate(depth_group_numeric = as.numeric(as.character(depth_group)))

# Pivot longer to make it tidy for ggplot
agg_long <- agg_results %>%
  pivot_longer(cols = c(J, k, morisita), names_to = "metric", values_to = "value")

# Plot
ggplot(agg_long, aes(x = depth_group_numeric, y = value)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  facet_wrap(~ metric + taxon, scales = "free") +  # Free x and y axes
  theme_minimal() +
  labs(
    title = "Intraspecific Aggregation Metrics by Depth",
    x = "Depth group (m)",
    y = "Aggregation measure"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)
  )

# _______________________________________________________________

                  # fish alone

df <- taxon_list[['fish']]

df_long <- df %>%
  pivot_longer(cols = 6:ncol(.), names_to = "species", values_to = "abundance") %>%
  filter(abundance > 0)  # exclude species not observed

df_grouped <- df_long %>%
  mutate(depth_group = as.character(depth_group)) %>%
  group_by(depth_group, species) %>%
  summarise(
    abund_vec = list(abundance),
    .groups = "drop"
  ) %>%
  mutate(
    J = map_dbl(abund_vec, calc_J),
    k = map_dbl(abund_vec, calc_k),
    morisita = map_dbl(abund_vec, calc_morisita),
    taxon = 'fish'
  )

df_grouped <- df_grouped %>%
  mutate(depth_group_numeric = as.numeric(as.character(depth_group)))

# Pivot longer to make it tidy for ggplot
df_grouped_long <- df_grouped %>%
  pivot_longer(cols = c(J, k, morisita), names_to = "metric", values_to = "value")

# Plot
ggplot(df_grouped_long, aes(x = depth_group_numeric, y = value)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "blue") +
  facet_wrap(~ metric, scales = "free") +  # Free x and y axes
  theme_minimal() +
  labs(
    title = "Intraspecific Aggregation Metrics by Depth",
    x = "Depth group (m)",
    y = "Aggregation measure"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)
  )




##############

# Conclusions:
# the results are consistent among the methods
# as well with the mobr analysis for the fish:
# aggregations are higher around 30 m and then go down.

# next step - choose a method and add it to the general graph:

# save J index to add it to the main graph:
j_results <- agg_long %>% filter(metric == "J")

# remove NA in the value column
j_results_clean <- j_results %>% as.data.frame() %>% filter(!is.na(value))%>%
  dplyr::select(-abund_vec)  

setwd(wd_processed_data)
write.csv(j_results_clean, "agg_j_results.csv")




























                          # model visualization
                        

library(dplyr)
library(tidyr)
library(ggplot2)
library(senlm)
library(purrr)

# wd:
wd_models <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed/models" 

setwd(wd_models)
one_model_a <- readRDS("one_model_a.rds")
one_model_r <- readRDS("one_model_r.rds")
one_model_e <- readRDS("one_model_e.rds")
one_model_agg <- readRDS("one_model_agg.rds")

subsets_taxons <- readRDS("subsets_taxons.rds")
subsets_taxons_clean <- readRDS("subsets_taxons_clean.rds")
# _____________________________________________________________

#   Template for plotting the raw data
dev.off() 
# fish (subset_taxons[[1]]) - abundance
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, 
     las = 1, xlab = "Depth (m)",  
     ylab = "abundance per sample", ylim = c(0,200))  
mygrey <- grey(level = 0.65, alpha = 0.4)  
points(subsets_taxons[[1]]$depth, subsets_taxons$fish$abundance, pch = 21, col = mygrey, bg = mygrey)

# fish - mixgaussian_negbin
plot(model_a_fish)
mygrey <- grey(level = 0.65, alpha = 0.4)  
points(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, 
       pch = 21, col = mygrey, bg = mygrey)

#   Alternative for plotting the Model - with the model name on top

# _____________________________________________________________

# abundance:

                                # mean function:

# 6 - fit the best model for visualization:
# mixgaussian_negbin
model_a_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[1]]$mean_fun, err_dist = one_model_a[[1]]$err_dist)
# mixgaussian_zinb
model_a_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[2]]$mean_fun, err_dist = one_model_a[[2]]$err_dist)
# gaussian_negbin
model_a_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[3]]$mean_fun, err_dist = one_model_a[[3]]$err_dist)
# gaussian_negbin
model_a_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[4]]$mean_fun, err_dist = one_model_a[[4]]$err_dist)

# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed

      # Fish - mixgaussian_negbin
    predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
    fitted_values <- predict(model_a_fish, predict.x)
    plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="Abundance", xlab = "Depth")
    lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

    # Sponges - mixgaussian_zinb
    predict.x <- seq(from=min(subsets_taxons[[2]]$depth), to = max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]])) 
    fitted_values <- predict(model_a_sponges, predict.x)
    plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="Abundance", xlab = "Depth")
    lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

    # Stony corals - gaussian_negbin
    predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
    fitted_values <- predict(model_a_stony, predict.x)
    plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="Abundance", xlab = "Depth")
    lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

    # Soft corals - gaussian_negbin
    predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
    fitted_values <- predict(model_a_soft, predict.x)
    plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="Abundance", xlab = "Depth")
    lines(predict.x, fitted_values, col = "#228875", lwd = 2) 

  # Reset layout to default
par(mfrow = c(1, 1))
 
# _____________________________________________________________
     
# richness:
     
                        # mean function:
     
# 6 - fit the best model for visualization:
# mixgaussian_zinb     
model_r_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[1]]$mean_fun, err_dist = one_model_r[[1]]$err_dist)
# sech_zinb
model_r_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[2]]$mean_fun, err_dist = one_model_r[[2]]$err_dist)
# gaussian_poisson
model_r_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[3]]$mean_fun, err_dist = one_model_r[[3]]$err_dist)
# gaussian_poisson
model_r_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[4]]$mean_fun, err_dist = one_model_r[[4]]$err_dist)

# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed

# Fish - mixgaussian_zinb
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_r_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="Richness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Sponges - sech_zinb
predict.x <- seq(from=min(subsets_taxons[[2]]$depth), to = max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]])) 
fitted_values <- predict(model_r_sponges, predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="Richness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Stony corals - gaussian_poisson
predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
fitted_values <- predict(model_r_stony, predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="Richness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Soft corals - gaussian_poisson
predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
fitted_values <- predict(model_r_soft, predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="Richness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#224188", lwd = 2) 

# Reset layout to default
par(mfrow = c(1, 1))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# A loop that does the same thing:

# Model fitting
model_r_list <- lapply(1:4, function(i) {
  senlm(data = subsets_taxons[[i]],
        xvar = "depth",
        yvar = "richness",
        mean_fun = one_model_r[[i]]$mean_fun,
        err_dist = one_model_r[[i]]$err_dist)
})



# Set up multi-panel layout and margins
par(mfrow = c(2, 2), mar = c(2, 2, 4, 1), oma = c(5, 5, 2, 1))  # reduced inner margins, added outer margins

for (i in 1:4) {
  data_i <- subsets_taxons[[i]]
  model_i <- model_r_list[[i]]
  predict.x <- seq(from = min(data_i$depth), to = max(data_i$depth), length.out = nrow(data_i))
  fitted_values <- predict(model_i, predict.x)
  
  plot(data_i$depth, data_i$richness, pch = 21, col = mygrey, bg = mygrey,
       main = titles[i], xaxt = "n", yaxt = "n")  # suppress x and y axis labels
  axis(1)  # re-add x axis
  axis(2)  # re-add y axis
  lines(predict.x, fitted_values, col = "#224188", lwd = 2)
}

# Add shared x and y axis labels outside the plotting area
mtext("Depth", side = 1, outer = TRUE, line = 3, cex = 1.2)
mtext("Richness", side = 2, outer = TRUE, line = 3, cex = 1.2)

# Reset layout
par(mfrow = c(1, 1))


# _____________________________________________________________

# evenness:

                        # mean function:

# 6 - fit the best model for visualization:
# mixgaussian_negbin
model_e_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[1]]$mean_fun, err_dist = one_model_e[[1]]$err_dist)
# mixgaussian_zinb
model_e_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[2]]$mean_fun, err_dist = one_model_e[[2]]$err_dist)
# gaussian_negbin
model_e_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[3]]$mean_fun, err_dist = one_model_e[[3]]$err_dist)
# gaussian_negbin
model_e_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[4]]$mean_fun, err_dist = one_model_e[[4]]$err_dist)

# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed

# Fish - mixgaussian_negbin
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_e_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Sponges - mixgaussian_zinb
predict.x <- seq(from=min(subsets_taxons[[2]]$depth), to = max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]])) 
fitted_values <- predict(model_e_sponges, predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Stony corals - gaussian_negbin
predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
fitted_values <- predict(model_e_stony, predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Soft corals - gaussian_negbin
predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
fitted_values <- predict(model_e_soft, predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2) 


################# If the previous chank runs, I can delete this######################

# # Reset layout to default
# par(mfrow = c(1, 1))
# # 6 - fit the best model for visualization:
# # gaussian_zig    
# model_e_fish <- senlm(data = subsets_taxons_clean[[1]], xvar = "depth", yvar = "sigma",  mean_fun = "gaussian", err_dist = "zig")
# # mixgaussian_ziig
# model_e_sponges <- senlm(data = subsets_taxons_clean[[2]], xvar = "depth", yvar = "sigma",  mean_fun = "mixgaussian", err_dist = "ziig")
# # mixgaussian_ziig
# model_e_stony <- senlm(data = subsets_taxons_clean[[3]], xvar = "depth", yvar = "sigma",  mean_fun = "mixgaussian", err_dist = "ziig")
# # gaussian_tweedie
# model_e_soft <- senlm(data = subsets_taxons_clean[[4]], xvar = "depth", yvar = "sigma",  mean_fun = "gaussian", err_dist = "tweedie")
# 
# # Set up a multi-panel layout: 4 rows and 2 columns
# par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed
# 
#       # Fish - mixgaussian_ziig
#     predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
#     fitted_values <- predict(model_e_fish, predict.x)
#     plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
#      main = "Fishes", ylab="Evenness", xlab = "Depth")
#     lines(predict.x, fitted_values, col = "#882247", lwd = 2)   
# 
#     # Sponges - mixgaussian_ziig
#     predict.x <- seq(from=min(subsets_taxons_clean[[2]]$depth), to = max(subsets_taxons_clean[[2]]$depth), length.out = nrow(subsets_taxons_clean[[2]])) 
#     fitted_values <- predict(model_e_sponges, predict.x)
#     plot(subsets_taxons_clean[[2]]$depth, subsets_taxons_clean[[2]]$sigma, pch = 21, col = mygrey, bg = mygrey,
#      main = "Sponges", ylab="Evenness", xlab = "Depth")
#     lines(predict.x, fitted_values, col = "#882247", lwd = 2)  
# 
#     # Stony corals - mixgaussian_ziig
#     predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
#     fitted_values <- predict(model_e_stony, predict.x)
#     plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$sigma, pch = 21, col = mygrey, bg = mygrey,
#      main = "Stony Corals", ylab="Evenness", xlab = "Depth")
#     lines(predict.x, fitted_values, col = "#882247", lwd = 2)  
# 
#     # Soft corals - gaussian_tweedie
#     predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
#     fitted_values <- predict(model_e_soft, predict.x)
#     plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$sigma, pch = 21, col = mygrey, bg = mygrey,
#     main = "Soft Corals", ylab="Evenness", xlab = "Depth")
#     lines(predict.x, fitted_values, col = "#882247", lwd = 2) 
# 
# 
# # Reset layout to default
# par(mfrow = c(1, 1))

# _____________________________________________________________

# aggregations:

                    # mean function:


# 6 - fit the best model for visualization:
# gaussian_zig    
model_e_fish <- senlm(data = subsets_taxons_clean[[1]], xvar = "depth", yvar = "sigma",  mean_fun = "gaussian", err_dist = "zig")
# mixgaussian_ziig
model_e_sponges <- senlm(data = subsets_taxons_clean[[2]], xvar = "depth", yvar = "sigma",  mean_fun = "mixgaussian", err_dist = "ziig")
# mixgaussian_ziig
model_e_stony <- senlm(data = subsets_taxons_clean[[3]], xvar = "depth", yvar = "sigma",  mean_fun = "mixgaussian", err_dist = "ziig")
# gaussian_tweedie
model_e_soft <- senlm(data = subsets_taxons_clean[[4]], xvar = "depth", yvar = "sigma",  mean_fun = "gaussian", err_dist = "tweedie")

# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed

# Fish - mixgaussian_ziig
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_e_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)   

# Sponges - mixgaussian_ziig
predict.x <- seq(from=min(subsets_taxons_clean[[2]]$depth), to = max(subsets_taxons_clean[[2]]$depth), length.out = nrow(subsets_taxons_clean[[2]])) 
fitted_values <- predict(model_e_sponges, predict.x)
plot(subsets_taxons_clean[[2]]$depth, subsets_taxons_clean[[2]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Stony corals - mixgaussian_ziig
predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
fitted_values <- predict(model_e_stony, predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Soft corals - gaussian_tweedie
predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
fitted_values <- predict(model_e_soft, predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2) 


# Reset layout to default
par(mfrow = c(1, 1))


# ______________________________________________________

# aggregations:


# 6 - fit the best model for visualization:
# summary(one_model_agg[[4]])
# gaussian_zig    
model_agg_fish <- senlm(data = subsets_taxons_agg[[1]], xvar = "depth_group", yvar = "value",  mean_fun = "gaussian", err_dist = "zig")
# mixgaussian_ziig
model_agg_sponges <- senlm(data = subsets_taxons_agg[[2]], xvar = "depth_group", yvar = "value",  mean_fun = "gaussian", err_dist = "gaussian")
# mixgaussian_ziig
model_agg_stony <- senlm(data = subsets_taxons_agg[[3]], xvar = "depth_group", yvar = "value",  mean_fun = "gaussian", err_dist = "tweedie")
# gaussian_tweedie
model_agg_soft <- senlm(data = subsets_taxons_agg[[4]], xvar = "depth_group", yvar = "value",  mean_fun = "gaussian", err_dist = "tweedie")

# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2))  # Adjust margins as needed

# Fish - 
predict.x <- seq(from=min(subsets_taxons_agg[[1]]$depth_group), to = max(subsets_taxons_agg[[1]]$depth_group), length.out = nrow(subsets_taxons_agg[[1]])) 
fitted_values <- predict(model_agg_fish, predict.x)
plot(subsets_taxons_agg[[1]]$depth_group, subsets_taxons_agg[[1]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="aggregations", xlab = "Depth")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)   

# Sponges -
predict.x <- seq(from=min(subsets_taxons_agg[[2]]$depth_group), to = max(subsets_taxons_agg[[2]]$depth_group), length.out = nrow(subsets_taxons_agg[[2]])) 
fitted_values <- predict(model_agg_sponges, predict.x)
plot(subsets_taxons_agg[[2]]$depth_group, subsets_taxons_agg[[2]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="aggregations", xlab = "Depth")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)  

# Stony corals - mixgaussian_ziig
predict.x <- seq(from=min(subsets_taxons_agg[[3]]$depth_group), to = max(subsets_taxons_agg[[3]]$depth_group), length.out = nrow(subsets_taxons_agg[[3]])) 
fitted_values <- predict(model_agg_stony, predict.x)
plot(subsets_taxons_agg[[3]]$depth_group, subsets_taxons_agg[[3]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="aggregations", xlab = "Depth")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)  

# Soft corals - 
predict.x <- seq(from=min(subsets_taxons_agg[[4]]$depth_group), to = max(subsets_taxons_agg[[4]]$depth_group), length.out = nrow(subsets_taxons_agg[[4]])) 
fitted_values <- predict(model_agg_soft, predict.x)
plot(subsets_taxons_agg[[4]]$depth_group, subsets_taxons_agg[[4]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="aggregations", xlab = "Depth")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2) 

# Reset layout to default
par(mfrow = c(1, 1))
























# ______________________________________________________

# ______________________________________________________

# creating a colorful version for richness

# Define models, names, and semi-transparent colors
models <- list(model_r_fish, model_r_sponges, model_r_stony, model_r_soft)
taxon_names <- c("Fishes", "Sponges", "Stony Corals", "Soft Corals")

# Base colors
base_colors <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a")
# Add transparency (alpha = 0.4)
point_colors <- scales::alpha(base_colors, 0.4)

# Set up 2x2 plot layout with shared axes
par(mfrow = c(2, 2), mar = c(3, 3, 2, 1), oma = c(4, 4, 2, 1))  # Adjust margins

# Plot loop
for (i in seq_along(subsets_taxons)) {
  data <- subsets_taxons[[i]]
  model <- models[[i]]
  taxon <- taxon_names[i]
  color <- point_colors[i]
  
  # Prediction
  predict.x <- seq(min(data$depth), max(data$depth), length.out = nrow(data))
  fitted_values <- predict(model, predict.x)
  
  # Plot
  plot(data$depth, data$richness,
       pch = 21, col = base_colors[i], bg = color,  # Border is solid, fill is transparent
       main = taxon, axes = FALSE, xlab = "", ylab = "",
       xlim = range(data$depth), ylim = range(data$richness, fitted_values)
  )
  lines(predict.x, fitted_values, col = "mygrey", lwd = 2)
  
  # Add axes only on outer plots
  if (i %% 2 == 1) axis(2)   # left side
  if (i > 2) axis(1)         # bottom row
  box()
}

# Shared axis titles
mtext("Depth (m)", side = 1, outer = TRUE, line = 2.5, cex = 1.2)
mtext("Richness per sample", side = 2, outer = TRUE, line = 2.5, cex = 1.2)


# ______________________________________________________

 

     
           # abundance and richness together:

# Set up an 8-panel layout: 4 rows and 2 columns
par(mfrow = c(4, 2), mar = c(4, 4, 2, 1))  # Adjust margins if needed

# Define a grey color for points
mygrey <- grey(level = 0.65, alpha = 0.4)

      # Fish - Richness
     plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, 
     las = 1, xlab = "Depth (m)", ylab = "Richness per sample", 
     ylim = c(0, 60), main = "Fish - Richness")
     points(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_r_fish, col = "blue", lwd = 2)  # Add trend line

    # Fish - Abundance
    plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, 
     las = 1, xlab = "Depth (m)", ylab = "Abundance per sample", 
     ylim = c(0, 200), main = "Fish - Abundance")
    points(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_a_fish, col = "blue", lwd = 2)  # Add trend line

    # Sponges - Richness
    plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, 
     las = 1, xlab = "Depth (m)", ylab = "Richness per sample", 
     ylim = c(0, 30), main = "Sponges - Richness")
    points(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_r_sponges, col = "blue", lwd = 2)  # Add trend line

    # Sponges - Abundance
    plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, 
     las = 1, xlab = "Depth (m)", ylab = "Abundance per sample", 
     ylim = c(0, 50), main = "Sponges - Abundance")
    points(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_a_sponges, col = "blue", lwd = 2)  # Add trend line

    # Stony Corals - Richness
    plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, 
     las = 1, xlab = "Depth (m)", ylab = "Richness per sample", 
     ylim = c(0, 80), main = "Stony Corals - Richness")
    points(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_r_stony, col = "blue", lwd = 2)  # Add trend line

    # Stony Corals - Abundance
    plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, 
     las = 1, xlab = "Depth (m)", ylab = "Abundance per sample", 
     ylim = c(0, 300), main = "Stony Corals - Abundance")
    points(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_a_stony, col = "blue", lwd = 2)  # Add trend line

    # Soft Corals - Richness
    plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, 
     las = 1, xlab = "Depth (m)", ylab = "Richness per sample", 
     ylim = c(0, 50), main = "Soft Corals - Richness")
    points(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_r_soft, col = "blue", lwd = 2)  # Add trend line

    # Soft Corals - Abundance
    plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, 
     las = 1, xlab = "Depth (m)", ylab = "Abundance per sample", 
     ylim = c(0, 300), main = "Soft Corals - Abundance")
    points(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, 
       pch = 21, col = mygrey, bg = mygrey)
    abline(model_a_soft, col = "blue", lwd = 2)  # Add trend line
  
# Reset layout to default
par(mfrow = c(1, 1))

     


# Extract fitted values and original data
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_r_fish, predict.x)
plot(depth_values, richness_values, pch = 21, col = mygrey, bg = mygrey,
     main = "Fish Model")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

depth_values <- subsets_taxons[[1]]$depth
richness_values <- subsets_taxons[[1]]$richness


# Create the plot with a custom title

lines(depth_values, fitted_values, col = "blue")
  
     
# Extract fitted values and original data
fitted_values <- predict(model_r_sponges, subsets_taxons[[2]]$depth)
depth_values <- subsets_taxons[[2]]$depth
richness_values <- subsets_taxons[[2]]$richness

# Create the plot with a custom title
plot(depth_values, richness_values, pch = 21, col = mygrey, bg = mygrey,
     main = "sponges Model")
lines(depth_values, fitted_values, col = "blue")     
     
     
     
     
     
     
     
     
     
     
# _____________________________________________________________


     

# JUST FISH - for Taiwan:
par(mfrow = c(2, 2), mar = c(5, 5, 4, 2)) # Adjust margins if needed)

# Fish - mixgaussian_negbin
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_a_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     ylab="Abundance", xlab = "Depth")
lines(predict.x, fitted_values, col = "#228875", lwd = 2)  


# Fish - mixgaussian_zinb
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_r_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, pch = 21, col = mygrey, bg = mygrey,
     ylab="Richness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  


# Fish - mixgaussian_ziig
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_e_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     ylab="Evenness", xlab = "Depth")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)    

# ______________________________________________________

dev.off() 
# Fish - mixgaussian_zinb
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_r_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, 
     pch = 16, col = mygrey, bg = mygrey, cex.lab = 1.5, cex.main = 1.5,
     main = "Fish Richness Across Depths", 
     ylab= expression("Species Richness (N)") , 
     xlab = expression("Depth (m)"))
lines(predict.x, fitted_values, col = "#882247", lwd = 3)  

     
     
     
     
     # doesnt work
     # 
     # 
     # predict_x_a <- purrr::map(taxas, ~ {
     #   seq(from=min(subsets_taxons[[.x]]$depth), to = max(subsets_taxons[[.x]]$depth), length.out = 100000)
     # })
     # 
     # predict.x_a <- seq(from=min(subsets_taxons[[1]]$depth), 
     #                    to = max(subsets_taxons[[1]]$depth), length.out = 100000)   
     # yfit.one_model_a <- predict(one_model_a[[1]], predict.x_a)  
     # lines(predict.x_a, yfit.one_model_a, col = "#00BFC4", lwd = 2) 
     # legend(x = 20, y = 80, legend = "mix.gauss + zinb", 
     #        lwd = 2, col = "#00BFC4", cex = 0.8)
     # 
     # 
     # predict.x_a <- seq(from=min(subsets_taxons[[1]]$depth), 
     #                    to = max(subsets_taxons[[1]]$depth), length.out = 100000)   
     # yfit.one_model_a <- predict(model_a_fish, predict.x)  
     # lines(predict.x, yfit.one_model_a, col = "#00BFC4", lwd = 2) 
     # legend(x = 20, y = 80, legend = "mix.gauss + zinb", 
     #        lwd = 2, col = "#00BFC4", cex = 0.8)



library(ggplot2)
library(dplyr)
library(purrr)

# Create a function to generate prediction data for each taxon
generate_predictions <- function(i, models) {
  predict.x <- seq(from = min(subsets_taxons[[i]]$depth), 
                   to = max(subsets_taxons[[i]]$depth), 
                   length.out = nrow(subsets_taxons[[i]]))
  
  fitted_values <- predict(models[[i]], predict.x)
  
  data.frame(
    depth = predict.x,
    richness = fitted_values,
    taxon = c("Fishes", "Sponges", "Stony Corals", "Soft Corals")[i]
  )
}

# Generate prediction data for all taxa
prediction_data <- map_df(1:4, generate_predictions, 
                          models = list(model_r_fish, model_r_sponges, model_r_stony, model_r_soft))

# Combine original data into a single data frame
original_data <- bind_rows(
  mutate(subsets_taxons[[1]], taxon = "Fishes"),
  mutate(subsets_taxons[[2]], taxon = "Sponges"),
  mutate(subsets_taxons[[3]], taxon = "Stony Corals"),
  mutate(subsets_taxons[[4]], taxon = "Soft Corals")
)

# Create the plot with faceting
ggplot() +
  geom_point(data = original_data, 
             aes(x = depth, y = richness), 
             shape = 21, color = mygrey, fill = mygrey) +
  geom_line(data = prediction_data, 
            aes(x = depth, y = richness), 
            color = "#224188", linewidth = 1) +
  facet_wrap(~ taxon, scales = "free", ncol = 4) +
  labs(y = "Richness", x = "") +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 11, face = "bold"),
    axis.title.y = element_text(face = "bold", size = 11),
    panel.spacing = unit(1, "lines")
  )

     # 
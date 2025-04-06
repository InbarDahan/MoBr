

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
best_models_aic_a <- readRDS("best_models_aic_a.rds")
best_models_aic_r <- readRDS("best_models_aic_r.rds")
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

# _____________________________________________________________


# evenness:

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

# visualization:

# # abundance:
# 
# # 
# predict_x_a <- purrr::map(taxas, ~ {
#   seq(from=min(subsets_taxons[[.x]]$depth), to = max(subsets_taxons[[.x]]$depth), length.out = 100000)
# })
# 
# 
#   
#      
     
     
     
     
     
     
     
     
     
     
     
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
     # 
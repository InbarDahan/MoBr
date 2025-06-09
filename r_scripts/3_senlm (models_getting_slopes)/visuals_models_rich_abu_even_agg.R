

                          # model visualization
                    # creating the plot for the 4 measures

library(dplyr)
library(tidyr)
library(ggplot2)
library(senlm)
library(purrr)

# wd:
wd_models <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed/models" 

setwd(wd_models)
# best model:
one_model_a <- readRDS("one_model_a.rds")
one_model_r <- readRDS("one_model_r.rds")
one_model_e <- readRDS("one_model_e.rds")
one_model_agg <- readRDS("one_model_agg.rds")

# 5 best models:
models_list_r <- readRDS("best_models_r.rds")
models_list_a <- readRDS("best_models_a.rds")
models_list_e <- readRDS("best_models_e.rds")

# data - clean if needed
subsets_taxons <- readRDS("subsets_taxons.rds")
subsets_taxons_even <- readRDS("subsets_taxons_even.rds") # without NA's and negative values
subsets_taxons_agg <- readRDS("subsets_taxons_agg.rds")
# _____________________________________________________________

      # Running the models for all of the diversity measures:

# richness models:
models_list_r[[2]] # - when the first models doesn't look the best (over fitting)

# mixgaussian_negbin     
model_r_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "richness", mean_fun = one_model_r[[1]]$richness$sech_negbin$model_info$mean_fun, err_dist = one_model_r[[1]]$richness$sech_negbin$model_info$err_dist)
# sech_zinb
model_r_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[2]]$richness$sech_zinb$model_info$mean_fun, err_dist = one_model_r[[2]]$richness$sech_zinb$model_info$err_dist)
# hofV_poisson
model_r_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "richness",   mean_fun = one_model_r[[3]]$richness$hofV_negbin$model_info$mean_fun, err_dist = one_model_r[[3]]$richness$hofV_negbin$model_info$err_dist)
#
model_r_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "richness",  mean_fun = one_model_r[[4]]$richness$gaussian_poisson$model_info$mean_fun, err_dist = one_model_r[[4]]$richness$gaussian_poisson$model_info$err_dist)
# - - - - -

# abundance models:

models_list_a[[4]]

# mixgaussian_zinbl
model_a_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[1]]$abundance$mixgaussian_zinbl$model_info$mean_fun, err_dist = one_model_a[[1]]$abundance$mixgaussian_zinbl$model_info$err_dist)
# sech_zinb
model_a_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[2]]$abundance$sech_zinb$model_info$mean_fun, err_dist = one_model_a[[2]]$abundance$sech_zinb$model_info$err_dist)
# gaussian_negbin
model_a_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[3]]$abundance$gaussian_negbin$model_info$mean_fun, err_dist = one_model_a[[3]]$abundance$gaussian_negbin$model_info$err_dist)
# gaussian_negbin
model_a_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[4]]$abundance$constant_negbin$model_info$mean_fun, err_dist = one_model_a[[4]]$abundance$constant_negbin$model_info$err_dist)

# - - - - -

# evenness models:

models_list_e[[2]]

# gaussian_zig
model_e_fish <- senlm(data = subsets_taxons_even[[1]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[1]]$model_info$mean_fun, err_dist = one_model_e[[1]]$model_info$err_dist)
# mixgaussian_ziigl
model_e_sponges <- senlm(data = subsets_taxons_even[[2]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[2]]$model_info$mean_fun, err_dist = one_model_e[[2]]$model_info$err_dist)
# mixgaussian_ziig
model_e_stony <- senlm(data = subsets_taxons_even[[3]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[3]]$model_info$mean_fun, err_dist = one_model_e[[3]]$model_info$err_dist)
# gaussian_tweedie
model_e_soft <- senlm(data = subsets_taxons_even[[4]], xvar = "depth", yvar = "sigma",  mean_fun = one_model_e[[4]]$model_info$mean_fun, err_dist = one_model_e[[4]]$model_info$err_dist)

# - - - - -

# aggregations models:

# gaussian_zig    
model_agg_fish <- senlm(data = subsets_taxons_agg[[1]], xvar = "depth_group", yvar = "value",  mean_fun = one_model_agg[[1]]$model_info$mean_fun, err_dist = one_model_agg[[1]]$model_info$err_dist)
# gaussian_gaussian
model_agg_sponges <- senlm(data = subsets_taxons_agg[[2]], xvar = "depth_group", yvar = "value", mean_fun = one_model_agg[[2]]$model_info$mean_fun, err_dist = one_model_agg[[2]]$model_info$err_dist)
# gaussian_tweedie
model_agg_stony <- senlm(data = subsets_taxons_agg[[3]], xvar = "depth_group", yvar = "value", mean_fun = one_model_agg[[3]]$model_info$mean_fun, err_dist = one_model_agg[[3]]$model_info$err_dist)
# gaussian_tweedie
model_agg_soft <- senlm(data = subsets_taxons_agg[[4]], xvar = "depth_group", yvar = "value", mean_fun = one_model_agg[[4]]$model_info$mean_fun, err_dist = one_model_agg[[4]]$model_info$err_dist)

# _____________________________________________________________

              # Plotting the united plot 

# Reset layout to default
par(mfrow = c(1, 1))


# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(4, 4), mar = c(1, 1, 2, 1.5))  # Adjust margins as needed
mygrey <- grey(level = 0.65, alpha = 0.3)  

# - - - - -

# Fish
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_r_fish, predict.x)   # for this function we need an 'senlm' object and not data frame and that's why we run the models in this script
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="Richness", xlab = "", cex.lab = 1.1, # Size of axis labels
     font.lab = 2) # Font: 1=plain, 2=bold, 3=italic, 4=bold italic
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Sponges
predict.x <- seq(from=min(subsets_taxons[[2]]$depth), to = max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]])) 
fitted_values <- predict(model_r_sponges, predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Sponges", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Stony corals 
predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
fitted_values <- predict(model_r_stony, predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Stony Corals", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#224188", lwd = 2)  

# Soft corals
predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
fitted_values <- predict(model_r_soft, predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Soft Corals", ylab="", xlab = "")
rect(xleft = 4, xright = 34, ybottom = par("usr")[3], ytop = par("usr")[4], 
     col = rgb(0.8, 0.8, 0.8, 0.3), border = NA)
lines(predict.x, fitted_values, col = "#224188", lwd = 2) 

# - - - - -

# abundance

# Fish 
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(model_a_fish, predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="Abundance", xlab = "", cex.lab = 1.1, font.lab = 2)
lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

# Sponges 
predict.x <- seq(from=min(subsets_taxons[[2]]$depth), to = max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]])) 
fitted_values <- predict(model_a_sponges, predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

# Stony corals
predict.x <- seq(from=min(subsets_taxons[[3]]$depth), to = max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]])) 
fitted_values <- predict(model_a_stony, predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#228875", lwd = 2)  

# Soft corals
predict.x <- seq(from=min(subsets_taxons[[4]]$depth), to = max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]])) 
fitted_values <- predict(model_a_soft, predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#228875", lwd = 2) 

# - - - - -

# evenness:

# Fish 
predict.x <- seq(from=min(subsets_taxons_even[[1]]$depth), to = max(subsets_taxons_even[[1]]$depth), length.out = nrow(subsets_taxons_even[[1]])) 
fitted_values <- predict(model_e_fish, predict.x)
plot(subsets_taxons_even[[1]]$depth, subsets_taxons_even[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="Evenness", xlab = "", cex.lab = 1.1, font.lab = 2)
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Sponges 
predict.x <- seq(from=min(subsets_taxons_even[[2]]$depth), to = max(subsets_taxons_even[[2]]$depth), length.out = nrow(subsets_taxons_even[[2]])) 
fitted_values <- predict(model_e_sponges, predict.x)
plot(subsets_taxons_even[[2]]$depth, subsets_taxons_even[[2]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Stony corals
predict.x <- seq(from=min(subsets_taxons_even[[3]]$depth), to = max(subsets_taxons_even[[3]]$depth), length.out = nrow(subsets_taxons_even[[3]])) 
fitted_values <- predict(model_e_stony, predict.x)
plot(subsets_taxons_even[[3]]$depth, subsets_taxons_even[[3]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#882247", lwd = 2)  

# Soft corals
predict.x <- seq(from=min(subsets_taxons_even[[4]]$depth), to = max(subsets_taxons_even[[4]]$depth), length.out = nrow(subsets_taxons_even[[4]])) 
fitted_values <- predict(model_e_soft, predict.x)
plot(subsets_taxons_even[[4]]$depth, subsets_taxons_even[[4]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#882247", lwd = 2) 

# - - - - -

# aggregations

# Fish 
predict.x <- seq(from=min(subsets_taxons_agg[[1]]$depth_group), to = max(subsets_taxons_agg[[1]]$depth_group), length.out = nrow(subsets_taxons_agg[[1]])) 
fitted_values <- predict(one_model_agg[[1]], predict.x)
plot(subsets_taxons_agg[[1]]$depth_group, subsets_taxons_agg[[1]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="aggregations", xlab = "", cex.lab = 1.1, font.lab = 2)
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)   

# Sponges
predict.x <- seq(from=min(subsets_taxons_agg[[2]]$depth_group), to = max(subsets_taxons_agg[[2]]$depth_group), length.out = nrow(subsets_taxons_agg[[2]])) 
fitted_values <- predict(one_model_agg[[2]], predict.x)
plot(subsets_taxons_agg[[2]]$depth_group, subsets_taxons_agg[[2]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)  

# Stony corals
predict.x <- seq(from=min(subsets_taxons_agg[[3]]$depth_group), to = max(subsets_taxons_agg[[3]]$depth_group), length.out = nrow(subsets_taxons_agg[[3]])) 
fitted_values <- predict(model_agg_stony, predict.x)
plot(subsets_taxons_agg[[3]]$depth_group, subsets_taxons_agg[[3]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2)  

# Soft corals 
predict.x <- seq(from=min(subsets_taxons_agg[[4]]$depth_group), to = max(subsets_taxons_agg[[4]]$depth_group), length.out = nrow(subsets_taxons_agg[[4]])) 
fitted_values <- predict(model_agg_soft, predict.x)
plot(subsets_taxons_agg[[4]]$depth_group, subsets_taxons_agg[[4]]$value, pch = 21, col = mygrey, bg = mygrey,
     main = "", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#CAFF70", lwd = 2) 

# - - - - -

# Add shared x and y axis labels outside the plotting area
mtext("Depth", side = 1, outer = TRUE, line = 3, cex = 1.2)
mtext("Diversity measure", side = 2, outer = TRUE, line = 3, cex = 1.2)

# Reset layout to default
par(mfrow = c(1, 1))


# ______________________________________________________

# ______________________________________________________


par(mfrow = c(4, 4), oma = c(4, 1, 3, 1))

# Define function to get rounded y-axis breaks
get_y_breaks <- function(values) {
  min_val <- floor(min(values, na.rm = TRUE))
  max_val <- ceiling(max(values, na.rm = TRUE))
  mid_val <- round((min_val + max_val) / 2)
  return(c(min_val, mid_val, max_val))
}

# ROW 1: RICHNESS
# Add row title
mtext("Richness", side = 2, outer = TRUE, line = -0.5, at = 0.875, cex = 1.2, font = 2, las = 1)

# Fishes
x_max <- ceiling(max(subsets_taxons[[1]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[1]]$depth), max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]]))
fitted_values <- predict(model_r_fish, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[1]]$richness)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Fishes", ylab="", xlab = "", font.lab = 2,      
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)), 
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# Sponges
x_max <- ceiling(max(subsets_taxons[[2]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[2]]$depth), max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]]))
fitted_values <- predict(model_r_sponges, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[2]]$richness)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Porifera", ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# Stony corals
x_max <- ceiling(max(subsets_taxons[[3]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[3]]$depth), max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]]))
fitted_values <- predict(model_r_stony, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[3]]$richness)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Scleractinia", ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# Soft corals
x_max <- ceiling(max(subsets_taxons[[4]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[4]]$depth), max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]]))
fitted_values <- predict(model_r_soft, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[4]]$richness)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, pch = 21, col = mygrey, bg = mygrey,
     main = "Octocorallia", ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
# Add shaded region for missing data (4-34m)
rect(4, par("usr")[3], 34, par("usr")[4], col = rgb(1, 0, 0, 0.2), border = NA)
lines(predict.x, fitted_values, col = "#084081", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# ROW 2: ABUNDANCE
# Add row title
mtext("Abundance", side = 2, outer = TRUE, line = -0.5, at = 0.625, cex = 1.2, font = 2, las = 1)

x_max <- ceiling(max(subsets_taxons[[1]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[1]]$depth), max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]]))
fitted_values <- predict(model_a_fish, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[1]]$abundance)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "", font.lab = 2,
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons[[2]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[2]]$depth), max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]]))
fitted_values <- predict(model_a_sponges, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[2]]$abundance)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons[[3]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[3]]$depth), max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]]))
fitted_values <- predict(model_a_stony, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[3]]$abundance)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons[[4]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons[[4]]$depth), max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]]))
fitted_values <- predict(model_a_soft, predict.x)
y_breaks <- get_y_breaks(subsets_taxons[[4]]$abundance)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
# Add shaded region for missing data (4-34m)
rect(4, par("usr")[3], 34, par("usr")[4], col = rgb(1, 0, 0, 0.2), border = NA)
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# ROW 3: EVENNESS
# Add row title
mtext("Evenness", side = 2, outer = TRUE, line = -0.5, at = 0.375, cex = 1.2, font = 2, las = 1)

x_max <- ceiling(max(subsets_taxons_even[[1]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons_even[[1]]$depth), max(subsets_taxons_even[[1]]$depth), length.out = nrow(subsets_taxons_even[[1]]))
fitted_values <- predict(model_e_fish, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_even[[1]]$sigma)
plot(subsets_taxons_even[[1]]$depth, subsets_taxons_even[[1]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "", font.lab = 2,
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_even[[2]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons_even[[2]]$depth), max(subsets_taxons_even[[2]]$depth), length.out = nrow(subsets_taxons_even[[2]]))
fitted_values <- predict(model_e_sponges, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_even[[2]]$sigma)
plot(subsets_taxons_even[[2]]$depth, subsets_taxons_even[[2]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_even[[3]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons_even[[3]]$depth), max(subsets_taxons_even[[3]]$depth), length.out = nrow(subsets_taxons_even[[3]]))
fitted_values <- predict(model_e_stony, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_even[[3]]$sigma)
plot(subsets_taxons_even[[3]]$depth, subsets_taxons_even[[3]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_even[[4]]$depth) / 10) * 10
predict.x <- seq(min(subsets_taxons_even[[4]]$depth), max(subsets_taxons_even[[4]]$depth), length.out = nrow(subsets_taxons_even[[4]]))
fitted_values <- predict(model_e_soft, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_even[[4]]$sigma)
plot(subsets_taxons_even[[4]]$depth, subsets_taxons_even[[4]]$sigma, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n")
# Add shaded region for missing data (4-34m)
rect(4, par("usr")[3], 34, par("usr")[4], col = rgb(1, 0, 0, 0.2), border = NA)
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)
axis(side = 2, at = y_breaks, las = 1)

# ROW 4: AGGREGATIONS
# Add row title
mtext("Aggregations", side = 2, outer = TRUE, line = -0.5, at = 0.125, cex = 1.2, font = 2, las = 1)

x_max <- ceiling(max(subsets_taxons_agg[[1]]$depth_group) / 10) * 10
predict.x <- seq(min(subsets_taxons_agg[[1]]$depth_group), max(subsets_taxons_agg[[1]]$depth_group), length.out = nrow(subsets_taxons_agg[[1]]))
fitted_values <- predict(model_agg_fish, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_agg[[1]]$value)
plot(subsets_taxons_agg[[1]]$depth_group, subsets_taxons_agg[[1]]$value, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "", font.lab = 2,
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", yaxt = "n")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)
axis(side = 1, at = seq(0, x_max, by = 10))
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_agg[[2]]$depth_group) / 10) * 10
predict.x <- seq(min(subsets_taxons_agg[[2]]$depth_group), max(subsets_taxons_agg[[2]]$depth_group), length.out = nrow(subsets_taxons_agg[[2]]))
fitted_values <- predict(model_agg_sponges, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_agg[[2]]$value)
plot(subsets_taxons_agg[[2]]$depth_group, subsets_taxons_agg[[2]]$value, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", yaxt = "n")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)
axis(side = 1, at = seq(0, x_max, by = 10))
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_agg[[3]]$depth_group) / 10) * 10
predict.x <- seq(min(subsets_taxons_agg[[3]]$depth_group), max(subsets_taxons_agg[[3]]$depth_group), length.out = nrow(subsets_taxons_agg[[3]]))
fitted_values <- predict(model_agg_stony, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_agg[[3]]$value)
plot(subsets_taxons_agg[[3]]$depth_group, subsets_taxons_agg[[3]]$value, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", yaxt = "n")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)
axis(side = 1, at = seq(0, x_max, by = 10))
axis(side = 2, at = y_breaks, las = 1)

x_max <- ceiling(max(subsets_taxons_agg[[4]]$depth_group) / 10) * 10
predict.x <- seq(min(subsets_taxons_agg[[4]]$depth_group), max(subsets_taxons_agg[[4]]$depth_group), length.out = nrow(subsets_taxons_agg[[4]]))
fitted_values <- predict(model_agg_soft, predict.x)
y_breaks <- get_y_breaks(subsets_taxons_agg[[4]]$value)
plot(subsets_taxons_agg[[4]]$depth_group, subsets_taxons_agg[[4]]$value, pch = 21, col = mygrey, bg = mygrey,
     ylab="", xlab = "",
     xlim = c(0, x_max), ylim = c(min(y_breaks), max(y_breaks)),
     xaxs = "i", yaxs = "i", yaxt = "n")
# Add shaded region for missing data (4-34m)
rect(4, par("usr")[3], 34, par("usr")[4], col = rgb(1, 0, 0, 0.2), border = NA)
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)
axis(side = 1, at = seq(0, x_max, by = 10))
axis(side = 2, at = y_breaks, las = 1)

# Add shared x-axis label
mtext("Depth (m)", side = 1, outer = TRUE, line = 3, cex = 1.5, font = 2)

# # Add overall title
# mtext("Marine Biodiversity Metrics Across Depth Gradients", side = 3, outer = TRUE, line = 1, cex = 1.8, font = 2)


# ______________________________________________________

# ______________________________________________________

##############################################################################
# doesn't work well, I need to fix it. 
# I am trying to:
# 1) I want the depth to have 10 m breaks, instead of 20
# 2) I want that in each individuals graph, in the y axis, I will have only the lowest value (or zero?), a middle value and the max value and not anything in between
# 3) I want to have the headers in between rows (richness, abundance, evenness and aggregations)
# 4) I want to add the confidence intervals for the model estimate


# ______________________________________________________

# ______________________________________________________



colors_warm <- c("#E31A1C", "#FB6A4A", "#FDAE6B", "#FFD700")

colors_ocean <- c("#084081", "#2B8CBE", "#7BCCC4", "#BAE4BC")

colors_earth <- c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3")

colors_sunset <- c("#762A83", "#9970AB", "#E7A6C7", "#FEE0B6")


"#FB6A4A", "#9970AB", "#E7A6C7", "#FEE0B6")"

c("#1F78B4", "#33A02C", "#E31A1C", "#FDBF6F")




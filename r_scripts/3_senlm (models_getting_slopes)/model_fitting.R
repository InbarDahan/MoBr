# for this script the adapted msenlm function "msenlm_safe is required (r script "msenlm_safe function")


                      # model fitting
                # fitting non linear curves to
     # communities response to environmental gradients

# devtools::install_github("PRIMER-e/senlm")
library(dplyr)
library(tidyr)
library(senlm)
library(purrr) # to skip models crushing
library(mgcv)
# _____________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_models <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed/models" 
# _____________________________________________________________

# read: 
setwd(wd_processed_data)
combined_data <- read.csv("combined_long.csv")
combined_data = combined_data %>% dplyr::select(-X)

agg <- read.csv("agg_morisita_results.csv")
agg <- agg %>% dplyr::select(-X)
# _____________________________________________________________

# check for zero inflation in my data:

# rich + abu + even
zero_inflation_summary <- combined_data %>%
  group_by(taxon, diversity_measure) %>%
  summarize(
    total_values = n(),
    zero_count = sum(value == 0),
    zero_proportion = zero_count / total_values
  ) %>%
  arrange(desc(zero_proportion))

#  - seems like only the sponges data is zero inflated (abundance and richness)
#  - but low counts are present across all taxon

# agg
# rich + abu + even
zero_inflation_agg <- agg %>%
  group_by(taxon) %>%
  summarize(
    total_values = n(),
    zero_count = sum(value == 0),
    zero_proportion = zero_count / total_values
  ) %>%
  arrange(desc(zero_proportion))

# - no zero inflation

# _____________________________________________________________

                    # preparing the data:

# transposing the diversity measure column into 3 independent columns
transformed_data <- combined_data %>% pivot_wider(names_from = diversity_measure, values_from = value)
transformed_data <- as.data.frame(transformed_data) 

# list of all taxa:
# abu + rich + even:
taxas <- transformed_data %>% pull(taxon) %>% unique() %>% as.list %>% purrr::set_names(.)

# ~ loop:
subsets_taxons <- purrr::map(taxas,
                               function(taxa) {
                                 transformed_data %>%
                                   dplyr::filter(taxon == taxa)
                               })
# # save the models output:
# setwd(wd_models)
# saveRDS(subsets_taxons, file = "subsets_taxons.rds")


subsets_taxons_agg <- purrr::map(taxas,
                             function(taxa) {
                               agg %>%
                                 dplyr::filter(taxon == taxa)
                             })

# # save the models output:
# setwd(wd_models)
# saveRDS(subsets_taxons_agg, file = "subsets_taxons_agg.rds")
# _____________________________________________________________

                  # setting models to fit
# models types:
# - - - - - - - - -
# count models - discreet data (abundance + richness) 
count_models <- set_models (mean_class = "main", err_class= "count", method = "crossed") # for abundance and richness
# abundance models - continuous rather then discreet (evenness)
abundance_models <- set_models (mean_class="main", err_class= c("abundance"), method = "crossed") # for sigm
# _____________________________________________________________

              # fitting models - getting the best one

# abundance:

# 1 - define and fit the models:
fitted_models_abu <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = count_models, data = data_subset, xvar = "depth", yvar = "abundance", conf.level=0.95)})

# 2 - Create a list of the top abundance models (based on AIC values) for each taxon:
best_models_a <- purrr::map(taxas, ~ {
               taxa_models <- fitted_models_abu[[.x]]$abundance
               aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
               names(aic_values) <- names(taxa_models)
               ordered_aic_values <- sort(aic_values)
               head(ordered_aic_values)})

# 3 - Extract the best model for each taxon:
one_model_a <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_abu[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))})

# save the models output:
setwd(wd_models)
saveRDS(best_models_a, "best_models_a.rds")
saveRDS(one_model_a, "one_model_a.rds")

# - - - - - - - - - - - - - - 

# richness:

# 1 - define and fit the models:
fitted_models_rich <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = count_models, data = data_subset, xvar = "depth", yvar = "richness", conf.level=0.95)})

# 2 - Create a list of the top richness models (based on AIC values) for each taxon:
best_models_r <- purrr::map(taxas, ~ {
               taxa_models <- fitted_models_rich[[.x]]$richness
               aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
               names(aic_values) <- names(taxa_models)
               ordered_aic_values <- sort(aic_values)
               head(ordered_aic_values)})

# 3 - get the best model
one_model_r <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_rich[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))})

# save the models output:
setwd(wd_models)
saveRDS(best_models_r, "best_models_r.rds")
saveRDS(one_model_r, "one_model_r.rds")
# - - - - - - - - - - - - - -

# evenness:

subsets_taxons_even <- purrr::map(subsets_taxons, function(data_subset) {
  data_subset %>% filter(sigma >= 0 & !is.na(sigma))}) 

# save the models output:
setwd(wd_models)
saveRDS(subsets_taxons_even, file = "subsets_taxons_even.rds")

# 1 - define and fit the models:
fitted_models_even <- purrr::map(subsets_taxons_even, function(data_subset) {
  msenlm_safe(models = abundance_models, data = data_subset, xvar = "depth", yvar = "sigma", conf.level=0.95)})

# 2 - Create a list of the top evenness models (based on AIC values) for each taxon:
best_models_e <- purrr::map(taxas, ~ {
  taxa_models <- fitted_models_even[[.x]]$sigma                                              # Pull all sigma models for current taxon
  valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"])) # Keep only valid models with AIC
  aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"]) # Extract and sort AIC values
  ordered_model_names <- names(sort(aic_values))
  valid_models[ordered_model_names[1:min(5, length(ordered_model_names))]]# Get the top 5 full model objects
})

# 3 - get the best model
one_model_e <- purrr::map(taxas, ~ {
  taxa_models <- fitted_models_even[[.x]]$sigma
  valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"]))
  aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"])
  ordered_model_names <- names(sort(aic_values))
  valid_models[[ordered_model_names[1]]] # Return just the best model (lowest AIC)
})

# save the models output:
setwd(wd_models)
saveRDS(best_models_e, "best_models_e.rds")
saveRDS(one_model_e, "one_model_e.rds")

# - - - - - - - - - - - - - -

# aggregations:

### *** ### 
# Solved problem - the senlm package do not deal with negative data and the j index has negative values. 
# Hence I choose the morisita index, where the lowest value is 0. 
### *** ### 

# 1 - define and fit the models:
fitted_models_agg <- purrr::map(subsets_taxons_agg, function(data_subset) {
  msenlm_safe(models = abundance_models, data = data_subset, xvar = "depth_group", yvar = "value", conf.level=0.95)})

# 2 - Create a list of the top agg models (based on AIC values) for each taxon:
best_models_agg <- purrr::map(taxas, ~ {
  taxa_models <- fitted_models_agg[[.x]]$value                                              # Pull all value models for current taxon
  valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"])) # Keep only valid models with AIC
  aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"]) # Extract and sort AIC values
  ordered_model_names <- names(sort(aic_values))
  valid_models[ordered_model_names[1:min(5, length(ordered_model_names))]]# Get the top 5 full model objects
})

# 3 - get the best model
one_model_agg <- purrr::map(taxas, ~ {
  taxa_models <- fitted_models_agg[[.x]]$value
  valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"]))
  aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"])
  ordered_model_names <- names(sort(aic_values))
  valid_models[[ordered_model_names[1]]] # Return just the best model (lowest AIC)
})
# 
# save the models output:
setwd(wd_models)
saveRDS(best_models_agg, "best_models_agg.rds")
saveRDS(one_model_agg, "one_model_agg.rds")

# _____________________________________________________________


              # Plotting the united plot 

# Reset layout to default
par(mfrow = c(1, 1))


# Set up a multi-panel layout: 4 rows and 2 columns
par(mfrow = c(4, 4), mar = c(1, 1, 2, 1.5))  # Adjust margins as needed
mygrey <- grey(level = 0.65, alpha = 0.4)  

########################
# here I didnt read the models. its because the agg doest work
########################

# - - - - -

# Fish
predict.x <- seq(from=min(subsets_taxons[[1]]$depth), to = max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]])) 
fitted_values <- predict(one_model_r, predict.x)   # for this function we need an 'senlm' object and not data frame and that's why we run the models in this script
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

# _____________________________________________________________




# after changes from Yoni:

# setting parameters for plot layout
par(mfrow = c(4, 4),
    oma = c(5, 5, 3, 1),   # external gaps
    mar = c(2, 2, 2, 1),   # internal gaps
    cex.lab = 1.5,         # title and axis size
    cex.axis = 1.3,        # size of axis numbers
    bty = "n")             # removing the outline

# - - - - - Richness - - - - -

# Fishes
predict.x <- seq(min(subsets_taxons[[1]]$depth), max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]]))
fitted_values <- predict(one_model_r[[1]], predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$richness, pch = 21, col = "black", bg = "black",
     main = "Fishes", ylab="Richness", xlab = "", font.lab = 2)
lines(predict.x, fitted_values, col = "#084081", lwd = 3)

# Sponges
predict.x <- seq(min(subsets_taxons[[2]]$depth), max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]]))
fitted_values <- predict(one_model_r[[2]], predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$richness, pch = 21, col = "black", bg = "black",
     main = "Sponges", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)

# Stony corals
predict.x <- seq(min(subsets_taxons[[3]]$depth), max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]]))
fitted_values <- predict(one_model_r[[3]], predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$richness, pch = 21, col = "black", bg = "black",
     main = "Stony Corals", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)

# Soft corals
predict.x <- seq(min(subsets_taxons[[5]]$depth), max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]]))
fitted_values <- predict(one_model_r[[4]], predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$richness, pch = 21, col = "black", bg = "black",
     main = "Soft Corals", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#084081", lwd = 3)

# - - - - - Abundance - - - - -

predict.x <- seq(min(subsets_taxons[[1]]$depth), max(subsets_taxons[[1]]$depth), length.out = nrow(subsets_taxons[[1]]))
fitted_values <- predict(one_model_a[[1]], predict.x)
plot(subsets_taxons[[1]]$depth, subsets_taxons[[1]]$abundance, pch = 21, col = "black", bg = "black",
     ylab="Abundance", xlab = "", font.lab = 2)
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)

predict.x <- seq(min(subsets_taxons[[2]]$depth), max(subsets_taxons[[2]]$depth), length.out = nrow(subsets_taxons[[2]]))
fitted_values <- predict(one_model_r[[2]], predict.x)
plot(subsets_taxons[[2]]$depth, subsets_taxons[[2]]$abundance, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)

predict.x <- seq(min(subsets_taxons[[3]]$depth), max(subsets_taxons[[3]]$depth), length.out = nrow(subsets_taxons[[3]]))
fitted_values <- predict(one_model_r[[3]], predict.x)
plot(subsets_taxons[[3]]$depth, subsets_taxons[[3]]$abundance, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)

predict.x <- seq(min(subsets_taxons[[4]]$depth), max(subsets_taxons[[4]]$depth), length.out = nrow(subsets_taxons[[4]]))
fitted_values <- predict(one_model_r[[4]], predict.x)
plot(subsets_taxons[[4]]$depth, subsets_taxons[[4]]$abundance, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#2B8CBE", lwd = 3)

# - - - - - Evenness - - - - -

predict.x <- seq(min(subsets_taxons_even[[1]]$depth), max(subsets_taxons_even[[1]]$depth), length.out = nrow(subsets_taxons_even[[1]]))
fitted_values <- predict(one_model_e[[1]], predict.x)
plot(subsets_taxons_even[[1]]$depth, subsets_taxons_even[[1]]$sigma, pch = 21, col = "black", bg = "black",
     ylab="Evenness", xlab = "", font.lab = 2)
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)

predict.x <- seq(min(subsets_taxons_even[[2]]$depth), max(subsets_taxons_even[[2]]$depth), length.out = nrow(subsets_taxons_even[[2]]))
fitted_values <- predict(one_model_e[[2]], predict.x)
plot(subsets_taxons_even[[2]]$depth, subsets_taxons_even[[2]]$sigma, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)

predict.x <- seq(min(subsets_taxons_even[[3]]$depth), max(subsets_taxons_even[[3]]$depth), length.out = nrow(subsets_taxons_even[[3]]))
fitted_values <- predict(one_model_e[[3]], predict.x)
plot(subsets_taxons_even[[3]]$depth, subsets_taxons_even[[3]]$sigma, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)

predict.x <- seq(min(subsets_taxons_even[[4]]$depth), max(subsets_taxons_even[[4]]$depth), length.out = nrow(subsets_taxons_even[[4]]))
fitted_values <- predict(one_model_e[[4]], predict.x)
plot(subsets_taxons_even[[4]]$depth, subsets_taxons_even[[4]]$sigma, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#7BCCC4", lwd = 3)

# - - - - - Aggregations - - - - -

predict.x <- seq(min(subsets_taxons_agg[[1]]$depth_group), max(subsets_taxons_agg[[1]]$depth_group), length.out = nrow(subsets_taxons_agg[[1]]))
fitted_values <- predict(one_model_agg[[1]], predict.x)
plot(subsets_taxons_agg[[1]]$depth_group, subsets_taxons_agg[[1]]$value, pch = 21, col = "black", bg = "black",
     ylab="Aggregations", xlab = "", font.lab = 2)
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)

predict.x <- seq(min(subsets_taxons_agg[[2]]$depth_group), max(subsets_taxons_agg[[2]]$depth_group), length.out = nrow(subsets_taxons_agg[[2]]))
fitted_values <- predict(one_model_agg[[2]], predict.x)
plot(subsets_taxons_agg[[2]]$depth_group, subsets_taxons_agg[[2]]$value, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)

predict.x <- seq(min(subsets_taxons_agg[[3]]$depth_group), max(subsets_taxons_agg[[3]]$depth_group), length.out = nrow(subsets_taxons_agg[[3]]))
fitted_values <- predict(one_model_agg[[3]], predict.x)
plot(subsets_taxons_agg[[3]]$depth_group, subsets_taxons_agg[[3]]$value, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)

predict.x <- seq(min(subsets_taxons_agg[[4]]$depth_group), max(subsets_taxons_agg[[4]]$depth_group), length.out = nrow(subsets_taxons_agg[[4]]))
fitted_values <- predict(one_model_agg[[4]], predict.x)
plot(subsets_taxons_agg[[4]]$depth_group, subsets_taxons_agg[[4]]$value, pch = 21, col = "black", bg = "black", ylab="", xlab = "")
lines(predict.x, fitted_values, col = "#BAE4BC", lwd = 3)

# Add shared x-axis label
mtext("Depth (m)", side = 1, outer = TRUE, line = 3, cex = 1.5, font = 2)

# Add overall title
mtext("Marine Biodiversity Metrics Across Depth Gradients", side = 3, outer = TRUE, line = 1, cex = 1.8, font = 2)

# Reset layout to default
par(mfrow = c(1, 1))









colors_warm <- c("#E31A1C", "#FB6A4A", "#FDAE6B", "#FFD700")

colors_ocean <- c("#084081", "#2B8CBE", "#7BCCC4", "#BAE4BC")

colors_earth <- c("#8C510A", "#BF812D", "#DFC27D", "#F6E8C3")

colors_sunset <- c("#762A83", "#9970AB", "#E7A6C7", "#FEE0B6")


"#FB6A4A", "#9970AB", "#E7A6C7", "#FEE0B6")"

c("#1F78B4", "#33A02C", "#E31A1C", "#FDBF6F")











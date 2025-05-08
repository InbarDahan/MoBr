# for this script the adapted msenlm function "msenlm_safe is required (r script "msenlm_safe function")


                      # model fitting
                # fitting non linear curves to
     # communities response to environmental gradients

# devtools::install_github("PRIMER-e/senlm")
library(dplyr)
library(tidyr)
library(senlm)
library(purrr)
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

agg <- read.csv("agg_j_results.csv")
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
# save the models output:
setwd(wd_models)
saveRDS(subsets_taxons, file = "subsets_taxons.rds")


subsets_taxons_agg <- purrr::map(taxas,
                             function(taxa) {
                               agg %>%
                                 dplyr::filter(taxon == taxa)
                             })

# save the models output:
setwd(wd_models)
saveRDS(subsets_taxons_agg, file = "subsets_taxons_agg.rds")
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

# 1
Pars_c <- create_default_par_list (count_models)

# 2 - define and fit the models:
fitted_models_abu <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = count_models, data = data_subset, xvar = "depth", yvar = "abundance", conf.level=0.95)})

# 3 - Create a list of the top abundance models (based on AIC values) for each taxon:
best_models_a <- purrr::map(taxas, ~ {
               taxa_models <- fitted_models_abu[[.x]]$abundance
               aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
               names(aic_values) <- names(taxa_models)
               ordered_aic_values <- sort(aic_values)
               head(ordered_aic_values)})

# 4 - Extract the best model for each taxon:
one_model_a <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_abu[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))})

# save the models output:
setwd(wd_models)
write.csv(best_models_a, file = "best_models_a.csv")
write.csv(one_model_a, file = "one_model_a.csv")

# - - - - - - - - - - - - - -

# richness:

# 1
Pars_c_2 <- create_default_par_list (count_models)

# 2 - define and fit the models:
fitted_models_rich <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = count_models, data = data_subset, xvar = "depth", yvar = "richness", conf.level=0.95)})

# 3 - Create a list of the top richness models (based on AIC values) for each taxon:
best_models_r <- purrr::map(taxas, ~ {
               taxa_models <- fitted_models_rich[[.x]]$richness
               aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
               names(aic_values) <- names(taxa_models)
               ordered_aic_values <- sort(aic_values)
               head(ordered_aic_values)})

# 4 - get the best model
one_model_r <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_rich[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))})


# save the models output:
setwd(wd_models)
write.csv(best_models_r, file = "best_models_r.csv")
write.csv(one_model_r, file = "one_model_r.csv")

# - - - - - - - - - - - - - -

# evenness:

subsets_taxons_clean <- purrr::map(subsets_taxons, function(data_subset) {
  data_subset %>% filter(sigma >= 0 & !is.na(sigma))}) 

# 1
Pars_a <- create_default_par_list (abundance_models)

# 2 - define and fit the models:
fitted_models_even <- purrr::map(subsets_taxons_clean, function(data_subset) {
  msenlm_safe(models = abundance_models, data = data_subset, xvar = "depth", yvar = "sigma", conf.level=0.95)})

# 3 - Create a list of the top evenness models (based on AIC values) for each taxon:
best_models_e <- purrr::map(taxas, ~ {
  taxa_models <- fitted_models_even[[.x]]$sigma                                              # Pull all sigma models for current taxon
  valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"])) # Keep only valid models with AIC
  aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"]) # Extract and sort AIC values
  ordered_model_names <- names(sort(aic_values))
  valid_models[ordered_model_names[1:min(5, length(ordered_model_names))]]# Get the top 5 full model objects
})

# 4 - get the best model
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
# Problem - the senlm package do not deal with negative data and the j index has negative values. 
# if I will is k or morisita than the regular code will work but for the j, an alternative
# modelling approach is needed
### *** ### 

fitted_models_agg <- purrr::map(subsets_taxons_agg, function(data_subset) {
  tryCatch(gam(value ~ s(depth_group), data = data_subset, family = gaussian()),
           error = function(e) NULL )})
# 
# abundance_models_agg <- set_models(mean_fun = c("gaussian"), err_dist = c("gaussian"))
# 
# # 2 - define and fit the models:
# fitted_models_agg <- purrr::map(subsets_taxons_agg, function(data_subset) {
#   msenlm_safe(models = abundance_models, data = data_subset, xvar = "depth", yvar = "value", conf.level=0.95)})
# 
# # 3 - Create a list of the top agg models (based on AIC values) for each taxon:
# best_models_agg <- purrr::map(taxas, ~ {
#   taxa_models <- fitted_models_agg[[.x]]$value                                              # Pull all value models for current taxon
#   valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"])) # Keep only valid models with AIC
#   aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"]) # Extract and sort AIC values
#   ordered_model_names <- names(sort(aic_values))
#   valid_models[ordered_model_names[1:min(5, length(ordered_model_names))]]# Get the top 5 full model objects
# })
# 
# # 4 - get the best model
# one_model_agg <- purrr::map(taxas, ~ {
#   taxa_models <- fitted_models_agg[[.x]]$value
#   valid_models <- purrr::keep(taxa_models, ~ !is.null(.x$IC["AIC"]) && !is.na(.x$IC["AIC"]))
#   aic_values <- purrr::map_dbl(valid_models, ~ .x$IC["AIC"])
#   ordered_model_names <- names(sort(aic_values))
#   valid_models[[ordered_model_names[1]]] # Return just the best model (lowest AIC)
# })
# 
# # save the models output:
# setwd(wd_models)
# saveRDS(best_models_agg, "best_models_agg.rds")
# saveRDS(one_model_agg, "one_model_agg.rds")




                      # model fitting
                # fitting non linear curves to
     # communities response to environmental gradients

# devtools::install_github("PRIMER-e/senlm")
library(dplyr)
library(tidyr)
library(senlm)
# _____________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
# _____________________________________________________________

# read: 
setwd(wd_processed_data)
combined_data <- read.csv("combined_long.csv")
combined_data = combined_data %>% dplyr::select(-X)
# _____________________________________________________________

# check for zero inflation in my data:
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
# _____________________________________________________________

                    # preparing the data:

# transposing the diversity measure column into 3 independent columns
transformed_data <- combined_data %>% pivot_wider(names_from = diversity_measure, values_from = value)
transformed_data <- as.data.frame(transformed_data) 

# list of all taxa:
taxas <- transformed_data %>% pull(taxon) %>% unique() %>% as.list %>% purrr::set_names(.)

# ~ loop:
subsets_taxons <- purrr::map(taxas,
                               function(taxa) {
                                 transformed_data %>%
                                   dplyr::filter(taxon == taxa)
                               })
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
               head(ordered_aic_values)
})

# save the models output:
setwd(wd_processed_data)

write.csv(best_models_a, file = "best_models_a.csv")
# 4 - Extract the best model for each taxon:
model_summaries_a <- purrr::map(taxas, ~ {
                   taxa_data <- fitted_models_abu[[.x]]
                   summary(msenlm.best(taxa_data, best="AICc"))
})

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
               head(ordered_aic_values)
})


# save the models output:
setwd(wd_processed_data)
write.csv(best_models_r, file = "best_models_r.csv")

# 4 - Extract the best model for each taxon:
model_summaries_r <- purrr::map(taxas, ~ {
                   taxa_data <- fitted_models_rich[[.x]]
                   summary(msenlm.best(taxa_data, best="AICc"))
})

# - - - - - - - - - - - - - -

# evenness:

# 1
Pars_a <- create_default_par_list (abundance_models)

# 2 - define and fit the models:
fitted_models_even <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = abundance_models, data = data_subset, xvar = "depth", yvar = "sigma", conf.level=0.95)})

# 3 - Create a list of the top evenness models (based on AIC values) for each taxon:
best_models_e <- purrr::map(taxas, ~ {
                 taxa_models <- fitted_models_even[[.x]]$sigma
                 aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
                 names(aic_values) <- names(taxa_models)
                 ordered_aic_values <- sort(aic_values)
                 head(ordered_aic_values)
})

# 4 - Extract the best model for each taxon:
model_summaries_e <- purrr::map(taxas, ~ {
                     taxa_data <- fitted_models_even[[.x]]
                     summary(msenlm.best(taxa_data, best="AICc"))
})




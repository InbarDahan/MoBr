

                   # model fitting
           # fitting non linear curves to
    # communities response to environmental gradients

# devtools::install_github("PRIMER-e/senlm")
library(dplyr)
library(tidyr)
library(senlm)
library(purrr)
# _____________________________________________________________

# wd:
wd_processed_data <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed" 
wd_models <- "C:/Users/inbar/OneDrive/desktop/r/chapter_2/MoBr/data/processed/models" 

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


# 3 - top models Aic list (no prameters):
best_models_aic_a <- purrr::map(taxas, ~ {     # apply on all the taxas in the list
  taxa_models <- fitted_models_abu[[.x]]$abundance   # pull the fitted models
  aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])  # pull their AIC values
  names(aic_values) <- names(taxa_models)  # pull thier names
  ordered_aic_values <- sort(aic_values)   # order from low to his values
  head(ordered_aic_values)                 # pull the head
})

# 4 - summary top models by ordered AIC (with parameters):
best_models_parameters_a <- purrr::map(taxas, ~ {     # apply on all the taxas in the list
  multi <- summary(fitted_models_abu[[1]])  # pull the fitted models
  multi[order(multi$AICc),c(4,5,6,24:28)]
})

# 5 - get the parameters of the best model:
one_model_a <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_abu[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))
})

# 6 - fit the best model for visualization:
fit_a_fish <- senlm(data = subsets_taxons[[1]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[1]]$mean_fun, err_dist = one_model_a[[1]]$err_dist)
fit_a_sponges <- senlm(data = subsets_taxons[[2]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[2]]$mean_fun, err_dist = one_model_a[[2]]$err_dist)
fit_a_stony <- senlm(data = subsets_taxons[[3]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[3]]$mean_fun, err_dist = one_model_a[[3]]$err_dist)
fit_a_soft <- senlm(data = subsets_taxons[[4]], xvar = "depth", yvar = "abundance",  mean_fun = one_model_a[[4]]$mean_fun, err_dist = one_model_a[[4]]$err_dist)


# save the models output:
setwd(wd_models)
write.csv(best_models_parameters_a, file = "best_models_parameters_a.csv")
write.csv(one_model_a, file = "one_model_a.csv")

# - - - - - - - - - - - - - -

# richness:

# 1
Pars_c <- create_default_par_list(count_models)

# 2 - define and fit the models:
fitted_models_rich <- purrr::map(subsets_taxons, function(data_subset) {
  msenlm(models = count_models, data = data_subset, xvar = "depth", yvar = "richness", conf.level=0.95)})


# 3 - top models Aic list (no prameters):
best_models_aic_r <- purrr::map(taxas, ~ {     # apply on all the taxas in the list
  taxa_models <- fitted_models_rich[[.x]]$richness   # pull the fitted models
  aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])  # pull their AIC values
  names(aic_values) <- names(taxa_models)  # pull thier names
  ordered_aic_values <- sort(aic_values)   # order from low to his values
  head(ordered_aic_values)                 # pull the head
})

# 4 - summary top models by ordered AIC (with parameters)
best_models_parameters_r <- purrr::map(taxas, ~ {     # apply on all the taxas in the list
  multi <- summary(fitted_models_rich[[1]])  # pull the fitted models
  multi[order(multi$AICc),c(4,5,6,24:28)]
})


# 5 - get the parameters of the best model
one_model_r <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_rich[[.x]]
  summary(msenlm.best(taxa_data, best="AICc"))
})


# save the models output:
setwd(wd_models)
write.csv(best_models_parameters_r, file = "best_models_parameters_r.csv")
write.csv(one_model_r, file = "one_model_r.csv")

# - - - - - - - - - - - - - -

# evenness:

####### DOESNT WORK - FIX

# filter out invalid rows where sigma is negative or NA:
subsets_taxons_clean <- purrr::map(subsets_taxons, function(data_subset) {
  data_subset %>% filter(sigma >= 0 & !is.na(sigma))}) 

# 1
Pars_a <- create_default_par_list (abundance_models)

# 2 - define and fit the models:
fitted_models_even <- purrr::map(subsets_taxons_clean, function(data_subset) {
  msenlm(models = abundance_models, data = data_subset, xvar = "depth", yvar = "sigma", conf.level=0.95)})

# 3 - Create a list of the top evenness models (based on AIC values) for each taxon:
best_models_e <- purrr::map(taxas, ~ {    
  taxa_models <- fitted_models_even[[.x]]$sigma
  aic_values <- purrr::map_dbl(taxa_models, ~.x$IC["AIC"])
  names(aic_values) <- names(taxa_models)
  ordered_aic_values <- sort(aic_values)
  head(ordered_aic_values)
})

setwd(wd_models)

# 4 - Extract the best model for each taxon:
model_summaries_e <- purrr::map(taxas, ~ {
  taxa_data <- fitted_models_even[[.x]]
  summary(msenlm.best(subsets_taxons_clean, best="AICc"))
})

# - - - - - - - - - - - - - -

###########################################################################
                        # solutions 
                #  evenness models doesnt work

# TRY sigma * 100 AND RUN AGAIN (to check if the problem is caused by low values)
# filter out invalid rows where sigma is negative or NA + multiply 100:
subsets_taxons_clean_100 <- purrr::map(subsets_taxons_clean, function(data_subset) {
  data_subset %>% filter(sigma >= 0 & !is.na(sigma)) %>%
  mutate(sigma = as.numeric(sigma) * 100)
  }) 

# run just for the fish data
fitted_models_even_fish <- msenlm(models = abundance_models, data = subsets_taxons_clean_100$fish, xvar = "depth", yvar = "sigma", conf.level=0.95)
fitted_models_even_sponges <- msenlm(models = abundance_models, data = subsets_taxons_clean_100$sponges, xvar = "depth", yvar = "sigma", conf.level=0.95)

# - - - - - - - - - - - - - -

# from Dor - loop that include jumping over models that do not run:
## ---> it runs but the output is an empty list - like it doesn't save the models.
fitted_models_even <- purrr::map(
  subsets_taxons_clean,
  function(data_subset) {
    tryCatch({
      msenlm(
        models = abundance_models,
        data = data_subset,
        xvar = "depth",
        yvar = "sigma",
        conf.level = 0.95
      )
    }, error = function(e) {
      message("Error encountered: ", e$message)
      return(NULL) # Return NULL or an appropriate fallback value
    })
  })

# - - - - - - - - - - - - - -

# Trying to run the models on individual data sets:
# fish clean data:
fitted_models_even_c_fish <- msenlm(models = abundance_models, data = subsets_taxons_clean$fish, xvar = "depth", yvar = "sigma", conf.level=0.95)
# regular data
fitted_models_even_fish <- msenlm(models = abundance_models, data = subsets_taxons$fish, xvar = "depth", yvar = "sigma", conf.level=0.95)
# fish clean data:
fitted_models_even_c_sponges <- msenlm(models = abundance_models, data = subsets_taxons_clean$sponges, xvar = "depth", yvar = "sigma", conf.level=0.95)
# regular data
fitted_models_even_sponges <- msenlm(models = abundance_models, data = subsets_taxons$sponges, xvar = "depth", yvar = "sigma", conf.level=0.95)

# - - - - - - - - - - - - - -

# EDU - data check
class(abundance_models)
class(subsets_taxons_clean)
length(subsets_taxons_clean)

a <- 10
b <- 3

total <- sum(a,b)

total <- try(sum(a,c), silent = T)


fitted_models_even <- purrr::map(subsets_taxons_clean, function(data_subset) {
  
  # Loop over each model combination in abundance_models
  purrr::map(abundance_models, function(model) {
    
    safely_msenlm <- purrr::safely(function(data) {
      
      msenlm(models = model, data = data, xvar = "depth", yvar = "sigma", conf.level = 0.95)
    })
    
    # Run the model with error handling
    result <- safely_msenlm(data_subset)
    
    # If the model fails, return NULL (or any other default value)
    if (!is.null(result$error)) {
      return(NULL)  # Replace NULL with a default value if desired
    }
    
    return(result$result)  # Return the fitted model result if successful
  })
})

###########################################################################
# DOESNT WORK BECAUSE ABUNDANCE MODELA IS A TABLE.

fitted_models_even <- list()  # Initialize an empty list to store results

for (i in 1:length(subsets_taxons_clean)) {
  data_subset <- subsets_taxons_clean[[i]]  # Get the current data subset for the taxon
  taxon_results <- list()  # List to store results for this particular data subset
  
  # for (j in 1:length(abundance_models)) {
    # model <- abundance_models[[j]]  # Get the current model
  for (j in 1:length(abundance_models)) {
    print(paste(abundance_models[[j, 1]], abundance_models[[j, 2]], abundance_models[[j, 3]]))
    
    # if (abundance_models[[j, 1]] == "" & abundance_models[[j, 2]] == "" &
    #     abundance_models[[j, 3]] == "") {
    #   next
    #}
    
    # Try to fit the model and capture any errors
    tryCatch({
      model_result <- msenlm(models = modeli, data = data_subset, xvar = "depth", yvar = "sigma", conf.level = 0.95)
      
      # Store the result if the model fitting is successful
      taxon_results[[j]] <- model_result
    }, error = function(e) {
      # If an error occurs, print a message and skip this model
      message(paste("Model failed for taxon ", i, " with model ", j, ": ", e$message))
      taxon_results[[j]] <- NULL  # Store NULL for the failed model
    })
  }
  
  # Store the results for the current data subset (taxon)
  fitted_models_even[[i]] <- taxon_results
}






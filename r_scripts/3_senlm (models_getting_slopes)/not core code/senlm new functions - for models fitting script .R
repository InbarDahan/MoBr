
# msenlm new function, that can skip unworking models:
msenlm_safe <- function(models = NULL, data = NULL, xvar = NULL, yvar = NULL, 
                        method = "crossed", conf.level = 0.95) 
{
  if (is.null(xvar) | is.null(yvar)) {
    stop("Must specify xvar and yvar!")
  }
  if (!is.character(xvar)) {
    stop("xvar must be character!")
  }
  if (length(xvar) > 1) {
    stop("Too many x variable names specified!")
  }
  if (any(is.na(match(xvar, names(data))))) {
    stop("xvar not valid!")
  }
  if (!is.character(yvar)) {
    stop("yvar must be character!")
  }
  if (any(is.na(match(yvar, names(data))))) {
    stop("Some yvar not valid!")
  }
  if ((method == "paired") & (nrow(models) != length(yvar))) {
    stop("Length of yvar and number of models must match if method='paired'!")
  }
  x <- data[, xvar]
  xname <- xvar
  y <- data[, yvar]
  yname <- yvar
  Fits <- vector(mode = "list", length = length(yvar))
  names(Fits) <- yvar
  if (method == "paired") {
    ModelNames <- rep(NA, length = nrow(models))
  }
  for (i in 1:length(Fits)) {
    if (method == "paired") {
      ModelFits <- vector(mode = "list", length = 1)
      ModelFits[[1]] <- tryCatch({
        senlm(model = models[i, ], data = data, 
              xvar = xvar, yvar = yvar[i], conf.level = conf.level)
      }, error = function(e) {
        warning(paste("Error in model", i, ":", e$message))
        NULL
      })
      names(ModelFits) <- if (!is.null(ModelFits[[1]])) ModelFits[[1]]$model else "Error"
    }
    if (method == "crossed") {
      ModelFits <- vector(mode = "list", length = nrow(models))
      ModelNames <- rep(NA, length = nrow(models))
      for (j in 1:length(ModelFits)) {
        Fit <- tryCatch({
          senlm(model = models[j, ], data = data, 
                xvar = xvar, yvar = yvar[i], conf.level = conf.level)
        }, error = function(e) {
          warning(paste("Error in model", j, ":", e$message))
          NULL
        })
        ModelFits[[j]] <- Fit
        ModelNames[j] <- if (!is.null(Fit)) Fit$model else paste("Error_", j, sep = "")
      }
      names(ModelFits) <- ModelNames
    }
    Fits[[i]] <- ModelFits
  }
  class(Fits) <- "msenlm"
  return(Fits)
}


# new msenlm_best function - that can fit models that didnt work in the list
mselnm_best_safe <- function(object, best = "AICc") {
  if (!is.null(best)) {
    GOF <- c("nll", "AIC", "AICc", "BIC")
    if (all(best != GOF)) {
      stop('best option must be equal to "nll", "AIC", "AICc", "BIC"!')
    }
  }
  
  BFits <- vector(mode = "list", length = length(object))
  names(BFits) <- names(object)
  
  for (i in seq_along(object)) {
    # Skip if current object is not a list or is NULL
    if (is.null(object[[i]]) || !is.list(object[[i]])) {
      warning(paste("Invalid object:", names(object)[i]))
      next
    }
    
    NModels <- length(object[[i]])
    ICMat <- as.data.frame(matrix(NA, ncol = 4, nrow = NModels))
    names(ICMat) <- c("nll", "AIC", "AICc", "BIC")
    
    for (j in seq_along(object[[i]])) {
      model <- object[[i]][[j]]
      
      # Check if model is a valid model object with IC and convergence
      if (!is.null(model) && 
          is.list(model) && 
          all(c("IC", "convergence", "model") %in% names(model)) &&
          model$convergence == 0) {
        ICMat[j, ] <- model$IC[grep("npar", names(model$IC), invert = TRUE)]
      }
    }
    
    # Remove NA values before finding minimum
    GOF <- ICMat[, best]
    valid_indices <- which(!is.na(GOF))
    
    if (length(valid_indices) > 0) {
      best_model_index <- valid_indices[which.min(GOF[valid_indices])]
      
      BFits[[i]][[1]] <- object[[i]][[best_model_index]]
      names(BFits[[i]]) <- BFits[[i]][[1]]$model
    } else {
      warning(paste("No valid models for object:", names(object)[i]))
      BFits[[i]] <- NULL
    }
  }
  
  return(BFits)
}


# mmil_model.R
#
# Methods for the S3 object "mmil_model"

#' Construct an MMIL model. 
#' 
#' @param model TO DO 
#' @param rho A scalar value between 0 and 1 indicating the probability that
#' that a cell is NOT associated with the outcome-of-interest, given that the
#' cell was sampled from a patient who DOES have the outcome of interest. In other
#' words, the expected number of non-disease-associated cells in each sick patient. 
#' @param zeta The probability that the inherited patient label of any given cell will be 1. 
#' @param num_cells TO DO 
#' @param hyperparameters TO DO 
#' @param num_cells_disease TO DO 
#' @param lls TO DO 
#'
#' @return An "mmil_model", an S3 object subclass of the original model class. 
#' Attributes include rho, zeta, num_cells, num_cells_disease, and hyperparameters. 
#'
#' @family mmil_model utilities
#'
new_mmil_model <- 
  function(
    model, 
    rho, 
    zeta, 
    num_cells, 
    num_cells_disease, 
    lls, 
    hyperparameters = list()
  ) {
    
    # check arguments
    if (!inherits(model, what = c("glmnet", "nnet", "glm", "model_fit"))) { 
      stop("model must be a glmnet, glm, nnet, or model_fit object.")
    }
    
    if (!all(c(length(rho) == 1, length(zeta) == 1))) { 
      stop("rho and zeta must be of length 1.")
    }
    
    if (!all(c(is.double(rho), is.double(zeta)))) { 
      stop("rho and zeta must be doubles.")
    }
    
    if (!all(c(length(num_cells) == 1, length(num_cells_disease) == 1))) { 
      stop("num_cells and num_cells_disease must be of length 1.")
    }
    
    if (!all(c(is.numeric(num_cells), is.numeric(num_cells_disease)))) { 
      stop("num_cells and num_cells_disease must be integers")
    }
    
    if (inherits(model, what = "glmnet")) { 
      type <- "glmnet"
    } else if (inherits(model, "nnet")) { 
      type <- "mlp"
    } else if (inherits(model, "glm")) { 
      type <- "glm"
    } else { 
      type <- "parsnip"
    }
    
    check_mmil_model_hyperparameters(
      type = type, 
      hyperparameters = hyperparameters
    )
    
    # assemble object if all arguments are okay
    mmil_model <- 
      structure(
        model, 
        rho = rho, 
        zeta = zeta, 
        num_cells = num_cells, 
        num_cells_disease = num_cells_disease, 
        lls = lls, 
        hyperparameters = hyperparameters, 
        class = c("mmil_model", class(model))
      )
    
    # return result
    return(mmil_model)
  }

validate_mmil_model <- function() { 
  NULL
} 

# -------------------------------------------------------------------------



#' Construct an MMIL model.
#' 
#' @param model TO DO 
#' @param rho A scalar value between 0 and 1 indicating the probability that
#' that a cell is NOT associated with the outcome-of-interest, given that the
#' cell was sampled from a patient who DOES have the outcome of interest. In other
#' words, the expected number of non-disease-associated cells in each sick patient. 
#' @param zeta The probability that the inherited patient label of any given cell will be 1. 
#' @param num_cells TO DO 
#' @param num_cells_disease TO DO 
#' @param hyperparameters TO DO 
#' @param lls TO DO 
#'
#' @return TO DO 
#' 
#' @export
#'
#' @family mmil_model utilities
#'
mmil_model <- 
  function(
    model, 
    rho = 0.5, 
    zeta = 0.5, 
    num_cells, 
    num_cells_disease, 
    lls = 0, 
    hyperparameters = list()
  ) {  
    result <- 
      new_mmil_model(
        model = model, 
        rho = rho, 
        zeta = zeta, 
        num_cells = num_cells, 
        num_cells_disease = num_cells_disease,
        lls = lls, 
        hyperparameters = hyperparameters
      )
    
    return(result)
    
  }

check_mmil_model_hyperparameters <- function(type, hyperparameters) { 
  if (type == "glmnet") { 
    if (!all(names(hyperparameters) %in% c("lambda", "alpha"))) { 
      stop("both lambda and alpha hyperparameters must be specified for a glmnet model")
    }
    if (!all(c(is.double(hyperparameters$lambda), is.double(hyperparameters$alpha)))) { 
      stop("lambda and alpha must be doubles.")
    }
  } else if (type == "mlp") { 
    if (!all(names(hyperparameters) %in% c("num_neurons", "decay"))) { 
      stop("both num_neurons and decay hyperparameters must be specified for a nnet model")
    }
    if (!all(
      c(is.double(hyperparameters$num_neurons), 
        is.double(hyperparameters$decay))
    )) { 
      stop("num_neurons must be an integer and decay must be a double.")
    }
  } 
  else {
    # placeholder to test tidymodels params
    NULL
  }

}

#' Make predictions for new data using an MMIL model. 
#'
#' @param object An mmil_model object.
#' @param ... Additional parameters to pass to the object subclass. 
#'
#' @return Predictions based on the subclass predict method
#' 
#' @export
#'
#' @examples
#' NULL
#' 
predict.mmil_model <- function(object, ...) { 
  base_predictions <- NextMethod(...)
  return(base_predictions)
}


#' Find an MMIL model's coefficients 
#'
#' @param object An MMIL model. 
#' 
#' @param ... Additional arguments to pass to the subclass-specific `coef` function. 
#'
#' @return TODO - Timothy
#' 
#' @export
#'
#' @examples
#' NULL
#' 
coef.mmil_model <- 
  function(object, ...) { 
    coefficients <- NextMethod(...)
    return(coefficients)
  }

#' Get the hyperparameters of an MMIL object
#'
#' @param mmil_model_object An MMIL model. 
#'
#' @return A list of model type-specific hyperparameters
#' 
#' @export
#'
#' @examples
#' NULL
mmil_get_hyperparameters <- 
  function(mmil_model_object) { 
    result <- attr(x = mmil_model_object, which = "hyperparameters")
    return(result)
  }


#' Get the training log-likelihood values of an MMIL object
#'
#' @param mmil_model_object An MMIL model. 
#'
#' @return A vector of log-likelihood values (one for each iteration)
#' 
#' @export
#'
#' @examples
#' NULL
mmil_get_log_likelihoods <- 
  function(mmil_model_object) { 
    result <- attr(x = mmil_model_object, which = "lls")
    return(result)
  }







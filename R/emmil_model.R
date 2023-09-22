# emmil_model.R
#
# Methods for the S3 object "emmil_model"

#' Construct an EMMIL model. 
#' 
#' @param model TO DO 
#' @param rho TO DO 
#' @param zeta TO DO 
#' @param num_cells TO DO 
#' @param hyperparameters TO DO 
#' @param num_cells_disease TO DO 
#' @param lls TO DO 
#'
#' @return An "emmil_model", an S3 object subclass of the original model class. 
#' Attributes include rho, zeta, num_cells, num_cells_disease, and hyperparameters. 
#'
#' @family emmil_model utilities
#'
new_emmil_model <- 
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
    
    check_emmil_model_hyperparameters(
      type = type, 
      hyperparameters = hyperparameters
    )
    
    # assemble object if all arguments are okay
    emmil_model <- 
      structure(
        model, 
        rho = rho, 
        zeta = zeta, 
        num_cells = num_cells, 
        num_cells_disease = num_cells_disease, 
        lls = lls, 
        hyperparameters = hyperparameters, 
        class = c("emmil_model", class(model))
      )
    
    # return result
    
    return(emmil_model)
  }

validate_emmil_model <- function() { 
  NULL
} 

# -------------------------------------------------------------------------



#' Construct an EMMIL model.
#' 
#' @param model TO DO 
#' @param rho TO DO 
#' @param zeta TO DO 
#' @param num_cells TO DO 
#' @param num_cells_disease TO DO 
#' @param hyperparameters TO DO 
#' @param lls TO DO 
#'
#' @return TO DO 
#' 
#' @export
#'
#' @family emmil_model utilities
#'
emmil_model <- 
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
      new_emmil_model(
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

check_emmil_model_hyperparameters <- function(type, hyperparameters) { 
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

#' Make predictions for new data using an EMMIL model. 
#'
#' @param object TO DO 
#' @param adjust TO DO 
#' @param ... TO DO 
#'
#' @return TO DO 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
predict.emmil_model <- 
  function(object, adjust = FALSE, ...) { 
    base_predictions <- NextMethod(...)
    if (adjust) { 
      z_placeholder <- rep(1, times = length(base_predictions))
      updated_predictions <- 
        update_y(
          predictions = base_predictions, 
          z = z_placeholder, 
          rho = attr(object, which = "rho"), 
          zeta = attr(object, which = "zeta")
        )
      return(updated_predictions)
    } else {
      return(base_predictions)
    }
  }


#' Find an EMMIL model's coefficients 
#'
#' @param object TO DO 
#' @param adjust TO DO 
#' @param ... TO DO 
#'
#' @return TO DO 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
coef.emmil_model <- 
  function(object, ...) { 
    coefficients <- NextMethod(...)
    return(coefficients)
  }

#' Get the hyperparameters of an EMMIL object
#'
#' @param emmil_model_object An EMMIL model. 
#'
#' @return A list of model type-specific hyperparameters
#' 
#' @export
#'
#' @examples
#' NULL
emmil_get_hyperparameters <- 
  function(emmil_model_object) { 
    result <- attr(x = emmil_model_object, which = "hyperparameters")
    return(result)
  }


#' Get the training log-likelihood values of an EMMIL object
#'
#' @param emmil_model_object An EMMIL model. 
#'
#' @return A vector of log-likelihood values (one for each iteration)
#' 
#' @export
#'
#' @examples
#' NULL
emmil_get_log_likelihoods <- 
  function(emmil_model_object) { 
    result <- attr(x = emmil_model_object, which = "lls")
    return(result)
  }

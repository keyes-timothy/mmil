# core_functions.R 
#
# This script contains functions for initializing, training, and visualizing 
# EMMIL models. 

# simulate data ----------------------------------------------------------------

#' Simulate data to demonstrate how to use the EMMIL model. 
#' 
#' This function simulates a rectangular expression matrix (X), an integer 
#' vector of observed, inherited patient labels (z) and an integer vector of
#' unobserved, ground truth cell labels (y) for an imaginary dataset. The 
#' expression matrix X is built such that each row represents a cell and each 
#' column represents a biological feature (i.e. gene or protein expression). All 
#' values in X are drawn from a normal Gaussian distribution, and expression vectors of 
#' cells for which y = 1 are shifted from those of cells for which y = 0 
#' using a user-specified numeric value. 
#'
#' @param num_cells An integer number of cells to simulate
#' 
#' @param num_features An integer number of features to simulate
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. This is the probability 
#' that a random person from your population-of-interest will have the condition-of-interest. 
#' 
#' @param shift An integer indicating the magnitude of the shift applied to cells 
#' whose ground truth labels are y = 1. The value of shift is added to the first 
#' `round(num_features / 2)` values of the feature vectors for all cells whose 
#' ground truth label is 1. 
#'
#' @return A list with 3 items: X (the expression matrix), z (an integer 
#' vector of inherited patient labels for each cell) and true_y (an integer vector 
#' indicating the "true" answer for whether or not a cell is associated with
#' the outcome-of-interest (1) or not (0)). 
#' 
#' @export
#'
#' @examples
#' NULL
emmil_simulate_data <- function(num_cells = 100, num_features = 10, rho = 0.5, zeta = 0.5, shift = 4) { 
  
  # simulate feature matrix and inherited patient labels for each cell
  X <- 
    matrix(stats::rnorm(num_cells * num_features), nrow = num_cells, ncol = num_features)
  z <- stats::rbinom(num_cells, 1, zeta)
  
  # simulate true association with the condition of interest ("Y") for each cell
  true_y <- rep(0, num_cells)
  true_y[z == 1] <- stats::rbinom(sum(z == 1), 1, (1 - rho)) # P(y = 1 | z = 1) = 1 - rho
  
  num_features_to_shift <- round(num_features / 2)
  X[true_y == 1, 1:num_features_to_shift] <- 
    X[true_y == 1, 1:num_features_to_shift] + shift   # Shift columns when y = 1 to distinguish the groups
  
  # return result
  result <- 
    list(X = X, z = z, true_y = true_y)
  
  return(result)
}

#' Simulate data to with a particular choice of beta, to measure model performance.
#'
#' @param num_cells An integer number of cells to simulate
#'
#' @param rho A scalar value between 0 and 1 indicating the probability that
#' that a cell is NOT associated with the outcome-of-interest, given that the
#' cell was sampled from a patient who DOES have the outcome of interest. In other
#' words, the expected number of non-disease associated cells in each sick patient.
#'
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an
#' inherited patient label will be 1. This is the probability
#' that a random person from your population-of-interest will have the condition-of-interest.
#'
#' @param beta True model parameters determining P(y = 1 | X).
#'
#' @return A list with 3 items: X (the expression matrix), z (a numeric
#' vector of inherited patient labels for each cell) and true_y (an integer vector
#' indicating the "true" answer for whether or not a cell is associated with
#' the outcome-of-interest (1) or not (0)).
#'
#' @export
#'
#' @examples
#' NULL
emmil_simulate_data_with_beta <- function(num_cells, beta, rho, zeta){
  mult = 100 # have to generate extra data to meet user parameters
  
  p = length(beta)
  
  # Goals
  npos = ceiling((1-rho) * zeta * num_cells)
  nneg = num_cells - npos
  
  X = matrix(rnorm(mult * num_cells * p), nrow=mult * num_cells, ncol=p)
  XX = X
  # XX[, 1] = XX[, 1] * XX[, 2]
  # XX[, 3] = XX[, 3]
  # XX[, 4] = XX[, 4]^2
  probs = 1/(1 + exp(-XX %*% beta))
  y = rbinom(mult * num_cells, 1, probs)
  
  keep_ix = sample(c(which(y == 1)[1:npos], which(y == 0)[1:nneg]))
  
  X = X[keep_ix, ]
  y = y[keep_ix]
  z = rep(0, num_cells)
  z[y == 1] = 1
  z[y == 0] = rbinom(sum(1-y), 1, rho*zeta/(rho * zeta + (1-zeta)))
  
  return(
    list(
      X = X,
      true_y = y,
      z = z
    )
  )
}



# model initialization ---------------------------------------------------------

#' A function that initializes y_0 (the y vector for the 0th iteration) for the EMMIL algorithm.
#' 
#' This function initializes the first vector of cell-level label predictions 
#' in an EMMIL model. All values of the vector will be between 0 and 1. 
#' Specifically, all entries for cells for which z = 0 (cells from samples without 
#' the condition-of-interest) will be 0, and entries for cells for which z = 1 
#' (cells from samples with the condition-of-interest) will be 1 - `rho`. 
#'
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease associated cells in each sick patient.
#'
#' @return A numeric vector of the same length as z representing the initial disease-
#' association probabilities for each cell in the dataset.
#' 
#' @export
#'
emmil_initialize_y <- function(z, rho){
  y <- rep(0, length(z))
  y[z == 1] <- 1 - rho
  return(y)
}


#' Calculate the case-control adjustment in an EMMIL model.
#' 
#' This function calculates the case-control adjustment used during EMMIL 
#' modeling fitting. In brief, this adjustment ____. TODO. ERIN
#'
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. This is the probability 
#' that a random cell from your population-of-interest will be sampled from a 
#' person who has the condition-of-interest.
#' 
#' @param num_cells An integer representing the total number of cells in the 
#' dataset. Only used if z is missing (otherwise, `length(z)` is used). TODO: Remove
#' 
#' @param num_cells_disease An integer representing the number of cells in the 
#' dataset from samples with the condition- or disease-of-interest. 
#' Only used if z is missing (otherwise, `sum(z)` is used). TODO: Remove
#'
#' @return A scalar value indicating the case control adjustment that must be 
#' added to our logit to adjust the intercept of our model. 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_calculate_case_control_adjustment <- 
  function(z, num_cells, num_cells_disease, rho, zeta) {
    # if (!missing(z)) { 
      num_cells <- length(z)
      num_cells_disease <- sum(z)
    # } else { 
    #   if (missing(num_cells) | missing(num_cells_disease)) { 
    #     stop("If z is not provided, then num_cells and num_cells_disease are required.")
    #   }
    # }
    
    case_control_intercept_adjustment <-
      -log((1 - rho) * num_cells_disease / (num_cells - (1 - rho) * num_cells_disease)) + 
      log((1 - rho) * zeta / (1 - (1 - rho) * zeta))
    
    return(case_control_intercept_adjustment)
  }

#' Calculate the sample-label adjustment in an EMMIL model.
#' 
#' This function calculates the sample-label adjustment used during EMMIL 
#' modeling fitting and prediction. In brief, this adjustment adjusts the probability
#' scores for cells sampled from sick patients (i.e. those for whom z = 1), for whom
#' we expect a higher degree of disease-association. TODO: ERIN (review)
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. This is the probability 
#' that a random cell from your population-of-interest will be sampled from a 
#' person who has the condition-of-interest.
#'
#' @return A scalar value indicating the case control adjustment that must be 
#' added to the predicted logit to effectively adjust the intercept of the model. 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_calculate_sample_label_adjustment <- function(rho, zeta) { 
  intercept_adjustment <- -log(rho * zeta/(1 - (1-rho) * zeta))
  
  return(intercept_adjustment)
  
  }


# fit models -------------------------------------------------------------------

#' Fit an EMMIL model using glmnet.
#' 
#' This function is a part of the high-level API for emmil. It fits an EMMIL model 
#' using `glmnet` to make instance label predictions at each iteration. 
#' 
#' @param X The numeric feature matrix of single-cell data. Each row is a cell, and each 
#' column is a feature. 
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param lambda A non-negative scalar value representing the total amount of 
#' regularization for the `glmnet` model. Sometimes called the "penalty" parameter 
#' for a LASSO or elastic net model. 
#' 
#' @param case_control_adjustment TODO
#' 
#' @param num_iterations An integer indicating the maximum
#' number of iterations to use when fitting the model.  
#' 
#' @param alpha A number between 0 and 1 indicating the proportion of L1 regularization
#' and L2 regularization in the model. A value of 1 (the default) indicates a pure 
#' LASSO model, whereas a value of 0 indicates a pure ridge regression model. Any 
#' values inbetween indicate an elastic net model. 
#' 
#' @param early_stopping_patience An integer number of iterations for which the improvement
#' in the log-likelihood must be lower than `early_stopping_tolerance` in a row before 
#' model training is stopped. Defaults to 3. Must be larger than 0. 
#' 
#' @param early_stopping_tolerance A scalar value indicating the minimum value by which the
#' training log-likelihood must improve to continue training. If the improvement in the 
#' log-likelihood is smaller than this value for `early_stopping_patience` iterations 
#' in a row, model training will be stopped before `num_iterations` iterations are fit.
#' Defaults to 1e-4. TODO: ERIN or TIM - Do we indicate that this value is for the 
#' normalized (i.e. divided by number of cells in the training set) log-likelihood?
#'
#' @return An emmil_model object that inherits from the \code{\link[glmnet]{glmnet}} class. 
#' 
#' @export
#' 
#' @importFrom glmnet glmnet
#'
#' @examples
#' NULL
#' 
emmil_fit_glmnet <- 
  function(
    X, 
    z, 
    rho, 
    zeta, 
    lambda, 
    alpha = 1, 
    case_control_adjustment = 0, 
    num_iterations = 20L, 
    early_stopping_patience = 3L, 
    early_stopping_tolerance = 1e-4
  ) { 
    
    if (early_stopping_patience < 1) { 
      stop("early_stopping_patience must be 1 or larger.")
    }
    
    # initialize model 
    y <- initialize_y(z, rho)
    y_list <- rep(list(y), num_iterations)
    lls <- rep(0, num_iterations)
    num_cells <- length(z)
    num_cells_disease <- sum(z)
    patience_counter <- 0
    
    # iterate EM steps
    for (i in 1:num_iterations) {
      
      # break if patience has been exhausted
      if (patience_counter >= early_stopping_patience) {
        lls <- lls[1:i-1]
        break
      }
      
      # maximization step 
      model <- glmnet::glmnet(x = X, y = cbind(1 - y, y), alpha = alpha, family = "binomial")
      probabilities <- stats::predict(model, newx = X, s = lambda, type = "response") 
      predictions <- logit(probabilities) # logit 
      lls[i] <- 
        loglik(
          z = z, 
          predictions = predictions, 
          rho = rho, 
          num_cells_disease = num_cells_disease, 
          num_cells = num_cells
        )
      
      # check patience 
      if (i > 1) { 
        ll_difference <- (lls[i])/num_cells - (lls[i-1])/num_cells
        if (ll_difference < early_stopping_tolerance) { 
          patience_counter <- patience_counter + 1
        } else { 
          patience_counter <- 0
        }
      }
      
      # expectation step 
      adjusted_predictions <- predictions + case_control_adjustment
      y <- update_y(adjusted_predictions, z, rho, zeta)
      y_list[i] <- list(y)
    }
    
    # return result
    result <- 
      new_emmil_model(
        model = model, 
        rho = rho, 
        zeta = zeta, 
        num_cells = num_cells, 
        num_cells_disease = num_cells_disease, 
        lls = lls, 
        hyperparamters <- list(lambda = lambda, alpha = alpha)
      )
    
    return(result)
  }



#' Fit an EMMIL model using a generalized linear model.
#'
#' @param X The numeric feature matrix of single-cell data. Each row is a cell, and each 
#' column is a feature. 
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param case_control_adjustment TODO
#' 
#' @param num_iterations An integer indicating the maximum
#' number of iterations to use when fitting the model.  
#' 
#' @param early_stopping_patience An integer number of iterations for which the improvement
#' in the log-likelihood must be lower than `early_stopping_tolerance` in a row before 
#' model training is stopped. Defaults to 3. Must be larger than 0. 
#' 
#' @param early_stopping_tolerance A scalar value indicating the minimum value by which the
#' training log-likelihood must improve to continue training. If the improvement in the 
#' log-likelihood is smaller than this value for `early_stopping_patience` iterations 
#' in a row, model training will be stopped before `num_iterations` iterations are fit.
#' Defaults to 1e-4. TODO: ERIN or TIM - Do we indicate that this value is for the 
#' normalized (i.e. divided by number of cells in the training set) log-likelihood?
#'
#' @return A list.
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_fit_glm <- function(
    X, 
    z, 
    rho, 
    zeta, 
    case_control_adjustment = 0,
    num_iterations = 20L, 
    early_stopping_patience = 3L, 
    early_stopping_tolerance = 1e-4
) { 
  
  if (early_stopping_patience < 1) { 
    stop("early_stopping_patience must be 1 or larger.")
  }
  
  # initialize model 
  y <- initialize_y(z, rho)
  y_list <- rep(list(y), num_iterations)
  lls <- rep(0, num_iterations)
  num_cells <- length(z)
  num_cells_disease <- sum(z)
  patience_counter <- 0
  
  X <- as.data.frame(X)
  X$y <- y 
  
  # iterate EM steps
  for(i in 1:num_iterations){
    
    # break if patience has been exhausted
    if (patience_counter >= early_stopping_patience) {
      lls <- lls[1:i-1]
      break
    }
    
    # maximization step 
    model <- stats::glm(y ~ ., data = X, family = "quasibinomial")
    probabilities <- stats::predict(model, newdata = X, type = "response") 
    predictions <- log(probabilities / (1 - probabilities))  
    lls[i] <- 
      loglik(
        z = z, 
        predictions = predictions, 
        rho = rho, 
        num_cells_disease = num_cells_disease, 
        num_cells = num_cells
      )
    
    # check patience 
    if (i > 1) { 
      ll_difference <- (lls[i])/num_cells - (lls[i-1])/num_cells
      if (ll_difference < early_stopping_tolerance) { 
        patience_counter <- patience_counter + 1
      } else { 
        patience_counter <- 0
      }
    }
    
    # expectation step 
    adjusted_predictions <- predictions + case_control_adjustment
    y <- update_y(adjusted_predictions, z, rho, zeta)
    y_list[i] <- list(y)
    X$y <- y
  }
  
  # return result
  result <- 
    new_emmil_model(
      model = model, 
      rho = rho, 
      zeta = zeta, 
      num_cells = num_cells, 
      num_cells_disease = num_cells_disease, 
      lls = lls, 
      hyperparamters <- list()
    )
  
  return(result)
}


#' Fit an EMMIL model using a multilayer perceptron.
#'
#' @param X The numeric feature matrix of single-cell data. Each row is a cell, and each 
#' column is a feature.  
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be from a patient with 
#' the condition-of-interest.
#' 
#' @param case_control_adjustment TO DO 
#' 
#' @param num_iterations An integer indicating the maximum
#' number of iterations to use when fitting the model.  
#'  
#' @param num_neurons TO DO 
#' 
#' @param decay TO DO 
#' 
#' @param early_stopping_patience An integer number of iterations for which the improvement
#' in the log-likelihood must be lower than `early_stopping_tolerance` in a row before 
#' model training is stopped. Defaults to 3. Must be larger than 0. 
#' 
#' @param early_stopping_tolerance A scalar value indicating the minimum value by which the
#' training log-likelihood must improve to continue training. If the improvement in the 
#' log-likelihood is smaller than this value for `early_stopping_patience` iterations 
#' in a row, model training will be stopped before `num_iterations` iterations are fit.
#' Defaults to 1e-4. TODO: ERIN or TIM - Do we indicate that this value is for the 
#' normalized (i.e. divided by number of cells in the training set) log-likelihood?
#' 
#' @param ... TO DO 
#'
#' @return TO DO 
#' 
#' @export
#' 
#' @importFrom nnet nnet
#' 
#' @examples
#' NULL
emmil_fit_mlp <- 
  function(
    X, 
    z, 
    rho, 
    zeta, 
    case_control_adjustment = 0, 
    num_iterations = 20L, 
    num_neurons, 
    decay, 
    early_stopping_patience = 3L, 
    early_stopping_tolerance = 1e-4,    
    ...
  ) { 
    
    if (early_stopping_patience < 1) { 
      stop("early_stopping_patience must be 1 or larger.")
    }
    
    # initialize model 
    y <- initialize_y(z, rho)
    lls <- rep(0, num_iterations)
    y_list <- rep(list(y), num_iterations)
    num_cells <- length(z)
    num_cells_disease <- sum(z)
    patience_counter <- 0
    
    # iterate EM steps
    for(i in 1:num_iterations){
      
      # break if patience has been exhausted
      if (patience_counter >= early_stopping_patience) {
        lls <- lls[1:i-1]
        break
      }
      
      # maximization step 
      model <- 
        nnet::nnet(x = X, y = y, size = num_neurons, decay = decay, trace = FALSE, ...)
      probabilities <- stats::predict(model, newdata = X, type = "raw") # probabilities
      # convert probabilities to logit
      predictions <- log(probabilities / (1 - probabilities))
      lls[i] <-
        loglik(
          z = z, 
          predictions = predictions, 
          rho = rho, 
          num_cells_disease = num_cells_disease, 
          num_cells = num_cells
        )
      
      # check patience 
      if (i > 1) { 
        ll_difference <- (lls[i])/num_cells - (lls[i-1])/num_cells
        if (ll_difference < early_stopping_tolerance) { 
          patience_counter <- patience_counter + 1
        } else { 
          patience_counter <- 0
        }
      }
      
      # expectation step 
      adjusted_predictions <- predictions + case_control_adjustment
      y <- update_y(adjusted_predictions, z, rho, zeta)
      y_list[i] <- list(y)
      
    }
    
    # return result
    result <- 
      new_emmil_model(
        model = model, 
        rho = rho, 
        zeta = zeta, 
        num_cells = num_cells, 
        num_cells_disease = num_cells_disease, 
        lls = lls, 
        hyperparamters <- list(num_neurons = num_neurons, decay = decay)
      )
  }



#' Title
#'
#' @param X TO DO
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be from a patient with 
#' the condition-of-interest.
#' 
#' @param case_control_adjustment TO DO
#' 
#' @param num_iterations TO DO
#' 
#' @param model_spec Any tidymodels model specification.
#'
#' @return TO DO
#' @export
#'
#' @examples
#' NULL
emmil_fit_tidymodels <- function(X, z, rho, zeta, case_control_adjustment = 0, num_iterations = 20L, model_spec) { 
  # initialize model 
  y <- initialize_y(z, rho)
  lls <- rep(0, num_iterations)
  num_cells <- length(z)
  num_cells_disease <- sum(z)
  
  # iterate EM steps
  for(i in 1:num_iterations){
    # maximization step 
    model <- 
      model_spec |> 
      parsnip::fit_xy(x = as.data.frame(X), y = y) #y = factor(y, levels = c(0, 1)))
    
    # find probabilities and then convert them to logit
    probabilities <- 
      stats::predict(model, new_data = as.data.frame(X))$.pred # logit 
    probabilities <- 
      pmax(probabilities, 1e-5) |> 
      pmin(1 - 1e-5)
    predictions <- log(probabilities / (1 - probabilities))
    
    lls[i] <- 
      loglik(
        z = z, 
        predictions = predictions, 
        rho = rho, 
        num_cells_disease = num_cells_disease, 
        num_cells = num_cells
      )
    
    # expectation step 
    adjusted_predictions <- predictions + case_control_adjustment
    y <- update_y(adjusted_predictions, z, rho, zeta)
  }
  
  # return result
  result <- 
    list(
      y = y, 
      lls = lls, 
      model = model
    )
  
  
  result <- 
    new_emmil_model(
      model = model, 
      rho = rho, 
      zeta = zeta, 
      num_cells = num_cells, 
      num_cells_disease = num_cells_disease, 
      lls = lls, 
      hyperparamters <- list()
    )
  
  return(result)
}


# adjust predictions using EMMIL's intercept adjustments -----------------------

#' Update the predictions of an EMMIL model using the case-control (i.e. sampling 
#' bias) adjustment
#'
#' @param base_predictions Unadjusted predictions (either logit or probability values) 
#' from an EMMIL model. 
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param prediction_type A string indicating how to interpret `base_predictions`. 
#' Valid options are "logit" (the default), which interprets `base_predictions` as a 
#' vector of logit values, and "prob", which interprets `base_predictions` as a 
#' vector of probability values. Note that `prediction_type` also specifies the 
#' output of this function (either adjusted logit or probability scores, respectively).
#'
#' @return If `prediction_type` == "logit", a vector of adjusted logit values. If 
#' `prediction_type` == "prob", a vector of adjusted probability values. 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_adjust_predictions_case_control <- 
  function(base_predictions, z, rho, zeta, prediction_type = c("logit", "prob")) { 
    case_control_intercept_adjustment <- 
      emmil_calculate_case_control_adjustment(z = z, rho = rho, zeta = zeta)
    
    if (prediction_type == "prob") { 
      base_predictions <- logit(base_predictions)
    }
    
    updated_predictions <- base_predictions + case_control_intercept_adjustment
    
    if (prediction_type == "prob") { 
      updated_predictions <- expit(updated_predictions)
    }
    
    return(updated_predictions)
  }


#' Update the predictions of an EMMIL model using the sample-label adjustment
#' 
#' @param base_predictions Unadjusted predictions (either logit or probability values) 
#' from an EMMIL model. 
#' 
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param zeta A scalar value between 0 and 1 indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param adjust_healthy_samples A boolean value indicating if all cells from healthy 
#' patients should have their predicted disease-association probability set to 0. 
#' Defaults to TRUE. 
#' 
#' @param prediction_type A string indicating how to interpret `base_predictions`. 
#' Valid options are "logit" (the default), which interprets `base_predictions` as a 
#' vector of logit values, and "prob", which interprets `base_predictions` as a 
#' vector of probability values. Note that `prediction_type` also specifies the 
#' output of this function (either adjusted logit or probability scores, respectively).
#'
#' @return If `prediction_type` == "logit", a vector of adjusted logit values. If 
#' `prediction_type` == "prob", a vector of adjusted probability values. 
#' 
#' @export
#'
#' @examples
#' NULL
emmil_adjust_predictions_sample_label <- 
  function(base_predictions, z, rho, zeta, adjust_healthy_samples = TRUE, prediction_type = c("logit", "prob")) { 
    sample_label_intercept_adjustment <- 
      emmil_calculate_sample_label_adjustment(rho = rho, zeta = zeta)
    
    if (prediction_type == "prob") { 
      base_predictions <- logit(base_predictions)
    }
    
    # update predictions
    updated_predictions <- base_predictions
    updated_predictions[z == 1] <- base_predictions[z == 1] + sample_label_intercept_adjustment
    
    if (adjust_healthy_samples) { 
      updated_predictions[z == 0] <- logit(0)
    }
  
    if (prediction_type == "prob") { 
      updated_predictions <- expit(updated_predictions)
    }
    
    return(updated_predictions)
  }

#' Update the predictions of an EMMIL model using the case-control (i.e. sampling 
#' bias) adjustment when a model was trained on a distribution with one set of rho
#' and zeta values, but applied to a distribution with a new set of rho and zeta values. 
#'
#' @param base_predictions Unadjusted predictions (either logit or probability values) 
#' from an EMMIL model. 
#' 
#' @param training_z A numeric vector of inherited patient labels for each cell 
#' from the dataset on which the model was originally trained. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param training_rho A scalar value between 0 and 1 for the dataset on which the 
#' original model was trained indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param training_zeta A scalar value between 0 and 1 for the dataset on which the 
#' original model was trained indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param test_rho A scalar value between 0 and 1 for the new dataset for which 
#' the model was used to make predictions indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease-associated cells in each sick patient.
#' 
#' @param test_zeta A scalar value between 0 and 1 for the new dataset for which 
#' the model was used to make predictions indicating the marginal probability that an 
#' inherited patient label will be 1. In other words, the probability 
#' that a random cell from your population-of-interest will be sampled from a patient with 
#' the condition-of-interest.
#' 
#' @param prediction_type A string indicating how to interpret `base_predictions`. 
#' Valid options are "logit" (the default), which interprets `base_predictions` as a 
#' vector of logit values, and "prob", which interprets `base_predictions` as a 
#' vector of probability values. Note that `prediction_type` also specifies the 
#' output of this function (either adjusted logit or probability scores, respectively).
#'
#' @return If `prediction_type` == "logit", a vector of adjusted logit values. If 
#' `prediction_type` == "prob", a vector of adjusted probability values. 
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_adjust_predictions_case_control_new_rho_and_zeta <- 
  function(
    base_predictions, 
    training_z, 
    training_rho, 
    training_zeta, 
    test_rho, 
    test_zeta, 
    prediction_type = c("logit", "prob")
  ) { 
    case_control_intercept_adjustment <- 
      # original model's case-control adjustment
      emmil_calculate_case_control_adjustment(
        z = training_z, 
        rho = training_rho, 
        zeta = training_zeta
      ) + 
      # new adjustment for the new values of rho and zeta
      log(
        ((1 - test_rho) * test_zeta) / 
          (1 - ((1 - test_rho) * test_zeta))
      ) - 
      log(
        ((1 - training_rho) * training_zeta) / 
          (1 - ((1 - training_rho) * training_zeta))
      )
    
    if (prediction_type == "prob") { 
      base_predictions <- logit(base_predictions)
    }
    
    updated_predictions <- base_predictions + case_control_intercept_adjustment
    
    if (prediction_type == "prob") { 
      updated_predictions <- expit(updated_predictions)
    }
    
    return(updated_predictions)
  }


# visualization ----------------------------------------------------------------

#' Plot the log-likelihood values over each iteration of EMMIL. 
#'
#' @param emmil_model_object An EMMIL model fit using the high-level API. 
#' @param lls A numeric vector of log-likelihood values. Will be ignored if 
#' emmil_model_object is provided. 
#' @param theme A ggplot2 theme. Defaults to ggplot2::theme_bw(). 
#' @param ... Optional additional arguments to pass to ggplot2::geom_point().
#'
#' @return A ggplot2 object. 
#' @export
#' 
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 theme_bw
#'
#' @examples
#' NULL
#' 
emmil_plot_lls <- function(emmil_model_object, lls, theme = ggplot2::theme_bw(), ...) { 
  
  if (!missing(emmil_model_object)) { 
  lls <- 
    emmil_model_object |> 
    emmil_get_log_likelihoods()
  } else if (missing(lls)) { 
    stop("Either emmil_model_object or lls must be provided.")
  }
  
  result <- 
    ggplot2::ggplot(
      ggplot2::aes(x = 1:length(lls), y = lls), 
      data = NULL
    ) + 
    ggplot2::geom_point(...) + 
    ggplot2::labs(
      x = "Iteration", 
      y = "Log-likelihood"
    ) + 
    theme
  
  return(result)
}




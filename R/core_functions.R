# helper functions 


####################################################
# An example on (easy) simulated data
####################################################

# simulate data ----------------------------------------------------------------

#' Simulate data to demonstrate how to use the EM-MIL model. 
#'
#' @param num_cells An integer number of cells to simulate
#' 
#' @param num_features An integer number of features to simulate
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
#' @param shift 
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
emmil_simulate_data <- function(num_cells = 100, num_features = 10, rho = 0.5, zeta = 0.5, shift = 4) { 
  
  # simulate feature matrix and inherited patient labels for each cell
  X <- 
    matrix(rnorm(num_cells * num_features), nrow = num_cells, ncol = num_features)
  z <- rbinom(num_cells, 1, zeta)
  
  # simulate true association with the condition of interest ("Y") for each cell
  true_y <- rep(0, num_cells)
  true_y[z == 1] <- rbinom(sum(z == 1), 1, (1 - rho)) # P(y = 1 | z = 1) = 1 - rho
  
  num_features_to_shift <- round(num_features / 2)
  X[true_y == 1, 1:num_features_to_shift] <- 
    X[true_y == 1, 1:num_features_to_shift] + shift   # Shift columns when y = 1 to distinguish the groups
  
  # return result
  result <- 
    list(X = X, z = z, true_y = true_y)
  
  return(result)
}


# model initialization ---------------------------------------------------------

#' A function that initializes y_0 (the y vector for the 0th iteration) for the EM-MIL algorithm.
#'
#' @param z A numeric vector of inherited patient labels for each cell. 1 if the 
#' cell was sampled from a patient with the disease and 0 otherwise. 
#' 
#' @param rho A scalar value between 0 and 1 indicating the probability that 
#' that a cell is NOT associated with the outcome-of-interest, given that the 
#' cell was sampled from a patient who DOES have the outcome of interest. In other 
#' words, the expected number of non-disease associated cells in each sick patient.
#'
#' @return A vector of the same length as z representing the initial disease-
#' association probabilities for each cell in the dataset.
#' 
#' @export
#'
emmil_initialize_y <- function(z, rho){
  y <- rep(0, length(z))
  y[z == 1] <- 1 - rho
  return(y)
}


#' Title
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
#' inherited patient label will be 1. In other words, something like the probabilty 
#' that a random person from your population-of-interest will have the condition-of-interest
#' (TODO: ASK ERIN ABOUT THIS). 
#' @param num_cells TO DO 
#' @param num_cells_disease TO DO  
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
    if (!missing(z)) { 
      num_cells <- length(z)
      num_cells_disease <- sum(z)
    } else { 
      if (missing(num_cells) | missing(num_cells_disease)) { 
        stop("If z is not provided, then num_cells and num_cells_disease are required.")
      }
    }
    
    case_control_intercept_adjustment <-
      -log((1 - rho) * num_cells_disease / (num_cells - (1 - rho) * num_cells_disease)) + 
      log((1 - rho) * zeta / (1 - (1 - rho) * zeta))
    
    return(case_control_intercept_adjustment)
  }


# fit models -------------------------------------------------------------------

#' Fit an EMMIL model using glmnet.
#'
#' @param X TODO
#' @param z TODO
#' @param rho TODO
#' @param zeta TODO
#' @param lambda TODO
#' @param case_control_adjustment TODO
#' @param num_iterations TODO
#' @param alpha TODO
#'
#' @return A list.
#' 
#' @export
#' 
#' @importFrom glmnet glmnet
#'
#' @examples
#' NULL
#' 
emmil_fit_glmnet <- function(X, z, rho, zeta, lambda, alpha = 1, case_control_adjustment = 0, num_iterations = 20L) { 
  
  # initialize model 
  y <- initialize_y(z, rho)
  y_list <- rep(list(y), num_iterations)
  lls <- rep(0, num_iterations)
  num_cells <- length(z)
  num_cells_disease <- sum(z)
      
  # iterate EM steps
  for (i in 1:num_iterations) {
    # maximization step 
    model <- glmnet::glmnet(x = X, y = cbind(1 - y, y), alpha = alpha, family = "binomial")
    probabilities <- predict(model, newx = X, s = lambda, type = "response") 
    predictions <- log(probabilities / (1 - probabilities)) # logit 
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
    y_list[i] <- list(y)
  }
  
  # return result
  result <- 
    list(
      y = y, 
      y_list = y_list, 
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
      hyperparamters <- list(lambda = lambda, alpha = alpha)
    )
  
  return(result)
}



#' Fit an EMMIL model using a generalized linear model.
#'
#' @param X TODO
#' @param z TODO
#' @param rho TODO
#' @param zeta TODO
#' @param lambda TODO
#' @param case_control_adjustment TODO
#' @param num_iterations TODO
#'
#' @return A list.
#' 
#' @export
#'
#' @examples
#' NULL
#' 
emmil_fit_glm <- function(X, z, rho, zeta, case_control_adjustment = 0, num_iterations = 20L) { 
  
  # initialize model 
  y <- initialize_y(z, rho)
  y_list <- rep(list(y), num_iterations)
  lls <- rep(0, num_iterations)
  num_cells <- length(z)
  num_cells_disease <- sum(z)
  
  X <- as.data.frame(X)
  X$y <- y 
  
  # iterate EM steps
  for(i in 1:num_iterations){
    # maximization step 
    model <- glm(y ~ ., data = X, family = "quasibinomial")
    probabilities <- predict(model, newdata = X, type = "response") 
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
    y_list[i] <- list(y)
    X$y <- y
  }
  
  # return result
  result <- 
    list(
      y = y, 
      y_list = y_list, 
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


#' Fit an EMMIL model using a multilayer perceptron.
#'
#' @param X TO DO 
#' @param z TO DO 
#' @param rho TO DO 
#' @param zeta TO DO 
#' @param case_control_adjustment TO DO 
#' @param num_iterations TO DO 
#' @param num_neurons TO DO 
#' @param decay TO DO 
#' @param ... TO DO 
#'
#' @return TO DO 
#' 
#' @export
#' 
#' @importFrom nnet nnet
#' 
#' 
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
    ...
  ) { 
    # initialize model 
    y <- initialize_y(z, rho)
    lls <- rep(0, num_iterations)
    y_list <- rep(list(y), num_iterations)
    num_cells <- length(z)
    num_cells_disease <- sum(z)
    
    # iterate EM steps
    for(i in 1:num_iterations){
      # maximization step 
      # model <- nnet::nnet(x = X, y = y, entropy = TRUE, trace = FALSE, ...)
      model <- 
        nnet::nnet(x = X, y = y, size = num_neurons, decay = decay, trace = FALSE, ...)
      probabilities <- predict(model, newdata = X, type = "raw") # probabilities
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
      
      # expectation step 
      adjusted_predictions <- predictions + case_control_adjustment
      y <- update_y(adjusted_predictions, z, rho, zeta)
      y_list[i] <- list(y)
      
    }
    
    # return result
    result <- 
      list(
        y = y, 
        lls = lls, 
        model = model,
        y_list = y_list
      )
    
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
#' @param z TO DO
#' @param rho TO DO
#' @param zeta TO DO
#' @param case_control_adjustment TO DO
#' @param num_iterations TO DO
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
      predict(model, new_data = as.data.frame(X))$.pred #, type = "prob")$.pred_1 # logit 
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


# visualization ----------------------------------------------------------------

#' Plot the log-likelihood values over each iteration of EMMIL. 
#'
#' @param lls A numeric vector of log-likelihood values. 
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
emmil_plot_lls <- function(lls, theme = ggplot2::theme_bw(), ...) { 
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

emmil_predict <- function() {
  return(NULL)
}




# EM-MIL helper functions

library(glmnet)

#' A function that initializes y_0 (the y vector for the 0th iteration) for the EM-MIL algorithm.
#'
#' @param z The vector of inherited parent labels for each cell. 
#' @param rho The probability that a cell from a diseased sample is not 
#' disease-associated (i.e. the probability that y = 0 when z = 1).
#'
#' @return A vector of the same length as z representing the initial disease-
#' association probabilities for each cell in the dataset.
#' 
#' @export
#'
initialize_y <- function(z, rho){
  y <- rep(0, length(z))
  y[z == 1] <- 1 - rho
  return(y)
}

#' Obtain a binary class probability from a logit
#'
#' @param logit A numeric vector of logit (log-odds) values
#'
#' @return A numeric vector of the same length as logit representing the probability
#' that each observation is of class 1. 
#' 
#'
p_from_logit <- function(logit) {
  prediction <- 1/(1 + exp(-logit))
  return(prediction)
}

# expectation step
# xb is the logit that comes out of the glmnet model 

#' After a maximization step is computed by the EM-MIL algorithm, update the
#' y-values (disease-association probabilities) for each cell in the next iteration 
#'
#' @param predictions The (logit) model predictions from the maximization step
#' @param z The inherited (parent) labels for each cell
#' @param rho The probability that a cell from a diseased sample is not 
#' disease-associated (i.e. the probability that y = 0 when z = 1). 
#' @param zeta The probability that the inherited patient label of any given cell will be 1. 
#' (Will this be complicated by the fact that the patients are not balanced in the dataset? 
#' I.e. patients will not have the same number of cells). 
#'
#' @return A numeric vector of updated probabilities that each cell is disease-associated.
#'
update_y <- function(predictions, z, rho, zeta) {
  intercept_adjustment <- -log(rho * zeta/(1 - (1-rho) * zeta))
  
  y <- rep(0, length(z))
  y[z == 1] <- p_from_logit(predictions[z == 1] + intercept_adjustment)
  return(y)
}

#' Compute the log-likelihood of the EM-MIL algorithm
#'
#' @param z The inherited (parent) labels for each cell
#' @param predictions The (logit) model predictions from the maximization step
#' @param num_cells_disease The number of cells in the (training) dataset that come 
#' from samples with the disease. 
#' @param num_cells The total number of cells in the (training) dataset.
#'
#' @return A numeric scalar indicating the log-likelihood for the current iteration.
#'
loglik <- function(z, predictions, num_cells_disease, num_cells){
  intercept_adjustment = ll = 0
  ll = 
    ll + 
    sum(
      log(
        exp(predictions)[z == 1] + 
          rho * num_cells_disease / (num_cells - num_cells_disease * (1 - rho))
      )
    )
  ll = 
    ll + 
    sum(z == 0) * log(1 - rho * num_cells_disease / (num_cells - num_cells_disease * (1 - rho)))
  ll = ll - sum(log(exp(predictions) + 1))
  return(ll)
}

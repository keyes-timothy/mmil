# EMMIL helper functions

#' A function that initializes y_0 (the y vector for the 0th iteration) for the EM-MIL algorithm.
#'
#' @param z The vector of inherited parent labels for each cell. 
#' 
#' @param rho The probability that a cell from a diseased sample is not 
#' disease-associated (i.e. the probability that y = 0 when z = 1).
#' 
#' @param true_y Optional. A vector of length length(z) containing the known
#' labels for each cell. 
#'
#' @return A vector of the same length as z representing the initial disease-
#' association probabilities for each cell in the dataset.
#' 
#' @export
#'
initialize_y <- function(z, rho, true_y = NULL){
  
  if(is.null(true_y)) { 
    y <- rep(0, length(z))
    y[z == 1] <- 1 - rho
    
  } else { 
    y <- true_y 
    y[is.na(y) & z == 0] <- 0
    y[is.na(y) & z == 1] <- 1 - rho
  }
  return(y)
}

#' Obtain a binary class probability from a logit
#' 
#' This is the expit function.
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

# an alias
expit <- p_from_logit


#' The logit function. 
#'
#' @param p A vector of numeric scores (generally, probabilities between 0 and 1).
#'
#' @return A vector of the same length as `p`, after logit transformation.
#' 
#'
#' @examples
#' NULL
#' 
logit <- function(p) {
  return(log(p / (1 - p)))
}

#' After a maximization step is computed by the EMMIL algorithm, update the
#' y-values (disease-association probabilities) for each cell in the next iteration 
#'
#' @param predictions The (logit) model predictions from the maximization step
#' @param z The inherited (parent) labels for each cell
#' @param rho The probability that a cell from a diseased sample is not 
#' disease-associated (i.e. the probability that y = 0 when z = 1). 
#' @param zeta The probability that the inherited patient label of any given cell will be 1. 
#' @param true_y Optional. A numeric vector of length length(z) containing the known
#' labels for each cell. 
#'
#' @return A numeric vector of updated probabilities that each cell is disease-associated.
#'
update_y <- function(predictions, z, rho, zeta, true_y = NULL) {
  # adjustment for cells sampled from from people with the disease
  intercept_adjustment <- 
    emmil_calculate_sample_label_adjustment(rho = rho, zeta = zeta)
  
  # update y
  if(is.null(true_y)) { 
    y <- rep(0, length(z))
    adjusted_logits <- predictions[z == 1] + intercept_adjustment
    y[z == 1] <- p_from_logit(adjusted_logits)
    
  } else { 
    y <- true_y 
    adjusted_logits <- predictions[is.na(y) & z == 1] + intercept_adjustment
    y[is.na(y) & z == 0] <- 0
    y[is.na(y) & z == 1] <- p_from_logit(adjusted_logits)
  }
  
  return(y)
}

#' Compute the log-likelihood of the EM-MIL algorithm
#'
#' @param z An integer vector containing the inherited (parent) labels for each cell
#' 
#' @param predictions A numeric vector containing the (logit) model predictions from the maximization step
#' 
#' @param num_cells_disease A scalar indicating the number of cells in the (training) dataset that come 
#' from samples with the disease. 
#' 
#' @param num_cells A scalar indicating the total number of cells in the (training) dataset.
#' 
#' @param rho The probability that a cell from a diseased sample is not 
#' disease-associated (i.e. the probability that y = 0 when z = 1).
#'
#' @return A numeric scalar indicating the log-likelihood for the current iteration.
#' 
#' @export
#'
emmil_loglik <- function(z, predictions, rho, num_cells_disease, num_cells) {
  ll <- 0
  ll <- 
    ll + 
    sum(
      log(
        exp(predictions)[z == 1] + 
          rho * num_cells_disease / (num_cells - num_cells_disease * (1 - rho))
      )
    )
  ll <-
    ll + 
    sum(z == 0) * log(1 - rho * num_cells_disease / (num_cells - num_cells_disease * (1 - rho)))
  ll = ll - sum(log(exp(predictions) + 1))
  return(ll)
}

# alias
loglik <- emmil_loglik

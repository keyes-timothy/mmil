# helper functions 

source(here::here("scripts", "emmil_helper_functions.R"))



####################################################
# An example on (easy) simulated data
####################################################
library(glmnet)
library(ggplot2)

# simulate data 
num_cells <- 100
num_features <- 10
rho <- 0.5
zeta <- 0.5
X <- matrix(rnorm(n * p), nrow = num_cells, ncol = num_features)
z <- rbinom(num_cells, 1, zeta)

true_y <- rep(0, num_cells)
true_y[z == 1] <- rbinom(sum(z == 1), 1, (1 - rho)) # P(y = 1 | z = 1) = 1 - rho
X[true_y == 1, 1:5] <- X[true_y == 1, 1:5] + 4    # Shift columns by 4 when y = 1 to distinguish the groups

# Here we address our sampling bias.
# This should be pretty close to 0, because we don't have any sampling bias!
# Including for completeness.
# Question from Tim: Can we estimate rho and zeta (but especially zeta) from the training data?
num_cells_disease <- sum(z) # number of cells that come from samples with the disease, i.e. with z = 1

case_control_intercept_adjustment <-
  -log((1 - rho) * num_cells_disease / (num_cells - (1 - rho) * num_cells_disease)) + 
  log((1 - rho) * zeta / (1 - (1 - rho) * zeta))

# Choose some lambda
lambda <- 0.005 # chosen at random

y <- initialize_y(z, rho)

num_iterations <- 20
lls <- rep(0, num_iterations)

for(i in 1:num_iterations){
  # maximization step 
  model <- glmnet(x = X, y = cbind(1 - y, y), family = 'binomial')
  predictions <- predict(model, newx = X, s = lambda) # logit 
  lls[i] = loglik(z, predictions, num_cells_disease, num_cells)
  
  # expectation step 
  adjusted_predictions <- predictions + case_control_intercept_adjustment
  y <- update_y(adjusted_predictions, z, rho, zeta)
}

# This should be increasing:
data.frame(
  iteration = 1:num_iterations, 
  log_likelihood = lls
) |> 
  ggplot() + 
  geom_point(aes(x = iteration, y = log_likelihood)) + 
  labs(x = "iteration", y = "log likelihood")

# This should differentiate between the groups.
# Let's only look at the instances where z == 1.
# (We know the labels when z == 0.)
ggplot() +
  geom_boxplot(aes(x = 0, y = y[true_y == 0 & z == 1])) +
  geom_boxplot(aes(x = 1, y = y[true_y == 1 & z == 1])) +
  labs(x = "true label", y = "predicted probability")


####################################################
# Helpers
####################################################
initialize.y <- function(z, rho){
  y = rep(0, length(z))
  y[z == 1] = 1-rho
  return(y)
}

p.from.logit <- function(x) 1/(1 + exp(-x))

update.y <- function(xb, z, rho, zeta){
  intercept.adjustment = -log(rho*zeta/(1 - (1-rho)*zeta))

  y = rep(0, length(z))
  y[z == 1] = p.from.logit(xb[z == 1] + intercept.adjustment)
  return(y)
}

loglik <- function(z, xb, n1, n){
  intercept.adjustment = ll = 0
  ll = ll + sum(log(exp(xb)[z == 1] + rho*n1/(n - n1*(1-rho))))
  ll = ll + sum(z == 0) * log(1 - rho*n1/(n - n1*(1-rho)))
  ll = ll - sum(log(exp(xb) + 1))
  return(ll)
}
####################################################
# An example on (easy) simulated data
####################################################
library(glmnet)
library(ggplot2)

n = 100
p = 10
rho = 0.5
zeta = 0.5
X  = matrix(rnorm(n * p), nrow = n, ncol=p)
z  = rbinom(n, 1, zeta)
n1 = sum(z)
n  = length(z)

true.y = rep(0, n)
true.y[z == 1] = rbinom(sum(z == 1), 1, (1-rho)) # P(y = 1 | z = 1) = 1-rho
X[true.y == 1, 1:5] = X[true.y == 1, 1:5] + 4    # Shift columns by 4 when y = 1

# Here we address our sampling bias.
# This should be pretty close to 0, because we don't have any sampling bias!
# Including for completeness.
case.control.intercept.adjustment =
  -log((1-rho) * n1/(n - (1-rho)*n1)) + log((1-rho) * zeta/(1 - (1-rho)*zeta))

# Choose some lambda
lambda = .005 # chosen at random

y = initialize.y(z, rho)

niter = 20
lls = rep(0, niter)
for(i in 1:niter){
  model = glmnet(X, cbind(1-y, y), family = 'binomial')
  xb    = predict(model, X, s = lambda)
  lls[i] = loglik(z, xb, n1, n)

  xb    = xb + case.control.intercept.adjustment
  y     = update.y(xb, z, rho, zeta)
}

# This should be increasing:
ggplot() + geom_point(aes(x=1:niter, y=lls)) + labs(x = "iteration", y = "log likelihood")

# This should differentiate between the groups.
# Let's only look at the instances where z == 1.
# (We know the labels when z == 0.)
ggplot() +
  geom_boxplot(aes(x=0, y=y[true.y == 0 & z == 1])) +
  geom_boxplot(aes(x=1, y=y[true.y == 1 & z == 1])) +
  labs(x = "true label", y = "predicted probability")


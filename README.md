# Mixture LASSO

Biomedical datasets often consist of cells sampled from sick and healthy people. We usually know which patients are sick, but we do not know which cells are diseased, and obtaining gold standard labels can be burdensome. Pinpointing disease-associated cells is important for understanding biology, diagnosing and monitoring disease, and discovering novel therapeutics. 

Here, we show how to use Mixture LASSO, an expectation maximization method to train and calibrate a cell-level classifier using patient-level labels and the assumption that healthy people have no diseased cells. Mixture LASSO is versatile: it can be applied to any classifier that uses the logit link function, and it can be used to perform Platt scaling to calibrate models trained by any method in this missing label setting. 

## How to install
The `devtools` package offers a function to install R packages hosted on Github:
install_github("keyes-timothy/emmil").

There are other methods to install R packages from Github, [here](https://cran.r-project.org/web/packages/githubinstall/vignettes/githubinstall.html) is a useful link on this topic.

## How to use
To use this package, we can first simulate data:

```
ntrain = 500;
data.train = emmil_simulate_data_with_beta(n = ntrain,         # number of samples to create
                                           beta = c(1, 2, 3),  # true model coefficients
                                           rho = 0.5,          # fraction of healthy cells in sick people
                                           zeta = 0.5)         # fraction of cells sampled from sick people   
```

If the cell labels were known, we could train a classifier as usual.
```
library(glmnet)
optimal.model <- glmnet(data.train$X, data.train$true_y, family = "binomial")
```

Without knowing the true cell labels, we can use Mixture LASSO:
```
# number of cells sampled from sick people:
n1 <- sum(data$z) 

# account for any bias in our (finite) sample:
case.control.adj <- -log((1-rho)*n1/(n_train - (1-rho)*n1))
case.control.adj <- case.control.adj + log((1-rho)*zeta/(1 - (1-rho)*zeta)) 

# train the Mixture LASSO model:
mixture.model    <- emmil_fit_glmnet(data.train$X, data.train$z, 
                                  rho = rho,
                                  zeta = zeta, 
                                  lambda = 1e-4, # a lambda must be specified
                                  case_control_adjustment = case.control.adj,
                                  num_iterations = 20)
```

The object `mixture.model` is an `emmil_model object` that inherits from `glmnet` (and can be used as any `glmnet` model).


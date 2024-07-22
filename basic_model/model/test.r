rm(list = ls())

setGeneric("init", function(object) {
    standardGeneric("init")
})

source("model/logistic_manuel.r")
source("utils.r")

library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(ggplot2)
library(glmnet)
library(SGL)
library(DMwR)
library(themis)
library(reshape2)

load("model/mymatrix.RData")
load("model/myy.RData")
# print(X)
# print(y)

set.seed(54)

# n_rows <- 1000
# n_var <- 10
# lambda <- 0.006250

# covariates <- matrix(rnorm(n_var * n_rows, mean = 2, sd = 3), nrow = n_rows, ncol = n_var)
# beta_zeros <- rep(0, n_var %/% 2)
# beta_random <- runif(n_var - n_var %/% 2, min = -1, max = 1)
# true_beta <- c(beta_zeros, beta_random)
# intercept <- 0.5
# prod_scal <- intercept + as.vector(covariates %*% true_beta)
# proba <- 1 / (1 + exp(-prod_scal))
# y <- rbinom(n_rows, 1, prob = proba)


# logistic_classic <- new("logistic_manuel",
#     beta_init = rep(0, n_var),
#     X_init = covariates, y = y, lambda = lambda, intercept_init = 0.0,
#     thresh_cycle = 1e-10, thresh_newton = 1e-9, cycle_ite_max = 1e5, n_step_newton_max = 100
# )
# logistic_classic <- run(logistic_classic)
# beta <- logistic_classic@beta_final
# intercept <- logistic_classic@intercept_final

Sys.sleep(2)

n_rows <- nrow(X)
n_var <- ncol(X)
lambda <- 0.006250
covariates <- X

# print(dim(covariates))
# print(length(y))
# print(n_rows)

logistic_classic <- new("logistic_manuel",
    beta_init = rep(0, n_var),
    X_init = covariates, y = y, lambda = lambda, intercept_init = 0.0,
    thresh_cycle = 1e-10, thresh_newton = 1e-9, cycle_ite_max = 1e5, n_step_newton_max = 100
)
logistic_classic <- run(logistic_classic)
beta <- logistic_classic@beta_final
intercept <- logistic_classic@intercept_final


# logistic_classic <- glmnet:::glmnet.fit(
#     x = covariates, y = y, family = binomial(), alpha = 1, lambda = lambda, intercept = TRUE, maxit = 1e7,
#     thresh = 1e-11, weights = rep(1, n_rows) / n_rows
# )
# beta <- as.numeric(logistic_classic$beta)
# intercept <- logistic_classic$a0




print(dim(covariates))


print(beta)

-crit_logistic(x = covariates, y = y, beta = beta, intercept = intercept, lambda = lambda) # doit être minimisé

library(glmnet)
setClass("Test", representation(x = "numeric", y = "numeric"))

test <- new("Test", x = 1, y = 2)
print(slotNames(test)[1])
X <- matrix(rnorm(100), 10, 10)
y <- rnorm(10)
fit <- cv.glmnet(X, y)
print(fit)

setClass(
    "logistic_manuel",
    representation(
        beta_init = "numeric", # OBLIGATOIRE
        X_init = "matrix", # OBLIGATOIRE
        y = "numeric", # OBLIGATOIRE
        lambda = "numeric", # OBLIGATOIRE
        intercept_init = "numeric", # OBLIGATOIRE
        thresh_cycle = "numeric", # OBLIGATOIRE
        thresh_newton = "numeric", # OBLIGATOIRE
        cycle_ite_max = "numeric", # OBLIGATOIRE
        n_step_newton_max = "numeric", # OBLIGATOIRE
        X = "matrix",
        beta = "numeric",
        intercept = "numeric",
        beta_memory = "numeric",
        intercept_memory = "numeric",
        weights = "numeric",
        proba = "numeric",
        N = "numeric",
        p = "numeric",
        beta_final = "numeric",
        intercept_final = "numeric",
        prod_scal = "numeric",
        z = "numeric",
        weighted_squares = "numeric",
        r = "numeric",
        support = "logical",
        delta_j = "numeric",
        max_delta = "numeric",
        non_zero = "logical",
        means = "numeric",
        stds = "numeric"
    ),
    prototype()
)


####### Les weights etaitent *N avec glmnet!!!!!

setMethod("init", "logistic_manuel", function(object) {
    # object@means <- colMeans(object@X_init)
    # object@stds <- apply(object@X_init, 2, sd)
    # object@X <- scale(object@X_init, center = object@means, scale = object@stds)
    object@X <- object@X_init
    object@beta <- object@beta_init
    object@intercept <- object@intercept_init
    object@N <- dim(object@X)[1]
    object@p <- dim(object@X)[2]
    object@beta_memory <- copy(object@beta)
    object@intercept_memory <- object@intercept
    object@support <- rep(TRUE, object@p) # Utile pour décider sur quelles variables boucler le cycle
    object@non_zero <- rep(TRUE, object@p) # Contient vraiment qui est zero dans beta
    return(object)
})

setGeneric("calc_weights_and_z", function(object) {
    standardGeneric("calc_weights_and_z")
})

setMethod("calc_weights_and_z", "logistic_manuel", function(object) {
    object@prod_scal <- copy(object@intercept + as.vector(object@X %*% object@beta))
    object@proba <- 1 / (1 + exp(-object@prod_scal))
    object@weights <- object@proba * (1 - object@proba) / object@N
    for (i in 1:object@N) {
        if (abs(object@proba[i]) < 1e-5) {
            object@proba[i] <- 0
            object@weights[i] <- 1e-5 / object@N
        }
        if (abs(object@proba[i] - 1) < 1e-5) {
            object@proba[i] <- 1
            object@weights[i] <- 1e-5 / object@N
        }
    }
    object@z <- object@prod_scal + (object@y - object@proba) / (object@weights * object@N)
    # print(1 / (object@weights[39] * object@N))
    # print(object@N)
    # print(object@weights[39])
    # print(object@z[39])
    # print((object@y - object@proba)[39])
    # Sys.sleep(10)
    object@r <- copy(object@z - object@prod_scal)
    object@weighted_squares <- as.vector(object@weights %*% object@X^2)
    object@weighted_squares <- c(object@weighted_squares, sum(object@weights))
    return(object)
})

setGeneric("update_descent", function(object, ...) {
    standardGeneric("update_descent")
})

setMethod("update_descent", "logistic_manuel", function(object, j) {
    somme_r_j <- sum(object@weights * object@r * object@X[, j])
    somme_beta_j <- object@weighted_squares[j] * object@beta[j]
    li_soft <- soft_thresh(somme_r_j + somme_beta_j, object@lambda)
    # if (j == 52) {
    #     print(paste("sum", somme_r_j + somme_beta_j))
    #     print(paste("beta", object@beta[j]))
    # }
    # print(paste(abs(somme_r_i + somme_beta_j), object@lambda))
    value <- li_soft$value
    non_zero_soft_thresh <- li_soft$non_zero_soft_thresh
    new_beta_j <- value / object@weighted_squares[j]
    if (object@non_zero[j] | non_zero_soft_thresh) {
        object@r <- object@r + as.vector(object@X[, j] * (object@beta[j] - new_beta_j))
    }
    object@delta_j <- object@weighted_squares[j] * (object@beta[j] - new_beta_j)^2
    object@beta[j] <- new_beta_j
    if (non_zero_soft_thresh) {
        object@support[j] <- TRUE
        object@non_zero[j] <- TRUE
        # print("oui")
    } else {
        object@support[j] <- FALSE
        object@non_zero[j] <- FALSE
    }
    return(object)
})

setGeneric("update_intercept_descent", function(object) {
    standardGeneric("update_intercept_descent")
})

setMethod("update_intercept_descent", "logistic_manuel", function(object) {
    v <- object@z - as.vector(object@X %*% object@beta)
    new_intercept <- sum(object@weights * v) / sum(object@weights)
    object@r <- object@r + object@intercept - new_intercept
    object@delta_j <- object@weighted_squares[object@p + 1] * (new_intercept - object@intercept)^2
    object@intercept <- new_intercept
    return(object)
})

setGeneric("one_cycle", function(object) {
    standardGeneric("one_cycle")
})

setMethod("one_cycle", "logistic_manuel", function(object) {
    object@max_delta <- 0
    object <- update_intercept_descent(object)
    object@max_delta <- max(object@max_delta, object@delta_j)
    for (j in 1:object@p) {
        if (object@support[j]) {
            # objective_local(object)
            if (j == 52) {
                # print("on teste:")
                # print(sum(object@z))
                # print(abs(sum(object@weights * object@X[, j] * (object@z - object@intercept))))
                # print(sum(object@weights)))
            }
            object <- update_descent(object, j)
            object@max_delta <- max(object@max_delta, object@delta_j)
        }
    }
    # print(object@beta[51])
    return(object)
})

setGeneric("one_Newton", function(object) {
    standardGeneric("one_Newton")
})

setMethod("one_Newton", "logistic_manuel", function(object) {
    object <- calc_weights_and_z(object)
    ite_cycle <- 0
    object@max_delta <- Inf
    former_support <- rep(FALSE, object@p)
    while ((ite_cycle == 0 | !all(former_support == object@support)) & ite_cycle < object@cycle_ite_max %/% object@p) {
        # print(object@max_delta > object@thresh_cycle & ite_cycle < object@cycle_ite_max)
        while (object@max_delta > object@thresh_cycle & ite_cycle < object@cycle_ite_max %/% object@p) {
            object <- one_cycle(object)
            ite_cycle <- ite_cycle + 1
            # print(paste("cycle:", ite_cycle))
            # if (ite_cycle > 2000) {
            #     print(object@beta)
            #     Sys.sleep(2)
            # }
            # print(object@max_delta > object@thresh_cycle & ite_cycle < object@cycle_ite_max)
        }
        # if (ite_cycle == object@cycle_ite_max) {
        #     print("Warning: Cycle itératif max atteint. Le résultat n'est pas convergent.")
        #     print(paste("le max delta vaut", object@max_delta))
        # }
        former_support <- copy(object@support)
        object@support <- rep(TRUE, object@p)
        object <- one_cycle(object)
        ite_cycle <- ite_cycle + 1
        # print(paste("cycle", ite_cycle))
    }
    return(object)
})

setGeneric("objective_local", function(object) {
    standardGeneric("objective_local")
})

setMethod("objective_local", "logistic_manuel", function(object) {
    print(paste("local", 0.5 * sum(object@weights *
        (object@z - object@intercept - as.vector(object@X %*% object@beta))^2) + object@lambda * norm(as.matrix(object@beta), type = "1")))
})

setGeneric("objective_global", function(object) {
    standardGeneric("objective_global")
})

setMethod("objective_global", "logistic_manuel", function(object) {
    astuce_exp_log <- ifelse(object@intercept + as.vector(object@X %*% object@beta) > 50,
        object@intercept + as.vector(object@X %*% object@beta), log(1 + exp(object@intercept + as.vector(object@X %*% object@beta)))
    )
    value <- -(1 / object@N) * sum(object@y * (object@intercept + as.vector(object@X %*% object@beta))
        + astuce_exp_log)
    +object@lambda * norm(as.matrix(object@beta), type = "1")
    if (abs(value > 50)) {
        print(paste("global ouch!!!!", value))
        X <- object@X
        y <- object@y
        save(X, file = "model/mymatrix.RData")
        save(y, file = "model/myy.RData")
    }
})


#### Reste à faire newton par dessus (nb d'itérations prédéterminé)

setGeneric("run", function(object) {
    standardGeneric("run")
})


setMethod("run", "logistic_manuel", function(object) {
    object <- init(object)
    ite_newton <- 0
    max_delta_newton <- Inf
    while (ite_newton < object@n_step_newton_max & max_delta_newton > object@thresh_newton) {
        # print(paste("ite_newt", ite_newton))
        objective_global(object)
        object <- one_Newton(object)
        ite_newton <- ite_newton + 1
        intercept_delta <- object@weighted_squares[object@p + 1] * (object@intercept - object@intercept_memory)^2
        delta_beta <- object@weighted_squares[1:object@p] * (object@beta - object@beta_memory)^2
        delta_global <- append(delta_beta, intercept_delta)
        max_delta_newton <- max(delta_global)
        if (max_delta_newton < object@thresh_newton) {
            # print("Totale convergence dans Newton")
            # print(paste("le max delta vaut", max_delta_newton))
        }
    }
    object@beta_final <- object@beta # * object@stds + object@means
    # print(object@beta_final - object@beta_init)
    object@intercept_final <- object@intercept
    # print("beta:")
    # print(object@beta)
    return(object)
})


# TOUS LES BETA SONT ANNULES: la suite demain

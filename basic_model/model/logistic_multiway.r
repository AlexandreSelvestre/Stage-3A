source("model/logistic_manuel.r")

li_caret_multiway <- list()

li_caret_multiway$library <- "glmnet"

li_caret_multiway$type <- "Classification"

li_caret_multiway$parameters <- data.frame(parameter = c("lambda", "R"), class = c("numeric", "integer"), label = c(
    "la valeur du paramètre lambda", "le rang R"
))

create_grid_multiway <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        lambda <- log(seq(exp(0), exp(0.05), length.out = len + 1)[2:(len + 1)])
        R <- 1:4
    } else {
        lambda <- log(runif(len, min = exp(0.001), max = exp(0.05)))
        R <- 1:4
    }
    data_frame_grid <- expand.grid(lambda = lambda, R = R)
    return(data_frame_grid)
}

better_create_grid_multiway <- function(x, y, len = NULL, search = "grid", lambda_min = 0.001, lambda_max = 0.05, R_min = 1, R_max = 4) {
    if (search == "grid") {
        lambda <- log(seq(exp(lambda_min), exp(lambda_max), length.out = len)[1:len])
        R <- R_min:R_max
    } else {
        lambda <- log(runif(len, min = exp(lambda_min), max = exp(lambda_max)))
        R <- R_min:R_max
    }
    data_frame_grid <- expand.grid(lambda = lambda, R = R)
    return(data_frame_grid)
}

li_caret_multiway$grid <- create_grid_multiway


fit_multiway <- function(x, y, wts, param, lev, last, weights_dict, classProbs, index, eps, ite_max, n_iter_per_reg, k_smote, do_smote, index_variable, is_binary) {
    li_norm <- renormalize_in_model_fit_index_mode(x, index_variable, is_binary)
    x <- li_norm$new_x
    classe_min <- names(which.min(table(y)))
    classe_maj <- setdiff(levels(y), classe_min)

    if (do_smote) {
        li <- apply_smote(x, y, k_smote)
        x <- li$x
        y <- li$y
    } else {
        li <- apply_boot(x, y)
        x <- li$x
        y <- li$y
    }

    # Attention aux indexes qui ont changé!!!

    # On prend en entrée un dataframe
    # Formalisme pour les indices: il faut indiquer qui est dans quel mode dans index. On met -1 pour les variables cliniques
    # On supposera que l'ordre est la même pour toutes les variables: pas de décalage. Le tester au début: checker que c'est bien le même nom de colonne

    weights <- numeric(length(y))
    weights[y == classe_maj] <- weights_dict[[classe_maj]]
    weights[y == classe_min] <- weights_dict[[classe_min]]

    y_numeric <- ifelse(y == classe_min, 1, 0) # CCK donne 1 CHC donne 0
    R <- param$R
    names_x <- colnames(x)
    names_x_multimodal <- names_x[index != -1]


    J <- 0
    continue <- TRUE
    count <- 1
    while (continue) {
        value <- index[count]
        if (value == 1) {
            J <- J + 1
            count <- count + 1
        }
        if (value > 1) {
            continue <- FALSE
        }
        if (value == -1) {
            count <- count + 1
        }
    }

    modes <- unique(index)
    modes <- modes[modes != -1]
    K <- length(modes)
    beta_K_0 <- rnorm(R * K, mean = 0, sd = 1)
    ### Formalisme utilisé pour les matrices : celui du papier


    ## Commencer par définir Z_J en sommant sur les modes. Attention au problème lignes colonnes inversées. Pas encore de variables cliniques
    mat_x_modes <- as.matrix(x)[, index != -1]
    mat_x_restant <- as.matrix(as.matrix(x)[, index == -1])
    n_restant <- ncol(mat_x_restant)
    ####### Petite boucle de l'algorithme: utile dans la grande #########
    petite_boucle <- function(Z_init, vec_Q, beta_previous, intercept_previous, n_iter_per_reg) {
        Z <- rbind(Z_init, t(mat_x_restant))
        vec_Q <- append(vec_Q, rep(1, n_restant))
        # print("taille vec:")
        # print(length(vec_Q))
        Q_inv <- diag(vec_Q)
        Q <- inverse_diag(Q_inv)
        if (any(is.na(Q))) {
            print("Q cata")
            print(diag(Q))
            print("sepa")
            print(vec_Q)
        }
        if (any(is.na(Z))) {
            print("Z cata")
        }

        ########### Avec GLMNET
        Z <- as.matrix(Z)
        logistic_classic <- glmnet:::glmnet.fit(
            x = t(Q %*% Z), y = y_numeric, family = binomial(), alpha = 1,
            weights = weights / dim(Z)[2], lambda = param$lambda, intercept = TRUE, maxit = 1e7,
            thresh = 1e-11
        )
        faux_beta <- as.numeric(logistic_classic$beta)
        beta <- Q %*% faux_beta
        intercept <- logistic_classic$a0
        # print(beta)

        ############ Avec ma version du GLMNET

        # # print(as.vector(Q_inv %*% beta_previous))


        # # logistic_classic <- new("logistic_manuel",
        # #     beta_init = as.vector(Q_inv %*% beta_previous),
        # #     X_init = t(Q %*% Z), y = y_numeric, lambda = param$lambda, intercept_init = intercept_previous,
        # #     thresh_cycle = 1e-11, thresh_newton = 1e-11, cycle_ite_max = 1e6, n_step_newton_max = n_iter_per_reg
        # # )


        # # mat <- t(Q %*% Z)
        # # save(mat, file = "model/mymatrix.RData")
        # # save(y_numeric, file = "model/myy.RData")
        # # print(y_numeric)
        # logistic_classic <- run(logistic_classic)
        # # print(dim(t(Q %*% Z)))
        # # Sys.sleep(2)
        # faux_beta <- logistic_classic@beta_final
        # beta <- Q %*% faux_beta
        # intercept <- logistic_classic@intercept_final


        return(list(beta = beta, intercept = intercept, Z = Z, Q = Q, Q_inv = Q_inv))
    }



    ####### Grande boucle de l'algorithme à itérer autant que nécessaire ##########
    grande_boucle <- function(beta_K, beta_J_complet_previous, beta_K_complet_previous, intercept_previous_K, intercept_previous_J, n_iter_per_reg) {
        #   Définir Z_J_init (sans les variables cliniques)
        li_Z_J_init <- list() # init car sans clinique
        for (r in 1:R) {
            Z_J_init <- matrix(0, nrow = dim(mat_x_modes)[1], ncol = J)
            for (k in 1:K) {
                bloc_k <- mat_x_modes[, ((k - 1) * J + 1):(k * J)]
                Z_J_init <- Z_J_init + beta_K[(r - 1) * K + k] * bloc_k
                li_Z_J_init[[r]] <- Z_J_init
            }
        }
        Z_J_init <- li_Z_J_init[[1]]
        if (R >= 2) {
            for (r in 2:R) {
                Z_J_init <- cbind(Z_J_init, li_Z_J_init[[r]])
            }
        }
        Z_J_init <- t(Z_J_init)

        ## Définir Q_J de transformation pour la régularisation (attention aux colonnes inversées)

        vec_diag <- c()
        for (r in 1:R) {
            norme_1_r <- norm(as.matrix(beta_K[((r - 1) * K + 1):(r * K)]), type = "1")
            diag_elem_r <- rep(norme_1_r, J)
            vec_diag <- append(vec_diag, diag_elem_r)
        }



        ## Appliquer la petite boucle

        li_petite_boucle <- petite_boucle(
            Z_init = Z_J_init, vec_Q = vec_diag, beta_previous = beta_J_complet_previous,
            intercept_previous = intercept_previous_J, n_iter_per_reg = n_iter_per_reg
        )
        beta_J <- li_petite_boucle$beta[1:(J * R)]
        if (n_restant > 0) {
            beta_autre_J <- li_petite_boucle$beta[((J * R) + 1):length(li_petite_boucle$beta)]
        } else {
            beta_autre_J <- c()
        }
        Z_J <- li_petite_boucle$Z
        intercept_J <- li_petite_boucle$intercept
        intercept_previous_J <- intercept_J
        Q_J_full <- li_petite_boucle$Q
        Q_J_inv <- li_petite_boucle$Q_inv
        # print(beta_J)
        ## Opérer maintenant sur Z_K et Q_K


        li_Z_K_init <- list() # init car sans clinique
        for (r in 1:R) {
            Z_K_init <- matrix(0, nrow = dim(mat_x_modes)[1], ncol = K)
            for (j in 1:J) {
                col_to_extract <- seq(j, by = J, length.out = K)
                bloc_j <- mat_x_modes[, col_to_extract]
                # print(beta_J[(r - 1) * J + j])
                Z_K_init <- Z_K_init + beta_J[(r - 1) * J + j] * bloc_j
                li_Z_K_init[[r]] <- Z_K_init
            }
        }
        Z_K_init <- li_Z_K_init[[1]]
        if (R >= 2) {
            for (r in 2:R) {
                Z_K_init <- cbind(Z_K_init, li_Z_K_init[[r]])
            }
        }
        Z_K_init <- t(Z_K_init)
        # print(Z_K_init)
        if (all(abs(Z_K_init) < 1e-10)) {
            print("Z_K est nul")
        }

        ## Définir Q_K de transformation pour la régularisation (attention aux colonnes inversées)

        vec_diag <- c()
        for (r in 1:R) {
            norme_1_r <- norm(as.matrix(beta_J[((r - 1) * J + 1):(r * J)]), type = "1")
            diag_elem_r <- rep(norme_1_r, K)
            vec_diag <- append(vec_diag, diag_elem_r)
        }
        ## Appliquer à nouveau la petite boucle
        # print(beta_K_complet_previous)
        # print(dim(Z_K_init))
        li_petite_boucle <- petite_boucle(
            Z_init = Z_K_init, vec_Q = vec_diag, beta_previous = beta_K_complet_previous,
            intercept_previous = intercept_previous_K, n_iter_per_reg = n_iter_per_reg
        )
        Z_K <- li_petite_boucle$Z
        beta_K <- li_petite_boucle$beta[1:(K * R)]
        if (norm(as.matrix(beta_K), type = "2") > 0.00001) {
            beta_K <- beta_K / norm(as.matrix(beta_K), type = "2")
        } # RENORMALISER
        intercept_K <- li_petite_boucle$intercept
        if (n_restant > 0) {
            beta_autre_K <- li_petite_boucle$beta[(K * R + 1):length(li_petite_boucle$beta)]
        } else {
            beta_autre_K <- c()
        }
        Q_K_full <- li_petite_boucle$Q
        Q_K_inv <- li_petite_boucle$Q_inv

        beta_J_complet <- c(beta_J, beta_autre_J)
        beta_K_complet <- c(beta_K, beta_autre_K)

        # Calculer les deux critères
        # print(beta_J_complet)
        # print(beta_K_complet)
        if (all(abs(Z_K_init) < 1e-10)) {
            # print(Q_K_full)
            # print(Q_K_inv)
            # print(beta_K_complet)
            # print("ANTIJU")
            # print(li_petite_boucle$beta)
        }

        crit_log_J <- crit_logistic(t(Q_J_full %*% Z_J), y_numeric, Q_J_inv %*% beta_J_complet, intercept_J, param$lambda)
        crit_log_K <- crit_logistic(t(Q_K_full %*% Z_K), y_numeric, Q_K_inv %*% beta_K_complet, intercept_K, param$lambda)
        # print(crit_log_K)
        rapport <- abs(crit_log_J - crit_log_K) / abs(crit_log_J)

        return(list(
            beta_J = beta_J, beta_K = beta_K, beta_autre = beta_autre_K, intercept_J = intercept_J, intercept_K = intercept_K,
            rapport = rapport, crit_log_J = crit_log_J, crit_log_K = crit_log_K, beta_J_complet = beta_J_complet, beta_K_complet = beta_K_complet
        ))
    }

    ##### Exécuter la grande boucle autant que les critères le demandent
    continue <- TRUE
    beta_K <- beta_K_0
    iteration <- 1
    memoire_crit_J <- 1
    memoire_crit_K <- 1
    # beta_J_complet_previous <- 0.000001 * rnorm(J * R + n_restant, mean = 0, sd = 1)
    # beta_K_complet_previous <- 0.000001 * rnorm(K * R + n_restant, mean = 0, sd = 1)
    beta_J_complet_previous <- rep(0, J * R + n_restant)
    beta_K_complet_previous <- rep(0, K * R + n_restant)
    intercept_previous_J <- 0.0
    intercept_previous_K <- 0.0
    debut_time <- Sys.time()
    while (continue) {
        # print(iteration)
        # print(iteration)
        li_grande_boucle <- grande_boucle(beta_K, beta_J_complet_previous, beta_K_complet_previous,
            intercept_previous_J = intercept_previous_J,
            intercept_previous_K = intercept_previous_K, n_iter_per_reg = n_iter_per_reg
        )
        beta_J <- li_grande_boucle$beta_J
        beta_K <- li_grande_boucle$beta_K
        beta_autre <- li_grande_boucle$beta_autre
        intercept_J <- li_grande_boucle$intercept_J
        intercept_K <- li_grande_boucle$intercept_K
        # print(beta_K)

        beta_J_complet_previous <- li_grande_boucle$beta_J_complet
        beta_K_complet_previous <- li_grande_boucle$beta_K_complet
        intercept_previous_J <- intercept_J
        intercept_previous_K <- intercept_K

        rapport <- li_grande_boucle$rapport
        crit_log_J <- li_grande_boucle$crit_log_J
        crit_log_K <- li_grande_boucle$crit_log_K
        # print(rapport)
        if (rapport < eps) {
            continue <- FALSE
        }
        if (iteration >= ite_max) {
            continue <- FALSE
            print(paste("Warning, pas de convergence. Le rapport vaut:", rapport))
        }
        if (abs(crit_log_J - memoire_crit_J) < 1e-10 & abs(crit_log_K - memoire_crit_K) < 1e-10) {
            continue <- FALSE
            print("Warning, LOOP!!")
            print(paste("Rapport:", rapport))
            print(paste("Crit_J:", crit_log_J))
            print(paste("Crit_K:", crit_log_K))
        }
        iteration <- iteration + 1
        memoire_crit_J <- crit_log_J
        memoire_crit_K <- crit_log_K
        delta_t <- Sys.time() - debut_time
        if (delta_t > 30) {
            continue <- FALSE
            print("Warning, trop long")
            print(paste("Rapport:", rapport))
            print(paste("Crit_J:", crit_log_J))
            print(paste("Crit_K:", crit_log_K))
        }
    }
    futur_fit <- list(beta_J = beta_J, beta_K = beta_K, beta_autre = beta_autre, intercept = intercept_K, R = R, J = J, K = K, lev = lev, li_norm = li_norm, classe_maj = classe_maj, classe_min = classe_min)
    futur_fit$beta_unfolded <- get_beta_multiway(futur_fit)
    return(futur_fit)
}

li_caret_multiway$fit <- fit_multiway


get_beta_multiway <- function(modelFit) {
    beta_J <- modelFit$beta_J
    beta_K <- modelFit$beta_K
    beta_autre <- modelFit$beta_autre
    R <- modelFit$R
    J <- modelFit$J
    K <- modelFit$K
    beta <- rep(0, J * K + length(beta_autre))
    for (r in 1:R) {
        beta_r <- rep(0, J * K)
        for (j in 1:J) {
            for (k in 1:K) {
                beta_r[(k - 1) * J + j] <- beta_J[(r - 1) * J + j] * beta_K[(r - 1) * K + k]
            }
        }
        beta[1:(J * K)] <- beta[1:(J * K)] + beta_r
    }
    if (length(beta_autre) > 0) {
        beta[(J * K + 1):length(beta)] <- beta_autre #### Attention: faiblesse, les variables cliniques sont obligatoirement à la fin
    }
    return(beta)
}

predict_multiway <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- get_beta_multiway(modelFit)
    intercept <- modelFit$intercept
    lev <- modelFit$lev
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta + intercept
    })
    proba <- 1 / (1 + exp(-value))
    predicted_labels <- ifelse(proba > 0.5, classe_min, classe_maj)
    return(predicted_labels)
}

li_caret_multiway$predict <- predict_multiway


li_caret_multiway$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- get_beta_multiway(modelFit)
    intercept <- modelFit$intercept
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta + intercept
    })
    proba_class_1 <- 1 / (1 + exp(-value))
    proba_class_0 <- 1 - proba_class_1
    str_min <- as.character(classe_min)
    str_maj <- as.character(classe_maj)
    return(setNames(data.frame(proba_class_0, proba_class_1), c(str_maj, str_min)))
}

li_caret_multiway$loop <- NULL


setMethod("train_method", "apply_model", function(object) {
    # x, y, wts, param, lev, last, weights, classProbs, index, R, eps, ite_max, beta_K_0

    li <- reorder_in_modes(object@train_cols[, object@col_x], index_mode = object@index_mode, index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary, name_mode = object@name_mode, name_variable = object@name_variable, name_bloc = object@name_bloc)
    object@train_cols[, object@col_x] <- li$x
    colnames(object@train_cols)[colnames(object@train_cols) %in% object@col_x] <- colnames(li$x)

    li <- reorder_in_modes(object@test_set[, object@col_x], index_mode = object@index_mode, index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary, name_mode = object@name_mode, name_variable = object@name_variable, name_bloc = object@name_bloc)
    object@test_set[, object@col_x] <- li$x ### suite...
    colnames(object@test_set)[colnames(object@test_set) %in% object@col_x] <- colnames(li$x)

    li <- reorder_in_modes(object@data_used[, object@col_x], index_mode = object@index_mode, index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary, name_mode = object@name_mode, name_variable = object@name_variable, name_bloc = object@name_bloc)
    object@data_used[, object@col_x] <- li$x ### suite...
    colnames(object@data_used)[colnames(object@data_used) %in% object@col_x] <- colnames(li$x)
    write_xlsx(object@data_used, "..\\data\\no.xlsx")

    object@col_x <- setdiff(names(object@data_used), c(object@info_cols$exclude_cols, object@name_y))


    object@index_variable <- li$index_variable
    object@name_variable <- li$name_variable
    object@index_mode <- li$index_mode
    object@name_mode <- li$name_mode
    object@index_bloc <- li$index_bloc
    object@name_bloc <- li$name_bloc
    object@is_binary <- li$is_binary


    grid <- better_create_grid_multiway(
        x = object@train_cols[, object@col_x], y = object@y_train, len = object@tuneLength,
        search = object@search, lambda_min = object@lambda_min, lambda_max = object@lambda_max, object@R_min, object@R_max
    )
    # print(object@index_variable)

    if (object@do_parallel) {
        numCores <- detectCores()
        cl <- makePSOCKcluster(numCores - 1)
        registerDoParallel(cl)
        clusterEvalQ(cl, {
            files <- list.files("./utils", full.names = TRUE, pattern = "\\.r$")
            for (file in files) {
                source(file)
            }
        })
        clusterExport(cl, varlist = c("get_beta_multiway"))
    }
    object@model <- caret::train(
        y = object@y_train, x = object@train_cols[, object@col_x], index = object@index_mode,
        method = li_caret_multiway, trControl = object@cv, metric = "AUC",
        tuneLength = object@tuneLength, weights_dict = object@weights, tuneGrid = grid,
        eps = object@eps, ite_max = object@ite_max, n_iter_per_reg = object@n_iter_per_reg,
        k_smote = object@k_smote, do_smote = object@do_smote, index_variable = object@index_variable, is_binary = object@is_binary
    )
    if (object@do_parallel) {
        stopCluster(cl)
    }
    return(object)
})

setMethod("get_results", "apply_model", function(object) {
    # x_train <- reorder_in_modes(object@train_cols[, object@col_x], index_mode = object@index_mode, index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary)$x
    # x_test <- reorder_in_modes(object@test_set[, object@col_x], index_mode = object@index_mode, index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary)$x
    object@predictions <- as.vector(predict(object@model, newdata = as.matrix(object@test_set[, object@col_x])))
    object@predictions_proba <- predict(object@model, newdata = as.matrix(object@test_set[, object@col_x]), type = "prob")
    object@predictions_train_proba <- predict(object@model, newdata = as.matrix(object@train_cols[, object@col_x]), type = "prob")
    return(object)
})

setMethod("importance_method", "apply_model", function(object) {
    object@beta_final <- object@model$finalModel$beta_unfolded
    best_beta_J <- object@model$finalModel$beta_J
    best_beta_K <- object@model$finalModel$beta_K
    best_beta_autre <- object@model$finalModel$beta_autre
    Variable_names_seul <- unique(object@name_variable[object@index > -0.5])
    Variables_temps_seul <- unique(object@name_mode[object@index > -0.5])
    Variable_names <- c()
    Variables_temps <- c()
    R <- object@model$bestTune$R
    for (r in 1:R) {
        variable_name_local <- paste0(Variable_names_seul, "_", r)
        Variable_names <- append(Variable_names, variable_name_local)
        variables_temps_local <- paste0(Variables_temps_seul, "_", r)
        Variables_temps <- append(Variables_temps, variables_temps_local)
    }
    variable_importance <- data.frame(Variable = Variable_names, Overall = abs(best_beta_J))
    # print(variable_importance)
    name_df <- paste0("variable_importance_", object@id_term)
    object@li_df_var_imp[[name_df]] <- variable_importance
    variable_importance <- variable_importance[order(-variable_importance$Overall), ]
    variable_importance <- subset(variable_importance, Overall > 0.0001)

    image <- ggplot2::ggplot(variable_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistic_multiway/importance_var_multi", "_", object@id_term, ".png"), image)

    time_importance <- data.frame(Variable = Variables_temps, Overall = abs(best_beta_K))
    name_df <- paste0("time_importance_", object@id_term)
    object@li_df_var_imp[[name_df]] <- time_importance
    time_importance <- time_importance[order(-time_importance$Overall), ]
    time_importance <- subset(time_importance, Overall > 0.0001)

    image <- ggplot2::ggplot(time_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistic_multiway/importance_time", "_", object@id_term, ".png"), image)

    if (length(best_beta_autre) > 0) {
        other_names <- colnames(object@train_cols[, object@col_x])[object@index_mode == -1]
        other_importance <- data.frame(Variable = other_names, Overall = abs(best_beta_autre))
        name_df <- paste0("other_importance_", object@id_term)
        object@li_df_var_imp[[name_df]] <- other_importance
        other_importance <- other_importance[order(-other_importance$Overall), ]
        other_importance <- subset(other_importance, Overall > 0.0001)

        image <- ggplot2::ggplot(other_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
            geom_bar(stat = "identity") +
            coord_flip() +
            theme_light() +
            xlab("Variable") +
            ylab("Importance") +
            ggtitle("Variable Importance")
        ggsave(paste0("plots/logistic_multiway/importance_other", "_", object@id_term, ".png"), image)
    }


    df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)
    object@li_box_plots[[object@id_term]] <- df_long

    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        geom_boxplot() +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0("plots/logistic_multiway/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)

    return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    return(data.frame(Variable = names(object@data_used[, object@col_x]), Overall = abs(object@model$finalModel$beta_unfolded)))

})

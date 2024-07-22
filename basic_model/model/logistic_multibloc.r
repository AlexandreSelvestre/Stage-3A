li_caret_multiway <- list()

li_caret_multiway$library <- "glmnet"

li_caret_multiway$type <- "Classification"

li_caret_multiway$parameters <- data.frame(parameter = c("lambda", "R"), class = c("numeric", "integer"), label = c(
    "la valeur du paramètre lambda", "le rang R"
))

create_grid_multiway <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        lambda <- log(seq(exp(0), exp(0.05), length.out = len + 1)[2:(len + 1)])
        R <- c(1, 2, 3, 4)
    } else {
        lambda <- log(runif(len, min = exp(0.001), max = exp(0.05)))
        R <- c(1, 2, 3, 4)
    }
    data_frame_grid <- expand.grid(lambda = lambda, R = R)
    return(data_frame_grid)
}

better_create_grid_multiway <- function(x, y, len = NULL, search = "grid", lambda_min = 0.001, lambda_max = 0.05, R_min = 1, R_max = 5, tune_R = 2) {
    if (search == "grid") {
        lambda <- seq(lambda_min, lambda_max, length.out = len)[1:len]
    } else {
        lambda <- runif(len, min = lambda_min, max = lambda_max)
    }
    R <- round(seq(R_min, R_max, length.out = tune_R))
    data_frame_grid <- expand.grid(lambda = lambda, R = R)
    return(data_frame_grid)
}

li_caret_multiway$grid <- create_grid_multiway



find_modes <- function(x, index) {
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
            stop("Valeur -1 anormale durant la recherche de modes!")
        }
    }
    modes <- unique(index)
    modes <- modes[modes > -0.5]
    K <- length(modes)

    return(list(J = J, K = K))
}



fit_multiway <- function(x, y, wts, param, lev, last, weights_dict, classProbs, index, index_bloc, eps, ite_max, n_iter_per_reg, k_smote, do_smote, index_variable, is_binary) {
    # ici, index est bien sûr index_mode
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
    # On prend en entrée un dataframe
    # Formalisme pour les indices: il faut indiquer qui est dans quel mode dans index. On met -1 pour les variables cliniques
    # On supposera que l'ordre est la même pour toutes les variables: pas de décalage. Le tester au début: checker que c'est bien le même nom de colonne
    # L'index_bloc associe à chaque variable multimodale le numéro de son bloc (shape, texture, first order). Il vaut -1 sur les variables cliniques
    weights <- numeric(length(y))
    weights[y == classe_maj] <- weights_dict[[classe_maj]]
    weights[y == classe_min] <- weights_dict[[classe_min]]

    y_numeric <- ifelse(y == classe_min, 1, 0) # CCK donne 1 CHC donne 0
    R <- param$R

    ## Commencer par définir Z_J en sommant sur les modes. Attention au problème lignes colonnes inversées. Pas encore de variables cliniques
    mat_x_modes <- as.matrix(x)[, index_bloc > -0.5]
    mat_x_restant <- as.matrix(x)[, index_bloc == -1]
    n_restant <- ncol(mat_x_restant)

    # Il y a L blocs
    if (all(index_bloc > -0.5)) {
        L <- length(unique(index_bloc))
    } else {
        L <- length(unique(index_bloc)) - 1
    }
    # print(L)
    li_x_multi_bloc <- list()
    col_num <- 0
    for (l in 1:L) {
        li_x_multi_bloc[[l]] <- as.matrix(x)[, index_bloc == l]
        col_num <- col_num + ncol(li_x_multi_bloc[[l]])
    }
    if (col_num != ncol(mat_x_modes)) {
        print(col_num)
        stop("Problème de correspondance des colonnes des blocs!")
    }

    ###### On passe aux modes
    li_dim <- list()
    for (l in 1:L) {
        x_bloc <- li_x_multi_bloc[[l]]
        index_bloc_local <- index[index_bloc == l]
        # print(index_bloc_local)
        # print(length(index_bloc_local))
        # print(dim(x_bloc))
        dim_modes <- unlist(find_modes(x_bloc, index_bloc_local))
        li_dim[[l]] <- dim_modes
    }

    ### Formalisme utilisé pour les matrices : celui du papier


    ####### Petite boucle de l'algorithme: utile dans la grande #########
    petite_boucle <- function(Z_init, vec_Q) {
        Z <- rbind(Z_init, t(mat_x_restant))
        vec_Q <- append(vec_Q, rep(1, n_restant))
        # print("taille vec:")
        # print(length(vec_Q))
        Q_inv <- diag(vec_Q)
        Q <- inverse_diag(Q_inv)
        if (any(is.na(Q))) {
            print(cat("Matrice Q", diag(Q)))
            print(cat("vec_Q en cause", vec_Q))
            stop("Q est nan")
        }
        if (any(is.na(Z))) {
            stop("Z est nan")
        }
        Z <- as.matrix(Z)
        logistic_classic <- glmnet:::glmnet.fit(
            x = t(Q %*% Z), y = y_numeric, family = binomial(), alpha = 1,
            weights = weights / dim(Z)[2], lambda = param$lambda, intercept = TRUE, maxit = 1e7,
            thresh = 1e-8
        )
        faux_beta <- as.numeric(logistic_classic$beta)
        beta <- Q %*% faux_beta
        intercept <- logistic_classic$a0
        return(list(beta = beta, intercept = intercept, Z = Z, Q = Q, Q_inv = Q_inv))
    }



    ####### Grande boucle de l'algorithme à itérer autant que nécessaire ##########
    grande_boucle <- function(li_beta_K) {
        # print("gdbc")
        #   Définir Z_J_init (sans les variables cliniques)
        create_Z_J_init <- function(l) {
            mat_x_bloc <- li_x_multi_bloc[[l]]
            beta_K <- li_beta_K[[l]]
            K <- li_dim[[l]][2]
            J <- li_dim[[l]][1]
            li_Z_J_init_bloc <- list() # init car sans clinique
            for (r in 1:R) {
                Z_J_init_bloc <- matrix(0, nrow = dim(mat_x_bloc)[1], ncol = J)
                for (k in 1:K) {
                    bloc_k <- mat_x_bloc[, ((k - 1) * J + 1):(k * J)]
                    Z_J_init_bloc <- Z_J_init_bloc + beta_K[(r - 1) * K + k] * bloc_k
                    li_Z_J_init_bloc[[r]] <- Z_J_init_bloc
                }
            }
            Z_J_init_bloc <- li_Z_J_init_bloc[[1]]
            if (R >= 2) {
                for (r in 2:R) {
                    Z_J_init_bloc <- cbind(Z_J_init_bloc, li_Z_J_init_bloc[[r]])
                }
            }
            Z_J_init_bloc <- t(Z_J_init_bloc)
            if (any(is.na(Z_J_init_bloc))) {
                # print(Z_J_init_bloc)
                # print(li_beta_K)
                stop(" Z_J_init contient des NA")
            }
            return(Z_J_init_bloc)
        }

        ## Assembler les Z_J_l en un grand Z_J_init
        Z_J_init <- create_Z_J_init(1)
        if (L > 1) {
            for (l in 2:L) {
                Z_J_init <- rbind(Z_J_init, create_Z_J_init(l))
            }
        }
        # print(Z_J_init)


        ## Définir Q_J de transformation pour la régularisation (attention aux colonnes inversées)
        create_vec_diag_J <- function(l) {
            beta_K <- li_beta_K[[l]]
            K <- li_dim[[l]][2]
            J <- li_dim[[l]][1]
            vec_diag_bloc <- c()
            for (r in 1:R) {
                norme_1_r <- norm(as.matrix(beta_K[((r - 1) * K + 1):(r * K)]), type = "1")
                diag_elem_r <- rep(norme_1_r, J)
                vec_diag_bloc <- append(vec_diag_bloc, diag_elem_r)
            }
            return(vec_diag_bloc)
        }

        vec_diag <- create_vec_diag_J(1)
        if (L > 1) {
            for (l in 2:L) {
                vec_diag <- append(vec_diag, create_vec_diag_J(l))
            }
        }

        li_petite_boucle <- petite_boucle(Z_init = Z_J_init, vec_Q = vec_diag)
        beta_J_complet <- li_petite_boucle$beta
        beta_J_full <- li_petite_boucle$beta[1:(length(li_petite_boucle$beta) - n_restant)]
        li_beta_J <- list()
        precedent <- 0
        for (l in 1:L) {
            li_beta_J[[l]] <- beta_J_full[(precedent + 1):(precedent + li_dim[[l]][1] * R)]
            precedent <- precedent + li_dim[[l]][1] * R
        }
        if (n_restant > 0) {
            beta_autre_J <- li_petite_boucle$beta[(length(li_petite_boucle$beta) - n_restant + 1):length(li_petite_boucle$beta)]
        } else {
            beta_autre_J <- c()
        }
        Z_J <- li_petite_boucle$Z
        intercept_J <- li_petite_boucle$intercept
        Q_J_full <- li_petite_boucle$Q
        Q_J_inv <- li_petite_boucle$Q_inv
        # print(beta_J_complet)



        #########################    Opérer maintenant sur Z_K et Q_K

        create_Z_K_init <- function(l) {
            mat_x_bloc <- li_x_multi_bloc[[l]]
            beta_J <- li_beta_J[[l]]
            # print(beta_J)
            J <- li_dim[[l]][1]
            K <- li_dim[[l]][2]
            li_Z_K_init_bloc <- list() # init car sans clinique
            for (r in 1:R) {
                Z_K_init_bloc_r <- matrix(0, nrow = dim(mat_x_bloc)[1], ncol = K)
                for (j in 1:J) {
                    col_to_extract <- seq(j, by = J, length.out = K)
                    bloc_j <- mat_x_bloc[, col_to_extract]
                    # print(beta_J[(r - 1) * J + j])
                    Z_K_init_bloc_r <- Z_K_init_bloc_r + beta_J[(r - 1) * J + j] * bloc_j
                    li_Z_K_init_bloc[[r]] <- Z_K_init_bloc_r
                }
                # if (r == 2) {
                #     print(Z_K_init_bloc_r)
                # }
            }
            Z_K_init_bloc <- li_Z_K_init_bloc[[1]]
            if (R >= 2) {
                for (r in 2:R) {
                    Z_K_init_bloc <- cbind(Z_K_init_bloc, li_Z_K_init_bloc[[r]])
                }
            }
            Z_K_init_bloc <- t(Z_K_init_bloc)
            # print("done")
            # print(Z_K_init_bloc)
            if (all(abs(Z_K_init_bloc) < 1e-10)) {
                # print(beta_J)
                print("Z_K est nul : tout sera écrasé à 0")
            }
            return(Z_K_init_bloc)
        }

        Z_K_init <- create_Z_K_init(1)
        if (L > 1) {
            for (l in 2:L) {
                Z_K_init <- rbind(Z_K_init, create_Z_K_init(l))
            }
        }

        ## Définir Q_K de transformation pour la régularisation (attention aux colonnes inversées)

        create_vec_diag_K <- function(l) {
            beta_J <- li_beta_J[[l]]
            J <- li_dim[[l]][1]
            K <- li_dim[[l]][2]
            vec_diag_bloc <- c()
            for (r in 1:R) {
                norme_1_r <- norm(as.matrix(beta_J[((r - 1) * J + 1):(r * J)]), type = "1")
                diag_elem_r <- rep(norme_1_r, K)
                vec_diag_bloc <- append(vec_diag_bloc, diag_elem_r)
            }
            return(vec_diag_bloc)
        }


        vec_diag <- create_vec_diag_K(1)
        if (L > 1) {
            for (l in 2:L) {
                vec_diag <- append(vec_diag, create_vec_diag_K(l))
            }
        }


        ### Appliquer à nouveau la petite boucle

        # print(beta_K_complet_previous)
        # print(dim(Z_K_init))
        li_petite_boucle <- petite_boucle(Z_init = Z_K_init, vec_Q = vec_diag)
        Z_K <- li_petite_boucle$Z
        beta_K_full <- li_petite_boucle$beta[1:(length(li_petite_boucle$beta) - n_restant)]
        li_beta_K <- list()
        precedent <- 0
        for (l in 1:L) {
            li_beta_K[[l]] <- beta_K_full[(precedent + 1):(precedent + li_dim[[l]][2] * R)]
            precedent <- precedent + li_dim[[l]][2] * R
        }

        ## Renormaliser
        for (l in 1:L) {
            beta_K <- li_beta_K[[l]]
            K <- li_dim[[l]][2]
            for (r in 1:R) {
                norm_beta_K_r <- norm(as.matrix(beta_K[((r - 1) * K + 1):(r * K)]), type = "2")
                if (norm_beta_K_r > 0) {
                    beta_K[((r - 1) * K + 1):(r * K)] <- beta_K[((r - 1) * K + 1):(r * K)] / norm_beta_K_r
                }
                li_beta_K[[l]] <- beta_K
            }
        }
        intercept_K <- li_petite_boucle$intercept
        if (n_restant > 0) {
            beta_autre_K <- li_petite_boucle$beta[(length(li_petite_boucle$beta) - n_restant + 1):length(li_petite_boucle$beta)]
        } else {
            beta_autre_K <- c()
        }
        Q_K_full <- li_petite_boucle$Q
        Q_K_inv <- li_petite_boucle$Q_inv
        beta_K_complet <- li_petite_boucle$beta

        # Calculer les deux critères
        # print(beta_J_complet)
        # print(beta_K_complet)
        crit_log_J <- crit_logistic(t(Q_J_full %*% Z_J), y_numeric, Q_J_inv %*% beta_J_complet, intercept_J, param$lambda)
        crit_log_K <- crit_logistic(t(Q_K_full %*% Z_K), y_numeric, Q_K_inv %*% beta_K_complet, intercept_K, param$lambda)
        rapport <- abs(crit_log_J - crit_log_K) / abs(crit_log_J)
        return(list(
            li_beta_J = li_beta_J, li_beta_K = li_beta_K, beta_autre = beta_autre_K, intercept_J = intercept_J, intercept_K = intercept_K,
            rapport = rapport, crit_log_J = crit_log_J, crit_log_K = crit_log_K, beta_J_complet = beta_J_complet, beta_K_complet = beta_K_complet
        ))
    }

    ##### Exécuter la grande boucle autant que les critères le demandent
    continue <- TRUE
    ###### Attention beta_J et beta_K n'existent plus, ce sont des listes!! (plus simple à manipuler avec les dimensions variables)
    li_beta_K <- list()
    for (l in 1:L) {
        li_beta_K[[l]] <- rnorm(li_dim[[l]][2] * R, mean = 0, sd = 1)
    }
    iteration <- 1
    memoire_crit_J <- 1
    memoire_crit_K <- 1
    debut_time <- Sys.time()
    while (continue) {
        # print(iteration)
        li_grande_boucle <- grande_boucle(li_beta_K)
        li_beta_J <- li_grande_boucle$li_beta_J
        li_beta_K <- li_grande_boucle$li_beta_K
        beta_autre <- li_grande_boucle$beta_autre
        intercept_J <- li_grande_boucle$intercept_J
        intercept_K <- li_grande_boucle$intercept_K
        # print(li_beta_K[[0]])

        rapport <- li_grande_boucle$rapport
        crit_log_J <- li_grande_boucle$crit_log_J
        crit_log_K <- li_grande_boucle$crit_log_K
        # print(crit_log_J)
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
    futur_fit <- list(
        li_beta_J = li_beta_J, li_beta_K = li_beta_K, beta_autre = beta_autre, intercept = intercept_K, R = R, li_x_multi_bloc = li_x_multi_bloc,
        li_dim = li_dim, lev = lev, index = index, index_bloc = index_bloc, li_norm = li_norm, classe_maj = classe_maj, classe_min = classe_min
    )

    futur_fit$beta_unfolded <- get_beta_full(futur_fit)
    return(futur_fit)
}

li_caret_multiway$fit <- fit_multiway


get_beta_bloc <- function(beta_J, beta_K, R, J, K) {
    beta <- rep(0, J * K)
    for (r in 1:R) {
        beta_r <- rep(0, J * K)
        for (j in 1:J) {
            for (k in 1:K) {
                beta_r[(k - 1) * J + j] <- beta_J[(r - 1) * J + j] * beta_K[(r - 1) * K + k]
            }
        }
        beta[1:(J * K)] <- beta[1:(J * K)] + beta_r
    }
    return(beta)
}


get_beta_full <- function(modelFit) {
    li_x_multi_bloc <- modelFit$li_x_multi_bloc
    li_dim <- modelFit$li_dim
    index <- modelFit$index
    index_bloc <- modelFit$index_bloc
    li_beta_J <- modelFit$li_beta_J
    li_beta_K <- modelFit$li_beta_K
    beta_autre <- modelFit$beta_autre
    L <- length(li_x_multi_bloc)
    R <- modelFit$R
    size_beta_modes <- sum(unlist(lapply(li_dim, function(x) {
        return(x[1] * x[2])
    }))) # taille de beta sans beta_autre
    beta_modes <- rep(NA, size_beta_modes)
    for (l in 1:L) {
        beta_K <- li_beta_K[[l]]
        beta_J <- li_beta_J[[l]]
        J <- li_dim[[l]][1]
        K <- li_dim[[l]][2]
        beta_bloc <- get_beta_bloc(beta_J, beta_K, R, J, K)
        beta_modes[index_bloc[index_bloc != -1] == l] <- beta_bloc
        # Attention ordre temps
    }
    beta_final <- c(beta_modes, beta_autre)
    if (any(is.na(beta_final))) {
        print(beta_final)
        stop("Beta final est NA")
    }
    return(beta_final)
}



predict_multiway <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- get_beta_full(modelFit)
    intercept <- modelFit$intercept
    lev <- modelFit$lev
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta + intercept
    })
    proba <- 1 / (1 + exp(-value))
    predicted_labels <- ifelse(proba > 0.5, classe_min, classe_maj)
    # print(predicted_labels)
    return(predicted_labels)
}

li_caret_multiway$predict <- predict_multiway


li_caret_multiway$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- get_beta_full(modelFit)
    modelFit$beta_unfolded <- beta
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
        search = object@search, lambda_min = object@lambda_min, lambda_max = object@lambda_max, R_min = object@R_min, R_max = object@R_max, tune_R = object@tune_R
    )

    # dossier <- "logs"
    # fichiers <- list.files(dossier, full.names = TRUE)
    # file.remove(fichiers)


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
        clusterExport(cl, varlist = c("get_beta_full", "get_beta_bloc", "find_modes"))
    }


    # if (object@do_parallel) {
    #     numCores <- detectCores()
    #     cl <- makePSOCKcluster(numCores - 1)
    #     registerDoParallel(cl)
    #     clusterEvalQ(cl, source("utils.r"))
    #     clusterExport(cl, varlist = c("get_beta_full", "get_beta_bloc", "find_modes"))
    # }

    object@model <- caret::train(
        y = object@y_train, x = object@train_cols[, object@col_x], index = object@index_mode,
        method = li_caret_multiway, trControl = object@cv, metric = "AUC",
        tuneLength = object@tuneLength, weights_dict = object@weights, tuneGrid = grid, eps = object@eps, ite_max = object@ite_max, n_iter_per_reg = object@n_iter_per_reg,
        index_bloc = object@index_bloc, k_smote = object@k_smote, do_smote = object@do_smote, index_variable = object@index_variable, is_binary = object@is_binary
    )
    if (object@do_parallel) {
        stopCluster(cl)
    }
    return(object)
})

setMethod("get_results", "apply_model", function(object) {
    object@predictions <- as.vector(predict(object@model, newdata = as.matrix(object@test_set[, object@col_x])))
    object@predictions_proba <- predict(object@model, newdata = as.matrix(object@test_set[, object@col_x]), type = "prob")
    object@predictions_train_proba <- predict(object@model, newdata = as.matrix(object@train_cols[, object@col_x]), type = "prob")
    return(object)
})

setMethod("importance_method", "apply_model", function(object) {
    object@beta_final <- object@model$finalModel$beta_unfolded
    # print(object@model$finalModel$li_beta_K)
    li_best_beta_J <- object@model$finalModel$li_beta_J
    li_best_beta_K <- object@model$finalModel$li_beta_K
    best_beta_autre <- object@model$finalModel$beta_autre
    L <- length(li_best_beta_J)
    Variables_temps_seul <- unique(object@name_mode[object@index_mode > -0.5])
    Variables_names_seul <- unique(object@name_variable[object@index_variable > -0.5])
    n_var_radio_differentes <- length(Variables_names_seul)
    Variables_names_long <- c()
    Variables_temps <- c()
    R <- object@model$bestTune$R
    for (r in 1:R) {
        Variables_names_long <- append(Variables_names_long, paste0(Variables_names_seul, "_R=", r))
        Variables_temps <- append(Variables_temps, paste0(Variables_temps_seul, "_R=", r))
    }
    # print(Variables_names_seul)
    # print(Variables_names_long)
    li_variables_names <- list()
    for (l in 1:L) {
        li_variables_names[[l]] <- c()
    }
    li_variables_names <- vector("list", L)
    count <- 0
    for (name in Variables_names_long) {
        count <- count %% n_var_radio_differentes + 1
        num_bloc <- object@index_bloc[count]
        # print(length(li_variables_names))
        # print(count)
        li_variables_names[[num_bloc]] <- append(li_variables_names[[num_bloc]], name)
        count <- count + 1
    }
    grand_vec_temps <- c()
    grand_var_temps <- c()
    for (l in 1:L) {
        grand_vec_temps <- append(grand_vec_temps, li_best_beta_K[[l]])
        grand_var_temps <- append(grand_var_temps, paste0(Variables_temps, "_l=", l))
    }

    time_importance <- data.frame(Variable = grand_var_temps, Overall = abs(grand_vec_temps))
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
    ggsave(paste0("plots/logistic_multibloc/importance_time", "_", object@id_term, ".png"), image)

    grand_vec_var <- c()
    grand_var_var <- c()
    for (l in 1:L) {
        grand_vec_var <- append(grand_vec_var, li_best_beta_J[[l]])
        grand_var_var <- append(grand_var_var, paste0(li_variables_names[[l]], "_l = ", l))
        # print(li_best_beta_J[[l]])
        # print(paste0(li_variables_names[[l]], "_l = ", l))
    }

    variable_importance <- data.frame(Variable = grand_var_var, Overall = abs(grand_vec_var))
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
    ggsave(paste0("plots/logistic_multibloc/importance_var_multi", "_", object@id_term, ".png"), image)



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
        ggsave(paste0("plots/logistic_multibloc/importance_other", "_", object@id_term, ".png"), image)
    }

    df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)
    object@li_box_plots[[object@id_term]] <- df_long

    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        geom_boxplot() +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0("plots/logistic_multibloc/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)

    return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    return(data.frame(Variable = names(object@data_used[, object@col_x]), Overall = abs(object@model$finalModel$beta_unfolded)))
})


# lambda opti pour n = 1000: 0.00135

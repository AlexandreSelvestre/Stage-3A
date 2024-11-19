setClass("logistic_multibloc",
    contains = "apply_model",
    slots = representation(
        lambda_min = "numeric",
        lambda_max = "numeric",
        ite_max = "ANY",
        eps = "numeric",
        tuneLength = "numeric",
        index_type = "character",
        tune_R = "integer",
        R_min = "integer",
        R_max = "integer",
        weights = "list",
        n_iter_per_reg = "integer",
        regression = "logical",
        same_R = "logical",
        li_R = "list",
        lambda = "numeric",
        max_temps = "numeric"
    ),
    prototype = prototype(
        index_mode = c(),
        ite_max = NULL
    )
)


li_caret_multibloc <- list()

li_caret_multibloc$library <- "glmnet"

li_caret_multibloc$type <- "Classification"

create_param_for_caret_bloc <- function(L) {
    vec_param <- sapply(1:L, function(l) {
        return(paste0("R_", l))
    })
    vec_param <- c(vec_param, "lambda")
    vec_label <- sapply(1:L, function(l) {
        return(paste0("le rang R du bloc ", l))
    })
    vec_label <- c(vec_label, "la valeur du paramètre lambda")
    return(data.frame(parameter = vec_param, class = rep("numeric", L + 1), label = vec_label))
}

# False grid dans utils_model_specific

better_create_grid_multibloc <- function(x, y, len = NULL, search = "grid", L, lambda_min = 0.001, lambda_max = 0.05, R_min = 1, R_max = 5, tune_R = 2, li_R = NULL, same_R) {
    ### Attention, si li_R est renseigné et same_R est TRUE alors il faut avoir autant de rangs proposés par bloc (ils seront sur 1 même ligne)
    ## tune_R est ignoré si on fournit li_R
    if (search == "grid") {
        lambda <- exp(log(10) * seq(log10(lambda_min), log10(lambda_max), length.out = len))
    } else {
        lambda <- exp(log(10) * runif(len, min = log10(lambda_min), max = log10(lambda_max)))
    }
    if (length(li_R) == 0) {
        li_R_data_frame <- lapply(1:L, function(l) round(seq(R_min, R_max, length.out = tune_R)))
    } else {
        li_R_data_frame <- li_R
    }
    if (same_R) {
        R_data_frame <- as.data.frame(li_R_data_frame)
    } else {
        R_data_frame <- expand.grid(li_R_data_frame)
    }
    li_grid <- lapply(seq_along(lambda), function(i) {
        df <- data.frame(R_data_frame, lambda = lambda[i])
        return(df)
    })
    data_frame_grid <- do.call(rbind, li_grid)
    vec_names <- sapply(1:L, function(l) {
        return(paste0("R_", l))
    })
    vecnames <- c(vec_names, "lambda")
    data_frame_grid <- setNames(data_frame_grid, vecnames)
    # print(as.data.frame(data_frame_grid))
    return(as.data.frame(data_frame_grid))
}


fit_multiway <- function(x, y, wts, param, lev, last, weights_dict, classProbs, li_index_modes, index_bloc, eps, ite_max, n_iter_per_reg, k_smote, sampling_choice, index_variable, is_binary, classe_1 = NULL, max_temps) {
    li_norm <- renormalize_in_model_fit_index_mode(x, index_variable, index_bloc, is_binary)


    x <- li_norm$new_x
    classe_min <- names(which.min(table(y)))
    classe_maj <- setdiff(levels(y), classe_min)

    if (sampling_choice == "smote") {
        li <- apply_smote(x, y, k_smote)
        x <- li$x
        y <- li$y
    }
    if (sampling_choice == "up") {
        li <- apply_boot(x, y)
        x <- li$x
        y <- li$y
    }


    weights <- numeric(length(y))
    weights[y == classe_maj] <- weights_dict[[classe_maj]]
    weights[y == classe_min] <- weights_dict[[classe_min]]

    y_numeric <- convert_y(y, classe_1)
    if (is.null(classe_1)) {
        classe_1 <- classe_min
    }
    classe_0 <- setdiff(levels(y), classe_1)

    # on signalera par un -1 dans chaque index_mode, l'absence de lien de la variable avec le mode en question (soit car numéro de mode trop élevé, soit car variable tabulaire)
    # Index bloc est fiable pour indiquer les variables tabulaires. index_variable aussi. Attention: index_variable a le droit d'être local (propre à chaque bloc) aussi bien que global (commun à toute la liste de tenseurs)


    ## On changera ici la terminologie: slice pour l'indexation au sein d'un mode. Désormais, on désignera par mode une voie entière, i.e. une direction de parcours des tenseurs. Certains tenseurs en ont plus que d'autres. (tout comme, pour un même numéro de mode, certains tenseurs auront plus de slices que d'autres: pas de réelle continuité physique pour un même numéro de mode entre plusieurs tenseurs)

    # li_index_modes contient tous les index_modes les uns à la suite des autres (un index_mode par mode). Au cours de la fonction, on lui adjoindra index_variable à la fin. Attention, ne pas le faire soi-même, c'est fait automatiquement

    ######### Initialisation spécifique au modèle multibloc ##########

    different_blocs <- sort(unique(index_bloc[index_bloc != -1]))
    L <- length(different_blocs)

    # if (length(li_index_modes) == 0) {
    #     li_index_modes$single_mode <- index_mode
    # }

    # li_index_modes$variable <- index_variable


    # le numéro de variable est une slice comme une autre... et ce sera le dernier. On fera une passe du dernier au premier mode pour la fonction boucle_mode

    M_max <- length(li_index_modes)

    for (l_num in different_blocs) {
        for (m in seq_len(M_max)) {
            index_m_l <- li_index_modes[[m]][index_bloc == l_num]
            if (length(index_m_l[index_m_l > -0.5]) > 0) {
                li_index_modes[[m]][index_bloc == l_num][index_m_l > -0.5] <- contiguous_ranks(index_m_l[index_m_l > -0.5])
            }
        }
    }


    li_different_slices_per_mode <- lapply(seq_len(M_max), function(m) {
        li_different_slices <- lapply(different_blocs, function(l_num) {
            index_slices <- li_index_modes[[m]]
            index_slice_local <- index_slices[index_bloc == l_num]
            index_slice_local <- index_slice_local[index_slice_local > -0.5]
            different_slices <- sort(unique(index_slice_local))
            return(different_slices)
        })
        return(li_different_slices)
    }) # Liste par mode de liste par bloc indiquant les numéros de slices pour ce mode et ce bloc. Si le bloc n'est pas concerné par le mode, la case vaut c()

    for (m in seq_along(li_different_slices_per_mode)) {
        names(li_different_slices_per_mode[[m]]) <- as.character(different_blocs)
    }

    mat_x_tens <- as.matrix(x)[, index_bloc > -0.5]
    mat_x_tab <- as.matrix(x)[, index_bloc == -1]
    n_tab <- ncol(mat_x_tab)

    current_li_R <- lapply(seq_len(L), function(l) {
        R <- param[[paste0("R_", l)]]
        return(R)
    }) # Liste des rang R de chaque bloc à cette étape de la cv
    names(current_li_R) <- as.character(different_blocs)
    li_dim <- lapply(different_blocs, function(l_num) {
        l_char <- as.character(l_num)
        vec_dim <- sapply(seq_len(M_max), function(m) {
            li_different_slices <- li_different_slices_per_mode[[m]][[l_char]]
            return(length(li_different_slices))
        })
        return(vec_dim)
    })
    names(li_dim) <- as.character(different_blocs)

    li_least_slices_mode_per_bloc <- lapply(different_blocs, function(l_num) {
        l_char <- as.character(l_num)
        indices_non_zero <- which(li_dim[[l_char]] > 0)
        vec_dim_contracted <- li_dim[[l_char]][indices_non_zero]
        simplest_mode_contracted <- which.min(vec_dim_contracted)
        simplest_mode <- indices_non_zero[simplest_mode_contracted]
        return(simplest_mode)
    })
    names(li_least_slices_mode_per_bloc) <- as.character(different_blocs) # util pour compléter les modes manquants dans certains blocs sans exploser les temps de calcul
    li_x_multi_bloc_pos <- list()
    col_num <- 0
    # Hyper important: on crée une liste avec pour chaque bloc un array représentant le tenseur avec les colonnes dans le bon ordre (lexicographique vs indices tensoriels)

    for (l_num in different_blocs) {
        l_char <- as.character(l_num)
        li_x_multi_bloc_pos[[l_char]] <- reorder_local(as.matrix(x)[, index_bloc == l_num], li_index_modes, li_dim[[l_char]], which(index_bloc == l_num))
        if (any(is.na(li_x_multi_bloc_pos[[l_char]]$mat))) {
            print(l_num)
            stop("Echec NA reorder")
        }
        col_num <- col_num + ncol(li_x_multi_bloc_pos[[l_char]]$mat) # Sert seulement pour le test
    }
    if (col_num != ncol(mat_x_tens)) {
        print(col_num)
        stop("Problème de correspondance des colonnes des blocs!")
    }

    # Rappel: les variables sont à la fin de li_dim: l'ordre des dimensions est celui des modes dans li_index_modes


    ### Fonction qui appelle glmnet.fit quand on lui donne les bonnes matrices Z et Q (Q est donné sous forme de vecteur diagonal). Cette fonction se charge de remultiplier par Q pour obtenir le beta final. Attention, elle ne rétablit PAS les changements d'intercept!

    function_fit <- function(Z, vec_Q) {
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
            x = Z %*% Q, y = y_numeric, family = binomial(), alpha = 1,
            weights = weights / dim(Z)[1], lambda = param$lambda, intercept = TRUE, maxit = 1e8,
            thresh = 1e-9
        )
        faux_beta <- as.numeric(logistic_classic$beta)
        intercept <- logistic_classic$a0
        converged <- logistic_classic$converged


        # logistic_classic <- penalized::penalized(response = y_numeric, penalized = Z %*%Q, lambda1 = param$lambda, lambda2 = 0, model = "logistic", epsilon = 10^(-8), maxiter = 1000, standardize = FALSE, trace = FALSE)
        # faux_beta <- logistic_classic@penalized
        # intercept <- logistic_classic@unpenalized
        # print(logistic_classic@converged)

        beta <- Q %*% faux_beta
        return(list(beta = beta, intercept = intercept, Z = Z, Q = Q, Q_inv = Q_inv, converged = converged))
    }

    construct_Z_m <- function(m_initial, li_beta_modes_per_bloc) {
        # Cette fonction effectue tout le travail de création de Z pour un m donné

        li_Z_m_tens <- lapply(different_blocs, function(l_num) {
            l_char <- as.character(l_num)
            basic_Z_m_l_construction <- function(m) {
                # Fonction de base pour construire Z_m_l sur un mode donné (qui n'est potentiellement pas le m de départ si le mode n'existe pas dans le bloc!)
                R <- current_li_R[[l_char]]
                different_slices <- li_different_slices_per_mode[[m]][[l_char]]
                potential_zeros <- li_dim[[l_char]][1:m]
                num_zeros <- length(potential_zeros[potential_zeros == 0])
                m_local <- m - num_zeros # position du mode m parmi les vrais modes du bloc l
                vec_dim_l <- li_dim[[l_char]]
                vec_dim_l_contracted <- vec_dim_l[vec_dim_l > 0] # sans les zéros
                M_l <- length(vec_dim_l_contracted)
                li_index_sum <- lapply(setdiff(seq_len(M_l), m_local), function(m_prime) {
                    return(seq_len(vec_dim_l_contracted[m_prime]))
                })
                cart_prod_index_sum <- as.matrix(expand.grid(li_index_sum))
                vec_base_mode <- li_x_multi_bloc_pos[[l_char]]$vec_base_mode
                x_l_ordered <- li_x_multi_bloc_pos[[l_char]]$mat
                Z_m_l <- matrix(0, nrow = nrow(mat_x_tens), ncol = R * length(different_slices))
                for (a in seq_len(nrow(cart_prod_index_sum))) {
                    index_col_normal_format <- sapply(
                        seq_len(vec_dim_l_contracted[m_local]),
                        function(k_m) {
                            num_col_tens_format <- append(as.vector(cart_prod_index_sum[a, ]), k_m, after = m_local - 1) - 1
                            num_col_normal_format <- sum(vec_base_mode * num_col_tens_format) + 1
                            return(num_col_normal_format)
                        }
                    )
                    bloc_m <- x_l_ordered[, index_col_normal_format]
                    li_beta_contracted <- li_beta_modes_per_bloc[[l_char]][vec_dim_l > 0]
                    for (r in seq_len(R)) {
                        index_sum_local <- append(as.vector(cart_prod_index_sum[a, ]), NA, after = m_local - 1) # k_m est ajouté comme un "dummy": ne sera pas considéré: on le met à NA (pour générer une erreur en cas de pb)
                        vec_prod <- sapply(setdiff(seq_len(M_l), m_local), function(m_prime) {
                            return(li_beta_contracted[[m_prime]][(r - 1) * vec_dim_l_contracted[m_prime] + index_sum_local[m_prime]])
                        })
                        prod_beta <- prod(vec_prod)
                        Z_m_l_r <- bloc_m * prod_beta
                        if (any(is.na(Z_m_l_r))) {
                            # print(Z_m_l_r)
                            # print(r)
                            # print(l_num)
                            # print(index_col_normal_format)
                            # print(x_l_ordered[,163])
                            # print(bloc_m)
                            stop("NA")
                        }
                        # print(length(different_slices))
                        # print(m)
                        # print(l_num)
                        # print(dim(Z_m_l_r))
                        # print(dim(Z_m_l))
                        Z_m_l[, ((r - 1) * length(different_slices) + 1):(r * length(different_slices))] <- Z_m_l[, ((r - 1) * length(different_slices) + 1):(r * length(different_slices))] + Z_m_l_r
                    }
                }
                return(Z_m_l)
            }


            if (li_dim[[l_char]][m_initial] == 0) {
                Z_m_l <- basic_Z_m_l_construction(li_least_slices_mode_per_bloc[[l_char]])
            } else {
                Z_m_l <- basic_Z_m_l_construction(m_initial)
            }
            return(Z_m_l)
        })
        Z_m_tens <- do.call(cbind, li_Z_m_tens)

        # Gérer maintenant les variables tabulaires
        Z_m <- cbind(Z_m_tens, mat_x_tab)
        return(Z_m)
    }

    construct_vec_Q_m <- function(m_initial, li_beta_modes_per_bloc) {
        # Construit le vecteur vec_Q_m (pour le mode m)
        li_veq_Q_m <- lapply(different_blocs, function(l_num) {
            l_char <- as.character(l_num)
            R <- current_li_R[[l_char]]


            basic_Q_m_l_r_construction <- function(m, r) {
                # Sous fonction de base pour construire Q_m_l_r
                K_m_l <- li_dim[[l_char]][m]
                index_beta_present <- which(li_dim[[l_char]] > 0)
                vec_prod <- sapply(index_beta_present, function(m_prime) {
                    if (m_prime == m) {
                        return(1)
                    } else {
                        beta_mprime_l_r <- li_beta_modes_per_bloc[[l_char]][[m_prime]][((r - 1) * li_dim[[l_char]][m_prime] + 1):(r * li_dim[[l_char]][m_prime])]
                        return(norm(as.matrix(beta_mprime_l_r), type = "1"))
                    }
                })
                product_norm_beta <- prod(vec_prod)
                vec_Q_m_l_r <- rep(product_norm_beta, K_m_l)
                return(vec_Q_m_l_r)
            }


            li_vec_Q_m_l <- lapply(seq_len(R), function(r) {
                if (li_dim[[l_char]][m_initial] == 0) {
                    vec_Q_m_l_r <- basic_Q_m_l_r_construction(li_least_slices_mode_per_bloc[[l_char]], r)
                } else {
                    vec_Q_m_l_r <- basic_Q_m_l_r_construction(m_initial, r)
                }
                return(vec_Q_m_l_r)
            })
            veq_Q_m_l <- do.call(c, li_vec_Q_m_l)
            return(veq_Q_m_l)
        })
        vec_Q_m <- do.call(c, li_veq_Q_m)
        # Considérer les variables tabulaires
        vec_Q_m <- c(vec_Q_m, rep(1, n_tab))
        return(vec_Q_m)
    }


    update_li_beta_modes_per_bloc <- function(li_beta_modes_per_bloc, m, new_beta) {
        # Réalise le bon update, en tenant compte de quel mode a réellement été updaté dans chaque bloc. Ne pas considérer ici l'intercept (qui sera géré au niveau global)
        li_true_mode <- lapply(different_blocs, function(l_num) {
            l_char <- as.character(l_num)
            if (li_dim[[l_char]][m] == 0) {
                return(li_least_slices_mode_per_bloc[[l_char]])
            } else {
                return(m)
            }
        })
        names(li_true_mode) <- as.character(different_blocs)
        offset <- 0
        for (l_num in different_blocs) {
            l_char <- as.character(l_num)
            true_mode <- li_true_mode[[l_char]]
            R <- current_li_R[[l_char]]
            K_m_l <- li_dim[[l_char]][true_mode]
            new_beta_l_true_mode <- new_beta[(offset + 1):(offset + K_m_l * R)]
            li_beta_modes_per_bloc[[l_char]][[true_mode]] <- new_beta_l_true_mode
            offset <- offset + K_m_l * R
        }
        li_beta_modes_per_bloc$tab <- new_beta[(offset + 1):(offset + n_tab)]
        return(li_beta_modes_per_bloc)
    }



    reconstruct_beta_multibloc <- function(li_beta_modes_per_bloc) {
        # Reconstruit le beta complet à partir des beta de chaque mode. La reconstruction remet beta dans l'ordre initial des colonnes du tableau d'entrée (pas d'ordre lexicographique ici)
        beta_unfolded <- rep(0, ncol(x))
        for (j in seq_len(ncol(x))) {
            l_num <- index_bloc[j]
            l_char <- as.character(l_num)
            if (l_num == -1) {
                beta_unfolded[j] <- li_beta_modes_per_bloc$tab[j - ncol(mat_x_tens)]
                ### Les tabulaires doivent être à la fin: seule contrainte!!!
            } else {
                R <- current_li_R[[l_char]]
                for (r in seq_len(R)) {
                    vec_beta_l_r <- sapply(seq_len(M_max), function(m) {
                        li_beta_m_l <- li_beta_modes_per_bloc[[l_char]][[m]]
                        slice_num <- li_index_modes[[m]][j]
                        if (slice_num == -1) {
                            return(1)
                        } else {
                            return(li_beta_m_l[(r - 1) * li_dim[[l_char]][m] + slice_num])
                        }
                    })
                    beta_unfolded[j] <- beta_unfolded[j] + prod(vec_beta_l_r)
                }
            }
        }
        return(beta_unfolded)
    }

    initialize_li_beta_modes_per_bloc <- function() {
        li_beta_modes_per_bloc <- lapply(different_blocs, function(l_num) {
            l_char <- as.character(l_num)
            lapply(seq_len(M_max), function(m) {
                if (li_dim[[l_char]][m] == 0) {
                    return(c())
                } else {
                    return(rnorm(li_dim[[l_char]][m] * current_li_R[[l_char]], mean = 0, sd = 1))
                }
            })
        })
        names(li_beta_modes_per_bloc) <- as.character(different_blocs)
        return(li_beta_modes_per_bloc)
    }




    li_beta_modes_per_bloc <- initialize_li_beta_modes_per_bloc()
    continue <- TRUE
    vec_crit_memoire <- rep(Inf, M_max)
    vec_crit <- rep(Inf, M_max)
    n_ite <- 0
    debut_time <- Sys.time()

    while (continue) {
        vec_crit_memoire <- data.table::copy(vec_crit)
        for (m in seq_len(M_max)) {
            Z_m <- construct_Z_m(m, li_beta_modes_per_bloc)
            vec_Q_m <- construct_vec_Q_m(m, li_beta_modes_per_bloc)
            li_fit <- function_fit(Z_m, vec_Q_m)
            new_beta <- li_fit$beta
            intercept <- li_fit$intercept
            li_beta_modes_per_bloc <- update_li_beta_modes_per_bloc(li_beta_modes_per_bloc, m, new_beta)
            beta_unfolded <- reconstruct_beta_multibloc(li_beta_modes_per_bloc)
            crit <- crit_logistic(as.matrix(x), y_numeric, beta_unfolded, intercept, param$lambda)
            vec_crit[m] <- crit

            Q <- li_fit$Q
            Q_inv <- li_fit$Q_inv
            converged <- li_fit$converged
            crit_compare <- crit_logistic(Z_m %*% Q, y_numeric, Q_inv %*% new_beta, intercept, param$lambda)
            print(converged)
            print(crit)
            print(paste("crit_compare", crit_compare))
        }
        n_ite <- n_ite + 1
        rapport <- abs(vec_crit[M_max] - vec_crit[1]) / abs(vec_crit[1])
        if (is.na(rapport)) {
            print("rapport NA")
            rapport <- 2 * eps
        }
        if (rapport < eps) {
            continue <- FALSE
        }
        difference <- sum(abs(vec_crit_memoire - vec_crit))
        if (is.na(difference)) {
            print("difference NA")
            difference <- M_max
        }
        if (difference < 1e-10 * M_max) {
            continue <- FALSE
            print("Warning, LOOP!!")
            print(paste("Le rapport vaut:", rapport))
        }
        for (m in seq_len(M_max - 1)) {
            if (vec_crit[m] > vec_crit[m + 1]) {
                continue <- FALSE
                print("Recroissance de la fonction de coût")
            }
        }

        delta_time <- as.numeric(Sys.time() - debut_time, units = "secs")
        if (delta_time > max_temps) {
            print("Warning, temps de calcul trop long")
            print(paste("Le rapport vaut:", rapport))
            continue <- FALSE
        }
        if (n_ite > n_iter_per_reg) {
            continue <- FALSE
            print("Warning, nombre d'itérations trop élevé")
            print(paste("Le rapport vaut:", rapport))
        }
    }

    li_fit <- list(
        li_beta_modes_per_bloc = li_beta_modes_per_bloc, beta_unfolded = beta_unfolded, intercept = intercept, current_li_R = current_li_R,
        li_dim = li_dim, lev = lev, index_bloc = index_bloc, li_norm = li_norm, classe_maj = classe_maj, classe_min = classe_min, classe_1 = classe_1, classe_0 = classe_0
    )
    return(li_fit)
}


li_caret_multibloc$fit <- fit_multiway


predict_multiway <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- modelFit$beta_unfolded
    intercept <- modelFit$intercept
    lev <- modelFit$lev
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta + intercept
    })
    proba <- 1 / (1 + exp(-value))
    predicted_labels <- ifelse(proba > 0.5, classe_1, classe_0)
    return(predicted_labels)
}

li_caret_multibloc$predict <- predict_multiway


li_caret_multibloc$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta <- modelFit$beta_unfolded
    modelFit$beta_unfolded <- beta
    intercept <- modelFit$intercept
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta + intercept
    })
    proba_class_1 <- 1 / (1 + exp(-value))
    proba_class_0 <- 1 - proba_class_1
    str_1 <- as.character(classe_1)
    str_0 <- as.character(classe_0)
    return(setNames(data.frame(proba_class_0, proba_class_1), c(str_0, str_1)))
}

li_caret_multibloc$loop <- NULL


setMethod("train_method", "apply_model", function(object) {
    # Créer la liste de paramètres de rang
    L <- length(unique(object@index_bloc[object@index_bloc > -0.5]))
    li_caret_multibloc$parameters <- create_param_for_caret_bloc(L)
    # Générer la fausse fonction qui crée la tune grid
    create_grid_multibloc <- def_false_grid_multibloc(L)
    li_caret_multibloc$grid <- create_grid_multibloc


    grid <- better_create_grid_multibloc(
        x = object@train_cols[, object@col_x], y = object@y_train, len = object@tuneLength,
        search = object@search, L = L, lambda_min = object@lambda_min, lambda_max = object@lambda_max, R_min = object@R_min, R_max = object@R_max, tune_R = object@tune_R, li_R = object@li_R, same_R = object@same_R
    )

    if (!object@use_li_index_modes) {
        object@li_index_modes$single_mode <- object@index_mode
        object@li_name_modes$single_mode <- object@name_mode
    }
    object@li_index_modes$variable <- object@index_variable
    object@li_name_modes$variable <- object@name_variable

    if (object@parallel$do) {
        numCores <- detectCores()
        if (object@parallel$forking) {
            cl <- makeForkCluster(object@parallel$n_process)
            registerDoParallel(cl)
        } else {
            cl <- makePSOCKcluster(object@parallel$n_process)
            registerDoParallel(cl)
            clusterEvalQ(cl, {
                files <- list.files("./utils", full.names = TRUE, pattern = "\\.r$")
                for (file in files) {
                    source(file)
                }
            })
            clusterExport(cl, varlist = c("get_beta_full", "get_beta_bloc", "li_caret_multibloc"))
            name_vec <- as.character(seq_len(object@parallel$n_process))
            clusterApply(cl, name_vec, function(name) assign("worker_name", name, envir = .GlobalEnv))
        }
    }

    object@model <- caret::train(
        y = object@y_train, x = object@train_cols[, object@col_x], li_index_modes = object@li_index_modes,
        method = li_caret_multibloc, trControl = object@cv, metric = "AUC",
        tuneLength = object@tuneLength, weights_dict = object@weights, tuneGrid = grid, eps = object@eps, ite_max = object@ite_max, n_iter_per_reg = object@n_iter_per_reg,
        index_bloc = object@index_bloc, k_smote = object@k_smote, sampling_choice = object@sampling, index_variable = object@index_variable, is_binary = object@is_binary, classe_1 = object@classe_1, max_temps = object@max_temps
    )
    if (object@parallel$do) {
        stopCluster(cl)
    }
    object@beta_final <- object@model$finalModel$beta_unfolded
    return(object)
})


setMethod("importance_method", "apply_model", function(object) {
    ###### VIEUX ET NON UTILISE: un seul rang########
    # #print(object@model$finalModel$li_beta_K)
    # li_best_beta_J <- object@model$finalModel$li_beta_J
    # li_best_beta_K <- object@model$finalModel$li_beta_K
    # best_beta_autre <- object@model$finalModel$beta_autre
    #     L <- length(li_best_beta_J)
    #     Variables_temps_seul <- unique(object@name_mode[object@index_mode > -0.5])
    #     Variables_names_seul <- unique(object@name_variable[object@index_variable > -0.5])
    #     n_var_radio_differentes <- length(Variables_names_seul)
    #     Variables_names_long <- c()
    #     Variables_temps <- c()
    #     R <- object@model$bestTune$R
    #     for (r in 1:R) {
    #         Variables_names_long <- append(Variables_names_long, paste0(Variables_names_seul, "_R=", r))
    #         Variables_temps <- append(Variables_temps, paste0(Variables_temps_seul, "_R=", r))
    #     }
    #     # print(Variables_names_seul)
    #     # print(Variables_names_long)
    #     li_variables_names <- list()
    #     for (l in 1:L) {
    #         li_variables_names[[l]] <- c()
    #     }
    #     li_variables_names <- vector("list", L)
    #     count <- 0
    #     for (name in Variables_names_long) {
    #         count <- count %% n_var_radio_differentes + 1
    #         num_bloc <- object@index_bloc[count]
    #         # print(length(li_variables_names))
    #         # print(count)
    #         li_variables_names[[num_bloc]] <- append(li_variables_names[[num_bloc]], name)
    #         count <- count + 1
    #     }
    #     grand_vec_temps <- c()
    #     grand_var_temps <- c()
    #     for (l in 1:L) {
    #         grand_vec_temps <- append(grand_vec_temps, li_best_beta_K[[l]])
    #         grand_var_temps <- append(grand_var_temps, paste0(Variables_temps, "_l=", l))
    #     }

    #     time_importance <- data.frame(Variable = grand_var_temps, Overall = abs(grand_vec_temps))
    #     name_df <- paste0("time_importance_", object@id_term)
    #     object@li_df_var_imp[[name_df]] <- time_importance
    #     time_importance <- time_importance[order(-time_importance$Overall), ]
    #     time_importance <- subset(time_importance, Overall > 0.0001)

    #     image <- ggplot2::ggplot(time_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
    #         geom_bar(stat = "identity") +
    #         coord_flip() +
    #         theme_light() +
    #         xlab("Variable") +
    #         ylab("Importance") +
    #         ggtitle("Variable Importance")
    #     ggsave(paste0("plots/logistic_multibloc/importance_time", "_", object@id_term, ".png"), image)

    #     grand_vec_var <- c()
    #     grand_var_var <- c()
    #     for (l in 1:L) {
    #         grand_vec_var <- append(grand_vec_var, li_best_beta_J[[l]])
    #         grand_var_var <- append(grand_var_var, paste0(li_variables_names[[l]], "_l = ", l))
    #         # print(li_best_beta_J[[l]])
    #         # print(paste0(li_variables_names[[l]], "_l = ", l))
    #     }

    #     variable_importance <- data.frame(Variable = grand_var_var, Overall = abs(grand_vec_var))
    #     name_df <- paste0("variable_importance_", object@id_term)
    #     object@li_df_var_imp[[name_df]] <- variable_importance
    #     variable_importance <- variable_importance[order(-variable_importance$Overall), ]
    #     variable_importance <- subset(variable_importance, Overall > 0.0001)

    #     image <- ggplot2::ggplot(variable_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
    #         geom_bar(stat = "identity") +
    #         coord_flip() +
    #         theme_light() +
    #         xlab("Variable") +
    #         ylab("Importance") +
    #         ggtitle("Variable Importance")
    #     ggsave(paste0("plots/logistic_multibloc/importance_var_multi", "_", object@id_term, ".png"), image)



    #     if (length(best_beta_autre) > 0) {
    #         other_names <- colnames(object@train_cols[, object@col_x])[object@index_mode == -1]
    #         other_importance <- data.frame(Variable = other_names, Overall = abs(best_beta_autre))
    #         name_df <- paste0("other_importance_", object@id_term)
    #         object@li_df_var_imp[[name_df]] <- other_importance
    #         other_importance <- other_importance[order(-other_importance$Overall), ]
    #         other_importance <- subset(other_importance, Overall > 0.0001)

    #         image <- ggplot2::ggplot(other_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
    #             geom_bar(stat = "identity") +
    #             coord_flip() +
    #             theme_light() +
    #             xlab("Variable") +
    #             ylab("Importance") +
    #             ggtitle("Variable Importance")
    #         ggsave(paste0("plots/logistic_multibloc/importance_other", "_", object@id_term, ".png"), image)
    #     }

    #     df_cv <- object@model$resample
    #     df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    #     df_long <- melt(df_cv)
    #     object@li_box_plots[[object@id_term]] <- df_long

    #     box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
    #         geom_boxplot() +
    #         stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
    #         theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    #     ggsave(paste0("plots/logistic_multibloc/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)

    #     return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    return(data.frame(Variable = names(object@data_used[, object@col_x]), Overall = abs(object@model$finalModel$beta_unfolded)))
})


# R_1 R_2 lambda
# 6   5   1  1e-04
# Setting levels: control = Classe_0, case = Classe_1
# Setting direction: controls > cases
# [1] "La valeur de l'AUC de test est de 1"
# [1] "La valeur de l'AUC de validation sur chaque fold est de 1"
# [2] "La valeur de l'AUC de validation sur chaque fold est de 1"
# [3] "La valeur de l'AUC de validation sur chaque fold est de 1"
# [4] "La valeur de l'AUC de validation sur chaque fold est de 1"
# [5] "La valeur de l'AUC de validation sur chaque fold est de 1"
# [1] "Ce qui donne une moyenne d'AUC de 1"
# Setting levels: control = Classe_0, case = Classe_1
# Setting direction: controls > cases
# [1] "La valeur de l'AUC de train est de 1"
# [1] "End of analysis phase"
# [1] "L'erreur moyenne de reconstruction du pictogramme est de 0.185689089440778"
# [1] "actuelle moyenne AUC test 1 ite: 1"
# [1] "actuelle moyenne AUC val 1 ite: 1"
# [1] "moyenne AUC test 1"
# [1] "moyenne AUC val 1"



#     ####### Grande boucle de l'algorithme à itérer autant que nécessaire ##########
#     # Elle sera itérée jusqu'à ce que le changement entre le premier et le dernier critère soit inférieur à un epsilon... (+ condit de loop et de temps d'attente)
#     grande_boucle <- function(li_beta_K) {
#         #   Définir Z_J_init (sans les variables cliniques)
#         li_Z_J_init <- lapply(different_blocs, function(l_num) {
#             l_char <- as.character(l_num)
#             mat_x_bloc <- li_x_multi_bloc[[l_char]]
#             beta_K <- li_beta_K[[l_char]]
#             K <- li_dim[[l_char]]$K
#             J <- li_dim[[l_char]]$J
#             R <- current_li_R[[l_char]]
#             li_Z_J_init_bloc <- list() # init car sans clinique du bloc l_num
#             for (r in seq_len(R)) {
#                 Z_J_init_bloc <- matrix(0, nrow = nrow(mat_x_bloc), ncol = J)
#                 for (k in li_different_modes[[l_char]]) {
#                     bloc_k <- mat_x_bloc[, index_mode[index_bloc == l_num] == k]
#                     ordre_variables_k <- order(index_variable[index_bloc == l_num][index_mode[index_bloc == l_num] == k])
#                     bloc_k <- bloc_k[, ordre_variables_k]
#                     index_k <- which(li_different_modes[[l_char]] == k)
#                     Z_J_init_bloc <- Z_J_init_bloc + beta_K[(r - 1) * K + index_k] * bloc_k
#                     li_Z_J_init_bloc[[r]] <- Z_J_init_bloc
#                 }
#             }

#             Z_J_init_bloc <- do.call(cbind, li_Z_J_init_bloc)
#             if (any(is.na(Z_J_init_bloc))) {
#                 stop(" Z_J_init contient des NA")
#             }
#             return(Z_J_init_bloc)
#         })

#         ## Assembler les Z_J_l en un grand Z_J_init
#         Z_J_init <- do.call(cbind, li_Z_J_init)


#         ## Définir Q_J de transformation pour la régularisation (attention aux colonnes inversées)
#         li_vec_diag <- lapply(different_blocs, function(l_num) {
#             l_char <- as.character(l_num)
#             beta_K <- li_beta_K[[l_char]]
#             K <- li_dim[[l_char]]$K
#             J <- li_dim[[l_char]]$J
#             vec_diag_bloc <- c()
#             R <- current_li_R[[l_char]]
#             for (r in 1:R) {
#                 norme_1_r <- norm(as.matrix(beta_K[((r - 1) * K + 1):(r * K)]), type = "1")
#                 diag_elem_r <- rep(norme_1_r, J)
#                 vec_diag_bloc <- append(vec_diag_bloc, diag_elem_r)
#             }
#             return(vec_diag_bloc)
#         })

#         vec_diag <- do.call(c, li_vec_diag)

#         li_petite_boucle <- petite_boucle(Z_init = Z_J_init, vec_Q = vec_diag)
#         beta_J_complet <- li_petite_boucle$beta
#         beta_J_full <- li_petite_boucle$beta[1:(length(li_petite_boucle$beta) - n_restant)]
#         li_beta_J <- list()
#         precedent <- 0
#         for (l_num in different_blocs) {
#             l_char <- as.character(l_num)
#             R <- current_li_R[[l_char]]
#             li_beta_J[[l_char]] <- beta_J_full[(precedent + 1):(precedent + li_dim[[l_char]]$J * R)]
#             precedent <- precedent + li_dim[[l_char]]$J * R
#         }
#         if (n_restant > 0) {
#             beta_autre_J <- li_petite_boucle$beta[(length(li_petite_boucle$beta) - n_restant + 1):length(li_petite_boucle$beta)]
#         } else {
#             beta_autre_J <- c()
#         }
#         Z_J <- li_petite_boucle$Z
#         intercept_J <- li_petite_boucle$intercept
#         Q_J_full <- li_petite_boucle$Q
#         Q_J_inv <- li_petite_boucle$Q_inv
#         # print(beta_J_complet)


#         #########################    Opérer maintenant sur Z_K et Q_K

#         li_Z_K_init <- lapply(different_blocs, function(l_num) {
#             l_char <- as.character(l_num)
#             mat_x_bloc <- li_x_multi_bloc[[l_char]]
#             beta_J <- li_beta_J[[l_char]]
#             J <- li_dim[[l_char]]$J
#             K <- li_dim[[l_char]]$K
#             R <- current_li_R[[l_char]]
#             li_Z_K_init_bloc <- list() # init car sans clinique
#             different_variables <- li_different_variables[[l_char]]
#             index_variable_local <- index_variable[index_bloc == l_num]
#             for (r in 1:R) {
#                 Z_K_init_bloc_r <- matrix(0, nrow = nrow(mat_x_bloc), ncol = K)
#                 for (j in different_variables) {
#                     col_to_extract <- which(index_variable_local == j)
#                     index_j <- which(different_variables == j)
#                     bloc_j <- mat_x_bloc[, col_to_extract]
#                     ordre_mode_j <- order(index_mode[index_bloc == l_num][index_variable_local == j])
#                     bloc_j <- bloc_j[, ordre_mode_j]
#                     Z_K_init_bloc_r <- Z_K_init_bloc_r + beta_J[(r - 1) * J + index_j] * bloc_j
#                 }
#                 li_Z_K_init_bloc[[r]] <- Z_K_init_bloc_r
#             }
#             Z_K_init_bloc <- do.call(cbind, li_Z_K_init_bloc)
#             if (all(abs(Z_K_init_bloc) < 1e-10)) {
#                 print("Z_K est nul : tout sera écrasé à 0")
#             }
#             return(Z_K_init_bloc)
#         })

#         Z_K_init <- do.call(cbind, li_Z_K_init)


#         li_vec_diag_K <- lapply(different_blocs, function(l_num) {
#             l_char <- as.character(l_num)
#             beta_J <- li_beta_J[[l_char]]
#             J <- li_dim[[l_char]]$J
#             K <- li_dim[[l_char]]$K
#             vec_diag_bloc <- c()
#             R <- current_li_R[[l_char]]
#             for (r in 1:R) {
#                 norme_1_r <- norm(as.matrix(beta_J[((r - 1) * J + 1):(r * J)]), type = "1")
#                 diag_elem_r <- rep(norme_1_r, K)
#                 vec_diag_bloc <- append(vec_diag_bloc, diag_elem_r)
#             }
#             return(vec_diag_bloc)
#         })


#         vec_diag <- do.call(c, li_vec_diag_K)


#         ### Appliquer à nouveau la petite boucle

#         li_petite_boucle <- petite_boucle(Z_init = Z_K_init, vec_Q = vec_diag)
#         Z_K <- li_petite_boucle$Z
#         beta_K_full <- li_petite_boucle$beta[1:(length(li_petite_boucle$beta) - n_restant)]
#         li_beta_K <- list()
#         precedent <- 0
#         for (l_num in different_blocs) {
#             l_char <- as.character(l_num)
#             R <- current_li_R[[l_char]]
#             li_beta_K[[l_char]] <- beta_K_full[(precedent + 1):(precedent + li_dim[[l_char]]$K * R)]
#             precedent <- precedent + li_dim[[l_char]]$K * R
#         }

#         ## Renormaliser
#         for (l_num in different_blocs) {
#             l_char <- as.character(l_num)
#             beta_K <- li_beta_K[[l_char]]
#             K <- li_dim[[l_char]]$K
#             R <- current_li_R[[l_char]]
#             for (r in 1:R) {
#                 norm_beta_K_r <- norm(as.matrix(beta_K[((r - 1) * K + 1):(r * K)]), type = "2")
#                 if (norm_beta_K_r > 0) {
#                     beta_K[((r - 1) * K + 1):(r * K)] <- beta_K[((r - 1) * K + 1):(r * K)] / norm_beta_K_r
#                 }
#                 li_beta_K[[l_char]] <- beta_K
#             }
#         }
#         intercept_K <- li_petite_boucle$intercept
#         if (n_restant > 0) {
#             beta_autre_K <- li_petite_boucle$beta[(length(li_petite_boucle$beta) - n_restant + 1):length(li_petite_boucle$beta)]
#         } else {
#             beta_autre_K <- c()
#         }
#         Q_K_full <- li_petite_boucle$Q
#         Q_K_inv <- li_petite_boucle$Q_inv
#         beta_K_complet <- li_petite_boucle$beta

#         # Calculer les deux critères
#         crit_log_J <- crit_logistic(Z_J %*% Q_J_full, y_numeric, Q_J_inv %*% beta_J_complet, intercept_J, param$lambda)
#         crit_log_K <- crit_logistic(Z_K %*% Q_K_full, y_numeric, Q_K_inv %*% beta_K_complet, intercept_K, param$lambda)
#         rapport <- abs(crit_log_J - crit_log_K) / abs(crit_log_J)
#         return(list(
#             li_beta_J = li_beta_J, li_beta_K = li_beta_K, beta_autre = beta_autre_K, intercept_J = intercept_J, intercept_K = intercept_K,
#             rapport = rapport, crit_log_J = crit_log_J, crit_log_K = crit_log_K, beta_J_complet = beta_J_complet, beta_K_complet = beta_K_complet
#         ))
#     }

#     ##### Exécuter la grande boucle autant que les critères le demandent
#     continue <- TRUE
#     ###### Attention beta_J et beta_K n'existent plus, ce sont des listes!! (plus simple à manipuler avec les dimensions variables)
#     li_beta_K <- list()
#     for (l_num in different_blocs) {
#         l_char <- as.character(l_num)
#         R <- current_li_R[[l_char]]
#         li_beta_K[[l_char]] <- rnorm(li_dim[[l_char]]$K * R, mean = 0, sd = 1)
#     }
#     iteration <- 1
#     memoire_crit_J <- 1
#     memoire_crit_K <- 1
#     debut_time <- Sys.time()
#     while (continue) {
#         li_grande_boucle <- grande_boucle(li_beta_K)
#         li_beta_J <- li_grande_boucle$li_beta_J
#         li_beta_K <- li_grande_boucle$li_beta_K
#         beta_autre <- li_grande_boucle$beta_autre
#         intercept_K <- li_grande_boucle$intercept_K

#         rapport <- li_grande_boucle$rapport
#         crit_log_J <- li_grande_boucle$crit_log_J
#         crit_log_K <- li_grande_boucle$crit_log_K
#         # print(crit_log_J)
#         # print(crit_log_K)

#         if (crit_log_J != -Inf & crit_log_K != -Inf) {
#             if (rapport < eps) {
#                 continue <- FALSE
#             }
#         }
#         if (iteration >= ite_max) {
#             continue <- FALSE
#             print(paste("Warning, pas de convergence. Le rapport vaut:", rapport))
#         }
#         if (crit_log_J != -Inf & crit_log_K != -Inf) {
#             if (abs(crit_log_J - memoire_crit_J) < 1e-10 & abs(crit_log_K - memoire_crit_K) < 1e-10) {
#                 continue <- FALSE
#                 print("Warning, LOOP!!")
#                 print(paste("Rapport:", rapport))
#                 print(paste("Crit_J:", crit_log_J))
#                 print(paste("Crit_K:", crit_log_K))
#             }
#         }
#         iteration <- iteration + 1
#         memoire_crit_J <- crit_log_J
#         memoire_crit_K <- crit_log_K
#         delta_t <- as.numeric(Sys.time() - debut_time, units = "secs")
#         if (delta_t > 30) {
#             continue <- FALSE
#             print("Warning, trop long")
#             print(paste("Rapport:", rapport))
#             print(paste("Crit_J:", crit_log_J))
#             print(paste("Crit_K:", crit_log_K))
#         }
#     }
#     futur_fit <- list(
#         li_beta_J = li_beta_J, li_beta_K = li_beta_K, beta_autre = beta_autre, intercept = intercept_K, current_li_R = current_li_R, li_x_multi_bloc = li_x_multi_bloc,
#         li_dim = li_dim, lev = lev, index_mode = index_mode, index_bloc = index_bloc, index_variable = index_variable, li_norm = li_norm, classe_maj = classe_maj, classe_min = classe_min, classe_1 = classe_1, classe_0 = classe_0
#     )
#     # print(li_beta_K)
#     # print(li_beta_J)
#     futur_fit$beta_unfolded <- get_beta_full(futur_fit)
#     # print(futur_fit$beta_unfolded)
#     return(futur_fit)
# }

library(jsonlite)
library(openxlsx)
library(readxl)
library(writexl)
source("utils/utils.r")
import_folder("./utils")
import_folder("./StepReg_modif/R")
# set.seed(1)

get_beta <- function(config_extrac) {
    if (config_extrac$is_random) {
        beta_random <- config_extrac$beta_random
        R <- beta_random$R
        li_d <- beta_random$d
        L <- length(li_d)
        beta_K_random <- beta_random$beta_K
        K <- beta_random$K
        beta_K <- list()
        for (l in 1:L) {
            beta_K[[l]] <- list()
            for (r in 1:R) {
                beta_K[[l]][[r]] <- rnorm(K, beta_K_random$mu, beta_K_random$sigma)
            }
        }
        beta_J_random <- beta_random$beta_J
        beta_J <- list()
        for (l in 1:L) {
            beta_J[[l]] <- list()
            mu <- beta_J_random$vec_mu[l]
            sigma <- beta_J_random$vec_sigma[l]
            for (r in 1:R) {
                beta_J[[l]][[r]] <- rnorm(li_d[l], mu, sigma)
            }
        }
        beta_tab_random <- beta_random$beta_tab
        beta_tab <- rnorm(beta_tab_random$M, beta_tab_random$mu, beta_tab_random$sigma)
    } else {
        beta_determ <- config_extrac$beta_determ
        R <- beta_determ$R
        beta_K_init <- beta_determ$beta_K
        beta_J_init <- beta_determ$beta_J
        beta_J <- list()
        beta_K <- list()
        for (l in 1:length(beta_J_init)) {
            beta_J[[l]] <- beta_J_init[[l]][1:R]
            beta_K[[l]] <- beta_K_init[[l]][1:R]
        }
        beta_tab <- beta_determ$beta_tab
    }
    beta_K <- lapply(beta_K, function(beta_K_l) {
        lapply(beta_K_l, function(beta_K_l_r) {
            beta_K_l_r / ifelse(norm(as.matrix(beta_K_l_r), type = "2") > 0.0001, norm(as.matrix(beta_K_l_r), type = "2") > 0.0001, 1)
        })
    })
    K <- length(beta_K[[1]][[1]])
    L <- length(beta_J)
    li_d <- c()
    for (l in 1:L) {
        li_d <- c(li_d, length(beta_J[[l]][[1]]))
    }
    J <- sum(li_d)
    M <- length(beta_tab)
    beta <- c()
    for (l in 1:L) {
        beta_l <- c()
        for (k in 1:K) {
            beta_k_l <- rep(0, li_d[l])
            for (r in 1:R) {
                beta_k_l <- beta_k_l + beta_K[[l]][[r]][[k]] * beta_J[[l]][[r]]
            }
            beta_l <- c(beta_l, beta_k_l)
        }
        beta <- c(beta, beta_l)
    }
    beta <- c(beta, beta_tab)
    beta_0 <- rnorm(1)
    return(list(beta = beta, K = K, L = L, J = J, M = M, R = R, li_d = li_d, beta_K = beta_K, beta_J = beta_J, beta_tab = beta_tab, beta_0 = beta_0))
}

find_moy_beta_bloc <- function(K, d_l, beta_l) {
    moy_beta <- sapply(1:d_l, function(j) {
        value <- 0
        for (k in 1:K) {
            value <- value + beta_l[(k - 1) * d_l + j]
        }
        return(value / K)
    })
    return(moy_beta)
}

find_moy_beta_tens <- function(K, L, li_d, beta) {
    moy_beta <- c()
    index_min <- 1
    index_max <- 0
    for (l in 1:L) {
        index_max <- index_max + li_d[l] * K
        moy_beta_l <- find_moy_beta_bloc(K, li_d[l], beta[index_min:index_max])
        moy_beta <- c(moy_beta, moy_beta_l)
        index_min <- index_max + 1
    }
    return(moy_beta)
}

# Séparer entre mean_x_J et mean_x_K: structure multivoie
# On va créer moyenne de x_J par rang et par classe (en se basant sur beta). Cette moyenne tient compte de la structure de blocs et de tenseurs.
create_mean_x_J <- function(K, L, li_d, beta, classe_num, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global) {
    moy_beta <- find_moy_beta_tens(K, L, li_d, beta)
    # sign_beta <- ifelse(moy_beta >= 0, 1, -1)
    sign_beta <- moy_beta / mean(abs(moy_beta))
    li_sign_beta <- list()
    index_min <- 1
    index_max <- 0
    for (l in 1:L) {
        index_max <- index_max + li_d[l]
        li_sign_beta[[l]] <- sign_beta[index_min:index_max]
        index_min <- index_max + 1
    }
    amplif_factor_bloc <- runif(L, min_amplif_bloc, max_amplif_bloc)
    li_amplif_factor_tens <- lapply(1:L, function(l) {
        rnorm(li_d[l], amplif_factor_bloc[l], sigma_amplif_global) # Donner une structure en blocs
    })
    val_class <- ifelse(classe_num == 1, 1, -1)
    li_amplif_factor_tens <- lapply(1:L, function(l) {
        return(li_sign_beta[[l]] * li_amplif_factor_tens[[l]] * val_class)
    }) # Tenir compte de la classe
    return(li_amplif_factor_tens) # On liste sur le bloc L
}

# Un mu_J et mu_K par rang: tenir compte du multi_rangs
create_mean_x_K <- function(K, L, sigma_amplif_K) {
    li_mu_K <- lapply(1:L, function(l) {
        rnorm(K, 1, sigma_amplif_K)
    }) # On liste à nouveau sur le bloc L
    return(li_mu_K)
}

create_mean_tab <- function(K, J, M, beta, classe_num, sigma_amplif_global) {
    moy_beta_tab <- beta[(K * J + 1):(K * J + M)]
    # sign_beta_tab <- ifelse(moy_beta_tab >= 0, 1, -1)
    sign_beta_tab <- moy_beta_tab / mean(abs(moy_beta_tab))
    amplif_factor_tab <- rnorm(M, 1, sigma_amplif_global)
    val_class <- ifelse(classe_num == 1, 1, -1)
    return(beta[(K * J + 1):(K * J + M)] * amplif_factor_tab * val_class)
}


create_Sigma_bloc <- function(K, d_l) {
    mat <- matrix(rnorm(K * d_l, 0, 1), nrow = K * d_l, ncol = K * d_l)
    Sigma_bloc <- mat %*% t(mat)
    Sigma_bloc <- 50 * Sigma_bloc / norm(as.matrix(Sigma_bloc), type = "2")
    # print(Sigma_bloc)
    return(Sigma_bloc)
}

generate_li_param <- function(n, J, K, L, M, li_d, R, beta, classe_num, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global, sigma_amplif_K) {
    li_mean_x_K <- lapply(1:R, function(r) {
        create_mean_x_K(K, L, sigma_amplif_K)
    })
    li_mean_x_J <- lapply(1:R, function(r) {
        create_mean_x_J(K, L, li_d, beta, classe_num, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global)
    })
    mean_tab <- create_mean_tab(K, J, M, beta, classe_num, sigma_amplif_global)
    li_elem_K <- lapply(1:R, function(r) {
        lapply(1:L, function(l) {
            Sigma <- create_Sigma_bloc(K, 1)
            elem_K <- MASS::mvrnorm(n, li_mean_x_K[[r]][[l]], Sigma)
        })
    })
    li_elem_J <- lapply(1:R, function(r) {
        lapply(1:L, function(l) {
            Sigma <- create_Sigma_bloc(1, li_d[l])
            elem_J <- MASS::mvrnorm(n, li_mean_x_J[[r]][[l]], Sigma)
        })
    })
    elem_tab <- MASS::mvrnorm(n, mean_tab, diag(1, M))
    return(list(li_elem_K = li_elem_K, li_elem_J = li_elem_J, elem_tab = elem_tab))
}

construct_X <- function(li_elem_K, li_elem_J, elem_tab) {
    n <- nrow(li_elem_K[[1]][[1]])
    L <- length(li_elem_K[[1]])
    R <- length(li_elem_K)
    X <- matrix(0, nrow = n, ncol = 0)
    li_elem_K_good_order <- vector("list", length = L)
    li_elem_J_good_order <- vector("list", length = L)
    li_elem_K_good_order <- lapply(li_elem_K_good_order, function(x) vector("list", length = R))
    li_elem_J_good_order <- lapply(li_elem_J_good_order, function(x) vector("list", length = R))
    for (r in 1:R) {
        for (l in 1:L) {
            li_elem_K_good_order[[l]][[r]] <- li_elem_K[[r]][[l]]
            li_elem_J_good_order[[l]][[r]] <- li_elem_J[[r]][[l]]
        }
    }
    for (l in 1:L) {
        li_elem_K_l <- li_elem_K_good_order[[l]]
        li_elem_J_l <- li_elem_J_good_order[[l]]
        li_product_kron <- lapply(1:R, function(r) {
            elem_K <- li_elem_K_l[[r]]
            elem_J <- li_elem_J_l[[r]]
            li_mat_per_lines <- lapply(1:n, function(i) {
                kronecker(elem_K[i, ], elem_J[i, ])
            })
            return(do.call(rbind, li_mat_per_lines))
        })
        X_l <- Reduce("+", li_product_kron)
        X <- cbind(X, X_l)
    }
    X <- cbind(X, elem_tab)
    return(X)
}



create_param_explanatory <- function(K, L, J, M, li_d, R, n, beta, prop, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global, sigma_noise, sigma_amplif_K) {
    # n : Nombre d'individus à créer
    # L: nombre de blocs
    # J: nombre de variables
    # K: nombre de voies (de temps)
    # M: nombre de variables tabulaires
    # Bilan : J*K variables tensorielles
    # li_d: vecteur des dimensions de chaque bloc
    # R: rang des tenseurs
    # n: beta: coefficient beta de regression
    # prop: proportion de la classe 1
    # mu_amplif_bloc: les espérances des facteurs d'amplification des variables de chaque bloc seront tirés selon une loi normale de moyenne mu_amplif_bloc[[l]]
    # sigma_amplif_bloc: les espérances des facteurs d'amplification des variables de chaque bloc seront tirés selon une loi normale d'écart type sigma_amplif_bloc[[l]]
    # sigma_amplif_global: les facteurs d'amplification seront tirés selon une loi normale d'écart type sigma_amplif_global (et espérance tirée selon la loi des deux précédents paramètres)
    # sigma_noise: écart-type du bruit
    # sigma_amplif_K: écart_type des facteurs d'amplification des voies

    n_1 <- round(n * prop)
    n_0 <- n - n_1
    li_param_0 <- generate_li_param(n_0, J, K, L, M, li_d, R, beta, 0, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global, sigma_amplif_K)
    li_param_1 <- generate_li_param(n_1, J, K, L, M, li_d, R, beta, 1, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global, sigma_amplif_K)
    X_0 <- construct_X(li_param_0$li_elem_K, li_param_0$li_elem_J, li_param_0$elem_tab)
    write_xlsx(as.data.frame(X_0), "..//data//data_test.xlsx")
    X_0 <- cbind(rep(0, n_0), X_0)
    X_1 <- construct_X(li_param_1$li_elem_K, li_param_1$li_elem_J, li_param_1$elem_tab)
    X_1 <- cbind(rep(1, n_1), X_1)
    X <- rbind(X_1, X_0)
    eps <- matrix(rnorm(n * (J * K + M), 0, sigma_noise), nrow = n, ncol = J * K + M)
    X[, 2:ncol(X)] <- X[, 2:ncol(X)] + eps
    X <- as.matrix(as.data.frame(X))
    colnames(X)[1] <- "theoritical_class"
    return(X)
}


apply_logistic <- function(mat_tot, beta, beta_0, index_variable, is_binary) {
    df_renorm <- renormalize_in_model_fit_index_mode(mat_tot, index_variable, is_binary)$new_x
    mat_renorm <- as.matrix(df_renorm)
    prod_scal <- as.vector(mat_renorm %*% beta)
    proba <- 1 / (1 + exp(-beta_0 - prod_scal))
    prediction <- rep(-1, length(proba))
    for (i in 1:length(proba)) {
        if (proba[i] > 0.5) {
            prediction[i] <- 1
        } else {
            prediction[i] <- 0
        }
    }
    return(prediction)
}



extract_all <- function(config_extrac, sysname) {
    current_seed <- .Random.seed
    set.seed(42)
    recreate_beta <- config_extrac$recreate_beta
    if (recreate_beta) {
        li_beta <- get_beta(config_extrac)
        wb <- createWorkbook()
        addWorksheet(wb, "main page")
        saveRDS(li_beta, file = "../data/RDS/li_beta.rds")
        # li_beta_reloaded <- readRDS(file = "../data/RDS/li_beta.rds")
    } else {
        li_beta <- readRDS(file = "../data/RDS/li_beta.rds")
    }
    beta <- li_beta$beta
    beta_0 <- li_beta$beta_0
    K <- li_beta$K
    L <- li_beta$L
    J <- li_beta$J
    M <- li_beta$M
    R <- li_beta$R
    li_d <- li_beta$li_d
    # beta_K <- li_beta$beta_K
    # beta_J <- li_beta$beta_J
    # beta_tab <- li_beta$beta_tab

    n <- config_extrac$create_var$n

    names_col <- rep("", K * J + M)
    index_mode <- rep(0, K * J + M)
    index_bloc <- rep(0, K * J + M)
    starting_bloc_add <- 0
    num_var_add <- 0
    index_variable <- rep(0, K * J + M)
    for (l in 1:L) {
        for (k in 1:K) {
            start_index <- starting_bloc_add + (k - 1) * li_d[l] + 1
            end_index <- starting_bloc_add + k * li_d[l]
            index_mode[start_index:end_index] <- k
            for (j in 1:li_d[l]) {
                index_variable[(starting_bloc_add + (k - 1) * li_d[l] + j)] <- num_var_add + j
                names_col[(starting_bloc_add + (k - 1) * li_d[l] + j)] <- paste("variable", num_var_add + j, "mode", k, "bloc", l)
            }
        }
        index_bloc[(starting_bloc_add + 1):(starting_bloc_add + K * li_d[l])] <- l
        starting_bloc_add <- starting_bloc_add + K * li_d[l]
        num_var_add <- num_var_add + li_d[l]
    }
    for (m in 1:M) {
        index_mode[K * J + m] <- -1
        index_bloc[K * J + m] <- -1
        index_variable[K * J + m] <- -1
        names_col[K * J + m] <- paste("variable clinique", m)
    }

    is_binary <- rep(FALSE, K * J + M)


    if (config_extrac$resample_explain) {
        min_amplif_bloc <- config_extrac$create_var$min_amplif_bloc
        max_amplif_bloc <- config_extrac$create_var$max_amplif_bloc
        sigma_amplif_global <- config_extrac$create_var$sigma_amplif_global
        sigma_amplif_K <- config_extrac$create_var$sigma_amplif_K
        sigma_noise <- config_extrac$create_var$sigma_noise
        prop <- config_extrac$create_var$prop
        mat_tot <- create_param_explanatory(K, L, J, M, li_d, R, n, beta, prop, min_amplif_bloc, max_amplif_bloc, sigma_amplif_global, sigma_noise, sigma_amplif_K)
    } else {
        if (sysname == "Linux") {
            mat_tot_names <- as.matrix(read.csv("..//data//data_used.csv"))
        } else {
            mat_tot_names <- as.matrix(read.csv("..\\data\\data_used.csv"))
        }
        mat_tot[, 1] <- ifelse(mat_tot_names[, 1] == "Classe_1", 1, 0)
        mat_tot[, 2] <- ifelse(mat_tot_names[, 2] == "Classe_1", 1, 0)
        mat_tot <- as.matrix(apply(mat_tot, 2, as.numeric))
    }
    if (config_extrac$resample_explain | recreate_beta) {
        beta_class <- apply_logistic(mat_tot[, 2:ncol(mat_tot)], beta, beta_0, index_variable, is_binary)
        mat_tot <- cbind(beta_class, mat_tot)
        mat_tot_names <- as.data.frame(mat_tot)
        mat_tot_names[, 1] <- ifelse(mat_tot[, 1] == 1, "Classe_1", "Classe_0")
        mat_tot_names[, 2] <- ifelse(mat_tot[, 2] == 1, "Classe_1", "Classe_0")
        colnames(mat_tot_names)[3:ncol(mat_tot_names)] <- names_col

        if (sysname == "Linux") {
            write.csv(mat_tot_names, "..//data//data_used.csv", row.names = FALSE)
            write_xlsx(mat_tot_names, "..//data//data_used.xlsx")
        } else {
            write.csv(mat_tot_names, "..\\data\\data_used.csv", row.names = FALSE)
            write_xlsx(mat_tot_names, "..\\data\\data_used.xlsx")
        }
    }

    beta_class <- mat_tot[, 1]
    theoritical_class <- mat_tot[, 2]
    .Random.seed <<- current_seed

    if (config_extrac$good_y == "structural") {
        exclude_cols <- c("beta_class")
        explained_col <- c("theoritical_class")
    } else {
        exclude_cols <- c("theoritical_class")
        explained_col <- c("beta_class")
    }
    info_cols <- list(exclude_cols = exclude_cols, explained_col = explained_col)
    name_mode <- as.character(index_mode)
    name_bloc <- as.character(index_bloc)
    name_variable <- as.character(index_variable)

    saveRDS(info_cols, file = "../data/RDS/info_cols.rds")
    saveRDS(is_binary, file = "../data/RDS/is_binary.rds")
    saveRDS(index_mode, file = "../data/RDS/index_mode.rds")
    saveRDS(index_bloc, file = "../data/RDS/index_bloc.rds")
    saveRDS(index_variable, file = "../data/RDS/index_variable.rds")

    saveRDS(name_mode, file = "../data/RDS/name_mode.rds")
    saveRDS(name_bloc, file = "../data/RDS/name_bloc.rds")
    saveRDS(name_variable, file = "../data/RDS/name_variable.rds")


    prop_1 <- sum(theoritical_class) / length(theoritical_class)
    prop_0 <- 1 - prop_1
    cat("Il y a une proportion de classe 1 égale à ", prop_1, "et une proportion de classe 0 égale à", prop_0, "\n")
    accuracy_loc <- sum(beta_class == theoritical_class) / length(theoritical_class)
    cat("l'accuracy de la concordance entre la création de X et ce que prédit beta est de", accuracy_loc, "\n")
    faux_pos <- sum(beta_class == 1 & theoritical_class == 0)
    faux_neg <- sum(beta_class == 0 & theoritical_class == 1)
    cat("Il y a ", faux_pos, "faux positifs et ", faux_neg, "faux négatifs\n")
}
# config_extrac <- read_json("configs/config_extrac_simul.json", simplifyVector = TRUE, simplifyMatrix = FALSE)
# extract_all(config_extrac, sys = Sys.info()["sysname"])

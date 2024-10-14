# find_sigma_mu <- function(mat) {
#     df <- as.data.frame(mat)
#     cols_to_norm <- c()
#     age_exist <- FALSE
#     for (name in colnames(df)) {
#         if (substr(name, nchar(name) - 3, nchar(name)) %in% c("PORT", "ART_", "VEIN", "TARD")) {
#             cols_to_norm <- c(cols_to_norm, name)
#         }
#         if (name == "Age_at_disease") {
#             age_exist <- TRUE
#         }
#     }
#     time_injs <- unique(substr(cols_to_norm, nchar(cols_to_norm) - 3, nchar(cols_to_norm)))
#     n_time_inj <- length(time_injs)
#     true_cols_to_norm <- unique(substr(cols_to_norm, 1, nchar(cols_to_norm) - 5))
#     n_variables <- length(true_cols_to_norm)
#     df_mu <- as.data.frame(matrix(0, nrow = nrow(mat), ncol = length(cols_to_norm)))
#     df_mu[0, ] <- 0
#     colnames(df_mu) <- cols_to_norm
#     df_sigma <- as.data.frame(matrix(0, nrow = nrow(mat), ncol = length(cols_to_norm)))
#     df_sigma[0, ] <- 0
#     colnames(df_sigma) <- cols_to_norm
#     for (name in true_cols_to_norm) {
#         good_cols_name <- c()
#         for (i in 1:n_time_inj) {
#             good_cols_name <- c(good_cols_name, paste0(name, "_", time_injs[i]))
#         }
#         good_cols <- df[, good_cols_name]
#         long_vec <- unlist(good_cols)
#         sigma <- sd(long_vec, na.rm = TRUE)
#         mu <- mean(long_vec, na.rm = TRUE)
#         df_mu[good_cols_name] <- mu
#         df_sigma[good_cols_name] <- sigma
#     }
#     if (age_exist) {
#         mu_disease <- mean(df$Age_at_disease, na.rm = TRUE)
#         sigma_disease <- sd(df$Age_at_disease, na.rm = TRUE)
#     } else {
#         mu_disease <- 0
#         sigma_disease <- -1
#     }
#     return(list(mu = df_mu, sigma = df_sigma, mu_disease = mu_disease, sigma_disease = sigma_disease))
# }

# renormalize_in_model_fit <- function(x) {
#     x <- as.data.frame(x)
#     li <- find_sigma_mu(x)
#     df_mu <- li$mu
#     df_sigma <- li$sigma
#     col_max_radio <- ncol(df_mu)
#     x[, 1:col_max_radio] <- (x[, 1:col_max_radio] - df_mu) / df_sigma
#     if (li$sigma_disease > -0.5) {
#         x$Age_at_disease <- (x$Age_at_disease - li$mu_disease) / li$sigma_disease
#         li_age <- list(do = TRUE, mu = li$mu_disease, sigma = li$sigma_disease)
#     } else {
#         li_age <- list(do = FALSE)
#     }
#     return(list(x = x, df_mu = df_mu, df_sigma = df_sigma, li_age = li_age))
# }

# renormalize_in_model_pred <- function(newdata, modelFit) {
#     li_norm <- modelFit$li_norm
#     df_mu_former <- li_norm$df_mu
#     df_sigma_former <- li_norm$df_sigma
#     li_age <- li_norm$li_age
#     n_line_new <- nrow(newdata)
#     mat_mu <- as.matrix(df_mu_former[rep(1, n_line_new), ])
#     mat_sigma <- as.matrix(df_sigma_former[rep(1, n_line_new), ])

#     newdata <- as.data.frame(newdata)
#     mat_mu <- as.matrix(find_sigma_mu(newdata)$mu)
#     mat_sigma <- as.matrix(find_sigma_mu(newdata)$sigma)
#     n_col_max_radio <- ncol(mat_mu)
#     if (li_age[["do"]]) {
#         mu <- li_age[["mu"]]
#         sigma <- li_age[["sigma"]]
#         newdata$Age_at_disease <- (newdata$Age_at_disease - mu) / sigma
#     }

#     newdata <- as.matrix(newdata)
#     newdata[, 1:n_col_max_radio] <- (newdata[, 1:n_col_max_radio] - mat_mu) / mat_sigma
#     return(newdata)
# }

get_bloc_vec_liver <- function(vec, exist_temps = TRUE, different_clin = FALSE) {
    name_bloc <- rep("a", length(vec))
    index_bloc <- rep(NA, length(vec))
    for (j in 1:length(vec)) {
        name_col <- vec[j]
        if (!substr(name_col, nchar(name_col) - 3, nchar(name_col)) %in% c("PORT", "ART_", "VEIN", "TARD") & exist_temps & !grepl("average_filter", name_col)) {
            if (different_clin) {
                name_bloc[j] <- paste0("clinical_", j)
                index_bloc[j] <- -1
            } else {
                if (grepl("shape", name_col)) {
                    name_bloc[j] <- "shape_univariate"
                    index_bloc[j] <- -1
                } else {
                    name_bloc[j] <- "clinical"
                    index_bloc[j] <- -1
                }
            }
        } else {
            classed <- FALSE
            if (grepl("firstorder", name_col)) {
                if (grepl("sain", name_col)) {
                    name_bloc[j] <- "first order sain"
                    index_bloc[j] <- 4
                    classed <- TRUE
                } else {
                    name_bloc[j] <- "first order"
                    index_bloc[j] <- 1
                    classed <- TRUE
                }
            }
            if (grepl("shape", name_col)) {
                name_bloc[j] <- "shape"
                index_bloc[j] <- 2
                classed <- TRUE
            }
            if (grepl("glcm", name_col) | grepl("gldm", name_col) | grepl("glrl", name_col) | grepl("glsz", name_col) | grepl("ngtd", name_col)) {
                if (grepl("sain", name_col)) {
                    name_bloc[j] <- "texture sain"
                    index_bloc[j] <- 5
                    classed <- TRUE
                } else {
                    name_bloc[j] <- "texture"
                    index_bloc[j] <- 3
                    classed <- TRUE
                }
            }
            if (classed == FALSE) {
                stop("Erreur de classification des blocs")
            }
        }
    }
    return(list(index_bloc = index_bloc, index_name = name_bloc))
}

get_mode_vec_liver <- function(x) {
    index_mode <- rep(0, ncol(x))
    index_name <- rep("a", ncol(x))
    for (j in 1:ncol(x)) {
        name_col <- names(x)[j]
        classified <- FALSE
        if (substr(name_col, nchar(name_col) - 3, nchar(name_col)) == "ART_") {
            index_mode[j] <- 1
            index_name[j] <- "ART_"
            classified <- TRUE
        }
        if (substr(name_col, nchar(name_col) - 3, nchar(name_col)) == "PORT") {
            index_mode[j] <- 2
            index_name[j] <- "PORT"
            classified <- TRUE
        }
        if (substr(name_col, nchar(name_col) - 3, nchar(name_col)) == "VEIN") {
            index_mode[j] <- 4
            index_name[j] <- "VEIN"
            classified <- TRUE
        }
        if (substr(name_col, nchar(name_col) - 3, nchar(name_col)) == "TARD") {
            index_mode[j] <- 3 # pour rester contigus, on sait jamais...
            index_name[j] <- "TARD"
            classified <- TRUE
        }
        if (!classified) {
            index_mode[j] <- -1
            index_name[j] <- "clinical"
        }
    }
    return(list(index_mode = index_mode, index_name = index_name))
}

get_variable_vec_liver <- function(x) {
    df_loc <- data.frame(name_cov = colnames(x))
    index_name <- substr(df_loc$name_cov, 0, nchar(df_loc$name_cov) - 5)
    index_clin <- which(!substr(df_loc$name_cov, nchar(df_loc$name_cov) - 3, nchar(df_loc$name_cov)) %in% c("PORT", "ART_", "VEIN", "TARD"))
    index_name[index_clin] <- colnames(x)[index_clin]
    index_name_no_clin <- index_name[-index_clin]
    index_variable <- rep(0, ncol(x))
    ordered_levels <- unique(index_name_no_clin)
    index_variable[-index_clin] <- as.integer(factor(index_name_no_clin, levels = ordered_levels))
    index_variable[index_clin] <- -1
    return(list(index_variable = index_variable, index_name = index_name))
}

get_variable_vec_liver_big <- function(x) {
    df_loc <- data.frame(name_cov = colnames(x))
    index_name <- sapply(df_loc$name_cov, function(name_cov) {
        true_name_cov <- sub("slice.*", "", name_cov) # avant le mot slice
        if (substr(true_name_cov, nchar(true_name_cov) - 3, nchar(true_name_cov)) %in% c("PORT", "ART_", "VEIN", "TARD")) {
            return(substr(true_name_cov, 0, nchar(true_name_cov) - 5))
        }
        return(true_name_cov)
    })


    index_clin <- which(!substr(df_loc$name_cov, nchar(df_loc$name_cov) - 3, nchar(df_loc$name_cov)) %in% c("PORT", "ART_", "VEIN", "TARD") & !grepl("average_filter", df_loc$name_cov))

    index_name[index_clin] <- colnames(x)[index_clin]
    index_name_no_clin <- index_name[-index_clin]
    index_variable <- rep(0, ncol(x))
    ordered_levels <- unique(index_name_no_clin)
    index_variable[-index_clin] <- as.integer(factor(index_name_no_clin, levels = ordered_levels))
    index_variable[index_clin] <- -1
    index_name <- unname(index_name)
    return(list(index_variable = index_variable, index_name = index_name))
}

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

def_false_grid_multibloc <- function(L) {
    create_grid_multibloc <- function(x, y, len = NULL, search = "grid") {
        if (search == "grid") {
            lambda <- log(seq(exp(0), exp(0.05), length.out = len + 1)[2:(len + 1)])
            li_R_data_frame <- lapply(1:L, function(l) {
                c(1, 2, 3, 4)
            })
            setNames(li_R_data_frame, paste0("R_", 1:L))
        } else {
            lambda <- log(runif(len, min = exp(0.001), max = exp(0.05)))
            li_R_data_frame <- lapply(1:L, function(l) {
                c(1, 2, 3, 4)
            })
            li_R_data_frame <- setNames(li_R_data_frame, paste0("R_", 1:L))
        }
        li_tot <- c(li_R_data_frame, list(lambda = lambda))
        data_frame_grid <- expand.grid(li_tot)
        return(data_frame_grid)
    }
}


better_create_grid_multibloc <- function(x, y, len = NULL, search = "grid", L, lambda_min = 0.001, lambda_max = 0.05, R_min = 1, R_max = 5, tune_R = 2, li_R = NULL, same_R) {
    ### Attention, si li_R est renseigné et same_R est TRUE alors il faut avoir autant de rangs proposés par bloc (ils seront sur 1 même ligne)
    ## tune_R est ignoré si on fournit li_R
    print(li_R$block_2)
    if (search == "grid") {
        lambda <- seq(lambda_min, lambda_max, length.out = len)[1:len]
    } else {
        lambda <- runif(len, min = lambda_min, max = lambda_max)
    }
    if (is.null(li_R)) {
        li_R_data_frame <- lapply(1:L, round(seq(R_min, R_max, length.out = tune_R)))
    } else {
        li_R_data_frame <- li_R
    }
    if (same_R) {
        R_data_frame <- as.data.frame(li_R_data_frame)
    } else {
        R_data_frame <- expand.grid(li_R_data_frame)
        print(R_data_frame)
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
    return(as.data.frame(data_frame_grid))
}

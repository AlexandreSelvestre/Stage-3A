find_sigma_mu <- function(mat) {
    df <- as.data.frame(mat)
    cols_to_norm <- c()
    age_exist <- FALSE
    for (name in colnames(df)) {
        if (substr(name, nchar(name) - 3, nchar(name)) %in% c("PORT", "ART_", "VEIN", "TARD")) {
            cols_to_norm <- c(cols_to_norm, name)
        }
        if (name == "Age_at_disease") {
            age_exist <- TRUE
        }
    }
    time_injs <- unique(substr(cols_to_norm, nchar(cols_to_norm) - 3, nchar(cols_to_norm)))
    n_time_inj <- length(time_injs)
    true_cols_to_norm <- unique(substr(cols_to_norm, 1, nchar(cols_to_norm) - 5))
    n_variables <- length(true_cols_to_norm)
    df_mu <- as.data.frame(matrix(0, nrow = nrow(mat), ncol = length(cols_to_norm)))
    df_mu[0, ] <- 0
    colnames(df_mu) <- cols_to_norm
    df_sigma <- as.data.frame(matrix(0, nrow = nrow(mat), ncol = length(cols_to_norm)))
    df_sigma[0, ] <- 0
    colnames(df_sigma) <- cols_to_norm
    for (name in true_cols_to_norm) {
        good_cols_name <- c()
        for (i in 1:n_time_inj) {
            good_cols_name <- c(good_cols_name, paste0(name, "_", time_injs[i]))
        }
        good_cols <- df[, good_cols_name]
        long_vec <- unlist(good_cols)
        sigma <- sd(long_vec, na.rm = TRUE)
        mu <- mean(long_vec, na.rm = TRUE)
        df_mu[good_cols_name] <- mu
        df_sigma[good_cols_name] <- sigma
    }
    if (age_exist) {
        mu_disease <- mean(df$Age_at_disease, na.rm = TRUE)
        sigma_disease <- sd(df$Age_at_disease, na.rm = TRUE)
    } else {
        mu_disease <- 0
        sigma_disease <- -1
    }
    return(list(mu = df_mu, sigma = df_sigma, mu_disease = mu_disease, sigma_disease = sigma_disease))
}

renormalize_in_model_fit <- function(x) {
    x <- as.data.frame(x)
    li <- find_sigma_mu(x)
    df_mu <- li$mu
    df_sigma <- li$sigma
    col_max_radio <- ncol(df_mu)
    x[, 1:col_max_radio] <- (x[, 1:col_max_radio] - df_mu) / df_sigma
    if (li$sigma_disease > -0.5) {
        x$Age_at_disease <- (x$Age_at_disease - li$mu_disease) / li$sigma_disease
        li_age <- list(do = TRUE, mu = li$mu_disease, sigma = li$sigma_disease)
    } else {
        li_age <- list(do = FALSE)
    }
    return(list(x = x, df_mu = df_mu, df_sigma = df_sigma, li_age = li_age))
}

renormalize_in_model_pred <- function(newdata, modelFit) {
    li_norm <- modelFit$li_norm
    df_mu_former <- li_norm$df_mu
    df_sigma_former <- li_norm$df_sigma
    li_age <- li_norm$li_age
    n_line_new <- nrow(newdata)
    mat_mu <- as.matrix(df_mu_former[rep(1, n_line_new), ])
    mat_sigma <- as.matrix(df_sigma_former[rep(1, n_line_new), ])

    newdata <- as.data.frame(newdata)
    mat_mu <- as.matrix(find_sigma_mu(newdata)$mu)
    mat_sigma <- as.matrix(find_sigma_mu(newdata)$sigma)
    n_col_max_radio <- ncol(mat_mu)
    if (li_age[["do"]]) {
        mu <- li_age[["mu"]]
        sigma <- li_age[["sigma"]]
        newdata$Age_at_disease <- (newdata$Age_at_disease - mu) / sigma
    }

    newdata <- as.matrix(newdata)
    newdata[, 1:n_col_max_radio] <- (newdata[, 1:n_col_max_radio] - mat_mu) / mat_sigma
    return(newdata)
}

get_bloc_vec_liver <- function(vec, exist_temps = TRUE, different_clin = FALSE) {
    name_bloc <- rep("a", length(vec))
    index_bloc <- rep(0, length(vec))
    for (j in 1:length(vec)) {
        name_col <- vec[j]
        if (!substr(name_col, nchar(name_col) - 3, nchar(name_col)) %in% c("PORT", "ART_", "VEIN", "TARD") & exist_temps) {
            if (different_clin) {
                name_bloc[j] <- paste0("clinical_", j)
                index_bloc[j] <- -1
            } else {
                name_bloc[j] <- "clinical"
                index_bloc[j] <- -1
            }
        } else {
            classed <- FALSE
            if (grepl("firstorder", name_col)) {
                name_bloc[j] <- "first order"
                index_bloc[j] <- 1
                classed <- TRUE
            }
            if (grepl("shape", name_col)) {
                name_bloc[j] <- "shape"
                index_bloc[j] <- 2
                classed <- TRUE
            }
            if (grepl("glcm", name_col) | grepl("gldm", name_col) | grepl("glrl", name_col) | grepl("glsz", name_col) | grepl("ngtd", name_col)) {
                name_bloc[j] <- "texture"
                index_bloc[j] <- 3
                classed <- TRUE
            }
            if (classed == FALSE) {
                name_bloc[j] <- "clinical"
                index_bloc[j] <- -1
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
            index_name <- "ART_"
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
    index_name_no_clin <- index_name[-index_clin]
    index_variable <- rep(0, ncol(x))
    ordered_levels <- unique(index_name_no_clin)
    index_variable[-index_clin] <- as.integer(factor(index_name_no_clin, levels = ordered_levels))
    index_variable[index_clin] <- -1
    return(list(index_variable = index_variable, index_name = index_name))
}

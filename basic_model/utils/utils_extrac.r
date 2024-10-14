normal <- function(x) {
    if (is.character(x)) {
        y <- as.numeric(x)
        if (all(is.na(y))) {
            return(x)
        }
        x <- y
    }
    if (is.numeric(x)) {
        return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    } else {
        return(x)
    }
}


convert_to_num <- function(x) {
    if (is.character(x)) {
        y <- as.numeric(x)
        if (all(is.na(y))) {
            return(x)
        }
        x <- y
    }
    if (is.numeric(x)) {
        return(x)
    } else {
        return(x)
    }
}

change_dates <- function(tableau) {
    tableau <- data.frame(lapply(tableau, function(col) {
        if (inherits(col, "POSIXct")) {
            return(as.character(col))
        } else {
            return(col)
        }
    }))
    return(tableau)
}

change_genre <- function(elem) {
    if (is.na(elem)) {
        return(NA_integer_)
    }
    if (elem == "F") {
        return(1)
    }
    if (elem == "M") {
        return(0)
    }
}

name_product_terms <- function(x, small_groups = FALSE) {
    if (grepl("self", x) | grepl("crossed", x)) {
        if (grepl("self", x)) {
            type_col <- "self"
            extension <- substr(x, nchar(x) - 3, nchar(x))
            name_col <- substr(x, 6, nchar(x) - 10)
        } else {
            type_col <- "crossed"
            extension <- substr(x, nchar(x) - 8, nchar(x))
            name_col <- substr(x, 9, nchar(x) - 10)
        }
    } else {
        type_col <- "single"
        extension <- substr(x, nchar(x) - 3, nchar(x))
        name_col <- substr(x, 1, nchar(x) - 5)
    }
    if (small_groups) {
        return(paste(type_col, name_col, extension, sep = "_"))
    } else {
        paste(type_col, extension, sep = "_")
    }
}

assainir <- function(df, config_extrac, no_multivariate_col, crushed_cols) {
    time_inj <- config_extrac$time_inj
    col_diagnos <- grepl("diagnostics", colnames(df))
    df <- df[, !col_diagnos]
    col_to_examine <- !colnames(df) %in% c("slice_num", "classe_name", "temps_inj", "patient_num")
    indices_not_examined <- which(col_to_examine == FALSE)
    vec_sd <- as.vector(apply(df[, col_to_examine], 2, function(x) sd(x, na.rm = TRUE)))
    col_const <- vec_sd < 10^-6
    j <- 1
    # Remettre des FALSE pour les colonnes non examinées aux bons endroits
    while (j <= ncol(df)) {
        if (j %in% indices_not_examined) {
            col_const <- append(col_const, FALSE, after = j - 1)
        }
        j <- j + 1
    }
    # print(col_const)
    df <- df[, !col_const]

    # Dégager les mixtes si demandé
    if (config_extrac$kill_mixtes) {
        df <- df[df$classe_name != "Mixtes", ]
    }
    df <- df[df$temps_inj %in% time_inj, ]
    col_to_concatenate <- setdiff(colnames(df), c(no_multivariate_col, crushed_cols))
    df[["key"]] <- paste0(df$classe_name, "_", df$patient_num)
    return(df)
}

make_in_lines <- function(df_to_join, config_extrac, no_multivariate_col, crushed_cols, special_name_col = "_") {
    time_inj <- config_extrac$time_inj
    col_to_concatenate <- setdiff(colnames(df_to_join), c(no_multivariate_col, crushed_cols))
    lines <- lapply(unique(df_to_join$key), function(key) {
        df_to_join_key <- df_to_join[df_to_join$key == key, ]
        vec_line <- c()
        if (nrow(df_to_join_key) != length(time_inj)) {
            # print(paste("Key", key, "not found in sain data"))
            return(NULL)
        } else {
            for (time in time_inj) {
                vec_line_time <- unname(unlist(df_to_join_key[df_to_join_key$temps_inj == time, col_to_concatenate]))
                vec_line <- c(vec_line, vec_line_time)
            }
            vec_line <- c(key, vec_line)
            return(vec_line)
        }
    })
    lines <- Filter(Negate(is.null), lines)
    df_to_join_in_lines <- as.data.frame(do.call(rbind, lines))
    new_col_names <- c()

    # Parcourir chaque élément de config_extrac$new_names_time_inj
    for (time_name in config_extrac$new_names_time_inj) {
        # Pour chaque élément de col_to_concatenate, créer un nouveau nom de colonne
        for (col_name in col_to_concatenate) {
            new_col_names <- c(new_col_names, paste0(col_name, special_name_col, time_name))
        }
    }
    colnames(df_to_join_in_lines) <- c("key", new_col_names)
    return(df_to_join_in_lines)
}

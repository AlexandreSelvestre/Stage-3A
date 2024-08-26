f1_macro_summary <- function(data, lev = NULL, model = NULL) {
    f1_scores <- sapply(levels(data$obs), function(l) {
        if (sum(data$pred == l) == 0) {
            return(0)
        } else {
            precision <- sum(data$pred == l & data$obs == l) / sum(data$pred == l)
            recall <- sum(data$pred == l & data$obs == l) / sum(data$obs == l)
            f1 <- 2 * precision * recall / (precision + recall)
            if (precision < 0.0001 & recall < 0.0001) {
                f1 <- 0
            }
            return(f1)
        }
    })
    # print(f1_scores)
    f1_macro <- mean(f1_scores, na.rm = TRUE)
    out <- c(F1 = f1_macro)
    names(out) <- c("F1")
    out
}

weighted_accuracy_summary <- function(data, lev = NULL, model = NULL) {
    acc_scores <- sapply(levels(data$obs), function(l) {
        accuracy <- sum(data$pred == l & data$obs == l) / sum(data$obs == l)
        return(accuracy)
    })

    # Définir les poids pour chaque classe
    weights <- ifelse(levels(data$obs) == levels(data$obs)[which.min(table(data$obs))], 1, 1) # Le premier contrôle le minoritaire

    # Calculer la moyenne pondérée des précisions
    weighted_acc <- weighted.mean(acc_scores, w = weights)

    out <- c(Accuracy = weighted_acc)
    names(out) <- c("Weighted Accuracy")
    out
}


roc_calculus <- function(data, lev = NULL, model = NULL) {
    # Identify the minority class

    auc <- unname(twoClassSummary(data, lev, model)["ROC"])

    # Return a named vector of summary statistics
    out <- c(AUC = auc)
    # print(out)
    # print(data.frame(true = data$obs, pred = data$CHC))
    return(out)
}

custom_stats <- function(data, lev = NULL, model = NULL) {
    roc_data <- roc_calculus(data, lev, model)
    f1_data <- f1_macro_summary(data, lev, model)
    weighted_acc_data <- weighted_accuracy_summary(data, lev, model)
    return(c(roc_data, f1_data, weighted_acc_data))
}
return(custom_stats)


inverse_diag <- function(mat) {
    if (is.matrix(mat)) {
        diag(1 / (diag(mat) + 1e-10))
    } else {
        stop("L'entrée doit être une matrice.")
    }
}


crit_logistic <- function(x, y, beta, intercept, lambda) {
    options(digits = 12)
    criteria <- 1 / length(y) * (sum(y * (intercept + as.vector(x %*% beta)) - log(1 + exp(intercept + as.vector(x %*% beta))))) - lambda * norm(as.matrix(beta), type = "1")
    if (is.na(sum(as.vector(x %*% beta)))) {
        # print(x)
        # print(beta)
    }
    return(criteria)
}

soft_thresh <- function(x, lambda) {
    if (abs(x) < lambda) {
        return(list(value = 0, non_zero_soft_thresh = FALSE))
    } else {
        return(list(value = sign(x) * (abs(x) - lambda), non_zero_soft_thresh = TRUE))
    }
}


unregister_dopar <- function() {
    env <- foreach:::.foreachGlobals
    rm(list = ls(name = env), pos = env)
}



aggregate_prop <- function(formula, data, ...) {
    # On moyenne les pourcentages de non nullité mais on somme les importances
    # On rescale ensuite les importances (elles sont relatives)
    df_grouped <- aggregate(formula, data = data, FUN = sum)
    col_grp <- as.character(formula)[3] # Sur quoi on groupe
    col_values <- as.character(formula)[2] # Overall
    formula_str <- paste("Percentage ~", col_grp)
    formula_obj <- as.formula(formula_str)
    df_grouped_perc <- aggregate(formula_obj, data = data, FUN = mean)
    df_grouped$Percentage <- df_grouped_perc$Percentage
    df_grouped[[col_values]] <- df_grouped[[col_values]] / sum(df_grouped[[col_values]])
    df_grouped <- subset(df_grouped, df_grouped[[col_values]] > 0.0001 | df_grouped$Percentage > 0.01)
    return(df_grouped)
}

melt_mine <- function(df) {
    new_df <- as.data.frame(matrix(nrow = 0, ncol = 2))
    colnames(new_df) <- c("variable", "value")
    for (name_col in colnames(df)) {
        local_df <- data.frame(variable = name_col, value = df[[name_col]])
        new_df <- rbind(new_df, local_df)
    }
    return(new_df)
}

calculer_quantile <- function(serie, nouveau_nombre) {
    serie_mise_a_jour <- c(serie, nouveau_nombre)
    serie_mise_a_jour <- sort(serie_mise_a_jour)
    index_nouveau_nombre <- match(nouveau_nombre, serie_mise_a_jour)
    quantile_nouveau_nombre <- index_nouveau_nombre / length(serie_mise_a_jour)
    return(quantile_nouveau_nombre)
}

my_quantile <- function(serie, proba) {
    ordered <- sort(serie)
    index_quant <- floor(proba * length(ordered))
    return(ordered[index_quant])
}


import_folder <- function(path) {
    files <- list()
    files <- append(files, list.files(path = path, pattern = "\\.r$", full.names = TRUE))
    files <- append(files, list.files(path = path, pattern = "\\.R$", full.names = TRUE))
    lapply(files, source)
}

glue_mats <- function(li_beta_matrix) {
    # Coller côte à côte des matrices de nombres de lignes variables en remplissant de zéros là où c'est nécessaire
    # Étape 1: Trouver le nombre maximal de lignes
    max_rows <- max(sapply(li_beta_matrix, nrow))

    # Étape 2: Ajuster chaque matrice pour qu'elle ait le nombre maximal de lignes
    adjusted_matrices <- lapply(li_beta_matrix, function(matrix) {
        if (nrow(matrix) < max_rows) {
            # Calculer le nombre de lignes à ajouter
            rows_to_add <- max_rows - nrow(matrix)
            # Créer une matrice de zéros avec le bon nombre de lignes et de colonnes
            zero_rows <- matrix(0, nrow = rows_to_add, ncol = ncol(matrix))
            # Ajouter les lignes de zéros à la matrice originale
            matrix <- rbind(matrix, zero_rows)
        }
        return(matrix)
    })

    # Étape 3: Combiner toutes les matrices ajustées côte à côte
    beta_matrix_combined <- do.call(cbind, adjusted_matrices)
    return(beta_matrix_combined)
}


silent_run <- function(func, ...) {
    # Créer un fichier temporaire pour capturer les sorties
    temp_file <- tempfile()

    # Ouvrir une connexion au fichier temporaire
    temp_conn <- file(temp_file, open = "wt")

    # Rediriger les messages de sortie standard vers la connexion au fichier temporaire
    sink(temp_conn)

    # Rediriger les messages d'avertissement vers la connexion au fichier temporaire
    sink(temp_conn, type = "message")

    # Exécuter la fonction en mode silencieux avec les arguments supplémentaires
    result <- func(...)

    # Réinitialiser les sorties à leurs destinations d'origine
    sink()
    sink(type = "message")

    # Fermer la connexion au fichier temporaire
    close(temp_conn)

    # Supprimer le fichier temporaire
    unlink(temp_file)

    return(result)
}

complete_orthonormal_basis <- function(vec) {
    # Normaliser le vecteur d'entrée
    vec <- vec / sqrt(sum(vec^2))

    # Créer une matrice en ajoutant des vecteurs aléatoires
    n <- length(vec)
    mat <- matrix(rnorm(n * n), n, n)
    mat[, 1] <- vec
    # Appliquer la décomposition QR
    qr_decomp <- qr(mat)
    Q <- qr.Q(qr_decomp)

    return(Q)
}

# vec <- c(1, 1, 2, 3)
# vec <- vec / sqrt(sum(vec^2))
# Q <- complete_orthonormal_basis(vec)
# print(Q)
# print(vec)

li_caret_simple <- list()

li_caret_simple$library <- "glmnet"

li_caret_simple$type <- "Classification"

li_caret_simple$parameters <- data.frame(parameter = c("lambda"), class = c("numeric"), label = c(
    "la valeur du paramètre lambda"
))

create_grid_simple <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        lambda <- log(seq(exp(0), exp(0.05), length.out = len + 1)[2:(len + 1)])
    } else {
        lambda <- log(runif(len, min = exp(0.001), max = exp(0.05)))
    }
    data_frame_grid <- expand.grid(lambda = lambda)
    return(data_frame_grid)
}

li_caret_simple$grid <- create_grid_simple

fit_simple <- function(x, y, wts, param, lev, last, weights, classProbs, classe_1 = NULL) {
    lambda <- param$lambda
    y <- factor(y, levels = lev)
    weights <- rep(1, length(y))
    regression <- glmnet:::glmnet.fit(
        x = as.matrix(x), y = y, family = binomial(), alpha = 1, lambda = lambda, intercept = TRUE, maxit = 10^7, weights / length(y)
    )
    # print(regression$converged)
    beta <- as.numeric(regression$beta)
    intercept <- regression$a0
    # print(beta)
    return(list(beta = beta, intercept = intercept, lev = lev, classe_0 = "CHC", classe_1 = "CCK"))
}

li_caret_simple$fit <- fit_simple

li_caret_simple$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # newdata est bien une matrice


    beta <- modelFit$beta
    intercept <- modelFit$intercept
    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0

    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    predicted_labels <- ifelse(proba > 0.5, classe_1, classe_0)
    return(predicted_labels)
}

li_caret_simple$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    beta <- modelFit$beta
    intercept <- modelFit$intercept

    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0
    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    str_1 <- as.character(classe_1)
    str_0 <- as.character(classe_0)
    return(setNames(data.frame(1 - proba, proba), c(str_0, str_1)))
}

li_caret_simple$loop <- NULL

logistic_importance_beta_intra <- function(results, li_imp, show_model) {
    print(show_model)
    if (TRUE) {
        print(paste("meilleur modèle: logistic avec lambda =", results$bestTune$lambda))
    }
    beta <- results$final_model$beta
    beta_global <- li_imp$beta_global
    beta_non_zero <- li_imp$beta_non_zero
    beta <- results$finalModel$beta
    li_imp$beta_global <- beta_global + abs(beta)
    non_zero_local <- ifelse(abs(beta) > 1e-5, 1, 0)
    li_imp$beta_non_zero <- beta_non_zero + non_zero_local
    return(li_imp)
}

logistic_importance_beta_extra <- function(li_imp, df) {
    beta_non_zero <- li_imp$beta_non_zero
    beta_global <- li_imp$beta_global
    beta_non_zero <- 100 * beta_non_zero / n_runs
    beta_global <- beta_global / n_runs
    names_df <- colnames(df[setdiff(colnames(df), c("Tumeur"))])
    names(beta_global) <- names_df
    names(beta_non_zero) <- names_df

    beta_df <- data.frame(
        name = gsub("[_.]", " ", names(beta_global)), # Remplacer les underscores et les points par des espaces
        value = beta_global,
        percentage = beta_non_zero
    )

    # Réorganiser les niveaux du facteur name en fonction des valeurs de beta_global
    beta_df$name <- factor(beta_df$name, levels = beta_df$name[order(beta_df$value, decreasing = FALSE)])

    # Créer le barplot avec les axes inversés et les pourcentages en vert fluo
    plot <- ggplot(beta_df, aes(x = name, y = value)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        geom_text(aes(label = paste0(percentage, "%")),
            color = "limegreen",
            hjust = -0.1
        ) +
        labs(title = "", x = "", y = "") +
        expand_limits(y = max(beta_df$value) * 1.1) + # Ajuster les limites de l'axe des y
        theme(axis.text.y = element_text(size = 12)) # Augmenter la taille du texte des noms des barres


    ggsave("plots/barplot_beta_global.png", plot)
}

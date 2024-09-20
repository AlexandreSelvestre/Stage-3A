setClass("logistic_grp",
    contains = "apply_model",
    slots = representation(
        lambda_min = "numeric",
        tuneLength = "numeric",
        index_type = "character",
        regression = "logical",
        lambda = "numeric",
        tuneGrid = "data.frame"
    )
)



li_caret_sgl <- list()

li_caret_sgl$library <- "gglasso"

li_caret_sgl$type <- "Classification"

li_caret_sgl$parameters <- data.frame(parameter = c("lambda"), class = c("numeric"), label = c(
    "la valeur du paramètre lambda"
))

create_grid_sgl <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        lambda <- seq(0, 1, length.out = len)
    } else {
        lambda <- runif(len, min = 0, max = 1)
    }
    data_frame_grid <- expand.grid(lambda = lambda)
    return(data_frame_grid)
}

li_caret_sgl$grid <- create_grid_sgl

fit_sgl <- function(x, y, wts, param, lev, last, weights, classProbs, index, k_smote, sampling_choice, index_variable, is_binary, index_bloc) {
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

    lambda <- param$lambda
    y_numeric <- ifelse(y == classe_min, 1, -1)
    fited <- gglasso::gglasso(
        as.matrix(x), y_numeric,
        group = index, loss = "logit", lambda = lambda
    )
    fited$classProbs <- classProbs
    fited$lev <- lev
    fited$x <- x
    fited$li_norm <- li_norm
    fited$beta <- as.vector(fited$beta)
    fited$intercept <- fited$b0
    fited$classe_min <- classe_min
    fited$classe_maj <- classe_maj
    return(fited)
}
li_caret_sgl$fit <- fit_sgl

li_caret_sgl$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta_fited <- as.vector(modelFit$beta)
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta_fited + modelFit$intercept
    })
    proba <- 1 / (1 + exp(-value))
    predicted_labels <- ifelse(proba > 0.5, classe_min, classe_maj)

    return(predicted_labels)
}


li_caret_sgl$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    beta_fited <- as.vector(modelFit$beta)
    value <- apply(newdata, 1, function(ligne) {
        ligne %*% beta_fited + modelFit$intercept
    })
    proba_class_1 <- 1 / (1 + exp(-value))
    proba_class_0 <- 1 - proba_class_1
    str_min <- as.character(classe_min)
    str_maj <- as.character(classe_maj)
    return(setNames(data.frame(proba_class_0, proba_class_1), c(str_maj, str_min)))
}

li_caret_sgl$loop <- NULL



setMethod("train_method", "apply_model", function(object) {
    # Créer la liste index utilisée pour les groupes (petits groupes)
    vec_names <- colnames(object@train_cols[, object@col_x])
    df_names <- data.frame(name_cov = vec_names)

    if (object@index_type == "var") {
        index <- object@index_variable
    }
    if (object@index_type == "mode") {
        index <- object@index_mode
    }
    if (object@index_type == "bloc") {
        index <- object@index_bloc
    }


    for (i in 1:length(index)) {
        if (index[i] < -0.5) {
            index[i] <- max(index) + 1
        }
    }
    index <- as.integer(as.factor(index))


    x <- copy(object@train_cols[, object@col_x])
    y <- copy(object@y_train)
    li_norm_loc <- renormalize_in_model_fit_index_mode(x, object@index_variable, object@is_binary)
    x <- li_norm_loc$new_x
    classe_min <- names(which.min(table(y)))
    classe_maj <- setdiff(levels(y), classe_min)



    if (object@sampling == "smote") {
        li <- apply_smote(x, y, k_smote)
        x <- li$x
        y <- li$y
    }
    if (object@sampling == "up") {
        li <- apply_boot(x, y)
        x <- li$x
        y <- li$y
    }

    ### Générer les alphas suboptimaux (vis-à-vis de la liste des lambda) de manière aléatoire ou non
    y_numeric <- ifelse(y == classe_min, 1, -1)
    seq_lambda <- gglasso::gglasso(
        as.matrix(x), y_numeric,
        group = index, loss = "logit",
        nlambda = object@tuneLength, lambda.factor = object@lambda_min
    )$lambda
    tuneGrid <- data.frame(lambda = seq_lambda)

    print("tuneGrid created")
    object@tuneGrid <- tuneGrid
    start_time <- Sys.time()
    if (object@parallel$do) {
        numCores <- detectCores()
        cl <- makePSOCKcluster(object@parallel$n_process)
        registerDoParallel(cl)
        clusterEvalQ(cl, {
            files <- list.files("./utils", full.names = TRUE, pattern = "\\.r$")
            for (file in files) {
                source(file)
            }
        })
    }
    object@model <- caret::train(
        y = object@y_train, x = as.matrix(object@train_cols[, object@col_x]),
        method = li_caret_sgl, trControl = object@cv, metric = "AUC", tuneGrid = tuneGrid, index = index, k_smote = object@k_smote, sampling_choice = object@sampling,
        index_variable = object@index_variable, is_binary = object@is_binary, index_bloc = object@index_bloc
    )
    if (object@parallel$do) {
        stopCluster(cl)
    }
    return(object)
})


setMethod("importance_method", "apply_model", function(object) {
    print("start imp")
    best_beta <- object@model$finalModel$beta
    object@beta_final <- best_beta
    variable_importance <- data.frame(Variable = colnames(object@train_cols[, object@col_x]), Overall = abs(best_beta))
    object@li_df_var_imp <- variable_importance
    # variable_importance <- variable_importance[order(-variable_importance$Overall), ]
    variable_importance$Group <- object@name_mode
    variable_importance$small_Group <- object@name_variable


    variable_importance_grouped <- aggregate(Overall ~ Group, data = variable_importance, FUN = mean)
    variable_importance_small_grouped <- aggregate(Overall ~ small_Group, data = variable_importance, FUN = mean)
    variable_importance <- subset(variable_importance, Overall > 0.0001)
    variable_importance_grouped <- subset(variable_importance_grouped, Overall > 0.0001)
    variable_importance_small_grouped <- subset(variable_importance_small_grouped, Overall > 0.0001)


    image <- ggplot2::ggplot(variable_importance, aes(x = reorder(Variable, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance") +
        theme(axis.text.y = element_text(size = 2)) # Adjust the size as needed
    ggsave(paste0("plots/logistic_grp/importance", "_", object@id_term, ".png"), image, width = 10, height = 15)


    image <- ggplot(variable_importance_grouped, aes(x = reorder(Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistic_grp/big_groups_importance", "_", object@id_term, ".png"), image)


    image <- ggplot(variable_importance_small_grouped, aes(x = reorder(small_Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistic_grp/small_groups_importance", "_", object@id_term, ".png"), image)


    df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)

    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        geom_boxplot() +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0("plots/logistic_grp/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)

    return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    return(object@li_df_var_imp)
})

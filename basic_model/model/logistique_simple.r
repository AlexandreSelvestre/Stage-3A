library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(DMwR)
library(themis)


setClass("logistique_simple",
    contains = "apply_model",
    slots = representation(
        lambda_min = "numeric",
        lambda_max = "numeric",
        tuneLength = "integer",
        weights = "list",
        regression = "logical",
        lambda = "numeric"
    )
)

# simple <- FALSE

li_caret_simple <- list()

li_caret_simple$library <- "varclust"

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

fit_simple <- function(x, y, wts, param, lev, last, weights_dict, classProbs, k_smote, sampling_choice, index_variable, index_bloc, is_binary, classe_1 = NULL, penalty_adapt) {
    li_norm <- renormalize_in_model_fit_index_mode(x, index_variable, index_bloc, is_binary)
    ######## THE GOOD LINE FOR NORMALIZATION
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



    # FAIRE APRES SMOTE
    lambda <- param$lambda
    weights <- numeric(length(y))
    weights[y == classe_maj] <- weights_dict[[classe_maj]]
    weights[y == classe_min] <- weights_dict[[classe_min]]
    # write_xlsx(cbind(y, x), "..\\data\\no.xlsx")
    y_numeric <- convert_y(y, classe_1)
    if (is.null(classe_1)) {
        classe_1 <- classe_min
    }
    classe_0 <- setdiff(levels(y), classe_1)

    if (penalty_adapt) {
        different_blocs <- unique(index_bloc[index_bloc > -0.5])
        li_pen_per_bloc <- lapply(different_blocs, function(l_num) {
            nb_var_bloc <- length(which(index_bloc == l_num))
            penalty_unscaled <- 1 / nb_var_bloc
            return(penalty_unscaled)
        })
        names(li_pen_per_bloc) <- as.character(different_blocs)
        penalty.factor <- sapply(seq_len(ncol(x)), function(j) {
            if (index_bloc[j] < -0.5) {
                return(1)
            } else {
                l_char <- as.character(index_bloc[j])
                return(li_pen_per_bloc[[l_char]])
            }
        })
    } else {
        penalty.factor <- rep(1, ncol(x))
    }

    regression <- glmnet:::glmnet.fit(
        x = as.matrix(x), y = y_numeric, family = binomial(), alpha = 1,
        weights = weights / length(y_numeric), lambda = lambda, intercept = TRUE, maxit = 10^7, penalty.factor = penalty.factor
    )
    print(regression$converged)

    beta <- as.numeric(regression$beta)
    intercept <- regression$a0
    return(list(beta = beta, intercept = intercept, lev = lev, li_norm = li_norm, classe_min = classe_min, classe_maj = classe_maj, classe_1 = classe_1, classe_0 = classe_0, penalty.factor = penalty.factor))
}

li_caret_simple$fit <- fit_simple

li_caret_simple$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # newdata est bien une matrice
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)

    beta <- modelFit$beta
    intercept <- modelFit$intercept

    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    predicted_labels <- ifelse(proba > 0.5, classe_1, classe_0)
    return(predicted_labels)
}

li_caret_simple$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    beta <- modelFit$beta
    intercept <- modelFit$intercept
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    classe_1 <- modelFit$classe_1
    classe_0 <- modelFit$classe_0
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)
    # print(newdata[(nrow(newdata) - 1):nrow(newdata), ])
    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    str_1 <- as.character(classe_1)
    str_0 <- as.character(classe_0)
    return(setNames(data.frame(1 - proba, proba), c(str_0, str_1)))
}

li_caret_simple$loop <- NULL



setMethod("train_method", "logistique_simple", function(object) {
    tuneGrid <- expand.grid(lambda = exp(log(10) * seq(log10(object@lambda_min), log10(object@lambda_max), length = object@tuneLength)))
    if (object@parallel$do) {
        numCores <- detectCores()
        cl <- makePSOCKcluster(object@parallel$n_process)
        # cl <- makePSOCKcluster(2)
        registerDoParallel(cl)
        clusterEvalQ(cl, {
            files <- list.files("./utils", full.names = TRUE, pattern = "\\.r$")
            for (file in files) {
                source(file)
            }
        })
    }
    print(dim(object@train_cols[, object@col_x]))
    object@model <- caret::train(
        y = object@y_train, x = object@train_cols[, object@col_x],
        method = li_caret_simple, trControl = object@cv, metric = "AUC",
        tuneLength = 8, tuneGrid = tuneGrid,
        weights_dict = object@weights, k_smote = object@k_smote, sampling_choice = object@sampling,
        index_variable = object@index_variable, index_bloc = object@index_bloc, is_binary = object@is_binary,
        classe_1 = object@classe_1, penalty_adapt = object@penalty_adapt
    )
    if (object@parallel$do) {
        stopCluster(cl)
    }
    print("done")
    # object@coeffi <- coef(object@model$finalModel, s = object@lambda)

    object@lambda <- object@model$bestTune$lambda
    return(object)
})



setMethod("importance_method", "apply_model", function(object) {
    object@beta_final <- object@model$finalModel$beta
    vec_importance <- abs(object@model$finalModel$beta) ##### SUITE
    variable_importance <- data.frame(Variable = object@col_x, Overall = vec_importance)
    object@li_df_var_imp <- variable_importance

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
        ggtitle("Variable Importance")
    # + theme(axis.text.y = element_text(size = 3)) # Adjust the size as needed
    ggsave(paste0("plots/logistique_simple/importance", "_", object@id_term, ".png"), image)

    image <- ggplot(variable_importance_grouped, aes(x = reorder(Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistique_simple/big_groups_importance", "_", object@id_term, ".png"), image)


    image <- ggplot(variable_importance_small_grouped, aes(x = reorder(small_Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    ggsave(paste0("plots/logistique_simple/small_groups_importance", "_", object@id_term, ".png"), image)
    #######################################

    image <- df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)
    # print(df_long)
    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0("plots/logistique_simple/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)
    return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    object@beta_final <- object@model$finalModel$beta
    ### afficher les pénalités par bloc
    # penalty.factor <- object@model$finalModel$penalty.factor
    # penalty_vec <- abs(object@beta_final) * penalty.factor
    # penalty_term_per_bloc <- lapply(unique(object@index_bloc), function(l_num) {
    #     return(sum(penalty_vec[object@index_bloc == l_num]))
    # })
    # names(penalty_term_per_bloc) <- as.character(unique(object@index_bloc))
    # print(penalty_term_per_bloc)
    ###
    vec_importance <- abs(object@model$finalModel$beta)
    variable_importance <- data.frame(Variable = object@col_x, Overall = vec_importance)
    object@li_df_var_imp <- variable_importance
    # print(tail(object@beta_final))
    return(object@li_df_var_imp)
})

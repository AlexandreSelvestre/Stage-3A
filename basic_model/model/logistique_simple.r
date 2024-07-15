library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(DMwR)
library(themis)


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

fit_simple <- function(x, y, wts, param, lev, last, weights_dict, classProbs, k_smote, do_smote, index_variable, is_binary) {
    li_norm <- renormalize_in_model_fit_index_mode(x, index_variable, is_binary)
    ######## THE GOOD LINE FOR NORMALIZATION
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



    # FAIRE APRES SMOTE
    lambda <- param$lambda
    weights <- numeric(length(y))
    weights[y == classe_maj] <- weights_dict[[classe_maj]]
    weights[y == classe_min] <- weights_dict[[classe_min]]
    # write_xlsx(cbind(y, x), "..\\data\\no.xlsx")

    y_numeric <- ifelse(y == classe_min, 1, 0)

    regression <- glmnet:::glmnet.fit(
        x = as.matrix(x), y = y_numeric, family = binomial(), alpha = 1,
        weights = weights / length(y_numeric), lambda = lambda, intercept = TRUE, maxit = 1e8
    )

    beta <- as.numeric(regression$beta)
    intercept <- regression$a0
    return(list(beta = beta, intercept = intercept, lev = lev, li_norm = li_norm, classe_min = classe_min, classe_maj = classe_maj))
}

li_caret_simple$fit <- fit_simple

li_caret_simple$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    # newdata est bien une matrice
    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)

    beta <- modelFit$beta
    intercept <- modelFit$intercept


    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    predicted_labels <- ifelse(proba > 0.5, classe_min, classe_maj)
    return(predicted_labels)
}

li_caret_simple$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    beta <- modelFit$beta
    intercept <- modelFit$intercept
    classe_min <- modelFit$classe_min
    classe_maj <- modelFit$classe_maj

    df_mu <- modelFit$li_norm$df_mu
    df_sigma <- modelFit$li_norm$df_sigma
    newdata <- renormalize_in_model_pred_index_mode(newdata, df_mu, df_sigma)



    proba <- 1 / (1 + exp(-as.vector(newdata %*% beta) - intercept))
    str_min <- as.character(classe_min)
    str_maj <- as.character(classe_maj)
    return(setNames(data.frame(1 - proba, proba), c(str_maj, str_min)))
}

li_caret_simple$loop <- NULL



setMethod("train_method", "apply_model", function(object) {
    tuneGrid <- expand.grid(lambda = exp(seq(log(object@lambda_min), log(object@lambda_max), length = object@tuneLength)))
    if (object@do_PCA) {
        object@model <- caret::train(
            y = object@train_cols$classe_name, x = as.matrix(object@train_cols[, object@col_x]),
            method = "glmnet", trControl = object@cv, metric = "AUC",
            tuneLength = 8, family = "binomial", tuneGrid = tuneGrid, preProcess = "pca",
            weights = weights, index = object@index_variable, is_binary = object@is_binary
        )
    } else {
        if (object@do_parallel) {
            numCores <- detectCores()
            cl <- makePSOCKcluster(numCores - 1)
            # cl <- makePSOCKcluster(2)
            registerDoParallel(cl)
            clusterEvalQ(cl, {
                files <- list.files("./utils", full.names = TRUE, pattern = "\\.r$")
                for (file in files) {
                    source(file)
                }
            })
        }
        object@model <- caret::train(
            y = object@y_train, x = object@train_cols[, object@col_x],
            method = li_caret_simple, trControl = object@cv, metric = "AUC",
            tuneLength = 8, tuneGrid = tuneGrid,
            weights_dict = object@weights, k_smote = object@k_smote, do_smote = object@do_smote,
            index_variable = object@index_variable, is_binary = object@is_binary
        )
        # x, y, wts, param, lev, last, weights_dict, classProbs, k_smote, do_smote, index, is_binary
        if (object@do_parallel) {
            stopCluster(cl)
        }
        print("done")
        # object@coeffi <- coef(object@model$finalModel, s = object@lambda)
    }
    object@lambda <- object@model$bestTune$lambda
    return(object)
})



setMethod("get_results", "apply_model", function(object) {
    # print(paste(c("test_set puis train", dim(object@test_set), dim(object@train_set)), collapse = "x"))
    if (object@do_PCA) {
        test_explic_pca <- predict(object@model$preProcess, object@test_set[, object@col_x])
        train_explic_pca <- predict(object@model$preProcess, object@train_cols[, object@col_x])
        object@predictions <- as.vector(predict(object@model$finalModel, newx = as.matrix(test_explic_pca), s = object@lambda, type = "class"))
        object@predictions_proba <- predict(object@model$finalModel, newx = as.matrix(test_explic_pca), type = "response", s = object@lambda)
        object@predictions_proba <- data.frame("1-probability" = 1 - object@predictions_proba, "probability" = object@predictions_proba) # standardiser le format
        object@predictions_train_proba <- predict(object@model$finalModel, newx = as.matrix(train_explic_pca), type = "response", s = object@lambda)
        object@predictions_train_proba <- data.frame("1-probability" = 1 - object@predictions_train_proba, "probability" = object@predictions_train_proba) # standardiser le format
    } else {
        object@predictions <- as.vector(predict(object@model, newdata = as.matrix(object@test_set[, object@col_x])))
        object@predictions_proba <- predict(object@model, newdata = as.matrix(object@test_set[, object@col_x]), type = "prob")
        object@predictions_train_proba <- predict(object@model, newdata = as.matrix(object@train_cols[, object@col_x]), type = "prob")
    }

    # print(object@predictions_proba)
    return(object)
})



setMethod("importance_method", "apply_model", function(object) {
    if (object@do_PCA == FALSE) {
        vec_importance <- abs(object@model$finalModel$beta) ##### SUITE
        variable_importance <- data.frame(Variable = object@col_x, Overall = vec_importance)
        object@li_df_var_imp <- variable_importance
        if (object@include_products == FALSE) {
            variable_importance$Group <- object@name_mode
            variable_importance$small_Group <- object@name_variable
        }

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
    }
    image <- df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)
    object@li_box_plots[[object@id_term]] <- df_long
    # print(df_long)
    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0("plots/logistique_simple/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)
    return(object)
})

setMethod("get_df_imp", "apply_model", function(object) {
    return(object@li_df_var_imp)
})

# [1] "actuelle moyenne AUC test 0.737847222222222 ite: 40"
# [1] "actuelle moyenne AUC val 0.801590909090909 ite: 40"
# [1] "actuelle moyenne Accuracy test 0.674652777777778 ite: 40"
# [1] "actuelle somme des confusion matrix ite: 40 :"
# [1] "actuelle moyenne du macro F1 test 0.655859353795209 ite: 40"
# [1] "actuelle moyenne du F1 CCK test 0.455854978355 ite: 40"
# [1] "actuelle moyenne du F1 CHC test 0.855863729235 ite: 40"




# [1] "Voilà la matrice de confusion sommée"
#           Reference
# Prediction CCK CHC
#        CCK  83 122
#        CHC  77 598

# [1] "voilà la matrice de confusion en pourcentage: ... % de CCK ont été bien classés vs ... % mal classés"

#     CCK                  CHC
# CCK "
# CHC " ##### A ecrire: il a sorti n'importe quoi


# [1] "La balanced accuracy sur l'échantillon total vaut: 0.645401987353207"
# le f1 score CCK sur l'échantillon entier vaut: 0.454794520548
# le f1 score CHC sur l'échantillon entier vaut: 0.857347670251
# [1] "Le f1 macro sur l'échantillon total vaut: 0.656071095399421"

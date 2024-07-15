library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(DMwR)
library(themis)

li_state <- list(n_current = 0, lambda_current = -1, var_in_clust = list(), li_pca = list(), li_directions = list(), vec_dim = list(), segmentation = list(), new_x = NULL)

li_caret_select <- list()

li_caret_select$library <- c("varclust", "stats", "glm2")

li_caret_select$type <- "Classification"

li_caret_select$parameters <- data.frame(parameter = c("lambda", "n_centro"), class = c("numeric", "integer"), label = c(
    "le lambda du lasso",
    "Le nombre de centroïdes utilisés dans les k means"
))

create_grid_select <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        n_centro <- seq(from = 5, to = 30, length.out = len)
        lambda <- seq(from = 0.01, to = 1, length.out = 20)
    } else {
        n_centro <- sample(5:30, len)
        lambda <- runif(20, 0.01, 1)
    }
    data_frame_grid <- expand.grid(n_centro = n_centro, lambda = lambda)
    return(data_frame_grid)
}

li_caret_select$grid <- create_grid_select

li_caret_select$fit <- function(x, y, wts, param, lev, last, weights_dict, classProbs, dim_max, k_smote, do_smote) {
    li_state$state <<- FALSE

    li_norm <- renormalize_in_model_fit(x) ######## THE GOOD LINE FOR NORMALIZATION
    x <- li_norm$x
    y_numeric <- ifelse(y == "CCK", 1, 0) # CCK donne 1 CHC donne 0
    n_centro <- param$n_centro
    col_to_cluster <- setdiff(colnames(x), c("Gender"))

    if (do_smote) {
        li <- apply_smote(x, y, k_smote)
        x <- li$x
        y <- li$y
    } else {
        li <- apply_boot(x, y)
        x <- li$x
        y <- li$y
    }


    y_numeric <- ifelse(y == "CCK", 1, 0) # CCK donne 1 CHC donne 0

    weights <- numeric(length(y))
    weights[y == "CHC"] <- weights_dict[["CHC"]]
    weights[y == "CCK"] <- weights_dict[["CCK"]]
    n_centro <- param$n_centro
    col_to_cluster <- setdiff(colnames(x), c("Gender"))

    # Si on change de nb de variables ou alors ni les variables ni le lambda n'ont changé (i.e nouveau fold), alors on recalcule les pca
    if (li_state$n_current != n_centro | abs(li_state$lambda_current - param$lambda) < 0.001) {
        ######### Kmeans itéré avec varclust.
        x_to_cluster <- x[, col_to_cluster]
        rownames(x_to_cluster) <- NULL
        mlcc.res <- varclust::mlcc.reps(as.matrix(x_to_cluster), numb.clusters = n_centro, numb.runs = 10, max.dim = dim_max, scale = FALSE)
        # print(as.matrix(x_to_cluster[1:5, 1:2]))
        segmentation <- mlcc.res$segmentation
        li_pca_false <- mlcc.res$basis
        # print(segmentation)
        li_directions <- list()

        ## Enregistrer les infos importantes du clustering
        var_in_clust <- list()
        vec_dim <- c()
        li_pca <- list()
        for (i in 1:n_centro) {
            var_in_clust[[i]] <- col_to_cluster[segmentation == i]
            # print(var_in_clust[[i]])
            vec_dim <- c(vec_dim, ncol(li_pca_false[[i]]))
            li_directions[[i]] <- stats::prcomp(x[, var_in_clust[[i]]], scale. = FALSE, center = FALSE, rank. = vec_dim[i])$rotation
            prod_scal <- as.matrix(x[, var_in_clust[[i]]]) %*% li_directions[[i]]
            li_pca[[i]] <- prod_scal
        }



        ## Créer le nouveau dataframe
        new_x <- as.data.frame(matrix(ncol = 0, nrow = nrow(x)))
        start_index <- 1
        for (i in 1:n_centro) {
            new_x <- cbind(new_x, li_pca[[i]])
            names <- c()
            for (j in 1:ncol(li_pca[[i]])) {
                names <- c(names, paste("Cluster_number", i, "component", j, sep = "_"))
            }
            colnames(new_x)[start_index:ncol(new_x)] <- names
            start_index <- start_index + vec_dim[[i]]
            # print(colnames(new_x))
        }
        new_x <- cbind(new_x, x[c("Gender")])
        # li_state <- list(n_current = 0, var_in_clust = list(), li_pca = list(), li_directions = list(), vec_dim = list(), segmentation = list(), new_x = NULL)
        li_state$n_current <<- n_centro
        li_state$var_in_clust <<- var_in_clust
        li_state$li_pca <<- li_pca
        li_state$li_directions <<- li_directions
        li_state$vec_dim <<- vec_dim
        li_state$segmentation <<- segmentation
        li_state$new_x <<- new_x
        li_state$lambda_current <<- param$lambda
    } else {
        var_in_clust <- li_state$var_in_clust
        li_pca <- li_state$li_pca
        li_directions <- li_state$li_directions
        vec_dim <- li_state$vec_dim
        segmentation <- li_state$segmentation
        new_x <- li_state$new_x
        li_state$lambda_current <<- param$lambda
    }



    df_y <- data.frame(y = y_numeric)
    x_formula <- cbind(new_x, df_y)

    ### Faire un modèle descendant sur les variables identifiées par varclust
    # print(colnames(x_formula))
    # print(x_formula[1:10, ])
    # print(x_formula)

    # cor_matrix <- cor(as.matrix(new_x))
    # print(cor_matrix)
    # png("correlation_heatmap.png")
    # heatmap(cor_matrix)
    # dev.off()


    ######### DESCENDANT

    # regression <- stepwise(y ~ ., x_formula, type = "logit", strategy = "backward", metric = "AIC", test_method_glm = "Rao", best_n = 1, n_ite = 100)

    # coefficients <- regression[[4]][["Estimate"]]
    # intercept <- as.numeric(coefficients[1])
    # variable_kept <- regression[[4]][["Variable"]][2:length(coefficients)]
    # beta_reduced <- coefficients[2:length(coefficients)]
    # beta <- rep(NA, ncol(new_x))
    # for (i in 1:ncol(new_x)) {
    #     name <- colnames(new_x)[i]
    #     if (name %in% variable_kept) {
    #         beta[i] <- as.numeric(beta_reduced[variable_kept == name])
    #     } else {
    #         beta[i] <- 0
    #     }
    # }

    # print(variable_kept)
    # gamma <- beta
    # names(gamma) <- colnames(new_x)
    # print(gamma)

    ########### SIMPLE FIT

    # regression <- stats::glm(y ~ ., family = binomial(link = "logit"), x_formula, control = glm.control(epsilon = 1e-8, maxit = 100, trace = FALSE))
    # beta <- regression$coefficients[2:length(regression$coefficients)]
    # intercept <- regression$coefficients[1]

    ############ LASSO

    lambda <- param$lambda
    regression <- glmnet:::glmnet.fit(
        x = as.matrix(new_x), y = y_numeric, family = binomial(), alpha = 1,
        weights = weights / length(y_numeric), lambda = lambda, intercept = TRUE, maxit = 1e7
    )
    beta <- as.numeric(regression$beta)
    intercept <- regression$a0
    # print(beta)

    # print(length(beta))
    # print(dim(new_x)[2])
    # print(beta)

    return(list(
        beta = beta, intercept = intercept, lev = lev, li_directions = li_directions, vec_dim = vec_dim, segmentation = segmentation,
        var_in_clust = var_in_clust, n_centro = n_centro, new_names = colnames(new_x), li_norm = li_norm
    ))
}

li_caret_select$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- renormalize_in_model_pred(newdata, modelFit) ######## THE GOOD LINE FOR NORMALIZATION
    beta <- modelFit$beta
    intercept <- modelFit$intercept
    li_directions <- modelFit$li_directions
    n_centro <- modelFit$n_centro
    # vec_dim <- modelFit$vec_dim
    # segmentation <- modelFit$segmentation
    var_in_clust <- modelFit$var_in_clust

    df_x <- data.frame(Gender = as.data.frame(newdata)$Gender)
    for (i in 1:n_centro) {
        directions <- li_directions[[i]]
        bloc_to_cluster <- newdata[, var_in_clust[[i]]]
        bloc_pca <- as.data.frame(as.matrix(bloc_to_cluster) %*% directions)
        df_x <- cbind(bloc_pca, df_x)
    }
    # print(df_x)
    # print(beta)
    # print(intercept)
    # print(dim(df_x))
    # print(length(beta))
    proba <- 1 / (1 + exp(-as.matrix(df_x) %*% beta - intercept))
    predicted_labels <- ifelse(proba > 0.5,
        ifelse(modelFit$lev[2] == "CCK", modelFit$lev[2], modelFit$lev[1]),
        ifelse(modelFit$lev[1] == "CHC", modelFit$lev[1], modelFit$lev[2])
    )
    # print(length(predicted_labels))
    # print(nrow(newdata))
    # print("done predict")
    return(predicted_labels)
}

li_caret_select$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- renormalize_in_model_pred(newdata, modelFit) ######## THE GOOD LINE FOR NORMALIZATION
    beta <- modelFit$beta
    intercept <- modelFit$intercept
    li_directions <- modelFit$li_directions
    n_centro <- modelFit$n_centro
    # vec_dim <- modelFit$vec_dim
    # segmentation <- modelFit$segmentation
    var_in_clust <- modelFit$var_in_clust

    df_x <- data.frame(Gender = as.data.frame(newdata)$Gender)
    for (i in 1:n_centro) {
        directions <- li_directions[[i]]
        bloc_to_cluster <- newdata[, var_in_clust[[i]]]
        bloc_pca <- as.data.frame(as.matrix(bloc_to_cluster) %*% directions)
        df_x <- cbind(bloc_pca, df_x)
    }
    # print(dim(df_x))
    # print(length(beta))
    proba <- 1 / (1 + exp(-as.matrix(df_x) %*% beta - intercept))
    # print(proba)
    # print(nrow(newdata))
    # print("done prob")
    return(data.frame("CHC" = 1 - proba, "CCK" = proba))
}

li_caret_select$loop <- NULL



setMethod("train_method", "apply_model", function(object) {
    vec_centro <- round(seq(from = object@n_centro_min, to = object@n_centro_max, length.out = object@tuneLength))
    # tuneGrid <- expand.grid(n_centro = vec_centro)

    # Utile si lasso. Metre n_lambda = 1 si pas de lasso (et le dégager dans le code source!!)
    vec_lambda <- seq(from = object@lambda_min, to = object@lambda_max, length.out = object@n_lambda)
    tuneGrid <- expand.grid(lambda = vec_lambda, n_centro = vec_centro)




    object@model <- caret::train(
        y = object@train_cols$classe_name, x = as.matrix(object@train_cols[, object@col_x]),
        method = li_caret_select, trControl = object@cv, metric = "AUC",
        tuneLength = 8, tuneGrid = tuneGrid,
        weights_dict = object@weights, dim_max = object@dim_max,
        k_smote = object@k_smote, do_smote = object@do_smote,
    )
    print("done")
    # object@coeffi <- coef(object@model$finalModel, s = object@lambda)
    object@n_centro <- object@model$bestTune$n_centro
    return(object)
})

setMethod("get_results", "apply_model", function(object) {
    object@predictions <- as.vector(predict(object@model, newdata = as.matrix(object@test_set[, object@col_x])))
    object@predictions_proba <- predict(object@model, newdata = as.matrix(object@test_set[, object@col_x]), type = "prob")
    object@predictions_train_proba <- predict(object@model, newdata = as.matrix(object@train_cols[, object@col_x]), type = "prob")
    return(object)
})

setMethod("importance_method", "apply_model", function(object) {
    print(object@model$finalModel$beta)
    var_in_clust <- object@model$finalModel$var_in_clust
    var_in_clust <- lapply(var_in_clust, sort)
    max_length <- max(sapply(var_in_clust, length))
    var_in_clust <- lapply(var_in_clust, function(x) {
        length(x) <- max_length
        return(x)
    })
    # Convertir la liste en un data.frame
    df <- do.call(cbind.data.frame, var_in_clust)
    writexl::write_xlsx(df, "..\\data\\variables\\logistic_select_var_in_clust.xlsx")
    return(object)
})

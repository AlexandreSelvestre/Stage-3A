get_perfs <- function(best_name_model, best_results, li_roc, vec_accu, df_test, stats_chosen, sum_confus_mat, reclassify_Mixtes, df_mixte, mixte_class_sum = NULL, mixte_proba_sum = NULL) {
    print(paste("Moyenne balanced accuracy validation", mean(best_results$resample$balanced_acc)))
    print(paste("Moyenne AUC validation", mean(best_results$resample$ROC)))
    x_test <- df_test[, setdiff(colnames(df_test), c("Tumeur"))]
    x_test <- data.frame(lapply(x_test, function(col) {
        return(factor(col, levels = c(0, 1)))
    }))
    predictions_proba <- predict(best_results, newdata = x_test, type = "prob")
    predictions_levels <- predict(best_results, newdata = x_test)
    predictions_proba$obs <- df_test$Tumeur
    balanced_acc <- balanced_accuracy(df_test$Tumeur, predictions_levels)
    roc_test <- AUC(predictions_proba, unique(df_test$Tumeur))
    li_roc <- c(li_roc, roc_test)
    vec_accu <- c(vec_accu, balanced_acc)
    sum_confus_mat <- sum_confus_mat + confus_mat(df_test$Tumeur, predictions_levels, unique(df_test$Tumeur))
    print(paste("current models chosen:"))
    print(stats_chosen)
    print(paste("current AUC", roc_test, "current balanced accuracy", balanced_acc))
    print(sum_confus_mat)
    print(paste("mean AUC", mean(li_roc), "mean balanced accuracy", mean(vec_accu)))
    print(paste("std AUC", sd(li_roc), "std balanced accuracy", sd(vec_accu)))
    if (reclassify_Mixtes) {
        x_mixte <- df_mixte[, setdiff(colnames(df_mixte), c("Tumeur"))]
        x_mixte <- data.frame(lapply(x_mixte, function(col) {
            return(factor(col, levels = c(0, 1)))
        }))
        predictions_proba_mixte <- predict(best_results, newdata = x_mixte, type = "prob")[["CCK"]]
        predictions_class_mixte <- predict(best_results, newdata = x_mixte)
        mixte_class_sum <- mixte_class_sum + c(sum(predictions_class_mixte == "CHC"), sum(predictions_class_mixte == "CCK"))
        mixte_proba_sum[1] <- mixte_proba_sum[1] + mean(predictions_proba_mixte)
        mixte_proba_sum[2] <- mixte_proba_sum[2] + 1
        mixte_proba_sum[3] <- mixte_proba_sum[1] / mixte_proba_sum[2]
        print(mixte_class_sum)
        print(paste("proba moyenne des Mixtes pour être CCK:", mixte_proba_sum[3]))
    } else {
        mixte_class_sum <- NULL
        mixte_proba_sum <- NULL
    }
    return(list(vec_accu = vec_accu, li_roc = li_roc, sum_confus_mat = sum_confus_mat, mixte_class_sum = mixte_class_sum, mixte_proba_sum = mixte_proba_sum))
}

one_run <- function(df, li_models, li_perfs, li_imp, k_folds, rep_folds, stats_chosen, metric, sum_confus_mat, numCores, verbose, reclassify_Mixtes, df_mixte) {
    index_train <- as.vector(createDataPartition(y = df$Tumeur, p = 0.75, list = FALSE))
    index_test <- setdiff(seq_len(nrow(df)), index_train)
    df_train <- df[index_train, ]
    df_test <- df[index_test, ]
    df_test <- na.omit(df_test)
    y <- df_train$Tumeur
    y <- as.factor(y)
    x <- df_train[, setdiff(colnames(df_train), c("Tumeur"))]
    x <- data.frame(lapply(x, function(col) {
        return(factor(col, levels = c(0, 1)))
    }))
    folds <- createMultiFolds(y = y, k = k_folds, times = rep_folds)
    # indexOut <- lapply(folds, function(x) {
    #     out <- seq_along(y)
    #     out <- setdiff(out, x)
    #     out <- sapply(out, function(i) {
    #         if (any(is.na(df_train[i, ]))) {
    #             return(-1)
    #         } else {
    #             return(i)
    #         }
    #     })
    #     out <- as.integer(out[out != -1])
    #     return(out)
    # })
    li_results <- list()
    vec_perfs <- c()
    for (name_model in names(li_models)) {
        show_model <- FALSE
        cl <- makePSOCKcluster(numCores)
        registerDoParallel(cl)
        clusterEvalQ(cl, {
            files <- c()
            for (file in files) {
                source(file)
            }
        })
        clusterExport(cl, varlist = c("balanced_accuracy", "AUC"))
        name_vec <- as.character(seq_len(numCores))
        clusterApply(cl, name_vec, function(name) assign("worker_name", name, envir = .GlobalEnv))

        if (any(is.na(x))) {
            capture.output({
                li <- missForest::missForest(xmis = cbind(data.frame(y = y), x), ntree = 1000, maxiter = 10, mtry = 1)
            })
            x <- li$ximp[2:ncol(li$ximp)]
            print(li$OOBerror)
        }

        li_results[[name_model]] <- rlang::exec(caret::train, x, y,
            method = li_models[[name_model]]$method,
            metric = metric, trControl = caret::trainControl(method = "repeatedcv", index = folds, indexOut = NULL, classProbs = TRUE, summaryFunction = caret_metrics, sampling = li_models[[name_model]]$sampling, verbose = verbose),
            tuneGrid = li_models[[name_model]]$grid, !!!li_models[[name_model]]$others
        )
        stopCluster(cl)
        vec_perfs <- c(vec_perfs, mean(li_results[[name_model]]$resample[[metric]]))
    }
    names(vec_perfs) <- names(li_models)
    max_perf <- max(vec_perfs)
    best_model_names <- names(vec_perfs[vec_perfs == max_perf])
    best_name_model <- sample(best_model_names, 1)
    best_results <- li_results[[best_name_model]]

    stats_chosen[[best_name_model]] <- stats_chosen[[best_name_model]] + 1
    for (name_model in names(li_models)) {
        if ("logistic" == name_model) {
            li_imp[[name_model]] <- logistic_importance_beta_intra(li_results[["logistic"]], li_imp[[name_model]], show_model = best_name_model == "logistic")
        }
        if ("rf" == name_model) {
            li_imp[[name_model]] <- rf_importance_intra(li_results[["rf"]], li_imp[[name_model]], show_model = best_name_model == "rf")
        }
        if (best_name_model == "xgbDART") {
            print(li_results[["xgbDART"]]$bestTune)
        }
        if (best_name_model == "xgbLinear") {
            print(li_results[["xgbLinear"]]$bestTune)
        }
    }
    li_roc <- li_perfs$li_roc
    vec_accu <- li_perfs$vec_accu
    mixte_class_sum <- li_perfs$mixte_class_sum
    mixte_proba_sum <- li_perfs$mixte_proba_sum
    li_perfs <- get_perfs(best_name_model, best_results, li_roc, vec_accu, df_test, stats_chosen, sum_confus_mat, reclassify_Mixtes, df_mixte, mixte_class_sum, mixte_proba_sum)
    sum_confus_mat <- li_perfs$sum_confus_mat

    return(list(li_imp = li_imp, li_perfs = li_perfs, stats_chosen = stats_chosen, sum_confus_mat = sum_confus_mat))
}

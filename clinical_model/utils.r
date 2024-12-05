balanced_accuracy <- function(true_values, predictions) {
    vec_classes_tots <- sapply(unique(true_values), function(classe) {
        return(sum(true_values == classe))
    })
    vec_trues <- sapply(unique(true_values), function(classe) {
        return(sum(true_values == classe & predictions == classe))
    })
    vec_accu <- vec_trues / vec_classes_tots
    balanced_acc <- mean(vec_accu)
    return(balanced_acc)
}


AUC <- function(data, lev) {
    deja_vu <- c()
    vec_roc <- c()
    # OVO method for ROC AUC
    for (classe_1 in lev) {
        deja_vu <- c(deja_vu, classe_1)
        for (classe_2 in setdiff(lev, deja_vu)) {
            # print(c(classe_1, classe_2))
            data_loc <- data[data$obs == classe_1 | data$obs == classe_2, ]
            lev_loc <- c(classe_1, classe_2)
            data_loc$obs <- factor(data_loc$obs, levels = lev_loc)
            vec_roc <- c(vec_roc, pROC::auc(pROC::roc(response = data_loc$obs, predictor = data_loc[[classe_2]])))
        }
    }
    ROC <- mean(vec_roc)
    return(ROC)
}

caret_metrics <- function(data, lev = NULL, model = NULL) {
    balanced_acc <- balanced_accuracy(data$obs, data$pred)
    acc <- c(balanced_acc = balanced_acc)
    ROC <- c(ROC = AUC(data, lev))
    out <- c(acc, ROC)
    # print(ROC)
    return(out)
}

confus_mat <- function(true_values, predictions, levels) {
    mat <- matrix(data = 0, nrow = length(levels), ncol = length(levels))
    for (true_class in seq_along(levels)) {
        for (prediction in seq_along(levels)) {
            mat[true_class, prediction] <- sum(true_values == levels[true_class] & predictions == levels[prediction])
        }
    }
    dimnames(mat) <- list("True Values" = as.character(levels), "Predictions" = as.character(levels))
    return(mat)
}

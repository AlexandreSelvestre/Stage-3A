library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(glmnet)
library(SGL)
library(doParallel)
library(jsonlite)
library(pROC)
library(ExPanDaR)
library(NbClust)
library(EMCluster)
library(magrittr)
library(parallel)
library(pracma)
library(mvnfast)
library(Rfast)
library(DMwR)
library(themis)
library(reshape2)
library(rlang)

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
    y <- ifelse(y == "CCK", 1, 0)
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

balanced_accuracy <- function(true_values, predictions) {
    true_0 <- sum(true_values == 0 & predictions == 0)
    true_1 <- sum(true_values == 1 & predictions == 1)
    true_0_classed_1 <- sum(true_values == 0 & predictions == 1)
    true_1_classed_0 <- sum(true_values == 1 & predictions == 0)
    accu_0 <- true_0 / (true_0 + true_0_classed_1)
    accu_0 <- true_1 / (true_1 + true_1_classed_0)
    balanced_acc <- (accu_0 + accu_0) / 2
    return(c(balanced = balanced_acc))
}

balanced_accuracy_caret <- function(data, lev = NULL, model = NULL) {
    # Créer une table de confusion
    confusion_matrix <- table(data$obs, data$pred) # les true labels sont en ligne au vu de l'ordre donné : obs puis pred

    # Calculer la sensibilité (recall pour la classe positive)
    sensitivity <- confusion_matrix[2, 2] / (confusion_matrix[2, 2] + confusion_matrix[2, 1])

    # Calculer la spécificité (recall pour la classe négative)
    specificity <- confusion_matrix[1, 1] / (confusion_matrix[1, 1] + confusion_matrix[1, 2])

    # Calculer la balanced accuracy
    balanced_acc <- (sensitivity + specificity) / 2

    # Retourner la balanced accuracy
    ROC <- twoClassSummary(data, lev, model)["ROC"]
    acc <- c(balanced_acc = balanced_acc)
    out <- c(acc, ROC)
    # print(ROC)
    return(out)
}




# 1 en sexe= Femme
path_data <- "../data"
df <- as.data.frame(readxl::read_xlsx(paste0(path_data, "/stat.xlsx")))
bad_cols <- c("LR-5", "LR-M", "rim_APHE")
df <- df[, setdiff(colnames(df), bad_cols)]
df <- na.omit(df)
df <- data.frame(lapply(df, as.numeric))
df <- df[df$Tumeur != 2, ]
df$Tumeur <- ifelse(df$Tumeur == 1, "CHC", "CCK")
write_xlsx(df, paste0(path_data, "/stat_analysed.xlsx"))
n_runs <- 50
li_roc <- c()
vec_accu <- c()
beta_global <- rep(0, ncol(df) - 1)
beta_non_zero <- rep(0, ncol(df) - 1)
for (i in seq_len(n_runs)) {
    index_train <- as.vector(createDataPartition(y = df$Tumeur, p = 0.75, list = FALSE))
    index_test <- setdiff(seq_len(nrow(df)), index_train)
    df_train <- df[index_train, ]
    df_test <- df[index_test, ]
    y <- df_train$Tumeur
    y <- factor(y, levels = c("CHC", "CCK"))
    # print(y)
    x <- as.matrix(df_train[, setdiff(colnames(df_train), c("Tumeur"))])
    result <- caret::train(x, y, li_caret_simple,
        metric = "ROC", trControl = caret::trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = balanced_accuracy_caret, sampling = "up", verbose = TRUE),
        tuneGrid = expand.grid(lambda = c(0, log(seq(exp(0.001), exp(1), length.out = 10))))
    )

    print(result$bestTune$lambda)
    print(paste("La valeur de l'AUC de validation sur chaque fold est de", result$resample$ROC))
    print(paste("Ce qui donne une moyenne d'AUC de", mean(result$resample$ROC)))
    beta <- result$finalModel$beta
    x_test <- as.matrix(df_test[, setdiff(colnames(df_test), c("Tumeur"))])
    y_test <- ifelse(df_test$Tumeur == "CHC", 0, 1)
    predictions_proba <- predict(result, newdata = as.matrix(x_test), type = "prob")[, 2]
    predictions_levels <- ifelse(predictions_proba > 0.5, 1, 0)
    balanced_acc <- balanced_accuracy(y_test, predictions_levels)
    # print(predictions_proba[, 1])
    roc_test <- pROC::auc(pROC::roc(response = y_test, predictor = predictions_proba))
    li_roc <- c(li_roc, roc_test)
    vec_accu <- c(vec_accu, balanced_acc)
    print(paste("current AUC", roc_test, "current balanced accuracy", balanced_acc))
    print(paste("mean AUC", mean(li_roc), "mean balanced accuracy", mean(vec_accu)))
    print(paste("std AUC", sd(li_roc), "std balanced accuracy", sd(vec_accu)))
    beta_global <- beta_global + abs(beta)
    non_zero_local <- ifelse(abs(beta) > 1e-5, 1, 0)
    beta_non_zero <- beta_non_zero + non_zero_local
}
beta_non_zero <- 100 * beta_non_zero / n_runs
beta_global <- beta_global / n_runs
names_df <- colnames(df[setdiff(colnames(df), c("Tumeur"))])
names(beta_global) <- names_df
names(beta_non_zero) <- names_df
print("beta_global")
print(beta_global)
print("beta_non_zero")
print(beta_non_zero)

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


ggsave("barplot_beta_global.png", plot)

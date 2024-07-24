library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(DMwR)
library(themis)
library(reshape2)
library(SGL)





setClass(
  "apply_model",
  representation(
    data_used = "data.frame", # OBLIGATOIRE     dataframe complet des données
    kill_mixtes = "logical", # OBLIGATOIRE        tuer les mixtes
    kill_mixtes_col = "ANY", # OBLIGATOIRE  colonnes à tuer
    weights = "list", # OBLIGATOIRE     poids des classes
    name_model = "character", # OBLIGATOIRE     nom du modèle
    k = "numeric", # OPTIONNEL                  nombre de folds
    p = "numeric", # OPTIONNEL                  proportion de données pour le training set
    calc_probs = "logical", # OPTIONNEL         calculer les probabilités dans la random forest
    show_logs = "logical", # OPTIONNEL         afficher les logs durant l'entraînement
    n_trees = "numeric", # OPTIONNEL           nombre d'arbres
    mtry = "numeric", # OPTIONNEL            nombre de variables à tester pour chaque split
    rep = "numeric", # OPTIONNEL             nombre de répétitions
    do_PCA = "logical", # OPTIONNEL      faire le preProcess ou non
    sampling = "ANY", # OPTIONNEL    faire le SMOTE/ bootstrap ou rien
    k_smote = "numeric", # OPTIONNEL         nombre de voisins pour le SMOTE
    include_products = "logical", # OPTIONNEL  inclure les produits des variables
    R = "numeric", # OPTIONNEL                Rang multiway
    eps = "numeric", # OPTIONNEL              epsilon proportion convergence multiway
    ite_max = "numeric", # OPTIONNEL          nombre d'itérations max multiway
    search = "character", # OPTIONNEL         type de recherche,
    index = "numeric", # OPTIONNEL            index des données
    pas = "numeric", # OPTIONNEL              pas pour la descente de gradient
    id_term = "character", # OPTIONNEL        nombre d'identifiants
    do_parallel = "logical", # OPTIONNEL     faire du parallélisme
    time_inj = "character",
    train_set = "data.frame",
    test_set = "data.frame",
    cv = "ANY",
    model = "ANY",
    col_x = "character",
    tuneGrid = "data.frame",
    tuneLength = "numeric",
    predictions = "ANY",
    predictions_proba = "ANY",
    train_cols = "data.frame",
    predictions_train_proba = "ANY",
    training_index = "numeric",
    lambda_min = "ANY",
    lambda_max = "numeric",
    lambda = "numeric",
    coeffi = "ANY",
    alpha_values = "numeric",
    n_inner_folds = "numeric",
    alpha_origin = "numeric",
    lambdas = "numeric",
    tuneGrid_type = "character",
    vec_alphas = "ANY",
    index_type = "character",
    dict_temps = "ANY",
    n_iter_per_reg = "numeric",
    li_df_var_imp = "ANY",
    li_box_plots = "ANY",
    df_measures = "ANY",
    index_bloc = "numeric",
    index_mode = "numeric",
    index_variable = "numeric",
    name_bloc = "character",
    name_mode = "character",
    name_variable = "character",
    is_binary = "logical",
    li_imp = "ANY",
    n_centro_min = "numeric",
    n_centro_max = "numeric",
    n_centro = "numeric",
    do_smote = "logical",
    do_boot = "logical",
    dim_max = "numeric",
    n_lambda = "numeric",
    confus_mat = "ANY",
    R_min = "integer",
    R_max = "integer",
    df_danger = "ANY",
    info_cols = "ANY",
    y_tot = "ANY",
    y_train = "ANY",
    y_test = "ANY",
    name_y = "character",
    class_maj_min = "character",
    analyse_data = "ANY",
    beta_final = "numeric",
    tune_R = "numeric",
    classe_1 = "character",
    li_R = "ANY",
    same_R = "logical"
  ),
  prototype(
    k = 5,
    p = 0.75,
    rep = 1,
    n_trees = 500,
    mtry = c(10, 12, 15), # Pour l'instant pas de limite à la croissance des arbres : maxnodes non donné
    calc_probs = TRUE,
    show_logs = TRUE,
    do_PCA = FALSE,
    sampling = NULL,
    name_model = "random_forest",
    id_term = "1"
  )
)


setGeneric("init", function(object) {
  standardGeneric("init")
})

setMethod("init", "apply_model", function(object) {
  object@name_y <- object@info_cols$explained_col
  object@col_x <- setdiff(names(object@data_used), c(object@info_cols$exclude_cols, object@name_y))
  object@data_used[[object@name_y]] <- as.factor(object@data_used[[object@name_y]])
  object@y_tot <- object@data_used[[object@name_y]]
  class_majoritaire <- names(which.max(table(object@y_tot)))
  class_minoritaire <- names(which.min(table(object@y_tot)))
  object@class_maj_min <- c(class_majoritaire, class_minoritaire)
  object@data_used[, object@col_x] <- na.roughfix(object@data_used[, object@col_x])

  object@index_mode <- readRDS(file = "../data/RDS/index_mode.rds")
  object@index_bloc <- readRDS(file = "../data/RDS/index_bloc.rds")
  object@index_variable <- readRDS(file = "../data/RDS/index_variable.rds")
  object@is_binary <- readRDS(file = "../data/RDS/is_binary.rds")
  object@name_mode <- readRDS(file = "../data/RDS/name_mode.rds")
  object@name_bloc <- readRDS(file = "../data/RDS/name_bloc.rds")
  object@name_variable <- readRDS(file = "../data/RDS/name_variable.rds")

  object@li_df_var_imp <- list()
  object@li_box_plots <- list()
  object@df_measures <- as.data.frame(matrix(ncol = 5, nrow = 0))
  # print(object@col_x)
  if (object@sampling == "smote") {
    object@do_smote <- TRUE
  } else {
    object@do_smote <- FALSE
  }
  if (object@sampling == "up") {
    object@do_boot <- TRUE
  } else {
    object@do_boot <- FALSE
  }

  return(object)
})



setGeneric("split_met", function(object, training_index, folds) {
  standardGeneric("split_met")
})

setMethod("split_met", "apply_model", function(object, training_index, folds) {
  liste <- sepa_data_set(
    data = object@data_used, training_index = training_index, folds = folds,
    show = object@show_logs, calc_probs = object@calc_probs,
    summary = super, sampling = object@sampling, k_smote = object@k_smote, search = object@search, rep = object@rep, y_tot = object@y_tot
  )
  object@train_cols <- liste[[1]]
  object@test_set <- liste[[2]]
  object@cv <- liste[[3]]
  object@training_index <- liste[[4]]
  object@y_train <- liste$y_train
  object@y_test <- liste$y_test
  # object@train_cols <- object@train_set
  cat("nombre de la classe majoritaire dans train", length(object@y_train[object@y_train == object@class_maj_min[1]]), "\n")
  cat("nombre de la classe minoritaire dans train", length(object@y_train[object@y_train == object@class_maj_min[2]]), "\n")
  cat("nombre de la classe majoritaire dans test", length(object@y_test[object@y_test == object@class_maj_min[1]]), "\n")
  cat("nombre de la classe minoritaire dans test", length(object@y_test[object@y_test == object@class_maj_min[2]]), "\n")

  cat("rapport nb_minoritaire/nb_majoritaire en training", length(object@y_train[object@y_train == object@class_maj_min[2]]) / length(object@y_train[object@y_train == object@class_maj_min[1]]), "\n")
  cat("rapport nb_minoritaire/nb_majoritaire  en testing", length(object@y_test[object@y_test == object@class_maj_min[2]]) / length(
    object@y_test[object@y_test == object@class_maj_min[1]]
  ), "\n")
  return(object)
})

setGeneric("analyse_data", function(object) {
  standardGeneric("analyse_data")
})

setGeneric("train_method", function(object) {
  standardGeneric("train_method")
})



setGeneric("get_results", function(object) {
  standardGeneric("get_results")
})

setGeneric("compare", function(object) {
  standardGeneric("compare")
})

setMethod("compare", "apply_model", function(object) {
  classe_name <- object@test_set$classe_name
  # print(object@predictions_proba)
  prediction <- object@predictions_proba[, 1]
  index_test <- setdiff(1:nrow(object@data_used), object@training_index)
  index_CCK_test <- index_test[classe_name == "CCK"]
  prediction_CHC <- prediction[classe_name == "CHC"]
  prediction_CCK <- prediction[classe_name == "CCK"]
  prob_lim_low <- 0.45
  prob_lim_high <- 0.75
  lim_bad_CHC <- my_quantile(prediction_CHC, prob_lim_low)
  lim_good_CHC <- my_quantile(prediction_CHC, prob_lim_high)
  object@df_danger <- data.frame(
    is_bad = rep(0, length(index_CCK_test)), score_bad = rep(0, length(index_CCK_test)),
    is_good = rep(0, length(index_CCK_test)), score_good = rep(0, length(index_CCK_test))
  )
  rownames(object@df_danger) <- index_CCK_test
  for (i in 1:length(prediction_CCK)) {
    if (prediction_CCK[i] < lim_bad_CHC) {
      print(calculer_quantile(prediction_CHC, prediction_CCK[i]))
      print(lim_bad_CHC)
      print(prediction_CCK[i])
      print(prediction_CHC)
      object@df_danger$is_bad[i] <- 1
      score_badness <- (prob_lim_low - calculer_quantile(prediction_CHC, prediction_CCK[i])) / prob_lim_low
      object@df_danger$score_bad[i] <- score_badness
    }
    if (prediction_CCK[i] > lim_good_CHC) {
      object@df_danger$is_good[i] <- 1
      score_goodness <- (calculer_quantile(prediction_CHC, prediction_CCK[i]) - prob_lim_high) / (1 - prob_lim_high)
      object@df_danger$score_good[i] <- score_goodness
    }
  }
  return(object)
})



setGeneric("importance_method", function(object) {
  standardGeneric("importance_method")
})


setGeneric("analyse_results", function(object) {
  standardGeneric("analyse_results")
})


setGeneric("get_df_imp", function(object) {
  standardGeneric("get_df_imp")
})


setMethod("analyse_results", "apply_model", function(object) {
  best_params <- object@model$bestTune

  importance_list <- list(
    random_forest = TRUE, logistique_simple = TRUE, logistic_grp = TRUE, logistic_multiway = FALSE, logistic_multibloc = FALSE,
    logistic_select = TRUE
  )
  do_importance <- importance_list[[object@name_model]]
  # object@predictions <- as.factor(object@predictions)
  object@predictions <- factor(object@predictions, levels = levels(object@test_set[[object@name_y]]))
  object@confus_mat <- confusionMatrix(object@predictions, object@test_set[[object@name_y]])
  # print(object@confus_mat$table)



  if (do_importance) {
    object <- importance_method(object)
  }

  print(paste(paste("Le nombre de données est", dim(object@data_used)[1], "dont", length(object@test_set[[object@name_y]]), "dans le testing dataset"), "et le jeu de donné de la grid search: "))
  print(best_params)
  if (object@do_PCA) {
    print(paste("La PCA a conservé le nombre suivant de composantes:", dim(predict(object@model$preProcess, object@test_set[, object@col_x]))[2]))
  }
  vec_roc_res <- c(0, 0, 0, 0, 0)

  pdf("plots/ROC_curves.pdf")
  par(mfrow = c(2, 1))
  roc_test <- pROC::roc(response = object@test_set[[object@name_y]], predictor = object@predictions_proba[, 1])
  ##### Changer et checker structure
  plot(roc_test, main = "ROC Curve test set")
  abline(a = 0, b = 1)
  value_auc <- pROC::auc(roc_test)
  print(paste("La valeur de l'AUC de test est de", value_auc))
  vec_roc_res[3] <- value_auc

  print(paste("La valeur de l'AUC de validation sur chaque fold est de", object@model$resample$AUC))
  # print(object@model$resample[, c("Resample", "R", "lambda", "AUC")])
  print(paste("Ce qui donne une moyenne d'AUC de", mean(object@model$resample$AUC)))
  vec_roc_res[2] <- mean(object@model$resample$AUC, na.rm = TRUE)

  roc_train <- pROC::roc(response = object@train_cols[[object@name_y]], predictor = object@predictions_train_proba[, 1])
  plot(roc_train, main = "ROC Curve train set")
  abline(a = 0, b = 1)
  value_auc_train <- pROC::auc(roc_train)
  print(paste("La valeur de l'AUC de train est de", value_auc_train))
  vec_roc_res[1] <- value_auc_train
  dev.off()

  vec_roc_res[4] <- f1_macro_summary(data.frame(obs = object@test_set[[object@name_y]], pred = object@predictions))
  vec_roc_res[5] <- weighted_accuracy_summary(data.frame(obs = object@test_set[[object@name_y]], pred = object@predictions))

  object@df_measures <- rbind(object@df_measures, vec_roc_res)
  colnames(object@df_measures) <- c("AUC_train", "AUC_val", "AUC_test", "F1", "Acc")
  return(object)
})

# Degager les négatifs: d'où ils viennent

# [1] "Les stats des individus CCK sont"
#    is_bad   score_bad is_good score_good
# 2       0  0.00000000      10  9.4736842
# 8       1 -0.05263158       7  4.1052632
# 14      0  0.00000000       7  4.6315789
# 17     10  3.15789474       3  0.7631579
# 18      4  1.36842105       5  2.6315789
# 19      0  0.00000000      11  9.5526316
# 20      0  0.00000000       8  5.3684211
# 22      0  0.00000000      11 10.6052632
# 29      1  0.21052632       7  5.8157895
# 41      1 -0.05263158       4  2.6842105
# 43      1  0.47368421       3  0.8947368
# 46      0  0.00000000       6  5.0789474
# 57      3 -0.02631579       6  4.5526316
# 64      3  1.02631579       1  0.6052632
# 77      9  3.34210526       0  0.0000000
# 84      0  0.00000000      14 11.1052632

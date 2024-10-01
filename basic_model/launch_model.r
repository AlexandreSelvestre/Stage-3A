## Classe relative à l'importation des données du modèle
setClass(
  "import_data",
  representation(
    data_used = "data.frame", #   input: dataframe complet des données (extraites)
    info_cols = "ANY", #          input: liste indiquant col expliqué et cols est à exclure
    index_bloc = "numeric", #     input: vecteur d'indices des blocs des explicatives
    index_mode = "ANY", #     input: vecteur d'indices des modes des explicatives
    index_variable = "numeric", # input: vecteur d'indices des variables des explicatives
    name_bloc = "character", #    input: vecteur des noms des blocs des explicatives
    name_mode = "character", #    input: vecteur des noms des modes des explicatives
    name_variable = "character", # input: vecteur des noms des variables des explicatives
    is_binary = "logical", #      input: vecteur logique des variables explicatives binaires
    path_data = "character", #         input: chemin du fichier de données
    name_y = "character", #       variable: nom de la variable expliquée
    class_maj_min = "character", # variable: c(nom_class_maj, nom_class_min)
    col_x = "character", #        variable: noms de colonnes de variables explicatives
    train_cols = "data.frame", #  variable: partie de data_used utilisée pour train
    test_set = "data.frame", #     variable: partie de data_used utilisée pour test,
    li_index_modes = "list",
    li_name_modes = "list",
    use_li_index_modes = "logical" # variable: indique si on utilise li_index_modes
  ),
  prototype(
    li_index_modes = list(),
    li_name_modes = list()
  )
)

## Classe relative à la structuration des variables contenant les données dans le modèle
# (hérite de la classe "import_data")
setClass(
  "data_structuration",
  contains = "import_data",
  representation(
    training_index = "numeric", #     input: vecteur des indices des lignes du train
    cv = "ANY", #                     variable: argument train_control dans caret
    predictions = "ANY", #            variable: vecteur des facteurs prédits sue le test
    predictions_proba = "ANY", #      variable: vecteur des probas prédites sur le test
    predictions_train_proba = "ANY", # variable: vecteur des probas prédites sur le train
    y_tot = "ANY", #                  variable: vecteur expliqué
    y_train = "ANY", #                variable: partie train de y_tot
    y_test = "ANY", #                 variable: partie test de y_tot
    beta_final = "numeric" #          variable : le beta obtenu à la fin du modèle
  )
)


setClass(
  "apply_model",
  contains = "data_structuration",
  representation(
    name_model = "character", #     input: nom du modèle
    k = "numeric", #                input: nombre de folds en cross valid
    p = "numeric", #                input: proportion de données dans le training set
    show_logs = "logical", #        input: afficher les logs durant l'entraînement
    rep = "numeric", #              input: nombre de répétitions des folds dans la cv
    sampling = "ANY", #             input: indique si faire smote / bootstrap (up)/rien
    k_smote = "numeric", #          input: nombre de voisins dans smote
    search = "character", #         input: espacer régulièrement ou aléat dans gridSearch
    id_term = "character", #        input: Identifiant du run
    parallel = "ANY", #             input: liste indiquant si oui et comment paralléliser
    model = "ANY", #                input: nom du modèle utilisé
    classe_1 = "character", #       input: Définit qui est la classe 1
    analyse_data = "list", #        input: liste sur post analyse data (faire clusters?)
    do_product = "logical", #       input: indique si réutiliser ou non le Sigma picto
    df_measures = "ANY", #          variable: dataframe des perfoemances
    li_df_var_imp = "ANY", #        variable: dataframe imp variables (pas tous les modèles)
    confus_mat = "ANY", #           variable: matrice de confusion
    df_danger = "ANY", #            variable: to refit: df scores sur lignes mauvaises en radio
    score_recons = "numeric", #      variable: erreur L1 de reconstruction du picto,
    penalty_adapt = "logical" #     variable: penalty adaptatif pour le modèle
  ),
  prototype(
    k = 5,
    p = 0.75,
    rep = 1,
    show_logs = TRUE,
    name_model = "logistique_simple",
    id_term = "1"
  )
)


setGeneric("init", function(object) {
  standardGeneric("init")
})

setMethod("init", "apply_model", function(object) {
  path_data <- object@path_data
  object@name_y <- object@info_cols$explained_col
  object@col_x <- setdiff(names(object@data_used), c(object@info_cols$exclude_cols, object@name_y))
  object@data_used[[object@name_y]] <- as.factor(object@data_used[[object@name_y]])
  object@y_tot <- object@data_used[[object@name_y]]
  class_majoritaire <- names(which.max(table(object@y_tot)))
  class_minoritaire <- names(which.min(table(object@y_tot)))
  object@class_maj_min <- c(class_majoritaire, class_minoritaire)
  object@data_used[, object@col_x] <- na.roughfix(object@data_used[, object@col_x])

  object@index_mode <- readRDS(file = paste0(path_data, "/RDS/index_mode.rds"))
  object@index_bloc <- readRDS(file = paste0(path_data, "/RDS/index_bloc.rds"))
  object@index_variable <- readRDS(file = paste0(path_data, "/RDS/index_variable.rds"))
  object@is_binary <- readRDS(file = paste0(path_data, "/RDS/is_binary.rds"))
  object@name_mode <- readRDS(file = paste0(path_data, "/RDS/name_mode.rds"))
  object@name_bloc <- readRDS(file = paste0(path_data, "/RDS/name_bloc.rds"))
  object@name_variable <- readRDS(file = paste0(path_data, "/RDS/name_variable.rds"))

  if (object@use_li_index_modes) {
    object@li_index_modes <- readRDS(file = paste0(path_data, "/RDS/li_index_modes.rds"))
    object@li_name_modes <- readRDS(file = paste0(path_data, "/RDS/li_name_modes.rds"))
  } else {
    object@li_index_modes <- list()
    object@li_name_modes <- list()
  }

  object@df_measures <- as.data.frame(matrix(ncol = 5, nrow = 0))


  return(object)
})



setGeneric("split_met", function(object, training_index, folds) {
  standardGeneric("split_met")
})

setMethod("split_met", "apply_model", function(object, training_index, folds) {
  liste <- sepa_data_set(
    data = object@data_used, training_index = training_index, folds = folds,
    show = object@show_logs, calc_probs = TRUE,
    summary = super, sampling = object@sampling, k_smote = object@k_smote, search = object@search, rep = object@rep, y_tot = object@y_tot
  )
  object@train_cols <- liste[[1]]
  object@test_set <- liste[[2]]
  object@cv <- liste[[3]]
  object@training_index <- liste[[4]]
  object@y_train <- liste$y_train
  object@y_test <- liste$y_test
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

setMethod("get_results", "apply_model", function(object) {
  object@predictions <- as.vector(predict(object@model, newdata = as.matrix(object@test_set[, object@col_x])))
  object@predictions_proba <- predict(object@model, newdata = as.matrix(object@test_set[, object@col_x]), type = "prob")
  object@predictions_train_proba <- predict(object@model, newdata = as.matrix(object@train_cols[, object@col_x]), type = "prob")
  return(object)
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
    random_forest = TRUE, logistique_simple = FALSE, logistic_grp = TRUE, logistic_multiway = FALSE, logistic_multibloc = FALSE,
    logistic_select = TRUE
  )
  do_importance <- importance_list[[object@name_model]]
  object@predictions <- factor(object@predictions, levels = levels(object@test_set[[object@name_y]]))
  object@confus_mat <- confusionMatrix(object@predictions, object@test_set[[object@name_y]])



  if (do_importance) {
    object <- importance_method(object)
  }

  print(paste(paste("Le nombre de données est", dim(object@data_used)[1], "dont", length(object@test_set[[object@name_y]]), "dans le testing dataset"), "et le jeu de donné de la grid search: "))
  print(best_params)
  vec_roc_res <- c(0, 0, 0, 0, 0)

  pdf(paste0(object@path_data, "/plots/ROC_curves.pdf"))
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

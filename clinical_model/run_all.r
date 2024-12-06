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
library(purrr)
library(xgboost)
library(plyr)
library(missForest)



####### AVANCER DISPERSION
########### paramétrage
name_models <- c("rf") # "logistic", "rf", "xgbDART", "xgbLinear"
n_runs <- 2
k_folds <- 10
rep_folds <- 1
metric <- "ROC" # "ROC" ou "balanced_acc"
numCores <- detectCores() - 1
analyze_data <- FALSE
dedoubler_shape <- TRUE
verbose <- TRUE
impute_data <- TRUE # impossible d'analyser les data en cas d'imputation: à cause des NA!!
kill_tumor <- "Mixte" # CHC CCK ou Mixte (ou "")
reclassify_Mixtes <- TRUE # ne sera considéré que si on a seulment étudié CCK et CHC

if (kill_tumor != "Mixte") {
    reclassify_Mixtes <- FALSE
} else {
    df_mixte <- NULL
}

########### Importations

if ("logistic" %in% name_models) {
    source("logistic.r")
}
if ("rf" %in% name_models) {
    source("rf.r")
}
source("import.r")

if (analyze_data) {
    source("analyze_data.r")
}

source("one_run.r")
source("utils.r")


############## Initialization
li_perfs <- list(li_roc = c(), vec_accu = c())
if (reclassify_Mixtes) {
    li_perfs$mixte_class_sum <- c("CHC" = 0, "CCK" = 0)
    li_perfs$mixte_proba_sum <- list()
} else {
    li_perfs$mixte_class_sum <- NULL
    li_perfs$mixte_proba_sum <- NULL
}
li_imp <- list()
li_models <- list()
if ("logistic" %in% name_models) {
    li_imp[["logistic"]] <- list()
    li_imp[["logistic"]]$beta_global <- rep(0, ncol(df) - 1)
    li_imp[["logistic"]]$beta_non_zero <- rep(0, ncol(df) - 1)
    li_models[["logistic"]] <- list(
        method = li_caret_simple, grid = expand.grid(lambda = c(0, log(seq(exp(0.001), exp(1), length.out = 10)))),
        others = list(), sampling = "up"
    )
}
if ("rf" %in% name_models) {
    li_models[["rf"]] <- list(method = li_rf, grid = expand.grid(mtry = c(1, 2, 3, 4, 5, 6)), others = list(ntree = 5000), sampling = NULL)
    li_imp[["rf"]] <- list()
}
if ("xgbDART" %in% name_models) {
    li_models[["xgbDART"]] <- list(
        method = "xgbDART", grid = expand.grid(
            nrounds = c(100, 200),
            max_depth = c(1, 2, 3, 5),
            eta = c(0.01, 0.05),
            gamma = c(0, 1, 5),
            subsample = c(0.5, 0.7),
            colsample_bytree = c(0.5, 0.7),
            rate_drop = c(0, 0.1, 0.2),
            skip_drop = c(0, 0.1, 0.2),
            min_child_weight = c(1, 2, 5)
        ),
        others = list(), sampling = "up"
    )
    li_imp[["xgbDART"]] <- list()
}

if ("xgbLinear" %in% name_models) {
    li_models[["xgbLinear"]] <- list(method = "xgbLinear", grid = expand.grid(
        nrounds = c(5, 10, 25, 50, 75, 100),
        lambda = 0,
        alpha = 0,
        eta = c(0.0001, 0.001)
    ), others = list(), sampling = "up")
    li_imp[["xgbLinear"]] <- list()
}
stats_chosen <- rep(0, length(name_models))
names(stats_chosen) <- name_models
levels_mat <- unique(df$Tumeur)
sum_confus_mat <- matrix(data = 0, nrow = length(levels_mat), ncol = length(levels_mat))

############## Runs
for (i in seq_len(n_runs)) {
    li_run <- one_run(df, li_models, li_perfs, li_imp, k_folds, rep_folds, stats_chosen, metric, sum_confus_mat, numCores, verbose, reclassify_Mixtes, df_mixte)
    li_imp <- li_run$li_imp
    li_perfs <- li_run$li_perfs
    stats_chosen <- li_run$stats_chosen
    sum_confus_mat <- li_run$sum_confus_mat
}

############## Analyzing runs
if ("logistic" %in% name_models) {
    logistic_importance_beta_extra(li_imp[["logistic"]], df)
}
if ("rf" %in% name_models) {
    rf_importance_extra(li_imp[["rf"]], df)
}
if (reclassify_Mixtes) {
    mixte_proba_sum <- li_perfs$mixte_proba_sum
    png(filename = "plots/mixte_proba_distribution.png", width = 800, height = 600)
    hist(Reduce("c", mixte_proba_sum), breaks = 101, main = "Répartition des éléments de mixte_proba_sum", xlab = "Valeurs", ylab = "Fréquence", col = "blue", border = "black")
    dev.off()
    print(paste("length:", length(mixte_proba_sum)))
    li_compos_per_compos <- lapply(mixte_proba_sum[[1]], function(elem) {
        return(c(elem))
    })
    # print(mixte_proba_sum)
    # print(li_compos_per_compos)
    if (length(mixte_proba_sum) > 1) {
        for (r in 2:length(mixte_proba_sum)) {
            for (j in 1:length(li_compos_per_compos)) {
                li_compos_per_compos[[j]] <- c(li_compos_per_compos[[j]], mixte_proba_sum[[r]][j])
            }
        }
    }
    vec_sd <- sapply(li_compos_per_compos, sd)
    names(vec_sd) <- df_mixte$patient_num
    vec_mean <- sapply(li_compos_per_compos, mean)
    names(vec_mean) <- df_mixte$patient_num
    print("proba par patient:")
    print(vec_mean)
    # print(sapply(li_compos_per_compos, mean))
    # print(vec_sd)
    print(paste("Ecart type moyen d'une même proba mixte entre runs", mean(vec_sd)))
}

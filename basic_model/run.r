rm(list = ls())

path_data <- "../data"

current_dir <- getwd()
if (current_dir == "/gpfs/users/selvestra/basic_model") {
    path_data <- "/gpfs/workdir/selvestra/data"
}

# library(Rmpi)
# library(doMPI)
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


sysname <- Sys.info()["sysname"]
set.seed(11) ## seed du modèle
seed_model <- .Random.seed
set.seed(7) # seed pour la cross validation
seed_cv <- .Random.seed
set.seed(2) ## seed pour la partition des données
seed_partition <- .Random.seed
.Random.seed <- seed_model

config_run <- config::get(file = "configs/config_run.yml", config = "my_config")

if (config_run$simulated_data) {
    config_model <- config::get(file = "configs/config_model.yml", config = "simu")
    if (config_run$simple_generation) {
        config_extrac <- config::get(file = "configs/extrac/config_simul_basic.yml", config = "my_config")
        name_config <- "simu"
        source("extrac/simul_basic.r")
    } else {
        config_extrac <- config::get(file = "configs/extrac/config_join_picto.yml", config = "my_config")
        name_config <- "simu"
        source("extrac/join_picto.r")
    }
} else {
    if (config_run$big_data) {
        config_model <- config::get(file = "configs/config_model.yml", config = "radio_big")
        config_extrac <- config::get(file = "configs/extrac/config_extrac_multi_slice.yml", config = "my_config")
        name_config <- "radio_big"
        source("extrac/extrac_multi_slice.r")
    } else {
        config_model <- config::get(file = "configs/config_model.yml", config = "radio")
        config_extrac <- config::get(file = "configs/extrac/config_extrac.yml", config = "my_config")
        name_config <- "radio"
        source("extrac/extrac.r")
    }
}

model_name <- config_model$name_model
li_model_configs <- list(
    "random_forest" = config::get(file = "configs/model/random_forest.yml", config = name_config),
    "logistique_simple" = config::get(file = "configs/model/logistique_simple.yml", config = name_config),
    "logistic_grp" = config::get(file = "configs/model/logistic_grp.yml", config = name_config),
    "logistic_multiway" = config::get(file = "configs/model/logistic_multiway.yml", config = name_config),
    "logistic_multibloc" = config::get(file = "configs/model/logistic_multibloc.yml", config = name_config),
    "logistic_select" = config::get(file = "configs/model/logistic_select.yml", config = name_config)
)

model_config <- li_model_configs[[model_name]]
config <- modifyList(config_model, config_extrac)
config <- modifyList(config, model_config)
config$path_data <- path_data
regression <- config$regression



path_model <- paste0("model/", model_name, ".r")
if (regression) {
    path_importance <- "run_importance/run_regression.r"
} else {
    stop("Not implemented yet")
}
path_plot <- paste0(path_data, "/plots/", model_name)




source("main.r")
source("./utils/utils.r")
source("launch_model.r")
source("analyse_data/find_outliers.r")
source("analyse_data/heatmap.r")

import_folder("./utils")
# import_folder("./StepReg_modif/R")

unregister_dopar()

source(path_model)
source(path_importance)

id_li <- 1:config_run$num_runs

sum_test <- 0
sum_val <- 0
ite <- 0
######### Ameliorer pour permettre simulated data
performance <- list(AUC_test = c(), AUC_val = c(), Acc = c(), F1_macro = c(), F1_CCK = c(), F1_CHC = c())
is_null <- list()
imp_li <- list()
li_confus <- list()
if (sysname == "Linux") {
    path <- paste0(path_data, "/data_used.csv")
} else {
    path <- "..\\data\\data_used.csv"
}


if (config_run$extrac_first) {
    extract_all(config, sysname)
}
data_used_local <- as.data.frame(read.csv(path))




# variables <- data_used_local[, setdiff(colnames(data_used_local), c("patient_num", "keys"))] ###### To change!!!!!
# index_CCK <- rownames(variables[variables$classe_name == "CCK", ])
# df_danger <- data.frame(is_bad = rep(0, length(index_CCK)), score_bad = rep(0, length(index_CCK)), is_good = rep(0, length(index_CCK)), score_good = rep(0, length(index_CCK)))
# rownames(df_danger) <- index_CCK
# vec_bad_lambda <- c()
# vec_good_lambda <- c()


if (config$minimal_information) {
    # Attention, en minimal info, on ne récupère pas la moyenne des résultats des pictos...
    for (id_term in id_li) {
        ite <- ite + 1
        list_execute <- execute(config, config_run, as.character(id_term), seed_cv, seed_partition, sysname)
        inference <- list_execute$inference
        seed_cv <- list_execute$seed_cv
        seed_partition <- list_execute$seed_partition
        test_result <- inference@df_measures[["AUC_test"]]
        sum_test <- test_result + sum_test
        val_result <- inference@df_measures[["AUC_val"]]
        sum_val <- val_result + sum_val
        print(paste("actuelle moyenne AUC test", sum_test / ite, "ite:", ite))
        print(paste("actuelle moyenne AUC val", sum_val / ite, "ite:", ite))
        unregister_dopar()
    }
    print(paste("moyenne AUC test", sum_test / length(id_li)))
    print(paste("moyenne AUC val", sum_val / length(id_li)))
} else {
    for (id_term in id_li) {
        ite <- ite + 1
        list_execute <- execute(config, config_run, as.character(id_term), seed_cv, seed_partition, sysname)
        inference <- list_execute$inference
        seed_cv <- list_execute$seed_cv
        seed_partition <- list_execute$seed_partition
        li_intermediaire <- run_imp_intra(inference, imp_li, performance, li_confus, is_null, ite)
        imp_li <- li_intermediaire$imp_li
        performance <- li_intermediaire$performance
        li_confus <- li_intermediaire$li_confus
        is_null <- li_intermediaire$is_null


        # inference <- compare(inference)
        # df_danger_loc <- inference@df_danger
        # df_danger[rownames(df_danger_loc), ] <- df_danger[rownames(df_danger_loc), ] + df_danger_loc
        # print("Les stats des individus CCK sont")
        # print(df_danger)
        # new_lambda <- inference@model$bestTune[[1]]
        # print(new_lambda)
        # if (performance$AUC_test[[length(performance$AUC_test)]] < 0.62) {
        #     vec_bad_lambda <- c(vec_bad_lambda, new_lambda)
        #     cat("moyenne des lambdas catastrophiques:", mean(vec_bad_lambda), "leur liste est", "\n")
        #     print(vec_bad_lambda)
        # }
        # if (performance$AUC_test[[length(performance$AUC_test)]] > 0.8) {
        #     vec_good_lambda <- c(vec_bad_lambda, new_lambda)
        #     cat("moyenne des lambdas très bons:", mean(vec_good_lambda), "leur liste est", "\n")
        #     print(vec_good_lambda)
        # }
        unregister_dopar()
    }

    run_imp_extra(imp_li, performance, li_confus, is_null, length(id_li), path_plot, inference)
}

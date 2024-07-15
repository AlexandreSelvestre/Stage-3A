rm(list = ls())

library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(ggplot2)
library(glmnet)
library(SGL)
library(doParallel)
library(jsonlite)
library(pROC)
library(ExPanDaR)


set.seed(8) ## seed maîtresse: celle des modèles
seed_model <- .Random.seed
set.seed(1)
seed_cv <- .Random.seed
.Random.seed <- seed_model

config_run <- config::get(file = "configs/config_run.yml", config = "my_config")
if (config_run$simulated_data) {
    config_model <- config::get(file = "configs/config_model.yml", config = "simu")
    config_extrac <- read_json("configs/config_extrac_simul.json", simplifyVector = TRUE, simplifyMatrix = FALSE)
    name_config <- "simu"
    source("extrac_simul.r")
} else {
    config_model <- config::get(file = "configs/config_model.yml", config = "radio")
    config_extrac <- config::get(file = "configs/config_extrac.yml", config = "my_config")
    name_config <- "radio"
    source("extrac.r")
}
model_name <- config_model$name_model

li_model_configs <- list(
    "random_forest" = config::get(file = "configs/random_forest.yml", config = name_config),
    "logistique_simple" = config::get(file = "configs/logistique_simple.yml", config = name_config),
    "logistic_grp" = config::get(file = "configs/logistic_grp.yml", config = name_config),
    "logistic_multiway" = config::get(file = "configs/logistic_multiway.yml", config = name_config),
    "logistic_multibloc" = config::get(file = "configs/logistic_multibloc.yml", config = name_config),
    "logistic_select" = config::get(file = "configs/logistic_select.yml", config = name_config)
)

model_config <- li_model_configs[[model_name]]
config <- modifyList(config_model, config_extrac)
config <- modifyList(config, model_config)
regression <- config$regression


path_model <- paste0("model/", model_name, ".r")
if (regression) {
    path_importance <- "run_importance/run_regression.r"
} else {
    stop("Not implemented yet")
}
path_plot <- paste0("plots/", model_name)

source("main.r")
source("./utils/utils.r")
source("launch_model.r")


import_folder("./utils")
import_folder("./StepReg_modif/R")

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
data_used_local <- as.data.frame(read.csv("..\\data\\data_used.csv"))
variables <- data_used_local[, setdiff(colnames(data_used_local), c("patient_num", "keys"))] ###### To change!!!!!
index_CCK <- rownames(variables[variables$classe_name == "CCK", ])
df_danger <- data.frame(is_bad = rep(0, length(index_CCK)), score_bad = rep(0, length(index_CCK)), is_good = rep(0, length(index_CCK)), score_good = rep(0, length(index_CCK)))
rownames(df_danger) <- index_CCK
vec_bad_lambda <- c()
vec_good_lambda <- c()

print("Start")

if (config$minimal_information) {
    for (id_term in id_li) {
        ite <- ite + 1
        list_execute <- execute(config, config_run, as.character(id_term), seed_cv)
        inference <- list_execute$inference
        seed_cv <- list_execute$seed_cv
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
        list_execute <- execute(config, config_run, as.character(id_term), seed_cv)
        inference <- list_execute$inference
        seed_cv <- list_execute$seed_cv
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

# Refaire le run long et en particulier les plots...


# Sur 40 runs sans pondération
# dataframe des dangers avec n_mauvais: 14 n_bons: 13
#    mauvais bons
# 2    0.375  0.500
# 8    0.625  0.500
# 14   0.125  0.500
# 17   1.000  0.250
# 18   0.750  0.250
# 19   0.250  0.500
# 20   0.250  0.375
# 22   0.500  0.375
# 29   0.625  0.250
# 41   0.250  0.500
# 43   0.375  0.125
# 46   0.000  0.875
# 57   0.250  0.250
# 64   0.375  0.000
# 77   0.875  0.250
# 84   0.375  1.000

# Sur 40 runs avec pondération:
# dataframe des dangers avec n_mauvais: 14 n_bons: 13
#       mauvais       bons
# 2  0.25627817 0.64628821
# 8  0.63103670 0.47161572
# 14 0.01030264 0.44978166
# 17 1.00000000 0.19213974
# 18 0.65743722 0.25764192
# 19 0.31036703 0.34061135
# 20 0.26207341 0.35371179
# 22 0.44365744 0.22270742
# 29 0.51835158 0.41048035
# 41 0.26207341 0.51528384
# 43 0.49774630 0.13973799
# 46 0.00000000 1.00000000
# 57 0.13329041 0.27947598
# 64 0.27237605 0.00000000
# 77 0.66773986 0.03930131
# 84 0.25627817 0.46288210


# git rm -r --cached .
# git add .
# git commit -m "Retiré les fichiers ignorés de l'index"

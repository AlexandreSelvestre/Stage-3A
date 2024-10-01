rm(list = ls())

path_data <- "../data"

current_dir <- getwd()
if (current_dir == "/gpfs/users/selvestra/basic_model") {
    # .libPaths("/gpfs/workdir/selvestra/R_packages")
    path_data <- "/gpfs/workdir/selvestra/data"
}

library(Rmpi)
library(doMPI)
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


# cl <- startMPIcluster(count = 2)
# registerDoMPI(cl)
# closeCluster(cl)
# x <- foreach(i = 1:3) %dopar% {
#     sqrt(i)
# }

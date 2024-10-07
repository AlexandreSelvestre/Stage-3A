# rm(list = ls())

# path_data <- "../data"

# current_dir <- getwd()
# if (current_dir == "/gpfs/users/selvestra/basic_model") {
#     # .libPaths("/gpfs/workdir/selvestra/R_packages")
#     path_data <- "/gpfs/workdir/selvestra/data"
# }

library(Rmpi)
library(doMPI)
# library(readxl)
# library(writexl)
# library(glue)
# library(data.table)
# library(caret)
# library(randomForest)
# library(ggplot2)
# library(glmnet)
# library(SGL)
# library(doParallel)
# library(jsonlite)
# library(pROC)
# library(ExPanDaR)
# library(NbClust)
# library(EMCluster)
# library(magrittr)
# library(parallel)
# library(pracma)
# library(mvnfast)
# library(Rfast)
# library(DMwR)
# library(themis)
# library(reshape2)
# library(rlang)

# srun --mpi=pmi2 Rscript toy.r
# srun --mpi=pmi2 ./mpi_info


print(paste("cpu", mpi.universe.size()))
print(paste("rank", mpi.comm.rank()))
print(paste("size MPI comm 1", mpi.comm.size(0)))
cl <- startMPIcluster(maxcores = 2, includemaster = TRUE)
registerDoMPI(cl)
results <- foreach(i = 1:4) %dopar% {
    res <- list()
    # res$sqrt <- sqrt(i)
    res$cpu <- paste("cpu", mpi.universe.size())
    res$rank <- paste("rank", mpi.comm.rank())
    res$size <- paste("size MPI comm 1", mpi.comm.size(0))
    res
}
# print(paste("cpu", mpi.universe.size()))
# print(paste("rank", mpi.comm.rank()))
# print(paste("size MPI comm 1", mpi.comm.size(0)))
closeCluster(cl)
print(results)
mpi.quit()

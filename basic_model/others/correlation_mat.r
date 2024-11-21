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
source("./utils/utils.r")
import_folder("./utils")

path_data <- "../data"

df_3D <- read_excel(paste0(path_data, "/data_radio_in_lines.xlsx"))
df_3D <- as.data.frame(df_3D)[setdiff(colnames(df_3D), c("keys", "patient_num", "classe_name"))]
df_3D <- df_3D[, grepl("gldm", colnames(df_3D))]
df_3D <- na.omit(df_3D)
index_variable <- get_variable_vec_liver(df_3D)$index_variable
index_bloc <- get_bloc_vec_liver(names(df_3D), exist_temps = TRUE, different_clin = FALSE)$index_bloc

df_3D <- renormalize_in_model_fit_index_mode(df_3D, index_variable, index_bloc, is_binary = NULL)$new_x

print(dim(df_3D))
cor_mat <- cor(df_3D)
print(dim(cor_mat))
# Transformer la matrice des corrélations en un format long
cor_mat_melted <- melt(cor_mat)
# Tracer la heatmap de la matrice des corrélations
p <- ggplot(data = cor_mat_melted, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
        low = "blue", high = "red", mid = "white",
        midpoint = 0, limit = c(-1, 1), space = "Lab",
        name = "Correlation"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(
        angle = 45, vjust = 1,
        size = 12, hjust = 1
    )) +
    coord_fixed()
# Enregistrer la heatmap dans un fichier .png
ggsave("./others/correlation_matrix.png", plot = p, width = 10, height = 8)

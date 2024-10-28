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
library(gglasso)

path_data <- "../data"

df <- read_excel(glue("{path_data}/radiomiques_global.xlsx"))
df$key <- paste0(df$classe_name, "_", df$patient_num)
df <- df[!df$key %in% c("CCK_16", "CHC_204"), ]
df <- df[df$classe_name == "Mixtes", ]


df <- df[, c("key", "temps_inj", colnames(df)[grepl("original", colnames(df))])]

df_art <- df[df$temps_inj == "ART", ]
colnames(df_art)[2:ncol(df_art)] <- paste0(colnames(df_art)[2:ncol(df_art)], "_ART")
df_port <- df[df$temps_inj == "PORT", ]
colnames(df_port)[2:ncol(df_art)] <- paste0(colnames(df_port)[2:ncol(df_art)], "_PORT")
df_vein <- df[df$temps_inj == "VEIN", ]
colnames(df_vein)[2:ncol(df_art)] <- paste0(colnames(df_vein)[2:ncol(df_art)], "_VEIN")
df_tard <- df[df$temps_inj == "TARD", ]
colnames(df_tard)[2:ncol(df_art)] <- paste0(colnames(df_tard)[2:ncol(df_art)], "_TARD")

df <- merge(df_art, df_port, by = "key")
# df <- merge(df, df_vein, by = "key")
df <- merge(df, df_tard, by = "key")



df <- na.omit(df)
write_xlsx(df, "../data/radiomiques_global_check.xlsx")
print(nrow(df))

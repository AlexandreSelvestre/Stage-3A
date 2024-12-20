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
library(NbClust)
library(EMCluster)
library(magrittr)
source("extrac/gene_x_scalar.r")

# Beta sera décomposé ligne par ligne pour rester cohérent
extract_all <- function(config_extrac, sysname) {
    path_data <- config_extrac$path_data
    # Récupérer les dossiers des pictogrammes d'intérêt
    vec_num_picto <- config_extrac$num_picto
    name_picto <- sapply(vec_num_picto, function(num_picto) {
        vec_path <- list.dirs(paste0(path_data, "/beta_picto"), recursive = FALSE, full.names = FALSE)
        good_dir <- vec_path[substr(vec_path, 1, length(as.character(num_picto))) == as.character(num_picto)]
        return(good_dir)
    })

    # Récupérer les matrices X et beta des pictogrammes sélectionnés
    li_beta_matrix <- lapply(name_picto, function(good_dir) {
        path <- paste0(path_data, "/beta_picto/", good_dir, "/Beta0.csv")
        beta_matrix <- as.matrix(read.csv(path, header = FALSE))
        return(beta_matrix)
    })

    # Mettre un même nombre de lignes (i.e. de modes K à chaque pictogramme, en complétant par des 0)
    max_row_beta <- max(sapply(li_beta_matrix, nrow))
    li_beta_matrix <- lapply(li_beta_matrix, function(beta_matrix) {
        if (nrow(beta_matrix) < max_row_beta) {
            rows_to_add <- max_row_beta - nrow(beta_matrix)
            zero_rows <- matrix(0, nrow = rows_to_add, ncol = ncol(beta_matrix))
            beta_matrix <- rbind(beta_matrix, zero_rows)
        }
        return(beta_matrix)
    })

    # Créer les noms de variables
    li_names <- lapply(1:length(name_picto), function(n) {
        beta_matrix <- li_beta_matrix[[n]]
        vec_names <- c()
        for (i in 1:nrow(beta_matrix)) {
            vec_names <- c(vec_names, sapply(1:ncol(beta_matrix), function(j) {
                paste0("bloc_", n, "_row_", i, "_col_", j)
            }))
        }
        return(vec_names)
    })
    # Mettre les noms de vecteurs en variables
    vec_names <- do.call(c, li_names)

    # Déplier beta
    li_beta_vec <- lapply(1:length(li_beta_matrix), function(n) {
        beta_matrix <- li_beta_matrix[[n]]
        beta_vec <- c(t(beta_matrix))
        return(beta_vec)
    })
    beta_vec <- do.call(c, li_beta_vec)

    # Nommer beta
    names(beta_vec) <- vec_names

    # Générer les données complexes
    print("start gene")
    X <- gene_x_scalar(config_extrac, beta_vec)
    print("end gene")

    colnames(X) <- vec_names





    # # renommer beta_vec
    # lapply(1:length(li_beta_vec), function(n) {
    #     names(li_beta_vec[[n]]) <- colnames(li_X[[n]])
    # })



    # coller les matrices beta
    beta_matrix <- glue_mats(li_beta_matrix)

    # Calculer le bruit puis y
    # proba <- 1 / (1 + exp(-X %*% beta_vec))
    # y <- rbinom(length(proba), size = 1, prob = proba)
    y <- c(rep("Classe_1", round(config_extrac$prop_class_1 * nrow(X))), rep("Classe_0", nrow(X) - round(config_extrac$prop_class_1 * nrow(X))))

    # Sauvegarder data_used
    data_used <- cbind(y, as.data.frame(X))
    colnames(data_used)[1] <- "beta_class"
    if (sysname == "Linux") {
        write.csv(data_used, paste0(path_data, "/data_used.csv"), row.names = FALSE)
        write_xlsx(data_used, paste0(path_data, "/data_used.xlsx"))
        write_xlsx(as.data.frame(beta_matrix), paste0(path_data, "/beta_picto/big_picto.xlsx"))
        create_heatmap(paste0(path_data, "/beta_picto/big_picto_", config_extrac$nom_spe, ".xlsx"), paste0(path_data, "/beta_picto/big_picto_heatmap_", config_extrac$nom_spe, ".png"), as.data.frame(beta_matrix))
        matrix_big_beta <- as.matrix(beta_matrix)
    } else {
        write.csv(data_used, "..\\data\\data_used.csv", row.names = FALSE)
        write_xlsx(data_used, "..\\data\\data_used.xlsx")
        write_xlsx(as.data.frame(beta_matrix), "..\\data\\beta_picto\\big_picto.xlsx")
    }

    # Créer les index utiles pour comprendre les données
    info_cols <- list(exclude_cols = c(), explained_col = c("beta_class"))
    is_binary <- rep(FALSE, ncol(X))
    index_bloc <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        return(rep(l, nrow(beta_matrix) * ncol(beta_matrix)))
    })))

    index_mode <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        vec_local <- unlist(lapply(1:nrow(beta_matrix), function(i) {
            return(rep(i, ncol(beta_matrix)))
        }))
        return(vec_local)
    })))
    value_max <- 0
    index_variable <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        if (l == 1) {
            n_col_previous <- 0
        } else {
            n_col_previous <- sum(sapply(1:(l - 1), function(i) {
                return(ncol(li_beta_matrix[[i]]))
            }))
        }
        vec_local <- unlist(lapply(1:nrow(beta_matrix), function(i) {
            return(1:ncol(beta_matrix) + n_col_previous)
        }))
    })))

    name_bloc <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        return(rep(paste0("bloc_", l), nrow(beta_matrix) * ncol(beta_matrix)))
    })))

    name_mode <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        vec_local <- unlist(lapply(1:nrow(beta_matrix), function(i) {
            return(rep(paste0("mode_", i), ncol(beta_matrix)))
        }))
        return(vec_local)
    })))

    name_variable <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        vec_local <- unlist(lapply(1:nrow(beta_matrix), function(i) {
            return(paste0("Variable_", 1:ncol(beta_matrix)))
        }))
    })))



    # print(index_variable)

    saveRDS(info_cols, file = paste0(path_data, "/RDS/info_cols.rds"))
    saveRDS(is_binary, file = paste0(path_data, "/RDS/is_binary.rds"))
    saveRDS(index_mode, file = paste0(path_data, "/RDS/index_mode.rds"))
    saveRDS(index_bloc, file = paste0(path_data, "/RDS/index_bloc.rds"))
    saveRDS(index_variable, file = paste0(path_data, "/RDS/index_variable.rds"))

    saveRDS(name_mode, file = paste0(path_data, "/RDS/name_mode.rds"))
    saveRDS(name_bloc, file = paste0(path_data, "/RDS/name_bloc.rds"))
    saveRDS(name_variable, file = paste0(path_data, "/RDS/name_variable.rds"))
    return(list(index_bloc = index_bloc, name_bloc = name_bloc, index_mode = index_mode, name_mode = name_mode, index_variable = index_variable, name_variable = name_variable, is_binary = is_binary, info_cols = info_cols, data_used = data_used, matrix_big_beta = matrix_big_beta))
}

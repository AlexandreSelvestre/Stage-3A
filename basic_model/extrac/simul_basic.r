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

### A ECRIRE: facile

# Beta sera décomposé ligne par ligne pour rester cohérent
extract_all <- function(config_extrac, sysname) {
    # Récupérer les dossiers des pictogrammes d'intérêt
    vec_num_picto <- config_extrac$num_picto
    name_picto <- sapply(vec_num_picto, function(num_picto) {
        vec_path <- list.dirs("../data/beta_picto", recursive = FALSE, full.names = FALSE)
        good_dir <- vec_path[substr(vec_path, 1, length(as.character(num_picto))) == as.character(num_picto)]
        return(good_dir)
    })

    # Récupérer les matrices X et beta des pictogrammes sélectionnés
    li_beta_matrix <- lapply(name_picto, function(good_dir) {
        path <- paste0("../data/beta_picto/", good_dir, "/Beta0.csv")
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

    # Générer les matrices X pour chaque pictogramme
    li_X <- mclapply(1:length(name_picto), function(n) {
        good_dir <- name_picto[n]
        beta_matrix <- li_beta_matrix[[n]]
        path <- paste0("../data/beta_picto/", good_dir, "/X.csv")
        X <- gene_x_scalar(config_extrac, beta_matrix)
        write.csv(X, path, row.names = FALSE)
        X <- as.matrix(read.csv(path))
        vec_names <- c()
        for (i in 1:nrow(beta_matrix)) {
            vec_names <- c(vec_names, sapply(1:ncol(beta_matrix), function(j) {
                paste0("bloc_", n, "_row_", i, "_col_", j)
            }))
        }
        colnames(X) <- vec_names
        return(X)
    }, mc.cores = detectCores() - 1)

    # Coller les X et déplier les beta
    X <- do.call(cbind, li_X)
    li_beta_vec <- lapply(1:length(li_beta_matrix), function(n) {
        beta_matrix <- li_beta_matrix[[n]]
        beta_vec <- c(t(beta_matrix))
        return(beta_vec)
    })

    # renommer beta_vec
    lapply(1:length(li_beta_vec), function(n) {
        names(li_beta_vec[[n]]) <- colnames(li_X[[n]])
    })

    # Coller les beta dépliés
    beta_vec <- do.call(c, li_beta_vec)

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
        write.csv(data_used, "..//data//data_used.csv", row.names = FALSE)
        write_xlsx(data_used, "..//data//data_used.xlsx")
        write_xlsx(as.data.frame(beta_matrix), "..//data//beta_picto//big_picto.xlsx")
        create_heatmap("..//data//beta_picto//big_picto.xlsx", "..//data//beta_picto//big_picto_heatmap.png")
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
    index_variable <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        vec_local <- unlist(lapply(1:nrow(beta_matrix), function(i) {
            return(1:ncol(beta_matrix))
        }))
    })))

    name_bloc <- unname(unlist(lapply(1:length(li_beta_matrix), function(l) {
        beta_matrix <- li_beta_matrix[[l]]
        return(rep("bloc l", nrow(beta_matrix) * ncol(beta_matrix)))
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


    saveRDS(info_cols, file = "../data/RDS/info_cols.rds")
    saveRDS(is_binary, file = "../data/RDS/is_binary.rds")
    saveRDS(index_mode, file = "../data/RDS/index_mode.rds")
    saveRDS(index_bloc, file = "../data/RDS/index_bloc.rds")
    saveRDS(index_variable, file = "../data/RDS/index_variable.rds")

    saveRDS(name_mode, file = "../data/RDS/name_mode.rds")
    saveRDS(name_bloc, file = "../data/RDS/name_bloc.rds")
    saveRDS(name_variable, file = "../data/RDS/name_variable.rds")
}

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

# source("./utils/utils.r")
# import_folder("./utils")

# config_extrac <- config::get(file = "configs/extrac/config_extrac_multi_slice.yml", config = "my_config")

extract_all <- function(config_extrac, sys_name) {
    brute_data <- read_excel("../data/multislice_excel.xlsx")
    brute_data <- as.data.frame(brute_data)


    ### Degager les colonnes fausses
    # Column1 est fantômatique
    exclure <- c("Column1")
    brute_data <- brute_data[, setdiff(colnames(brute_data), exclure)]

    # Les colonnes diagnostics sont inutiles
    col_diagnos <- substr(colnames(brute_data), 1, 11) == "diagnostics"
    brute_data <- brute_data[, !col_diagnos]

    # Toute colonne constante ne sert à rien
    vec_sd <- as.vector(apply(brute_data[, 5:ncol(brute_data)], 2, sd))
    col_const <- vec_sd < 10^-6
    col_const <- c(rep(FALSE, 4), col_const)
    brute_data <- brute_data[, !col_const]

    # Dégager les mixtes si demandé
    if (config_extrac$kill_mixtes) {
        brute_data <- brute_data[brute_data$classe_name != "Mixtes", ]
    }

    ### Donner des entiers aux slice_num
    brute_data$slice_num <- as.integer(brute_data$slice_num)

    ### refaire l'indexation pour avoir un indexe différent entre CHC, CCK et Mixtes
    li_index_type <- list() # liste par type de tumeur des index des patients
    for (i in seq_len(nrow(brute_data))) {
        type <- brute_data$classe_name[i]
        if (!type %in% names(li_index_type)) {
            li_index_type[[type]] <- c()
        }
        patient_num <- brute_data$patient_num[i]
        if (!patient_num %in% li_index_type[[type]]) {
            li_index_type[[type]] <- c(li_index_type[[type]], patient_num)
        }
    }

    li_index_ajouter_per_type <- list() # liste par type de tumeur des index à ajouter aux patient_num
    ajouter <- 0
    for (type in names(li_index_type)) {
        li_index_ajouter_per_type[[type]] <- ajouter
        ajouter <- ajouter + max(li_index_type[[type]])
    }

    # Noucelle organisation: dans l'ordre des types
    brute_data$key <- paste0(brute_data$classe_name, "_", brute_data$patient_num) # conserver l'ancienne indexation
    for (i in seq_len(nrow(brute_data))) {
        type <- brute_data$classe_name[i]
        former_patient_num <- brute_data$patient_num[i]
        new_num <- former_patient_num + li_index_ajouter_per_type[[type]]
        brute_data$patient_num[i] <- new_num
        type <- substr(type, 1, 3)
        key <- paste0(type, "_", former_patient_num)
    }


    ###### Degager toutes les lignes avec le mauvais temps considéré
    time_inj <- config_extrac$time_inj
    brute_data <- brute_data[brute_data$temps_inj %in% time_inj, ]

    ### Noter toutes les associations (patient_num, slice_num qui ne sont pas complètes...).
    ### On retirera les patients pour qui cela arrive à plus du tiers des slices. On le dégage aussi si cela se produit sur trois slices consécutives.
    patients_trop_incomplet <- c() # patient_num integer des patients qui ont trop de slices dont l'un des temps saute

    li_slice_per_time <- list() # liste par patient puis par slice num contenant le vecteur des temps présents
    li_num_slice_per_time <- list() # ici les slices seront juste un vecteur de numéros par patients (pas de vecteur de temps), les numéros seront des entiers (vs. des characters dans la liste précédente)



    for (i in seq_len(nrow(brute_data))) {
        patient_num <- as.character(brute_data$patient_num[i])
        time <- brute_data$temps_inj[i]
        slice_num <- as.character(brute_data$slice_num[i])
        if (!patient_num %in% names(li_slice_per_time)) {
            li_slice_per_time[[patient_num]] <- list()
            li_num_slice_per_time[[patient_num]] <- c()
        }
        if (!slice_num %in% names(li_slice_per_time[[patient_num]])) {
            li_slice_per_time[[patient_num]][[slice_num]] <- c()
            li_num_slice_per_time[[patient_num]] <- c(li_num_slice_per_time[[patient_num]], as.integer(slice_num))
        }
        if (!time %in% li_slice_per_time[[patient_num]][[slice_num]]) {
            li_slice_per_time[[patient_num]][[slice_num]] <- c(li_slice_per_time[[patient_num]][[slice_num]], time)
        }
    }

    ### Remettre en ordre
    for (patient_num in names(li_num_slice_per_time)) {
        li_num_slice_per_time[[patient_num]] <- sort(li_num_slice_per_time[[patient_num]])
    }
    li_distance_to_full_slice_per_slice_per_patient <- list() # indiquera la taille du trou pour les slices manquantes (trou avant et trou après)

    brute_data <- brute_data[order(brute_data$patient_num, brute_data$slice_num), ] # Ordonner

    ### Récupérer les bonnes slices
    nb_reshaped <- 0 # nb d'individus dont on a changé le nombre de slices de tumeur ( en dégageant début/fin)
    expected_size <- length(time_inj)
    patient_vide <- list()
    for (patient_num in names(li_slice_per_time)) {
        count_absent <- 0
        true_size <- 0
        consecutif <- 0
        start_tumeur <- FALSE
        exist_full_slice <- FALSE
        max_consecutif_for_this_slice <- 0
        before_end <- 0
        pos_i_to_remove <- c() # positions des "dummy slices" (sans tumeur: début et fin): à retirer absolument
        li_distance_to_full_slice_per_slice_per_patient[[patient_num]] <- list()
        dist_avant <- +Inf
        previous_slice_nums <- c()
        for (i in seq_along(li_num_slice_per_time[[patient_num]])) {
            slice_num <- li_num_slice_per_time[[patient_num]][i]
            slice_num <- as.character(slice_num)
            if (length(li_slice_per_time[[patient_num]][[slice_num]]) != expected_size) {
                if (config_extrac$strong_exclusion) {
                    condit <- length(li_slice_per_time[[patient_num]][[slice_num]]) > 0
                } else {
                    condit <- FALSE
                }
                if (condit) {
                    start_tumeur <- TRUE # alors la tumeur a démarré
                    before_end <- 0 # On se considère dans la tumeur: on en est donc pas encore sorti
                } else {
                    before_end <- before_end + 1 # On n'est pas dans la tumeur. Si on y revient pas après, on le saura
                    if (!start_tumeur) {
                        pos_i_to_remove <- c(pos_i_to_remove, i) # On est dans une "dummy slice" d'avant tumeur
                    }
                }

                if (start_tumeur) { # Si la tumeur avait démarré on est dans une slice manquante
                    count_absent <- count_absent + 1
                    consecutif <- consecutif + 1
                    true_size <- true_size + 1
                    li_distance_to_full_slice_per_slice_per_patient[[patient_num]][[slice_num]] <- c(dist_avant, Inf)
                    dist_avant <- dist_avant + 1
                    previous_slice_nums <- c(previous_slice_nums, slice_num)
                }
            } else { # On est dans une slice complète
                exist_full_slice <- TRUE
                li_distance_to_full_slice_per_slice_per_patient[[patient_num]][[slice_num]] <- c(0, 0)
                dist_avant <- 1 # distance à la tumeur si la prochaine slice est manquante
                if (length(previous_slice_nums) > 0) {
                    dist_apres <- length(previous_slice_nums)
                    for (num_previous_slice in previous_slice_nums) {
                        li_distance_to_full_slice_per_slice_per_patient[[patient_num]][[num_previous_slice]][2] <- dist_apres
                        dist_apres <- dist_apres - 1
                    }
                }
                previous_slice_nums <- c()
                # On met tout à jur se sachant en slice complète
                max_consecutif_for_this_slice <- max(max_consecutif_for_this_slice, consecutif)
                consecutif <- 0
                start_tumeur <- TRUE
                true_size <- true_size + 1
                before_end <- 0
            }
        }
        consecutif <- consecutif - before_end
        max_consecutif_for_this_slice <- max(max_consecutif_for_this_slice, consecutif)
        count_absent <- count_absent - before_end
        true_size <- true_size - before_end
        seq_end <- (length(li_num_slice_per_time[[patient_num]]) - before_end + 1):length(li_num_slice_per_time[[patient_num]])
        if (before_end > 0) { # Dégager les dummy tranches de la fin
            pos_i_to_remove <- c(pos_i_to_remove, seq_end)
        }
        no_suppress_slice_num <- FALSE
        if (!start_tumeur) { # La tumeur n'a jamais démarré
            count_absent <- length(li_num_slice_per_time[[patient_num]])
            max_consecutif_for_this_slice <- length(li_num_slice_per_time[[patient_num]])
            true_size <- length(li_num_slice_per_time[[patient_num]])
            no_suppress_slice_num <- TRUE
        }
        patient_vide[[patient_num]] <- !exist_full_slice
        max_tolerable <- true_size * config_extrac$max_prop_acceptable
        max_consecutif_acceptable <- true_size * config_extrac$max_prop_consecutif_acceptable


        # print(glue("Patient {patient_num} a {max_consecutif_for_this_slice} slices conseq incomplets sur {max_consecutif_acceptable} slices autorisées. Et true size: {true_size}, before_end {before_end}"))
        if (count_absent > max_tolerable || max_consecutif_for_this_slice > max_consecutif_acceptable) {
            patients_trop_incomplet <- c(patients_trop_incomplet, as.integer(patient_num))
        }
        if (no_suppress_slice_num) {
            pos_i_to_remove <- c()
        }

        if (length(pos_i_to_remove) > 0) {
            li_num_slice_per_time[[patient_num]] <- li_num_slice_per_time[[patient_num]][-pos_i_to_remove]
            nb_reshaped <- nb_reshaped + 1
        }
        if (length(li_num_slice_per_time[[patient_num]]) != true_size) {
            print(paste("Patient", patient_num, "n'a pas le bon nombre de slices"))
            print(li_num_slice_per_time[[patient_num]])
            stop("Erreur: patient n'a pas le bon nombre de true slices ")
        }
    }
    # print(patient_vide)

    ### Print le nombre de patients dans chaque catégorie avant puis après exclusion (ainsi que leur somme)
    # print(length(unique(brute_data[brute_data$classe_name == "CHC", ]$patient_num)))
    # print(length(unique(brute_data[brute_data$classe_name == "CCK", ]$patient_num)))
    # print(length(unique(brute_data$patient_num)))
    # brute_data <- brute_data[!brute_data$patient_num %in% patients_trop_incomplet, ]
    # # print(dim(brute_data))
    # print(length(unique(brute_data[brute_data$classe_name == "CHC", ]$patient_num)))
    # print(length(unique(brute_data[brute_data$classe_name == "CCK", ]$patient_num)))
    # print(length(unique(brute_data$patient_num)))

    # Définir les colonnes à exclure
    exclude_cols <- c("patient_num", "key") # non explicatives
    explained_col <- c("classe_name")
    crushed_cols <- c("slice_num", "temps_inj") # Disparaissent purement et simplement quand on réorganiuse les colonnes
    no_multivariate_col <- c(exclude_cols, explained_col)
    col_to_concatenate <- setdiff(colnames(brute_data), c(no_multivariate_col, crushed_cols))

    ### Extraire seulement les slices d'intérêt
    nb_to_extract <- config_extrac$n_slices
    li_patient <- lapply(unique(brute_data$patient_num), function(patient_num) {
        brute_data[brute_data$patient_num == patient_num, ]
        true_vec_slices <- li_num_slice_per_time[[as.character(patient_num)]]
        true_nb_slices <- length(true_vec_slices)
        if (patient_vide[[as.character(patient_num)]]) {
            return(matrix(NA, nrow = 0, ncol = 0))
        }
        vec_indices_to_extrac <- sapply(seq_len(nb_to_extract), function(i) {
            frac_to_extract <- i / nb_to_extract
            index_to_extract_float <- frac_to_extract * (true_nb_slices - 1) + 1
            # print(patient_num)
            # print(true_nb_slices)
            # print(li_distance_to_full_slice_per_slice_per_patient[[as.character(patient_num)]])
            distances_to_index_float <- data.table::copy(li_distance_to_full_slice_per_slice_per_patient[[as.character(patient_num)]][[round(index_to_extract_float)]]) # le c(avant, après) de la slice qu'on voudrait si elle existait...
            erreur <- round(index_to_extract_float) - index_to_extract_float
            distances_to_index_float_calculus <- data.table::copy(distances_to_index_float)
            distances_to_index_float_calculus[1] <- distances_to_index_float_calculus[1] + erreur
            distances_to_index_float_calculus[2] <- distances_to_index_float_calculus[2] - erreur
            avant_ou_apres <- which.min(distances_to_index_float_calculus)
            distances_to_index_float[1] <- -distances_to_index_float[1] # valeur algébrique du shift
            chosen_shift <- distances_to_index_float[avant_ou_apres]
            chosen_slice_position <- round(index_to_extract_float) + chosen_shift
            chosen_slice_num <- true_vec_slices[chosen_slice_position]
            temps_for_chosen_slice <- li_slice_per_time[[as.character(patient_num)]][[as.character(chosen_slice_num)]]
            if (length(temps_for_chosen_slice) != expected_size) {
                # print("start")
                # print(index_to_extract_float)
                # print(distances_to_index_float)
                # print(chosen_shift)
                # print(li_distance_to_full_slice_per_slice_per_patient[[as.character(patient_num)]])
                print(paste("Slice", chosen_slice_num, "patient", patient_num, "n'a pas le bon nombre de temps"))
                print(temps_for_chosen_slice)
                stop("Erreur: slice choisie n'a pas le bon nombre de temps")
            }
            return(chosen_slice_num)
        })

        df_patient <- brute_data[brute_data$patient_num == patient_num, ]
        new_df_patient <- as.data.frame(matrix(NA, nrow = 1, ncol = length(no_multivariate_col) + length(time_inj) * length(col_to_concatenate) * nb_to_extract))

        # Créer les noms de colonnes
        vec_names <- rep(col_to_concatenate, length(time_inj) * nb_to_extract)
        for (i in seq_along(time_inj)) {
            for (j in seq_len(nb_to_extract)) {
                borne_basse <- (i - 1) * length(col_to_concatenate) * nb_to_extract + (j - 1) * length(col_to_concatenate) + 1
                borne_haute <- (i - 1) * length(col_to_concatenate) * nb_to_extract + j * length(col_to_concatenate)
                vec_names[borne_basse:borne_haute] <- paste0(vec_names[borne_basse:borne_haute], "_slice_", j, "_", config_extrac$new_names_time_inj[i])
            }
        }

        vec_names <- c(no_multivariate_col, vec_names)
        colnames(new_df_patient) <- vec_names
        if (length(time_inj) != length(config_extrac$new_names_time_inj)) {
            stop("Erreur: pas le bon nombre de nouveaux noms de temps")
        }

        row <- c()
        for (time in time_inj) {
            for (j in seq_along(vec_indices_to_extrac)) {
                slice_num <- vec_indices_to_extrac[j]
                df_patient_slice <- df_patient[df_patient$slice_num == slice_num, ]
                row_time_slice <- df_patient_slice[df_patient_slice$temps_inj == time, col_to_concatenate]
                row_time_slice <- unname(unlist(row_time_slice))
                row <- c(row, row_time_slice)
            }
        }
        row <- c(unname(unlist(df_patient_slice[1, no_multivariate_col])), row)
        new_df_patient[1, ] <- row
        if (any(is.na(new_df_patient))) {
            stop("Erreur: il reste des NA dans le nouveau df")
        }
        return(new_df_patient)
    })
    names(li_patient) <- unique(as.character(brute_data$patient_num))
    data_radio_in_lines <- do.call(rbind, li_patient)
    ### On récupère la taille de la tumeur
    df_size <- as.data.frame(matrix(NA, nrow = nrow(data_radio_in_lines), ncol = 2))
    i <- 1
    for (patient_num in names(li_patient)) {
        if (ncol(li_patient[[patient_num]]) > 0) {
            true_size <- length(li_num_slice_per_time[[patient_num]])
            df_size[i, 1] <- data_radio_in_lines$key[i]
            df_size[i, 2] <- true_size
            i <- i + 1
        }
    }
    colnames(df_size) <- c("key", "size_tumeur")


    data_radio_in_lines <- data.frame(lapply(data_radio_in_lines, convert_to_num))
    if (config_extrac$kill_outliers == TRUE) {
        non_concerned_outliers <- c("patient_num", "classe_name", "key", "Gender")
        col_outli <- setdiff(names(data_radio_in_lines), non_concerned_outliers)
        to_be_treated <- data.frame(lapply(data_radio_in_lines[, col_outli], as.numeric))
        ancient <- copy(data_radio_in_lines)
        data_radio_in_lines[, col_outli] <- treat_outliers(to_be_treated, percentile = 0.02)
        print(paste("Outliers killed: ", sum(ancient[, col_outli] != data_radio_in_lines[, col_outli], na.rm = TRUE), "out of", nrow(ancient[, col_outli]) * ncol(ancient[, col_outli])))
    }

    data_patients <- read_excel("../data/Descriptif_patients.xlsx")
    data_patients <- change_dates(data_patients)
    data_patients <- data_patients[config_extrac$patients_col]
    data_patients$Gender <- unname(sapply(data_patients$Gender, change_genre))
    data_patients$key <- paste0(data_patients$classe_name, "_", data_patients$patient_num)

    data_radio_in_lines <- data_radio_in_lines[, setdiff(colnames(data_radio_in_lines), c("patient_num", "classe_name"))] # Pour éviter les doublons de colonnes


    data_used <- merge(data_radio_in_lines, data_patients, by = c("key"))
    data_used <- na.omit(data_used)

    ### Remettre en ordre les clonnes de data_used (les non explicatives au début)
    col_names <- colnames(data_used)
    two_col_key <- c("patient_num", "classe_name")
    other_cols <- col_names[!col_names %in% two_col_key]
    new_col_order <- c(two_col_key, other_cols)
    data_used <- data_used[, new_col_order]
    data_used <- merge(data_used, df_size, by = c("key"))

    data_used <- data_used[, setdiff(colnames(data_used), c("key"))]

    write_xlsx(data_used, "..//data//data_used.xlsx")
    write.csv(data_used, "..//data//data_used.csv", row.names = FALSE)

    ### Créer les autres infos indispensables, notamment pour du traitement multiway multibloc

    exclude_cols <- c("patient_num")
    explained_col <- c("classe_name")
    info_cols <- list(exclude_cols = exclude_cols, explained_col = explained_col)
    saveRDS(info_cols, file = "../data/RDS/info_cols.rds")

    is_binary <- rep(FALSE, ncol(data_used))
    if ("Gender" %in% names(data_used)) {
        is_binary[which(names(data_used) == "Gender")] <- TRUE
    }

    train_cols <- setdiff(names(data_used), c(exclude_cols, explained_col))

    li_index_bloc <- get_bloc_vec_liver(names(data_used[, train_cols]), exist_temps = TRUE, different_clin = FALSE)
    index_bloc <- li_index_bloc$index_bloc
    name_bloc <- li_index_bloc$index_name

    li_index_variable <- get_variable_vec_liver_big(data_used[, train_cols])
    index_variable <- li_index_variable$index_variable
    name_variable <- li_index_variable$index_name

    li_index_modes <- list()
    li_name_modes <- list()

    res_index_temps <- get_mode_vec_liver(data_used[, train_cols])
    li_index_modes$mode_temps <- res_index_temps$index_mode
    li_name_modes$mode_temps <- res_index_temps$index_name

    index_slice <- rep(-1, ncol(data_used[, train_cols]))
    name_slice <- rep("clinique", ncol(data_used[, train_cols]))
    # 2792 expected
    nb_multiway <- length(col_to_concatenate)
    for (i in seq_along(time_inj)) {
        for (j in seq_len(nb_to_extract)) {
            debut <- (i - 1) * nb_to_extract * nb_multiway + (j - 1) * nb_multiway + 1
            fin <- (i - 1) * nb_to_extract * nb_multiway + j * nb_multiway
            index_slice[debut:fin] <- j
            name_slice[debut:fin] <- paste0("slice_", j)
        }
    }

    li_index_modes$mode_slice <- index_slice
    li_name_modes$mode_slice <- name_slice

    saveRDS(index_bloc, file = "../data/RDS/index_bloc.rds")
    saveRDS(name_bloc, file = "../data/RDS/name_bloc.rds")

    saveRDS(li_index_modes, file = "../data/RDS/li_index_modes.rds")
    saveRDS(li_name_modes, file = "../data/RDS/li_name_modes.rds")

    saveRDS(index_variable, file = "../data/RDS/index_variable.rds")
    saveRDS(name_variable, file = "../data/RDS/name_variable.rds")

    saveRDS(is_binary, file = "../data/RDS/is_binary.rds")


    print("Fin extraction multislice")
}

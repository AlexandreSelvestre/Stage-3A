extract_all <- function(config_extrac, sysname) {
    ##### Importer les données ########
    # radiomic_global <- "/radiomiques_global.xlsx"
    radiomic_global <- "/global_excel.xlsx"
    path_data <- config_extrac$path_data
    if (sysname == "Linux") {
        data_radio <- read_excel(paste0(path_data, "/radiomiques_global.xlsx"))
        data_patients <- read_excel(paste0(path_data, "/Descriptif_patients.xlsx"))
    } else {
        data_radio <- read_excel("..\\data\\radiomiques_global.xlsx")
        data_patients <- read_excel("..\\data\\Descriptif_patients.xlsx")
    }

    #### Changer les dates en string ########
    data_radio <- change_dates(data_radio)
    data_patients <- change_dates(data_patients)

    #### Garder seulement les colonnes pertinentes pour la prédiction #######

    data_patients <- data_patients[config_extrac$patients_col]
    data_radio <- data_radio[, setdiff(names(data_radio), config_extrac$kill_col)]

    col_names_diagnos <- c()
    for (col_name in names(data_radio)) {
        if (substr(col_name, 1, 11) == "diagnostics") {
            col_names_diagnos <- append(col_names_diagnos, col_name)
        }
    }
    data_radio <- data_radio[, setdiff(names(data_radio), col_names_diagnos)]

    if (config_extrac$first_order_only) {
        keep_cols <- c("classe_name", "temps_inj", "patient_num")
        for (col_name in names(data_radio)) {
            if (grepl("firstorder", col_name)) {
                keep_cols <- append(keep_cols, col_name)
            }
        }
        data_radio <- data_radio[, keep_cols]
    }

    if (config_extrac$shape_only) {
        keep_cols <- c("classe_name", "temps_inj", "patient_num")
        for (col_name in names(data_radio)) {
            if (grepl("shape", col_name)) {
                keep_cols <- append(keep_cols, col_name)
            }
        }
        data_radio <- data_radio[, keep_cols]
    }

    if (config_extrac$texture_only) {
        keep_cols <- c("classe_name", "temps_inj", "patient_num")
        for (col_name in names(data_radio)) {
            if (grepl("glcm", col_name) | grepl("gldm", col_name) | grepl("glrl", col_name) | grepl("glsz", col_name) | grepl("ngtd", col_name)) {
                keep_cols <- append(keep_cols, col_name)
            }
        }
        data_radio <- data_radio[, keep_cols]
    }


    if (config_extrac$other_only) {
        kill_cols <- c()
        for (col_name in names(data_radio)) {
            if (grepl("firstorder", col_name)) {
                kill_cols <- append(kill_cols, col_name)
            }
        }
        data_radio <- data_radio[, setdiff(names(data_radio), kill_cols)]
    }

    #### Tuer les outliers dans les radiomiques ##########

    if (config_extrac$kill_outliers == TRUE) {
        non_concerned_outliers <- c("patient_num", "classe_name", "temps_inj")
        col_outli <- setdiff(names(data_radio), non_concerned_outliers)
        to_be_treated <- data.frame(lapply(data_radio[, col_outli], as.numeric))
        data_radio <- data.frame(lapply(data_radio, convert_to_num))
        # write_xlsx(data_radio, "..\\data\\sauver.xlsx")
        ancient <- copy(data_radio)
        data_radio[, col_outli] <- treat_outliers(to_be_treated, percentile = 0.01)
        print(paste("Outliers killed: ", sum(ancient != data_radio, na.rm = TRUE)))
        # write_xlsx(data_radio, "..\\data\\sauver_post.xlsx")
    }


    data_patients$Gender <- unname(sapply(data_patients$Gender, change_genre))

    former_col_names <- colnames(data_radio)
    unique_col_names <- c("keys", "patient_num", "classe_name")
    multiple_former_col_names <- former_col_names[setdiff(former_col_names, unique_col_names)]
    multiple_former_col_names <- setdiff(former_col_names, c("classe_name, temps_inj, patient_num"))


    new_col_names <- unique_col_names
    for (temps_inj in config_extrac$time_inj) {
        for (col_name in multiple_former_col_names) {
            new_col_names <- append(new_col_names, glue("{col_name}_{temps_inj}"))
        }
    }

    dict_temps <- list()
    for (i in 1:length(config_extrac$time_inj)) {
        dict_temps[[config_extrac$time_inj[i]]] <- i - 1
    }
    # Dictionnaire des indices des placements dans la nouvelle tables

    ## Créer la ligne "de base" qui définit les types des colonnes

    li_init <- data.frame(matrix(ncol = length(new_col_names), nrow = 1))
    colnames(li_init) <- new_col_names
    data_radio_in_lines <- data.frame(matrix(ncol = length(new_col_names), nrow = 0))
    colnames(data_radio_in_lines) <- new_col_names

    ## Créer les keys des lignes
    li_patient <- unique(data_radio[["patient_num"]])
    if (config_extrac$kill_mixtes) {
        types <- c("CHC", "CCK")
    } else {
        types <- c("CHC", "CCK", "Mixtes")
    }
    for (patient in li_patient) {
        for (type in types) {
            if (any(data_radio$classe_name == type & data_radio$patient_num == patient)) {
                key <- paste(patient, type, sep = ",")
                new_li <- copy(li_init)
                new_li[[1]] <- key
                new_li[[2]] <- patient
                new_li[[3]] <- type
                data_radio_in_lines <- rbind(data_radio_in_lines, new_li)
            }
        }
    }

    # write_xlsx(data_radio_in_lines, "..\\data\\test.xlsx")


    ## Remplir les lignes
    # data_radio_in_lines <- read_excel("..\\data\\test.xlsx")

    for (i in 1:dim(data_radio)[1]) {
        key <- paste(data_radio[["patient_num"]][i], data_radio[["classe_name"]][i], sep = ",")
        new_row <- data_radio_in_lines[data_radio_in_lines$keys == key, ]
        short_row <- data_radio[i, ]
        name_time_inj <- short_row[["temps_inj"]][1]
        if (name_time_inj == "ART") {
            name_time_inj <- "ART_"
        }
        if (name_time_inj %in% names(dict_temps)) {
            index_temps <- dict_temps[[name_time_inj]]
            start_index_row <- 4 + index_temps * length(multiple_former_col_names)
            new_row[1, start_index_row:(start_index_row + length(multiple_former_col_names) - 1)] <- short_row[, setdiff(colnames(short_row), c("classe_name", "temps_inj", "patient_num"))]
            data_radio_in_lines[data_radio_in_lines$keys == key, ] <- new_row
        }
    }


    # write_xlsx(data_radio_in_lines, "..\\data\\radio_in_lines.xlsx")
    data_used <- merge(data_radio_in_lines, data_patients, by = c("patient_num", "classe_name"))


    if (config_extrac$kill_empty == TRUE) {
        data_used <- na.omit(data_used)
    }

    # unique_col_names <- c("keys", "patient_num", "classe_name")


    if (sysname == "Linux") {
        write_xlsx(data_used, paste0(path_data, "/data_used.xlsx"))
        write.csv(data_used, paste0(path_data, "/data_used.csv"), row.names = FALSE)
    } else {
        write_xlsx(data_used, "..\\data\\data_used.xlsx")
        write.csv(data_used, "..\\data\\data_used.csv", row.names = FALSE)
    }


    exclude_cols <- c("patient_num", "keys")
    explained_col <- c("classe_name")
    info_cols <- list(exclude_cols = exclude_cols, explained_col = explained_col)
    saveRDS(info_cols, file = paste0(path_data, "/RDS/info_cols.rds"))

    # Créer les indices de modes, de variable et de binaire!
    # Ne pas présupposer l'existence de colonnes potentiellement inexistantes...

    is_binary <- rep(FALSE, ncol(data_used))
    if ("Gender" %in% names(data_used)) {
        is_binary[which(names(data_used) == "Gender")] <- TRUE
    }

    train_cols <- setdiff(names(data_used), c(exclude_cols, explained_col))

    li_index_bloc <- get_bloc_vec_liver(names(data_used[, train_cols]), exist_temps = TRUE, different_clin = FALSE)
    index_bloc <- li_index_bloc$index_bloc
    name_bloc <- li_index_bloc$index_name


    li_index_mode <- get_mode_vec_liver(data_used[, train_cols])
    index_mode <- li_index_mode$index_mode
    name_mode <- li_index_mode$index_name


    li_index_variable <- get_variable_vec_liver(data_used[, train_cols])
    index_variable <- li_index_variable$index_variable
    name_variable <- li_index_variable$index_name


    saveRDS(index_bloc, file = paste0(path_data, "/RDS/index_bloc.rds"))
    saveRDS(name_bloc, file = paste0(path_data, "/RDS/name_bloc.rds"))

    saveRDS(index_mode, file = paste0(path_data, "/RDS/index_mode.rds"))
    saveRDS(name_mode, file = paste0(path_data, "/RDS/name_mode.rds"))

    saveRDS(index_variable, file = paste0(path_data, "/RDS/index_variable.rds"))
    saveRDS(name_variable, file = paste0(path_data, "/RDS/name_variable.rds"))

    saveRDS(is_binary, file = paste0(path_data, "/RDS/is_binary.rds"))
}

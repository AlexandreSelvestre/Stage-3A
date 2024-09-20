library(readxl)
library(writexl)
library(glue)
library(data.table)
library(ExPanDaR)
source("utils.r")

########## Ici on extrait le multislices et on le préprocessera quand on en aura le temps ##############

extract_all <- function(config_extrac) {
    data_radio <- read_excel("..\\data\\multislice_excel.xlsx")
    data_patients <- read_excel("..\\data\\Descriptif_patients.xlsx")

    ###### Changer les dates en string ########
    data_radio <- change_dates(data_radio)
    data_patients <- change_dates(data_patients)

    ###### Garder seulement les colonnes pertinentes pour la prédiction #######
    data_patients <- data_patients[config_extrac$patients_col]
    data_radio <- data_radio[, setdiff(names(data_radio), config_extrac$radio_col)]

    col_names_diagnos <- c()
    for (col_name in names(data_radio)) {
        if (substr(col_name, 1, 11) == "diagnostics") {
            col_names_diagnos <- append(col_names_diagnos, col_name)
        }
    }
    data_radio <- data_radio[, setdiff(names(data_radio), col_names_diagnos)]
}

# Si premier ordre seul

if (config_extrac$first_order_only) {
    keep_cols <- c("classe_name", "temps_inj", "patient_num")
    for (col_name in names(data_radio)) {
        if (grepl("firstorder", col_name)) {
            keep_cols <- append(keep_cols, col_name)
        }
    }
    data_radio <- data_radio[, keep_cols]
}

# Si shape seul

if (config_extrac$shape_only) {
    keep_cols <- c("classe_name", "temps_inj", "patient_num")
    for (col_name in names(data_radio)) {
        if (grepl("shape", col_name)) {
            keep_cols <- append(keep_cols, col_name)
        }
    }
    data_radio <- data_radio[, keep_cols]
}


# Si texture seul

if (config_extrac$texture_only) {
    keep_cols <- c("classe_name", "temps_inj", "patient_num")
    for (col_name in names(data_radio)) {
        if (grepl("glcm", col_name) | grepl("gldm", col_name) | grepl("glrl", col_name) | grepl("glsz", col_name) | grepl("ngtd", col_name)) {
            keep_cols <- append(keep_cols, col_name)
        }
    }
    data_radio <- data_radio[, keep_cols]
}

# Si tout sauf first order

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
    # print(sapply(to_be_treated, class))
    data_radio <- data.frame(lapply(data_radio, convert_to_num))
    # write_xlsx(data_radio, "..\\data\\sauver.xlsx")
    ancient <- copy(data_radio)
    data_radio[, col_outli] <- treat_outliers(to_be_treated, percentile = 0.01)
    print(paste("Outliers killed: ", sum(ancient != data_radio, na.rm = TRUE)))
    # write_xlsx(data_radio, "..\\data\\sauver_post.xlsx")

    write_xlsx(data_radio, "..\\data\\sauver_post.xlsx")


    data_patients_simple <- data.frame(lapply(data_patients[, setdiff(names(data_patients), c(
        "patient_num"
    ))], normal))
    data_patients[, setdiff(names(data_patients), "patient_num")] <- data_patients_simple

    data_patients$Gender <- unname(sapply(data_patients$Gender, change_genre))


    ##On normalise vraiment comme ça? C'est étrange...
    data_radio_simple <- data.frame(lapply(data_radio[, setdiff(names(data_radio), c(
        "patient_num"
    ))], normal))
    data_radio[, setdiff(names(data_radio), "patient_num")] <- data_radio_simple
}

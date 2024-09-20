library(readxl)
library(writexl)
source("utils/utils_extrac.r")

original <- read_excel("../data/radiomiques_global.xlsx")
mine <- read_excel("../data/radiomic_global_refit.xlsx")

col_names_diagnos <- c()
for (col_name in names(original)) {
    if (substr(col_name, 1, 11) == "diagnostics") {
        col_names_diagnos <- append(col_names_diagnos, col_name)
    }
}

original <- original[, setdiff(names(original), col_names_diagnos)]

col_names_diagnos <- c()
for (col_name in names(mine)) {
    if (substr(col_name, 1, 11) == "diagnostics") {
        col_names_diagnos <- append(col_names_diagnos, col_name)
    }
}

mine <- mine[, setdiff(names(mine), col_names_diagnos)]
original <- original[, setdiff(names(original), "original_glcm_JointAverage")]

# ordre alphabetique colonnes
mine <- mine[, order(names(mine))]
original <- original[, order(names(original))]

# Gérer les types des colonnes (sinon mess up with order)
mine[] <- lapply(mine, convert_to_num)
original[] <- lapply(original, convert_to_num)

# ordre alphabetique lignes
mine <- mine[order(mine$classe_name, mine$patient_num, mine$temps_inj), ]
original <- original[order(original$classe_name, original$patient_num, original$temps_inj), ]

# Checker différences sur character
mine_char <- mine[, sapply(mine, is.character)]
original_char <- original[, sapply(original, is.character)]
mine_char <- as.matrix(mine_char)
original_char <- as.matrix(original_char)
num_differences <- sum(mine_char != original_char)

# Checker différences sur num
mine_num <- mine[, sapply(mine, is.numeric)]
original_num <- original[, sapply(original, is.numeric)]
mine_num <- as.matrix(mine_num)
original_num <- as.matrix(original_num)
val_difference <- sum(abs(mine_num - original_num))
print(val_difference)

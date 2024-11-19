import_folder <- function(path) {
    files <- list()
    files <- append(files, list.files(path = path, pattern = "\\.r$", full.names = TRUE))
    files <- append(files, list.files(path = path, pattern = "\\.R$", full.names = TRUE))
    # source("./utils/utils_extrac.r")
    # source("./utils/utils_math.r")
    # source("./utils/utils_model.r")
    # source("./utils/utils_model_medic_specific.r")
    # source("./utils/utils.r")
    li <- lapply(files, source)
}

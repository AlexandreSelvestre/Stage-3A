import_folder <- function(path) {
    files <- list()
    files <- append(files, list.files(path = path, pattern = "\\.r$", full.names = TRUE))
    files <- append(files, list.files(path = path, pattern = "\\.R$", full.names = TRUE))
    lapply(files, source)
}

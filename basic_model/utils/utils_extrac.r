normal <- function(x) {
    if (is.character(x)) {
        y <- as.numeric(x)
        if (all(is.na(y))) {
            return(x)
        }
        x <- y
    }
    if (is.numeric(x)) {
        return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
    } else {
        return(x)
    }
}


convert_to_num <- function(x) {
    if (is.character(x)) {
        y <- as.numeric(x)
        if (all(is.na(y))) {
            return(x)
        }
        x <- y
    }
    if (is.numeric(x)) {
        return(x)
    } else {
        return(x)
    }
}

change_dates <- function(tableau) {
    tableau <- data.frame(lapply(tableau, function(col) {
        if (inherits(col, "POSIXct")) {
            return(as.character(col))
        } else {
            return(col)
        }
    }))
    return(tableau)
}

change_genre <- function(elem) {
    if (is.na(elem)) {
        return(NA_integer_)
    }
    if (elem == "F") {
        return(1)
    }
    if (elem == "M") {
        return(0)
    }
}

name_product_terms <- function(x, small_groups = FALSE) {
    if (grepl("self", x) | grepl("crossed", x)) {
        if (grepl("self", x)) {
            type_col <- "self"
            extension <- substr(x, nchar(x) - 3, nchar(x))
            name_col <- substr(x, 6, nchar(x) - 10)
        } else {
            type_col <- "crossed"
            extension <- substr(x, nchar(x) - 8, nchar(x))
            name_col <- substr(x, 9, nchar(x) - 10)
        }
    } else {
        type_col <- "single"
        extension <- substr(x, nchar(x) - 3, nchar(x))
        name_col <- substr(x, 1, nchar(x) - 5)
    }
    if (small_groups) {
        return(paste(type_col, name_col, extension, sep = "_"))
    } else {
        paste(type_col, extension, sep = "_")
    }
}

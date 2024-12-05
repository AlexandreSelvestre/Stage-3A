path_data <- "../data"
df <- as.data.frame(readxl::read_xlsx(paste0(path_data, "/stat.xlsx")))
bad_cols <- c("LR-5", "LR-M", "rim_APHE")
df <- df[, setdiff(colnames(df), bad_cols)]
if (!impute_data) {
    df <- na.omit(df)
}
if (dedoubler_shape) {
    df$shape_2 <- as.integer(df$shape == 2)
    df$shape_1 <- as.integer(df$shape == 1)
    df <- df[, setdiff(colnames(df), "shape")]
}
df <- data.frame(lapply(df, as.numeric))
if (kill_tumor == "Mixte") {
    if (reclassify_Mixtes) {
        df_mixte <- df[df$Tumeur == 2, ]
        df_mixte$Tumeur <- "Mixte"
        df_mixte <- na.omit(df_mixte)
    }
    df <- df[df$Tumeur != 2, ]
}
if (kill_tumor == "CCK") {
    df <- df[df$Tumeur != 3, ]
}
if (kill_tumor == "CHC") {
    df <- df[df$Tumeur != 1, ]
}
df$Tumeur <- sapply(df$Tumeur, function(i) {
    if (i == 1) {
        return("CHC")
    }
    if (i == 3) {
        return("CCK")
    }
    if (i == 2) {
        return("Mixte")
    }
})
write_xlsx(df, paste0(path_data, "/stat_analysed.xlsx"))

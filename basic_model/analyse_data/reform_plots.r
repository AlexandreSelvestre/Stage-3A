setGeneric("reform_beta", function(object) {
    standardGeneric("reform_beta")
})

setMethod("reform_beta", "apply_model", function(object) {
    K <- length(unique(object@index_mode[object@index_mode > -0.5]))
    L <- length(unique(object@index_bloc[object@index_bloc > -0.5]))
    li_d <- sapply(1:L, function(l) {
        var_bloc_l <- object@index_variable[object@index_bloc == l]
        d_l <- length(unique(var_bloc_l))
        return(d_l)
    })
    # print(object@beta_final)
    J <- sum(li_d)
    li_mat <- lapply(1:L, function(l) {
        vec_beta_loc <- object@beta_final[object@index_bloc == l]
        mat_beta_loc <- matrix(vec_beta_loc, nrow = K, ncol = li_d[l], byrow = TRUE)
        return(mat_beta_loc)
    })
    mat_tot <- do.call(cbind, li_mat)
    path_xlsx <- paste0("../data/beta_picto/big_picto_reconstructed_", object@name_model, ".xlsx")
    path_heat <- paste0("../data/beta_picto/heatmap_", object@name_model, ".png")
    write_xlsx(as.data.frame(mat_tot), path_xlsx)
    create_heatmap(path_xlsx, path_heat)
    object <- compare_mat_beta(object, path_xlsx, path_heat)
    return(object)
})

setGeneric("compare_mat_beta", function(object, ...) {
    standardGeneric("compare_mat_beta")
})

setMethod("compare_mat_beta", "apply_model", function(object, path_xlsx, path_heat) {
    mat_reconstructed <- as.matrix(read_excel(path_xlsx))
    mat_origin <- as.matrix(silent_run(read_excel, "..//data//beta_picto//big_picto.xlsx"))
    mat_diff <- unname(mat_reconstructed) - unname(mat_origin) # OUILLE en données réelles
    erreur <- mean(abs(mat_diff))
    object@score_recons <- erreur
    print(paste("L'erreur moyenne de reconstruction du pictogramme est de", erreur))
    return(object)
})

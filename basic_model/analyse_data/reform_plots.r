setGeneric("reform_beta", function(object) {
    standardGeneric("reform_beta")
})

setMethod("reform_beta", "apply_model", function(object) {
    path_data <- object@path_data
    K <- length(unique(object@index_mode[object@index_mode > -0.5]))
    L <- length(unique(object@index_bloc[object@index_bloc > -0.5]))
    li_d <- sapply(1:L, function(l) {
        var_bloc_l <- object@index_variable[object@index_bloc == l]
        d_l <- length(unique(var_bloc_l))
        return(d_l)
    })
    object@beta_final <- object@model$finalModel$beta
    J <- sum(li_d)
    li_mat <- lapply(1:L, function(l) {
        vec_beta_loc <- object@beta_final[object@index_bloc == l]
        # print(li_d[l])
        # print(K)
        # print(length(vec_beta_loc))
        # print(vec_beta_loc)
        mat_beta_loc <- matrix(vec_beta_loc, nrow = K, ncol = li_d[l], byrow = TRUE)
        return(mat_beta_loc)
    })
    mat_tot <- do.call(cbind, li_mat)
    path_xlsx <- paste0(path_data, "/beta_picto/big_picto_reconstructed_", object@name_model, "_", object@nom_spe, ".xlsx")
    path_heat <- paste0(path_data, "/beta_picto/heatmap_", object@name_model, "_", object@nom_spe, ".png")
    write_xlsx(as.data.frame(mat_tot), path_xlsx)
    create_heatmap(path_xlsx, path_heat, as.data.frame(mat_tot))
    object <- compare_mat_beta(object, path_xlsx, path_heat, mat_tot)
    return(object)
})

setGeneric("compare_mat_beta", function(object, ...) {
    standardGeneric("compare_mat_beta")
})

setMethod("compare_mat_beta", "apply_model", function(object, path_xlsx, path_heat, mat_tot) {
    mat_reconstructed <- mat_tot
    mat_origin <- object@li_extrac$matrix_big_beta
    # mat_origin <- as.matrix(silent_run(read_excel, paste0(path_data, "/beta_picto/big_picto.xlsx")))
    mat_diff <- unname(mat_reconstructed) - unname(mat_origin) # OUILLE en données réelles
    erreur <- mean(abs(mat_diff))
    object@score_recons <- erreur
    print(paste("L'erreur moyenne de reconstruction du pictogramme est de", erreur))
    return(object)
})

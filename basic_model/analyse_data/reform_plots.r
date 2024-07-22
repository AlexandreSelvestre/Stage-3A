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
    J <- sum(li_d)
    li_mat <- lapply(1:L, function(l) {
        vec_beta_loc <- object@beta_final[object@index_bloc == l]
        mat_beta_loc <- matrix(vec_beta_loc, nrow = K, ncol = li_d[l], byrow = TRUE)
        return(mat_beta_loc)
    })
    mat_tot <- do.call(cbind, li_mat)
    write_xlsx(as.data.frame(mat_tot), paste0("../data/beta_picto/big_picto_reconstructed", object@name_model, ".xlsx"))
    return(object)
})

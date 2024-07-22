gene_x_scalar <- function(config_extrac, beta_matrix) {
    n_sample <- config_extrac$n_sample
    dist_intra <- config_extrac$dist_intra_class # distance typiqe intra_classe
    ecart_class <- config_extrac$ecart_class # ecart type
    prop <- config_extrac$prop_class_1 # proportion de l'échantillon
    beta_vec <- c(t(beta_matrix))
    mu_0 <- rep(0, length(beta_vec))
    mu_1 <- ecart_class * beta_vec / norm(as.matrix(beta_vec), type = "2")
    sigma <- dist_intra / sqrt(length(beta_vec))
    X_1 <- mvrnorm(n = round(prop * n_sample), mu_1, sigma * diag(length(beta_vec)), tol = 1e-06, empirical = FALSE)
    X_0 <- mvrnorm(n = n_sample - round(prop * n_sample), mu_0, sigma * diag(length(beta_vec)), tol = 1e-06, empirical = FALSE)
    return(rbind(X_1, X_0))
}

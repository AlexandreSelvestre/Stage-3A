gene_x_scalar <- function(config_extrac, beta_matrix) {
    n_sample <- config_extrac$n_sample
    dist_sepa <- config_extrac$dist_sepa # distance typique suivant beta dans une même classe
    dist_non_sepa <- config_extrac$dist_non_sepa # distance typique ortho à beta dans une même classe
    ecart_class <- config_extrac$ecart_class # distance entre les 2 classes
    prop <- config_extrac$prop_class_1 # proportion de l'échantillon
    beta_vec <- c(t(beta_matrix))
    mu_0 <- rep(0, length(beta_vec))
    mu_1 <- ecart_class * beta_vec / norm(as.matrix(beta_vec), type = "2")
    sigma_non_sepa <- dist_non_sepa
    sigma_sepa <- dist_sepa
    P <- complete_orthonormal_basis(beta_vec)
    vec_diag <- c(sigma_sepa^2, rep(sigma_non_sepa^2, length(beta_vec) - 1))
    D <- diag(vec_diag)
    Sigma <- P %*% D %*% t(P)
    X_1 <- mvrnorm(n = round(prop * n_sample), mu_1, Sigma, tol = 1e-06, empirical = FALSE)
    X_0 <- mvrnorm(n = n_sample - round(prop * n_sample), mu_0, Sigma, tol = 1e-06, empirical = FALSE)
    return(rbind(X_1, X_0))
}

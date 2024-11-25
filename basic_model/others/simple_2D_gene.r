library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(glmnet)
library(SGL)
library(doParallel)
library(jsonlite)
library(pROC)
library(ExPanDaR)
library(NbClust)
library(EMCluster)
library(magrittr)
library(parallel)
library(pracma)
library(mvnfast)
library(Rfast)
library(DMwR)
library(themis)
library(reshape2)
library(rlang)

source("./utils/utils.r")
import_folder("./utils")

beta <- c(-2, 1)
n_1 <- 1000
n_2 <- 1000
mu_1 <- -beta / (2 * norm(as.matrix(beta), type = "2"))
mu_2 <- beta / (2 * norm(as.matrix(beta), type = "2"))
D <- diag(c(0.03, 0.09))
P <- complete_orthonormal_basis(beta)
inside <- mat.mult(P, D)
Sigma <- Tcrossprod(inside, P)
X_1 <- rmvn(n_1, mu_1, Sigma)
X_2 <- rmvn(n_2, mu_2, Sigma)
X <- rbind(X_1, X_2)
y <- c(rep(0, n_1), rep(1, n_2))

df <- data.frame(X1 = X[, 1], X2 = X[, 2], y = as.factor(y))

# Tracer les points avec ggplot2
p <- ggplot(df, aes(x = X1, y = X2, color = y)) +
    geom_point() +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = "Explanatory variables", x = "X1", y = "X2") +
    theme_classic() + # Utiliser un thème avec un fond blanc et des axes gradués
    theme(
        axis.line = element_line(color = "black"),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_line(color = "grey90"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
    ) +
    coord_fixed(ratio = 1) + # Assurer un rapport d'aspect 1:1 pour un repère orthonormé
    scale_x_continuous(expand = c(0, 0), limits = c(min(df$X1), max(df$X1))) +
    scale_y_continuous(expand = c(0, 0), limits = c(min(df$X2), max(df$X2))) +
    geom_hline(yintercept = 0, color = "black") + # Ajouter l'axe horizontal
    geom_vline(xintercept = 0, color = "black") + # Ajouter l'axe vertical
    geom_abline(intercept = 0, slope = 2, color = "green", size = 1.5) # Ajouter la droite



# Afficher le graphique
ggsave("./others/2D.png", p, width = 10, height = 20)

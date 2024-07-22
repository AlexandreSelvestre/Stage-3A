setMethod("analyse_data", "apply_model", function(object) {
    variables <- object@data_used[, setdiff(names(object@data_used), object@info_cols$exclude_cols)]
    variables[[object@name_y]] <- 2 - as.integer(as.factor(variables[[object@name_y]]))
    variables <- variables[order(variables[[object@name_y]]), ]
    variables[, setdiff(colnames(variables), c(object@name_y))] <- renormalize_in_model_fit_index_mode(variables[, setdiff(colnames(variables), c(object@name_y))], object@index_variable, object@is_binary)$new_x
    print(dim(variables))
    if (any(is.na(as.matrix(variables)))) {
        cols_with_na <- colnames(variables)[colSums(is.na(variables)) > 0]

        # Afficher les noms des colonnes avec des NA
        print(cols_with_na)
        stop("Les données contiennent des valeurs NA ou infinies.")
    }

    var_to_clust <- variables[, setdiff(colnames(variables), c(object@name_y))]

    if (object@analyse_data$EM) {
        li_clust <- list()
        emobj <- init.EM(
            x = var_to_clust, nclass = 2, lab = NULL, EMC = .EMC,
            stable.solution = TRUE, min.n = NULL, min.n.iter = 10,
            method = "em.EM"
        )
        print("en cours")
        li_clust$Best.partition <- assign.class(var_to_clust,
            emobj = emobj, pi = NULL, Mu = NULL, LTSigma = NULL,
            lab = NULL, return.all = TRUE
        )$class
        print("fait accompli")
    } else {
        li_clust <- NbClust(data = var_to_clust, distance = "euclidean", min.nc = 2, max.nc = 2, method = "ward.D2", index = "kl")
        cat("Meilleur nombre de clusters", li_clust$Best.nc)
        print("Meilleure partition")
        print(li_clust$Best.partition)
        print("Comparaison aux indices de type:")
        print(variables[[object@name_y]]) # 1  = CCK, 0 = CHC
    }
    ## Attention à 46 et 22 chez les CCK

    pca_rotation <- prcomp(variables, center = FALSE, scale. = FALSE, rank. = 10)$rotation
    data_pca <- as.matrix(variables) %*% pca_rotation[, 3:4]


    # pca_rotation <- prcomp(variables, center = FALSE, scale. = FALSE)
    # library(factoextra)
    # fviz_pca_biplot(pca_rotation, habillage = variables$classe_name, geom.ind = "point", pointsize = 2, repel = TRUE)

    # Créer un dataframe pour la visualisation contenant les deux premières composantes principales
    pca_data <- data.frame(PC1 = data_pca[, 1], PC2 = data_pca[, 2], Classe = variables[[object@name_y]])
    pca_data$number <- paste(rownames(variables), unname(li_clust$Best.partition))
    # print(pca_data$number)
    print(li_clust$Best.partition)
    # Visualiser les deux premières composantes principales avec ggplot2
    # Détermination des limites des axes
    x_limits <- range(pca_data$PC1) * 1.1 # Ajout de marge
    y_limits <- range(pca_data$PC2) * 1.1 # Ajout de marge

    # Création du plot PCA
    image <- ggplot(pca_data, aes(x = PC1, y = PC2, color = as.factor(Classe))) +
        geom_point() +
        geom_text(aes(label = number), hjust = -0.3, vjust = 0.5, size = 2) + # Ajout des numéros des points
        labs(
            title = "PCA - 2D Representation of Data",
            x = "First Principal Component",
            y = "Second Principal Component",
            color = "Classe"
        ) +
        scale_color_discrete(name = "Classe") +
        scale_x_continuous(limits = x_limits) + # Ajustement des limites des axes x
        scale_y_continuous(limits = y_limits) + # Ajustement des limites des axes y
        theme_classic() + # Utilisation d'un thème classique pour assurer la visibilité des axes
        theme(
            axis.line = element_line(color = "black"), # Assurez-vous que les lignes des axes sont visibles
            axis.text = element_text(color = "black"), # Couleur du texte des axes
            axis.title = element_text(color = "black"), # Couleur du titre des axes
            plot.background = element_rect(fill = "white"), # Fond blanc
            panel.background = element_rect(fill = "white"), # Fond blanc pour le panel
            panel.grid.major = element_line(color = "grey", linetype = "dotted"), # Lignes de grille majeures
            panel.grid.minor = element_line(color = "lightgrey", linetype = "dotted") # Lignes de grille mineures
        )

    # Utilisation de ggsave pour sauvegarder l'image
    ggsave("./plots/clusters.png", image, width = 8, height = 6, dpi = 300)
})

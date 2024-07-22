png("chemin/vers/le/fichier/heatmap.png", width = 800, height = 600)
mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
couleurs <- colorRampPalette(c("blue", "white", "red"))(n = 299)
heatmap(matrice, Colv = NA, Rowv = NA, col = couleurs, scale = "none")

library(writexl)
library(readxl)

# Fonction pour créer une heatmap avec une palette de couleurs personnalisée
create_heatmap <- function(file_path, output_png, mat) {
    mat <- as.matrix(mat)
    png(output_png, width = 800, height = 600)
    # mat <- as.matrix(silent_run(read_excel, file_path))
    # write_xlsx(as.data.frame(mat), paste0(file_path, "_reform_oui.xlsx"))
    mat <- mat[nrow(mat):1, ] # Inverser l'ordre des lignes OUILLLLLLLE

    # Définir les couleurs et les breaks
    couleurs <- colorRampPalette(c("blue", "white", "red"))(n = 299)
    max_val <- max(abs(mat), na.rm = TRUE)
    breaks <- seq(-max_val, max_val, length.out = 300)

    heatmap(mat, Colv = NA, Rowv = NA, col = couleurs, breaks = breaks, scale = "none")
    dev.off()
}

# Appeler la fonction pour chaque heatmap
# create_heatmap("../data/beta_picto/big_picto_reconstructedlogistic_multibloc.xlsx", "../data/beta_picto/heatmap_multi.png")
# create_heatmap("../data/beta_picto/big_picto_reconstructedlogistique_simple.xlsx", "../data/beta_picto/heatmap_simple.png")
# create_heatmap("../data/beta_picto/big_picto.xlsx", "../data/beta_picto/heatmap_true.png")

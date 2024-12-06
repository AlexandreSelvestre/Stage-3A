# library(pheatmap)
# test <- matrix(rnorm(200), 20, 10)
# test[1:10, seq(1, 10, 2)] <- test[1:10, seq(1, 10, 2)] + 3
# test[11:20, seq(2, 10, 2)] <- test[11:20, seq(2, 10, 2)] + 2
# test[15:20, seq(2, 10, 2)] <- test[15:20, seq(2, 10, 2)] + 4
# colnames(test) <- paste("Test", 1:10, sep = "")
# rownames(test) <- paste("Gene", 1:20, sep = "")
# png("heatmap_test.png", width = 600, height = 1000)
# annotation_col <- data.frame(
#     CellType = factor(rep(c("CT1", "CT2"), 5)),
#     Time = 1:5
# )
# rownames(annotation_col) <- paste("Test", 1:10, sep = "")
# pheatmap(test, annotation_col = annotation_col)
# dev.off()




# df <- data.frame(a = c(1, 2), b = c(1, NA))
# out <- c(1, 2)
# out <- sapply(out, function(i) {
#     if (any(is.na(df[i, ]))) {
#         return(-1)
#     } else {
#         return(i)
#     }
# })
# print(out[out != -1])




# Exemple de matrice binaire
# binary_matrix <- matrix(c(0, 1, 0, 1, 1, 0, 1, 0), nrow = 4, ncol = 2)

# # Conversion de la matrice binaire en matrice de facteurs
# factor_matrix <- data.frame(lapply(as.data.frame(binary_matrix), function(col) {
#     factor(col, levels = c(0, 1))
# }))
# print(factor_matrix)
# # Vérification que chaque colonne est un facteur
# print(is.factor(factor_matrix[, 1]))
# print(is.factor(factor_matrix[, 2]))

li <- list(a = c(0, 1), b = c(2, 3))
# print(Reduce("c", li))
# print(unlist(li))
li$c <- c(0, 1)
print(li)

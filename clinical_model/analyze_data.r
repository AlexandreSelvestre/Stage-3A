library(factoextra)

levels_name <- paste0(as.character(unique(df$Tumeur)), collapse = "", sep = "_")
# print(df[, setdiff(colnames(df), c("Tumeur"))])
fit <- prcomp(df[, setdiff(colnames(df), c("Tumeur"))], scale = TRUE)
plot <- fviz_pca_ind(fit,
    geom.ind = "point", # show points only (nbut not "text")
    habillage = as.factor(df$Tumeur), # color by groups
    palette = c("#00AFBB", "#E7B800"),
    addEllipses = TRUE, # Concentration ellipses
    legend.title = "Tumeur",
    repel = TRUE, # Avoid text overlapping,
) + theme_bw()

print(paste0("plots/", levels_name, "pca_elipses.png"))
ggsave(paste0("plots/", levels_name, "pca_elipses.png"), plot)
plot <- fviz_pca_biplot(fit, habillage = as.factor(df$Tumeur)) + theme_bw()
ggsave(paste0("plots/", levels_name, "pca_decomp.png"), plot)
library(pheatmap)
png(paste0("plots/", levels_name, "heatmap.png"), width = 600, height = 1000)
annotation <- data.frame(Tumeur = ifelse(df$Tumeur == unique(df$Tumeur)[1], 0, 1))
# print(annotation)
rownames(annotation) <- rownames(df)
pheatmap(df[, setdiff(colnames(df), c("Tumeur"))], annotation_row = annotation) # scale = "column"
dev.off()

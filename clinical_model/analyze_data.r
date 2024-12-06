library(factoextra)

levels_name <- paste0(as.character(unique(df$Tumeur)), collapse = "", sep = "_")
# print(df[, setdiff(colnames(df), c("Tumeur"))])
if (kill_tumor != "") {
    rownames(df) <- df$patient_num
    fit <- prcomp(df[, setdiff(colnames(df), c("Tumeur", "patient_num"))], scale = TRUE)
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
    pheatmap(df[, setdiff(colnames(df), c("Tumeur", "patient_num"))], annotation_row = annotation) # scale = "column"
    dev.off()
}


li_lda <- MASS::lda(Tumeur ~ ., data = df[, setdiff(colnames(df), c("patient_num"))])
# Extraire les scores des individus sur les deux premières composantes discriminantes
lda_scores <- predict(li_lda)$x

# Créer un data frame avec les scores, les classes des individus et les numéros de patient
lda_df <- data.frame(lda_scores, Tumeur = df$Tumeur, patient_num = df$patient_num)
if (ncol(lda_scores) == 1) {
    lda_df$LD2 <- lda_df$LD1
}
plot <- ggplot(lda_df, aes(x = LD1, y = LD2, color = Tumeur)) +
    geom_point(size = 3) +
    geom_text(aes(label = patient_num), vjust = -1, hjust = 1) +
    labs(title = "LDA Plot", x = "LD1", y = "LD2") +
    theme_bw() +
    theme(legend.title = element_text(size = 10), legend.text = element_text(size = 8))

ggsave(paste0("plots/", levels_name, "lda_plot.png"), plot)

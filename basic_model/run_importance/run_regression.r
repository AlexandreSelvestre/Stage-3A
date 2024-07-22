run_imp_intra <- function(inference, imp_li, performance, li_confus, is_null, ite) {
    df_imp <- get_df_imp(inference)
    imp_li <- append(imp_li, list(df_imp))
    test_result <- inference@df_measures[["AUC_test"]]
    val_result <- inference@df_measures[["AUC_val"]]
    F1_macro_result <- inference@df_measures[["F1"]]
    Acc_result <- inference@df_measures[["Acc"]]

    performance$AUC_test <- append(performance$AUC_test, test_result)
    performance$AUC_val <- append(performance$AUC_val, val_result)
    performance$F1_macro <- append(performance$F1_macro, F1_macro_result)
    performance$Acc <- append(performance$Acc, Acc_result)
    li_confus <- append(li_confus, list(inference@confus_mat$table))
    confus_mat <- inference@confus_mat$table
    F1_CCK <- 2 * confus_mat[1, 1] / (sum(confus_mat[1, ]) + sum(confus_mat[, 1]))
    F1_CHC <- 2 * confus_mat[2, 2] / (sum(confus_mat[2, ]) + sum(confus_mat[, 2]))
    performance$F1_CCK <- append(performance$F1_CCK, F1_CCK)
    performance$F1_CHC <- append(performance$F1_CHC, F1_CHC)

    print(paste("actuelle moyenne AUC test", mean(performance$AUC_test), "ite:", ite))
    print(paste("actuelle moyenne AUC val", mean(performance$AUC_val), "ite:", ite))
    print(paste("actuelle moyenne Accuracy test", mean(performance$Acc), "ite:", ite))
    print(paste("actuelle somme des confusion matrix", "ite:", ite, ":"))
    print(paste("actuelle moyenne du macro F1 test", mean(performance$F1_macro), "ite:", ite))
    print(paste("actuelle moyenne du F1", colnames(confus_mat)[1], "test", mean(performance$F1_CCK), "ite:", ite))
    print(paste("actuelle moyenne du F1", colnames(confus_mat)[2], "test", mean(performance$F1_CHC), "ite:", ite))
    print(Reduce("+", li_confus))

    li_return <- list(imp_li = imp_li, performance = performance, li_confus = li_confus, is_null = is_null)
    return(li_return)
}



run_imp_extra <- function(imp_li, performance, li_confus, is_null, n_samples, path_plot, inference) {
    imp_sum <- imp_li[[1]]
    if (length(imp_li) > 1) {
        for (i in 2:length(imp_li)) {
            imp_sum$Overall <- imp_sum$Overall + imp_li[[i]]$Overall
        }
    }
    non_null <- data.frame(Variable = imp_li[[1]]$Variable, Overall = 0)
    for (df in imp_li) {
        vec_beta <- df$Overall
        vec_non_null <- vec_beta > 0.0001
        non_null$Overall <- non_null$Overall + vec_non_null
    }
    imp_average <- copy(imp_sum)
    imp_average$Overall <- imp_average$Overall / n_samples
    # imp_average <- subset(imp_average, Overall > 0.0001)
    imp_average$Overall <- imp_average$Overall / sum(imp_average$Overall)
    imp_average$Percentage <- 100 * non_null$Overall / n_samples

    ending_name <- "beta_value"
    ##### Utiliser l'inference pour les index et refaire les plots en général
    plot_global(imp_average, path_plot, ending_name)

    # ending_name <- "non_null"
    # plot_global(non_null, path_plot, ending_name, renorm = FALSE)
    performance_long <- melt_mine(as.data.frame(performance)[, setdiff(names(performance), "AUC_val")])
    box_plots_stats <- ggplot(performance_long, aes(x = variable, y = value)) +
        stat_summary(fun = median, geom = "point", shape = 20, size = 3, color = "red") +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    ggsave(paste0(path_plot, "/global_plot_stats.png"), box_plots_stats)

    mat_sum <- Reduce("+", li_confus)
    print("Voilà la matrice de confusion sommée")
    print(mat_sum)
    new_mat <- matrix(0, nrow = 2, ncol = 2)
    new_mat[1, 1] <- mat_sum[1, 1] / sum(mat_sum[, 1])
    new_mat[1, 2] <- mat_sum[1, 2] / sum(mat_sum[, 2])
    new_mat[2, 1] <- mat_sum[2, 1] / sum(mat_sum[, 1])
    new_mat[2, 2] <- mat_sum[2, 2] / sum(mat_sum[, 2])
    colnames(new_mat) <- colnames(mat_sum)
    rownames(new_mat) <- rownames(mat_sum)
    print("voilà la matrice de confusion en pourcentage: ... % de CCK ont été bien classés vs ... % mal classés")
    print(new_mat)

    final_accuracy_macro <- mean(c(mat_sum[1, 1] / sum(mat_sum[1, ]), mat_sum[2, 2] / sum(mat_sum[2, ])))
    print(paste("La balanced accuracy sur l'échantillon total vaut:", final_accuracy_macro))
    f_1_1 <- 2 * mat_sum[1, 1] / (sum(mat_sum[1, ]) + sum(mat_sum[, 1]))
    cat("le f1 score", colnames(mat_sum)[1], "sur l'échantillon entier vaut:", f_1_1, "\n")
    f_1_2 <- 2 * mat_sum[2, 2] / (sum(mat_sum[2, ]) + sum(mat_sum[, 2]))
    cat("le f1 score", colnames(mat_sum)[2], "sur l'échantillon entier vaut:", f_1_2, "\n")
    final_f1_macro <- mean(c(f_1_1, f_1_2))
    print(paste("Le f1 macro sur l'échantillon total vaut:", final_f1_macro))
}

library(readxl)
library(writexl)
library(glue)
library(data.table)
library(caret)
library(randomForest)
library(ggplot2)
library(DMwR)
library(themis)




setMethod("train_method", "apply_model", function(object) {
    weight_vec <- unlist(object@weights)
    # tuneGrid <- expand.grid(wcl = 1, mtry = object@mtry)
    tuneGrid <- expand.grid(mtry = object@mtry)
    if (object@do_PCA) {
        object@model <- caret::train(
            y = object@train_cols$classe_name, x = object@train_cols[, object@col_x],
            method = "rf", trControl = object@cv, n_trees = object@n_trees, na.action = na.roughfix, metric = "Weighted Accuracy",
            tuneLength = 8, tuneGrid = tuneGrid, preProcess = "pca", classwt = weight_vec
        )
    } else {
        object@model <- caret::train(
            y = object@train_cols$classe_name, x = object@train_cols[, object@col_x],
            method = "rf", trControl = object@cv, n_trees = object@n_trees, na.action = na.roughfix, metric = "AUC",
            tuneLength = 8, tuneGrid = tuneGrid, classwt = weight_vec
        )
    }
    return(object)
})


setMethod("get_results", "apply_model", function(object) {
    # print(paste(c("test_set puis train", dim(object@test_set), dim(object@train_set)), collapse = "x"))
    if (object@do_PCA) {
        test_explic_pca <- predict(object@model$preProcess, object@test_set[, object@col_x])
        train_explic_pca <- predict(object@model$preProcess, object@train_cols[, object@col_x])
        object@predictions <- predict(object@model$finalModel, newdata = test_explic_pca)
        object@predictions_proba <- predict(object@model$finalModel, newdata = test_explic_pca, type = "prob")
        object@predictions_train_proba <- predict(object@model$finalModel, newdata = train_explic_pca, type = "prob")
    } else {
        object@predictions <- predict(object@model$finalModel, newdata = object@test_set[, object@col_x])
        object@predictions_proba <- predict(object@model$finalModel, newdata = object@test_set[, object@col_x], type = "prob")
        object@predictions_train_proba <- predict(object@model$finalModel, newdata = object@train_cols[, object@col_x], type = "prob")
    }

    lignes_erreur_true_false <- object@test_set$classe_name != object@predictions
    lignes_erreur_testing <- object@test_set[lignes_erreur_true_false, c("keys", "original_firstorder_10Percentile_PORT", "classe_name")]
    lignes_erreur_testing <- lignes_erreur_testing[lignes_erreur_testing$classe_name == "CCK", ]
    print("erreurs globales:")
    print(lignes_erreur_testing)
    lignes_reussite_true_false <- object@test_set$classe_name == object@predictions
    lignes_reussite_testing <- object@test_set[lignes_reussite_true_false, c("keys", "original_firstorder_10Percentile_PORT", "classe_name")]
    lignes_reussite_testing <- lignes_reussite_testing[lignes_reussite_testing$classe_name == "CCK", ]
    print("réussites globales:")
    print(lignes_reussite_testing)


    # print(object@predictions_proba)
    return(object)
})



setMethod("importance_method", "apply_model", function(object) {
    puiss <- importance(object@model$finalModel)
    puiss_df <- data.frame(Variable = row.names(puiss), Importance = puiss[, 1])
    puiss_df <- puiss_df[order(-puiss_df$Importance), ]
    name_df <- paste("variable_importance", object@id_term, sep = "_")
    object@li_df_var_imp[[name_df]] <- puiss_df

    importance_graph <- ggplot(puiss_df, aes(x = reorder(Variable, Importance), y = Importance)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance") +
        theme(axis.text.y = element_text(size = 2)) # Adjust the size as needed
    print(importance_graph)

    # ggsave("plots/random_forest/importance_plot_split.png", importance_graph, width = 10, height = 20)
    ggsave(paste0("plots/random_forest/importance", "_", object@id_term, ".png"), importance_graph, width = 10, height = 20)


    puiss_df$Group <- substr(puiss_df$Variable, nchar(puiss_df$Variable) - 3, nchar(puiss_df$Variable))
    puiss_df_grouped <- aggregate(Importance ~ Group, data = puiss_df, FUN = mean)
    puiss_df_grouped <- puiss_df_grouped[order(-puiss_df_grouped$Importance), ]
    name_df <- paste("variable_importance_grouped_big", object@id_term, sep = "_")
    object@li_df_var_imp[[name_df]] <- puiss_df_grouped


    importance_graph_grouped <- ggplot(puiss_df_grouped, aes(x = reorder(Group, Importance), y = Importance)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance")
    print(importance_graph_grouped)

    # ggsave("plots/random_forest/importance_plot_grouped_big.png", importance_graph_grouped)
    ggsave(paste0("plots/random_forest/importance_big_groups", "_", object@id_term, ".png"), importance_graph_grouped)

    puiss_df$Group_small <- substr(puiss_df$Variable, 0, nchar(puiss_df$Variable) - 5)
    puiss_df_grouped_small <- aggregate(Importance ~ Group_small, data = puiss_df, FUN = mean)
    puiss_df_grouped_small <- puiss_df_grouped_small[order(-puiss_df_grouped_small$Importance), ]
    name_df <- paste("variable_importance_grouped_small", object@id_term, sep = "_")
    object@li_df_var_imp[[name_df]] <- puiss_df_grouped_small

    importance_graph_grouped <- ggplot(puiss_df_grouped_small, aes(x = reorder(Group_small, Importance), y = Importance)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Variable Importance") +
        theme(axis.text.y = element_text(size = 3)) # Adjust the size as needed
    print(importance_graph_grouped)

    # ggsave("plots/random_forest/importance_plot_grouped_small.png", importance_graph_grouped, width = 10, height = 20)
    ggsave(paste0("plots/random_forest/importance_small_groups", "_", object@id_term, ".png"), importance_graph_grouped, width = 10, height = 20)

    df_cv <- object@model$resample
    df_cv <- df_cv[, setdiff(names(df_cv), "Resample")]
    df_long <- melt(df_cv)
    object@li_box_plots[[object@id_term]] <- df_long
    box_plots_stats <- ggplot(df_long, aes(x = variable, y = value)) +
        geom_boxplot() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for readability
    # ggsave("plots/random_forest/box_plots_stats.png", box_plots_stats)
    ggsave(paste0("plots/random_forest/box_plots_stats", "_", object@id_term, ".png"), box_plots_stats)
    return(object)
})

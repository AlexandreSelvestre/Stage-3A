sepa_data_set <- function(data_used, summary, training_index, rep, show = TRUE, calc_probs = TRUE, sampling = NULL, k_smote = 5, search = "random", folds, y_tot) {
    # Pour bien se souvenir de comment gérer smote dans caret, même si on ne l'utilise plus
    # smotest <- list(
    #     name = "SMOTE with more neighbors!",
    #     func = function(x, y) {
    #         library(DMwR)
    #         dat <- if (is.data.frame(x)) x else as.data.frame(x)
    #         dat$.y <- y

    #         # Calculer le pourcentage perc.over
    #         class_majoritaire <- names(which.max(table(y)))
    #         class_minoritaire <- names(which.min(table(y)))
    #         perc.over <- (sum(y == class_majoritaire) / sum(y == class_minoritaire)) * 100
    #         perc.under <- 100 + 100 / floor(perc.over / 100) + 1

    #         dat <- SMOTE(.y ~ ., data = dat, k = k_smote, perc.over = perc.over)
    #         list(
    #             x = dat[, !grepl(".y", colnames(dat), fixed = TRUE)],
    #             y = dat$.y
    #         )
    #     },
    #     first = TRUE
    # )

    # if (sampling == "smote") {
    #     sampling <- smotest
    # }
    train_set <- data_used[training_index, ]
    test_set <- data_used[-training_index, ]
    y_train <- y_tot[training_index]
    y_test <- y_tot[-training_index]
    rownames(train_set) <- NULL
    cv <- trainControl(
        method = "repeatedcv", verboseIter = show, classProbs = calc_probs, repeats = rep,
        summaryFunction = custom_stats, search = search, sampling = NULL, returnResamp = "final",
        savePredictions = "final",
        index = folds, timingSamps = length(folds)
    )
    # set.seed(1)
    return(list(train_set = train_set, test_set = test_set, cv = cv, index = training_index, y_train = y_train, y_test = y_test))
}



apply_smote <- function(x, y, k = 3) {
    library(DMwR)
    dat <- if (is.data.frame(x)) x else as.data.frame(x)
    dat$.y <- y
    classe_majoritaire <- names(which.max(table(y)))
    classe_minoritaire <- names(which.min(table(y)))
    perc.over <- (sum(y == classe_majoritaire) / sum(y == classe_minoritaire)) * 100
    perc.under <- 100 + 100 / floor(perc.over / 100) + 1
    dat <- SMOTE(.y ~ ., data = dat, k = k, perc.over = perc.over, perc.under = perc.under)
    # print(sum(dat$.y == classe_majoritaire) / sum(dat$.y == classe_minoritaire))
    li <- list(x = dat[, setdiff(names(dat), ".y")], y = dat$.y)
    return(li)
}


apply_boot <- function(x, y) {
    x <- if (is.data.frame(x)) x else as.data.frame(x)
    classe_majoritaire <- names(which.max(table(y)))
    classe_minoritaire <- names(which.min(table(y)))
    difference <- sum(y == classe_majoritaire) - sum(y == classe_minoritaire)
    for (i in 1:difference) {
        index_minoraires <- which(y == classe_minoritaire)
        index_to_add <- sample(index_minoraires, 1)
        y <- c(y, y[index_to_add])
        x <- rbind(x, x[index_to_add, ])
    }

    # print(sum(y == classe_majoritaire) / sum(y == classe_minoritaire))

    return(list(x = x, y = y))
}



find_sigma_mu_index_mode <- function(x, index, is_binary) {
    # index_variable nécessaire: c'est index
    x <- as.matrix(x)
    # l'indice -1 est attribué aux variables tabulaires dans index
    # is_binary contient des 1 au niveau des variables à ne pas normaliser (binaires), même taille que index même si seule premiere occurence (en mode) regardée
    index_non_neg <- index[index >= 0]
    M <- length(index[index == -1])
    K <- sum(index_non_neg == index_non_neg[1])
    J <- length(unique(index_non_neg))
    df_mu <- as.data.frame(matrix(0, nrow = nrow(x), ncol = ncol(x)))
    df_sigma <- as.data.frame(matrix(1, nrow = nrow(x), ncol = ncol(x)))
    li_dico <- lapply(1:J, function(x) c())
    li_dico$tab <- c()
    # print(li_dico)
    for (a in 1:ncol(x)) {
        if (index[a] > -0.5) {
            # print(index[a])
            li_dico[[index[a]]] <- c(li_dico[[index[a]]], a)
        } else {
            li_dico$tab <- c(li_dico$tab, a)
        }
    }

    for (j in 1:J) {
        if (!is_binary[j]) {
            full_col <- c()
            for (col_num in li_dico[[j]]) {
                col <- x[, col_num]
                full_col <- c(full_col, col)
            }
            mu <- mean(full_col, na.rm = TRUE)
            sigma <- sd(full_col, na.rm = TRUE)
            if (sigma == 0) {
                sigma <- 1
                print(paste("Attention, la variable", colnames(x)[col_num], "a une variance nulle"))
            }
            for (col_num in li_dico[[j]]) {
                df_mu[, col_num] <- mu
                df_sigma[, col_num] <- sigma
            }
        }
    }
    for (col_num in li_dico$tab) {
        if (!is_binary[col_num]) {
            df_mu[, col_num] <- mean(x[, col_num], na.rm = TRUE)
            sigma <- sd(x[, col_num], na.rm = TRUE)
            if (sigma == 0) {
                sigma <- 1
                print(paste("Attention, la variable", colnames(x)[col_num], "a une variance nulle"))
            }
            df_sigma[, col_num] <- sigma
        }
    }
    x <- as.data.frame(x)
    colnames(df_mu) <- colnames(x)
    colnames(df_sigma) <- colnames(x)
    return(list(mu = df_mu, sigma = df_sigma))
}

renormalize_in_model_fit_index_mode <- function(x, index, is_binary = NULL) {
    # index_variable nécessaire
    if (is.null(is_binary)) {
        is_binary <- rep(FALSE, ncol(x))
    }
    x <- as.data.frame(x)
    li <- find_sigma_mu_index_mode(x, index, is_binary)
    df_mu <- li$mu
    df_sigma <- li$sigma
    new_x <- data.table::copy(x)
    # log_file <- paste0("./log/fit_log_", Sys.getpid(), ".txt")
    # cat("Début de fit_simple\n", file = log_file, append = TRUE)
    # cat("Valeur de beta: ", dim(x), dim(df_mu), dim(df_sigma), "\n", file = log_file, append = TRUE)
    new_x <- (new_x - df_mu) / df_sigma
    return(list(new_x = new_x, df_mu = df_mu, df_sigma = df_sigma))
}

renormalize_in_model_pred_index_mode <- function(newdata, df_mu, df_sigma) {
    newdata <- as.data.frame(newdata)
    n_lines <- nrow(newdata)
    newdata <- (newdata - df_mu[1:n_lines, ]) / df_sigma[1:n_lines, ]
    newdata <- as.matrix(newdata)
    return(newdata)
}

reorder_in_modes <- function(x, index_mode, index_variable, index_bloc, name_mode, name_variable, name_bloc, is_binary) {
    # Dans certains modèles qui suivent, on suppose que les numéros de variables sont croissant dans l'ordre des blocs (les variables de numéro maximal sont associées aux blocs de numéoro maximal). On évitera donc le cas contraire (même si cette fonction particulière s'en sortirait bien... )
    # L'application de cette fonction assure que les indices de sortie sont tous bien jolis et contigus... Ca sauve certains modèles un peu mal construits. Et on arrange tout sur le modèle du data used des radiomiques (par mode mais pas par bloc). En fait le nouvel arrangement de x se fout des blocs (de l'importance de l'ordre des variables si on veut faire du multibloc ensuite: ce dernier modèle ayant absolument besoin de cet ordonnancement.)
    different_variables <- order(unique(index_variable[index_variable > -0.5]))
    different_blocs <- order(unique(index_bloc[index_bloc > -0.5]))
    x <- as.data.frame(x)
    li_modes <- list()
    count_tab <- 1
    new_index_mode <- rep(-1, length(index_mode))
    new_index_variable <- rep(-1, length(index_variable))
    new_index_bloc <- rep(-1, length(index_bloc))
    new_is_binary <- rep(FALSE, length(is_binary))
    new_name_mode <- rep("", length(index_mode))
    new_name_variable <- rep("", length(index_variable))
    new_name_bloc <- rep("", length(index_bloc))
    K <- length(unique(index_mode[index_mode > -0.5]))
    J <- length(unique(index_variable[index_variable > -0.5]))
    for (i in 1:length(index_mode)) {
        if (!as.character(index_mode[i]) %in% names(li_modes)) {
            li_modes[[as.character(index_mode[i])]] <- as.data.frame(matrix(0, ncol = length(index_mode[index_mode == index_mode[i]]), nrow = nrow(x)))
        }
        num_var <- which(different_variables == index_variable[i])
        if (index_mode[i] > -0.5) {
            li_modes[[as.character(index_mode[i])]][, num_var] <- x[, i]
            colnames(li_modes[[as.character(index_mode[i])]])[num_var] <- colnames(x)[i]
            # cat(index_mode[i], colnames(x)[i], "\n")
            # print(colnames(li_modes[[as.character(index_mode[i])]])[count_tab])
        } else {
            li_modes[[as.character(index_mode[i])]][, count_tab] <- x[, i]
            colnames(li_modes[[as.character(index_mode[i])]])[count_tab] <- colnames(x)[i]
            count_tab <- count_tab + 1
        }
    }
    li_modes <- li_modes[order(as.numeric(names(li_modes)))]
    new_x <- as.data.frame(matrix(0, nrow = nrow(x), ncol = 0))
    for (df in li_modes[which(as.numeric(names(li_modes)) > -0.5)]) {
        new_x <- cbind(new_x, df)
    }
    # print(colnames(new_x))
    new_x <- cbind(new_x, li_modes[["-1"]])
    for (k in 1:K) {
        for (j in (1:J)) {
            new_index_mode[(k - 1) * J + j] <- k
            new_index_variable[(k - 1) * J + j] <- j
            where_mode <- which(index_mode == order(unique(index_mode[index_mode > -0.5]))[k])
            where_variable <- which(index_variable == different_variables[j])
            previous_index <- intersect(where_mode, where_variable)
            if (length(previous_index) > 1) {
                print(previous_index)
                stop("Indexes non identifiables")
            }
            new_index_bloc[(k - 1) * J + j] <- which(different_blocs == index_bloc[previous_index])
            new_is_binary[(k - 1) * J + j] <- is_binary[previous_index]
            new_name_mode[(k - 1) * J + j] <- name_mode[previous_index]
            new_name_variable[(k - 1) * J + j] <- name_variable[previous_index]
            new_name_bloc[(k - 1) * J + j] <- name_bloc[previous_index]
            # cat(previous_index, j, index_variable[previous_index], "\n")
        }
    }
    return(list(x = new_x, index_mode = new_index_mode, index_variable = new_index_variable, index_bloc = new_index_bloc, is_binary = new_is_binary, name_mode = new_name_mode, name_variable = new_name_variable, name_bloc = new_name_bloc))
}


plot_global <- function(imp_average, path_plot, ending_name, inference) {
    imp_plot <- subset(imp_average, Overall > 0.0001 | Percentage > 0.01)
    image <- ggplot2::ggplot(imp_plot, aes(x = reorder(Variable, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(Percentage, 1), "%"), y = 0), hjust = -0.5, color = "green") +
        geom_point(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, color = "Percentage of non zero beta coefficient for this variable"), size = 5) + # point invisible pour légende
        scale_color_manual(name = "", values = c("Percentage of non zero beta coefficient for this variable" = "green")) +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Relative variable importance") +
        theme(
            axis.text.y = element_text(size = 5), legend.position = "bottom",
            legend.text = element_text(size = 12)
        )
    ggsave(paste0(path_plot, "/global_all", "_", ending_name, ".png"), image, width = 10, height = 20)

    variable_importance <- imp_average
    variable_importance$bloc <- inference@name_bloc
    variable_importance_bloc <- aggregate_prop(Overall ~ bloc, data = variable_importance)

    image <- ggplot(variable_importance_bloc, aes(x = reorder(bloc, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(Percentage, 1), "%"), y = 0), hjust = -0.5, color = "green") +
        geom_point(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, color = "Percentage of non zero beta coefficient for this variable"), size = 5) + # point invisible pour légende
        scale_color_manual(name = "", values = c("Percentage of non zero beta coefficient for this variable" = "green")) +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Grouped relative variable importance of each bloc") +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 12)
        )
    ggsave(paste0(path_plot, "/global_blocs", "_", ending_name, ".png"), image) # bloc par bloc

    variable_importance <- variable_importance[variable_importance$bloc != "clinical", ]
    variable_importance$Group <- substr(variable_importance$Variable, nchar(variable_importance$Variable) - 3, nchar(variable_importance$Variable))
    variable_importance$small_Group <- substr(variable_importance$Variable, 0, nchar(variable_importance$Variable) - 5)
    if (renorm) {
        variable_importance_grouped <- aggregate_prop(Overall ~ Group, data = variable_importance, FUN = mean)
        variable_importance_small_grouped <- aggregate_prop(Overall ~ small_Group, data = variable_importance, FUN = mean)
    } else {
        variable_importance_grouped <- aggregate(Overall ~ Group, data = variable_importance, FUN = mean)
        variable_importance_small_grouped <- aggregate(Overall ~ small_Group, data = variable_importance, FUN = mean)
    }

    image <- ggplot(variable_importance_grouped, aes(x = reorder(Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(Percentage, 1), "%"), y = 0), hjust = -0.5, color = "green") +
        geom_point(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, color = "Percentage of non zero beta coefficient for this variable"), size = 5) + # point invisible pour légende
        scale_color_manual(name = "", values = c("Percentage of non zero beta coefficient for this variable" = "green")) +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Grouped relative variable importance of each time") +
        theme(
            legend.position = "bottom",
            legend.text = element_text(size = 12)
        )
    ggsave(paste0(path_plot, "/global_big_groups", "_", ending_name, ".png"), image) # temps par temps


    image <- ggplot(variable_importance_small_grouped, aes(x = reorder(small_Group, Overall), y = Overall)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = paste0(round(Percentage, 1), "%"), y = 0), hjust = -0.5, color = "green") +
        geom_point(data = data.frame(x = Inf, y = Inf), aes(x = x, y = y, color = "Percentage of non zero beta coefficient for this variable"), size = 5) + # point invisible pour légende
        scale_color_manual(name = "", values = c("Percentage of non zero beta coefficient for this variable" = "green")) +
        coord_flip() +
        theme_light() +
        xlab("Variable") +
        ylab("Importance") +
        ggtitle("Grouped relative variable importance of each quantity of interest") +
        theme(
            axis.text.y = element_text(size = 6), legend.position = "bottom",
            legend.text = element_text(size = 12)
        )

    ggsave(paste0(path_plot, "/global_small_groups", "_", ending_name, ".png"), image, width = 10, height = 20) # variable par variable
}
# library(readxl)
# library(writexl)
# data_used <- as.data.frame(read.csv("..\\data\\data_used.csv"))
# index_mode <- readRDS("..\\data\\RDS\\index_mode.rds")
# index_variable <- readRDS("..\\data\\RDS\\index_variable.rds")
# index_bloc <- readRDS("..\\data\\RDS\\index_bloc.rds")
# is_binary <- readRDS("..\\data\\RDS\\is_binary.rds")
# li <- reorder_in_modes(data_used[, 3:ncol(data_used)], index_mode, index_variable, index_bloc, is_binary)
# df <- li$x
# new_index_bloc <- li$index_bloc
# new_is_binary <- li$is_binary
# new_index_variable <- li$index_variable
# print("les binaires:")
# print(new_is_binary)
# print("les blocs:")
# print(new_index_bloc)
# print("les variables")
# print(new_index_variable)
# write_xlsx(df, "..\\data\\no.xlsx")

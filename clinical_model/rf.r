li_rf <- list()

li_rf$library <- c("randomForest", "missForest")

li_rf$type <- "Classification"

li_rf$parameters <- data.frame(parameter = c("mtry"), class = c("numeric"), label = c(
    "le nombre de variables à tirer"
))

create_grid_rf <- function(x, y, len = NULL, search = "grid") {
    if (search == "grid") {
        mtry <- seq(1, ncol(x), length.out = len)
    } else {
        mtry <- sample(1:ncol(x), len, replace = TRUE)
    }
    data_frame_grid <- expand.grid(mtry = mtry)
    return(data_frame_grid)
}

li_rf$grid <- create_grid_rf

li_rf$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ntree = 5000) {
    x <- data.frame(lapply(as.data.frame(x), function(col) {
        return(factor(col, levels = c(0, 1)))
    }))
    if (any(is.na(x))) {
        capture.output({
            li <- missForest::missForest(xmis = cbind(data.frame(y = y), x), ntree = 1000, maxiter = 10, mtry = 1)
        })
        x <- li$ximp[2:ncol(li$ximp)]
        print(li$OOBerror)
    }
    vec_classes <- sapply(lev, function(classe) {
        return(sum(y == classe))
    })
    names(vec_classes) <- as.character(lev)
    # class_min <- lev[which.min(vec_classes)]
    class_max <- lev[which.max(vec_classes)]
    for (classe in setdiff(lev, class_max)) {
        difference <- vec_classes[class_max] - vec_classes[classe]
        if (difference > 0) {
            additional <- x[sample(which(y == classe), difference, replace = TRUE), ]
            x <- rbind(x, additional)
            y <- append(y, factor(rep(classe, difference), levels = lev))
        }
    }
    mtry <- param$mtry
    res <- randomForest::randomForest(
        x = x, y = y, mtry = mtry, ntree = ntree, importance = TRUE
    )
    return(res)
}

li_rf$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- data.frame(lapply(as.data.frame(newdata), function(col) {
        return(factor(col, levels = c(0, 1)))
    }))
    return(predict(modelFit, newdata))
}

li_rf$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL) {
    newdata <- data.frame(lapply(as.data.frame(newdata), function(col) {
        return(factor(col, levels = c(0, 1)))
    }))
    return(predict(modelFit, newdata, type = "prob"))
}

li_rf$loop <- NULL




rf_importance_intra <- function(results, li_imp_rf, show_model) {
    if (show_model) {
        print(paste("meilleur modèle: rf avec mtry =", results$bestTune$mtry))
    }
    li_imp_rf[[length(li_imp_rf) + 1]] <- as.vector(results$finalModel$importance)
    return(li_imp_rf)
}

rf_importance_extra <- function(li_imp_rf, df) {
    vec_global_sum <- Reduce("+", li_imp_rf)
    vec_global_sum <- vec_global_sum / nrow(df)
    names(vec_global_sum) <- colnames(df[setdiff(colnames(df), c("Tumeur"))])

    df_res <- data.frame(Variable = gsub("[_.]", " ", names(vec_global_sum)), Importance = as.numeric(vec_global_sum))
    df_res$Variable <- factor(df_res$Variable, levels = df_res$Variable[order(df_res$Importance)])
    # Créer le barplot et sauvegarder en PNG

    plot <- ggplot(df_res, aes(x = Variable, y = Importance)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        theme(axis.text.y = element_text(size = 12)) +
        labs(title = "Importance des Variables", x = "", y = "")
    ggsave("plots/barplot_rf.png", plot)
}

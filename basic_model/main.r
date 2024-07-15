execute <- function(config, config_run, id_term, seed_cv) {
    if (config_run$do_extract) {
        extract_all(config)
    }
    config$id_term <- id_term
    data_used <- as.data.frame(read.csv("..\\data\\data_used.csv", check.names = FALSE))
    print(paste(c("dimensions data used", dim(data_used)))) 
    info_cols <- readRDS(file = "../data/RDS/info_cols.rds")
    config$data_used <- data_used
    config$info_cols <- info_cols
    slot_names <- slotNames(getClass("apply_model"))
    config_names <- names(config)
    arguments <- config[intersect(slot_names, config_names)]
    inference <- do.call("new", args = c("apply_model", arguments))

    inference <- init(inference)

    # Gérer la seed et la séparation train/ test + la séparation en folds
    seed_model <- .Random.seed
    .Random.seed <<- seed_cv
    training_index <- as.vector(createDataPartition(y = data_used[[info_cols$explained_col]], p = config$p, list = FALSE))
    training_set <- data_used[training_index, ]
    folds <- createMultiFolds(y = data_used[[info_cols$explained_col]][training_index], k = config$k, times = config$rep)
    global_number_folds <- list()
    for (i in 1:length(folds)) {
        fold <- folds[[i]]
        global_number_folds[[i]] <- setdiff(training_index, training_index[fold])
    }
    seed_cv <- .Random.seed
    .Random.seed <<- seed_model

    inference <- split_met(inference, training_index, folds)

    inference <- train_method(inference)

    print("End of learning phase")

    inference <- get_results(inference)

    print("End of prediction phase")

    inference <- analyse_results(inference)
    return(list(inference = inference, seed_cv = seed_cv))
}

# 18,19,43,77
# 20,22,41,77
# 8, 29, 77, 84
# 2, 17, 77, 84
# 2, 17, 22, 29

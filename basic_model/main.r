execute <- function(config, config_run, id_term, seed_cv, seed_partition, sysname) {
    if (config_run$do_extract & as.numeric(id_term) > 1) {
        extract_all(config, sysname)
    }
    config$id_term <- id_term
    path_data <- config$path_data
    if (sysname == "Linux") {
        path <- paste0(path_data, "/data_used.csv")
    } else {
        path <- "..\\data\\data_used.csv"
    }
    use_li_index_modes <- FALSE
    if (config_run$big_data) {
        use_li_index_modes <- TRUE
    }
    data_used <- as.data.frame(read.csv(path, check.names = FALSE))
    print(paste(c("dimensions data used", dim(data_used))))
    info_cols <- readRDS(file = paste0(path_data, "/RDS/info_cols.rds"))
    config$data_used <- data_used
    config$info_cols <- info_cols
    config$use_li_index_modes <- use_li_index_modes
    name_model <- config$name_model
    slot_names <- slotNames(getClass(name_model))
    config_names <- names(config)
    arguments <- config[intersect(slot_names, config_names)]
    inference <- do.call("new", args = c(name_model, arguments))
    inference <- init(inference)
    if (config$analyse_data$do & as.numeric(id_term) == 1) {
        analyse_data(inference)
        print("data analysed")
    }

    # Gérer la seed et la séparation train/ test + la séparation en folds
    seed_model <- .Random.seed
    .Random.seed <<- seed_partition
    if (!config_run$keep_partition | as.numeric(id_term) == 1) {
        training_index <- as.vector(createDataPartition(y = data_used[[info_cols$explained_col]], p = config$p, list = FALSE))
        saveRDS(training_index, file = paste0(path_data, "/RDS/training_index.rds"))
    } else {
        training_index <- readRDS(file = paste0(path_data, "/RDS/training_index.rds"))
    }
    seed_partition <- .Random.seed


    .Random.seed <<- seed_cv
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

    print("End of analysis phase")

    if (config_run$simulated_data) {
        source("./analyse_data/reform_plots.r")
        inference <- reform_beta(inference)
    }
    return(list(inference = inference, seed_cv = seed_cv, seed_partition = seed_partition))
}

# 18,19,43,77
# 20,22,41,77
# 8, 29, 77, 84
# 2, 17, 77, 84
# 2, 17, 22, 29

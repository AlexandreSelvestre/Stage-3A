# options(error = recover)
# ls(envir = .GlobalEnv): astuce!
config <- config::get(file = "config.yml", config = "my_config")
folder <- config$folder
full_name <- paste0(folder, "/conf.yml")
print(full_name)
sub_config <- config::get(file = full_name)
config <- modifyList(config, sub_config)
print(config$value)

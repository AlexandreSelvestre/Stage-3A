library(doParallel)
library(foreach)
library(parallel)
library(itertools)

threads <- 2
n <- 200

name_vec <- sapply(1:threads, function(x) paste0("worker_", x))
cl <- makePSOCKcluster(threads)
registerDoParallel(cl)
# Checker si on change l'iterator de place??
clusterExport(cl, "name_vec")

iterator <- isplitIndices(n, chunks = threads)
clusterApply(cl, name_vec, function(name) assign("worker_name", name, envir = .GlobalEnv)) # Global env indispensable car le local est détruit après chaque appel aux workers
res <- foreach(i = iterator, .combine = "c") %dopar% {
    cat("I do", sum(i), "and 'i' was", i, file = paste0(worker_name, ".txt"), append = TRUE)
    return(sum(i))
}

clusterEvalQ(cl, {
    name_file <- paste0(worker_name, "_archaique", ".txt")
    file.create(name_file)
    lines <- c(
        "first",
        "JUJU"
    )
    writeLines(lines, name_file)
})

stopCluster(cl)


# Le sauveur dans caret:
# # Debut:
# # Lister tous les fichiers .txt contenant "archaique" dans leur nom
# files_to_remove <- list.files(pattern = "archaique.*\\.txt$")

# # Supprimer les fichiers listés
# file.remove(files_to_remove)


### Dans le fit:

# name_file <- paste0(worker_name, "_archaique", ".txt")
# file.create(name_file)
# lines <- c(
#     "first",
#     "JUJU"
# )
# writeLines(lines, name_file)

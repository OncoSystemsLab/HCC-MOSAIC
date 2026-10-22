hcc_input_file <- function(name) {
    path <- Sys.getenv(name, unset = "")
    if (!nzchar(path)) stop(paste("Set the input-file environment variable:", name))
    if (!file.exists(path)) stop(paste("Input file not found for:", name))
    normalizePath(path, winslash = "/", mustWork = TRUE)
}

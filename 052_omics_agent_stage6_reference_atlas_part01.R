if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "10_reference_atlas", "37_accept_v12_rc4_release.R"))

qc <- accept_v12_rc4_release(PROJECT_DIR)

Sys.setenv(HCC_PROJECT_DIR = PROJECT_DIR)

local_lib <- Sys.getenv("HCC_R_LIBRARY", file.path(getwd(), ".library"))

if (dir.exists(local_lib)) {
    .libPaths(c(local_lib, .libPaths()))
}


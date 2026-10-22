if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

MOSAIC_DIR <- Sys.getenv("HCC_MOSAIC_ROOT", unset = "HCC_MOSAIC_PublicResource")

MOSAIC_DIR <- normalizePath(MOSAIC_DIR, winslash = "/", mustWork = TRUE)

source(file.path(MOSAIC_DIR, "tools", "cptac159_common.R"), encoding = "UTF-8")

hcc_cptac159_build_assets(MOSAIC_DIR)


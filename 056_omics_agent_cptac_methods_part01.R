if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

status <- system2(file.path(R.home("bin"), "Rscript"), c("--vanilla", shQuote(file.path(MOSAIC_DIR, "tools", 
    "update_cptac159.R")), "--step=validate"))

stopifnot(status == 0L)


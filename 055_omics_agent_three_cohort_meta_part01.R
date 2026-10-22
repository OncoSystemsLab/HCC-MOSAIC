if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

local({
    previous_wd <- setwd(MOSAIC_DIR)
    on.exit(setwd(previous_wd), add = TRUE)
    source("tools/build_three_cohort_meta.R", local = TRUE, encoding = "UTF-8")
})


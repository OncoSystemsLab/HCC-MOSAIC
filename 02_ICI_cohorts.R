# Prepare ICI expression inputs and response labels.
source("R/project.R", encoding = "UTF-8")
require_stages("01")
source(file.path(CFG$root, "R/03_ici.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
finish_stage("02", artifact("stage02_ici_pc_candidate_bundle.rds"))

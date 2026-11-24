# Align the HCC and ICI cohorts to a shared gene background.
source("R/project.R", encoding = "UTF-8")
require_stages("02")
source(file.path(CFG$root, "R/04_ssl.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
finish_stage("03", c(artifact("stage03_unlabeled_raw_bundle.rds"), artifact("stage03_ici_aligned_with_mask_bundle.rds")))

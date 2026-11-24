# Create the SSL split, preprocess expression, and audit samples.
source("R/project.R", encoding = "UTF-8")
require_stages("03")
source(file.path(CFG$root, "R/05_preprocess.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
source(file.path(CFG$root, "R/05_pca.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
finish_stage("04", artifact("stage04_SSL_prepared_candidate_bundle.rds"))

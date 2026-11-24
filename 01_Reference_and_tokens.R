# Build the GENCODE reference and audit the token panel.
source("R/project.R", encoding = "UTF-8")
require_stages("00")
source(file.path(CFG$root, "R/01_reference.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
source(file.path(CFG$root, "R/02_token_reference.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
finish_stage("01", c(artifact("stage01_reference_bundle.rds"), artifact("stage01_token_audit_bundle.rds")))

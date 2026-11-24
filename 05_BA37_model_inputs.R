# Export rank features, encoder inputs, and observation masks.
source("R/project.R", encoding = "UTF-8")
require_stages("04")
source(file.path(CFG$root, "R/06_token_audit.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
source(file.path(CFG$root, "R/06_export.R"), encoding = "UTF-8", local = new.env(parent = globalenv()))
finish_stage("05", file.path(out_dir("06_inputs"), c("export_manifest.tsv", "token_hierarchy.tsv", "ssl_gene_z.float32.bin")))

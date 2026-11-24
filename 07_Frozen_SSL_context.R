# Extract frozen SSL representations.
source("R/project.R", encoding = "UTF-8")
require_stages("06")
context <- out_dir("07_context"); dir.create(context, recursive = TRUE, showWarnings = FALSE)
run_py("extract_context.py", c("--data", out_dir("06_inputs"),
  "--model", file.path(CFG$root, "python/ssl_model.py"),
  "--initialization", out_dir("07_SSL/model_smoke/SSL37_UNTRAINED_init.pt"),
  "--checkpoint", out_dir("07_SSL/pretraining/SSL37_best.pt"),
  "--training", out_dir("07_SSL/pretraining"), "--output", context),
  file.path(context, "console.txt"), "BA37_CONTEXT_EXTRACTION_COMPLETED")
finish_stage("07", file.path(context, c("context_provenance.json", "context_files_manifest.tsv")))

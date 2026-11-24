# Pretrain the SSL encoder.
source("R/project.R", encoding = "UTF-8")
require_stages("05")
data <- out_dir("06_inputs")
model <- file.path(CFG$root, "python/ssl_model.py")
smoke <- out_dir("07_SSL/model_smoke"); training <- out_dir("07_SSL/pretraining")
dir.create(smoke, recursive = TRUE, showWarnings = FALSE)
dir.create(training, recursive = TRUE, showWarnings = FALSE)
run_py("ssl_model.py", c("--data", data, "--output", smoke), file.path(smoke, "console.txt"), "SSL37_MODEL_CHECK: PASS")
run_py("train_ssl.py", c("--data", data, "--model", model,
  "--initialization", file.path(smoke, "SSL37_UNTRAINED_init.pt"), "--output", training),
  file.path(training, "console.txt"), "SSL37_PRETRAINING_COMPLETED")
finish_stage("06", file.path(training, c("SSL37_best.pt", "training_summary.json")))

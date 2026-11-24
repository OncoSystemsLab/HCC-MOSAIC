# Refit the model, score 616 tumors, and compute attribution.
source("R/project.R", encoding = "UTF-8")
require_stages("08")
source(file.path(CFG$root, "R/freeze_audit.R"), encoding = "UTF-8")
FULL <- out_dir("09_final/01_all67_fixed_model")
SCORES <- out_dir("09_final/02_scores616")
ATTR <- out_dir("09_final/03_attribution")
for (d in c(FULL, SCORES, ATTR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
run_py("refit_all67.py", c("--data", out_dir("06_inputs"), "--context", out_dir("07_context"),
  "--loco", out_dir("08_LOCO"), "--fusion", file.path(CFG$root, "python/train_early.py"),
  "--loader", file.path(CFG$root, "python/load_inputs.py"), "--output", FULL), file.path(FULL, "console.txt"), "P009_BA37_ALL67_REFIT_COMPLETED")
run_py("score_616.py", c("--model", FULL, "--data", out_dir("06_inputs"), "--context", out_dir("07_context"),
  "--output", SCORES), file.path(SCORES, "console.txt"), "FIXED_P009_BA37_616_SCORING_COMPLETED")
run_py("attribution.py", c("--model", FULL, "--data", out_dir("06_inputs"), "--context", out_dir("07_context"),
  "--scores", file.path(SCORES, "P009_BA37_fixed_model_scores_616.tsv"), "--output", ATTR),
  file.path(ATTR, "console.txt"), "HIERARCHICAL_ATTRIBUTION_COMPLETED")
saveRDS(analysis_config(), out_dir("09_final/BA37_analysis_config.rds"))
code_sources <- c(list.files(file.path(CFG$root), pattern = "[.]R$|[.]Rproj$", full.names = TRUE),
  list.files(file.path(CFG$root, "R"), pattern = "[.]R$|[.]tsv$", full.names = TRUE),
  list.files(file.path(CFG$root, "python"), pattern = "[.]py$", full.names = TRUE),
  list.files(file.path(CFG$root, "environment"), full.names = TRUE))
relative <- substring(gsub("\\\\", "/", code_sources), nchar(CFG$root) + 2L)
targets <- file.path(out_dir("09_final/00_frozen_development/code_and_environment"), relative)
for (d in unique(dirname(targets))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
if (!all(file.copy(code_sources, targets, overwrite = FALSE)) ||
    any(unname(tools::md5sum(code_sources)) != unname(tools::md5sum(targets)))) stop("\u4EE3\u7801\u5FEB\u7167\u590D\u5236\u5931\u8D25\u3002")
write.table(data.frame(path = relative, md5 = unname(tools::md5sum(code_sources))),
  out_dir("09_final/00_frozen_development/code_environment_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
finish_stage("09", c(file.path(SCORES, "P009_BA37_fixed_model_scores_616.tsv"), file.path(ATTR, "token_importance_summary.tsv")))

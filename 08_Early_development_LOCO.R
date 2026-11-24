# Evaluate early fusion by development LOCO.
source("R/project.R", encoding = "UTF-8")
require_stages("07")
loco <- out_dir("08_LOCO"); panels <- file.path(loco, "fixed_panel")
dir.create(panels, recursive = TRUE, showWarnings = FALSE)
h <- read.delim(file.path(out_dir("06_inputs"), "token_hierarchy.tsv"), fileEncoding = "UTF-8")
stopifnot(nrow(h) == 37L)
roster <- data.frame(panel_id = "P009_BA37", size = 37L, token_name = h$token_name,
  concept_id = h$concept_id, module_id = h$module_id)
write.table(roster, file.path(panels, "candidate_panel_membership.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
run_py("train_early.py", c("--data", out_dir("06_inputs"), "--context", out_dir("07_context"),
  "--panels", panels, "--loader", file.path(CFG$root, "python/load_inputs.py"), "--output", loco),
  file.path(loco, "console.txt"), "FIXED_PANEL_SSL_FUSION_COMPLETED")
print(read.delim(file.path(loco, "fusion_method_summary.tsv")), row.names = FALSE)
finish_stage("08", file.path(loco, c("patient_predictions.tsv", "outer_cohort_metrics.tsv")))

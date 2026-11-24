source("Config.R", encoding = "UTF-8")
STATE_FILE <- file.path(PROJECT_ROOT, "run_state.rds")
if (!file.exists(STATE_FILE)) stop("\u8BF7\u9996\u5148\u8FD0\u884C00_Environment_and_inputs.R\u3002")
CFG <- readRDS(STATE_FILE)
if (unname(tools::md5sum(file.path(PROJECT_ROOT, "HCC_37token.xlsx"))) != CFG$token_md5 ||
    unname(tools::md5sum(file.path(PROJECT_ROOT, "gencode.v50.annotation.gtf.gz"))) != CFG$gtf_md5) {
  stop("37-token\u8868\u6216GENCODE\u5DF2\u6539\u53D8\uFF0C\u8BF7\u4ECE00\u5EFA\u7ACB\u65B0\u8FD0\u884C\uFF0C\u4E0D\u80FD\u6DF7\u7528\u65E7\u4EA7\u7269\u3002")
}
input_records <- readRDS(file.path(CFG$run_root, "00_inventory/stage02_input_inventory.rds"))$input_manifest
if (any(!file.exists(input_records$path)) || any(unname(tools::md5sum(input_records$path)) != input_records$md5)) {
  stop("\u672C\u6B21\u8FD0\u884C\u7684\u961F\u5217\u6570\u636E\u526F\u672C\u7F3A\u5931\u6216\u53D1\u751F\u53D8\u5316\uFF0C\u8BF7\u6838\u67E5\u8F93\u5165\u5E76\u4ECE00\u91CD\u5EFA\u3002")
}
out_dir <- function(key) file.path(CFG$run_root, key)
fig_dir <- function(key = "") out_dir(file.path("09_final/04_figures", key))
artifact <- function(filename) {
  tab <- read.delim(file.path(CFG$root, "R/artifact_registry.tsv"), stringsAsFactors = FALSE)
  hit <- match(filename, tab$filename)
  if (is.na(hit)) stop("\u672A\u6CE8\u518C\u7684\u4E2D\u95F4\u6587\u4EF6\uFF1A", filename)
  out_dir(file.path(tab$directory[hit], filename))
}
analysis_config <- function() list(model_id = "P009_BA37", root = out_dir("09_final"),
  inputs = out_dir("06_inputs"), loco = out_dir("08_LOCO"),
  full_model = out_dir("09_final/01_all67_fixed_model"),
  scores = out_dir("09_final/02_scores616"), attribution = out_dir("09_final/03_attribution"))
run_py <- function(script, args, log, marker = NULL) {
  dir.create(dirname(log), recursive = TRUE, showWarnings = FALSE)
  args <- vapply(args, shQuote, character(1))
  status <- system2(CFG$python, c("-I", shQuote(file.path(CFG$root, "python", script)), args), stdout = log, stderr = log)
  output <- readLines(log, warn = FALSE)
  cat(paste(tail(output, 30), collapse = "\n"), "\n")
  if (status != 0L || (!is.null(marker) && !any(output == marker))) stop("\u8FD0\u884C\u5931\u8D25\uFF1B\u5B8C\u6574\u65E5\u5FD7\uFF1A", log)
  invisible(output)
}
finish_stage <- function(name, expected) {
  if (any(!file.exists(expected))) stop("\u672C\u9636\u6BB5\u4EA7\u7269\u4E0D\u5B8C\u6574\uFF1A", paste(expected[!file.exists(expected)], collapse = ", "))
  CFG$completed <- unique(c(CFG$completed, name))
  saveRDS(CFG, STATE_FILE)
  writeLines(capture.output(sessionInfo()), file.path(CFG$root, "environment", paste0(name, "_R_sessionInfo.txt")))
  cat("\n\u9636\u6BB5\u5B8C\u6210\uFF1A", name, "\n\u672C\u6B21\u8FD0\u884C\uFF1A", CFG$run_root, "\n")
}
require_stages <- function(names) {
  if (!all(names %in% CFG$completed)) stop("\u8BF7\u5148\u5B8C\u6210\u9636\u6BB5\uFF1A", paste(setdiff(names, CFG$completed), collapse = ", "))
}

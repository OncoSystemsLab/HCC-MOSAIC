# Validate inputs and initialize a run.
source("Config.R", encoding = "UTF-8")
needed <- c("readxl", "ggplot2", "pROC", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("\u8BF7\u5148\u5B89\u88C5R\u5305\uFF1A", paste(missing, collapse = ", "))
if (!nzchar(PYTHON_EXE)) stop("Activate an environment with Python available on PATH.")
inputs <- file.path(PROJECT_ROOT, c("HCC_37token.xlsx", "gencode.v50.annotation.gtf.gz"))
stopifnot(all(file.exists(inputs)))
for (d in c("data", "outputs", "environment", "logs")) dir.create(file.path(PROJECT_ROOT, d), recursive = TRUE, showWarnings = FALSE)
missing_inputs <- RAW_INPUTS$path[!file.exists(RAW_INPUTS$path)]
if (length(missing_inputs)) stop("Missing required cohort files:\n", paste(missing_inputs, collapse = "\n"))
RAW_INPUTS$md5 <- unname(tools::md5sum(RAW_INPUTS$path))
RAW_INPUTS$size_bytes <- file.info(RAW_INPUTS$path)$size
write.table(RAW_INPUTS, file.path(PROJECT_ROOT, "data/source_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_root <- file.path(PROJECT_ROOT, "outputs", run_id)
if (dir.exists(run_root)) stop("\u8FD0\u884C\u76EE\u5F55\u5DF2\u5B58\u5728\uFF0C\u8BF7\u7A0D\u540E\u91CD\u8BD5\u3002")
dir.create(run_root, recursive = TRUE)
CFG <- list(root = PROJECT_ROOT, run_id = run_id, run_root = run_root,
  python = PYTHON_EXE, protocol = PROTOCOL, completed = "00",
  token_md5 = unname(tools::md5sum(inputs[1])), gtf_md5 = unname(tools::md5sum(inputs[2])))
saveRDS(CFG, file.path(PROJECT_ROOT, "run_state.rds"))
inventory <- file.path(run_root, "00_inventory"); dir.create(inventory)
input_manifest <- RAW_INPUTS; input_manifest$exists <- TRUE
input_manifest$provenance <- "versioned local input copy"
saveRDS(list(input_manifest = input_manifest), file.path(inventory, "stage02_input_inventory.rds"))
write.table(input_manifest, file.path(inventory, "source_manifest.tsv"), sep = "\t", row.names = FALSE, quote = FALSE, fileEncoding = "UTF-8")
env_log <- file.path(PROJECT_ROOT, "environment/python_check.txt")
status <- system2(PYTHON_EXE, c("-I", shQuote(file.path(PROJECT_ROOT, "python/check_environment.py"))), stdout = env_log, stderr = env_log)
cat(paste(readLines(env_log, warn = FALSE), collapse = "\n"), "\n")
if (status != 0L) stop("Python\u73AF\u5883\u68C0\u67E5\u5931\u8D25\uFF0C\u89C1\uFF1A", env_log)
system2(PYTHON_EXE, c("-I", "-m", "pip", "freeze"), stdout = file.path(PROJECT_ROOT, "environment/requirements.lock.txt"))
conda <- Sys.which("conda")
if (nzchar(conda) && nzchar(CONDA_ENVIRONMENT)) system2(conda, c("env", "export", "--no-builds", "-n", shQuote(CONDA_ENVIRONMENT)), stdout = file.path(PROJECT_ROOT, "environment/conda_environment.yml"))
writeLines(capture.output(sessionInfo()), file.path(PROJECT_ROOT, "environment/R_sessionInfo.txt"))
write.table(as.data.frame(installed.packages()[, c("Package", "Version")]), file.path(PROJECT_ROOT, "environment/R_packages.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
saveRDS(PROTOCOL, file.path(run_root, "training_protocol.rds"))
code <- c(list.files(PROJECT_ROOT, pattern = "[.]R$", full.names = TRUE),
  list.files(file.path(PROJECT_ROOT, "R"), pattern = "[.]R$|[.]tsv$", full.names = TRUE),
  list.files(file.path(PROJECT_ROOT, "python"), pattern = "[.]py$", full.names = TRUE))
write.table(data.frame(path = code, md5 = unname(tools::md5sum(code))), file.path(run_root, "code_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, fileEncoding = "UTF-8")
writeLines(capture.output(str(PROTOCOL)), file.path(run_root, "training_protocol.txt"))
cat("\n00\u5B8C\u6210\uFF1B\u6570\u636E\u526F\u672C\u5DF2\u6838\u67E5\u3002\n\u8FD0\u884C\u76EE\u5F55\uFF1A", run_root, "\n\u8BF7\u4F9D\u6B21\u8FD0\u884C01\u300102\u2026\u2026\u3002\n")

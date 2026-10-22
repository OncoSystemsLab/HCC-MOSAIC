if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

source_files <- c(Stage2_Agent = file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v1.R"), Functional_Scorer = file.path(PROJECT_DIR, 
    "07_agent", "HCC_FunctionalScorer_v1.R"), Association_Engine = file.path(PROJECT_DIR, "07_agent", 
    "HCC_AssociationEngine_v1.R"))

stopifnot(all(file.exists(source_files)))

extract_functions <- function(path, label) {
    expressions <- parse(file = path, keep.source = FALSE)
    found <- list()
    for (expr in expressions) {
        if (!is.call(expr)) 
            next
        operator <- as.character(expr[[1]])[1]
        if (!operator %in% c("<-", "=")) 
            next
        lhs <- expr[[2]]
        rhs <- expr[[3]]
        if (!is.symbol(lhs)) 
            next
        if (!is.call(rhs)) 
            next
        if (!identical(rhs[[1]], as.name("function"))) 
            next
        function_name <- as.character(lhs)
        arguments <- names(as.list(rhs[[2]]))
        found[[length(found) + 1L]] <- data.frame(Engine = label, File = basename(path), Function = function_name, 
            Arguments = paste(arguments, collapse = ", "), stringsAsFactors = FALSE)
    }
    if (length(found) == 0L) {
        return(data.frame(Engine = character(), File = character(), Function = character(), Arguments = character()))
    }
    do.call(rbind, found)
}

inventory <- do.call(rbind, Map(extract_functions, path = unname(source_files), label = names(source_files)))

rownames(inventory) <- NULL

stopifnot(nrow(inventory) > 0L, !anyDuplicated(paste(inventory$Engine, inventory$Function)))

required_functions <- c("score_HCC_modules_fixed_v1", "make_HCC_covariates_v1", "hcc_fit_association_v1")

missing_functions <- setdiff(required_functions, inventory$Function)

stopifnot(length(missing_functions) == 0L)

asset_paths <- c(BA4_Model = file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"), 
    Molecular_Config = file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"), 
    Molecular_Reference = file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"), 
    Functional_Reference = file.path(PROJECT_DIR, "00_config", "Stage3_3C_module_reference_69_v1.rds"), 
    Stage3_Cold_Start_Log = file.path(PROJECT_DIR, "08_reports", "Stage3_3C3_cold_start_log.txt"))

asset_inventory <- data.frame(Asset = names(asset_paths), Path = unname(asset_paths), Exists = file.exists(asset_paths), 
    stringsAsFactors = FALSE)

stopifnot(all(asset_inventory$Exists))

stage2_path <- source_files[["Stage2_Agent"]]

stage2_lines <- readLines(stage2_path, warn = FALSE)

target_functions <- c("run_HCC_OmicsAgent", "read_HCC_TPM_single", "predict_HCC_fingerprint_v1", "render_HCC_fingerprint_v1_1", 
    "render_HCC_fingerprint_v1_2")

target_inventory <- inventory[inventory$Function %in% target_functions, , drop = FALSE]

snippet_lines <- character()

for (fn in target_functions) {
    pattern <- paste0("^\\s*", fn, "\\s*<-\\s*function\\s*\\(")
    hits <- grep(pattern, stage2_lines, perl = TRUE)
    if (length(hits) == 0L) 
        next
    for (start in hits) {
        end <- min(length(stage2_lines), start + 19L)
        snippet_lines <- c(snippet_lines, paste0("===== ", fn, " ====="), sprintf("%04d: %s", seq.int(start, 
            end), stage2_lines[start:end]), "")
    }
}

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(inventory, file.path(report_dir, "Stage4_1A_function_inventory.csv"), row.names = FALSE)

write.csv(asset_inventory, file.path(report_dir, "Stage4_1A_asset_inventory.csv"), row.names = FALSE)

writeLines(snippet_lines, file.path(report_dir, "Stage4_1A_Stage2_entrypoint_source.txt"))

stopifnot(file.exists(file.path(report_dir, "Stage4_1A_function_inventory.csv")), nrow(inventory) >= 
    3L)

PROJECT_DIR <- "HCC_OmicsAgent"

run_HCC_OmicsAgent_v2 <- function(mode = c("patient_rna", "patient_multiomics", "cohort_analysis"), input, 
    project_dir, sample_uid = NULL, mutation_gene = "CTNNB1") {
    mode <- match.arg(mode)
    stopifnot(is.character(project_dir), length(project_dir) == 1L, dir.exists(project_dir))
    load_engine <- function(filename) {
        path <- file.path(project_dir, "07_agent", filename)
        if (!file.exists(path)) {
            stop("Missing engine file: ", path)
        }
        env <- new.env(parent = .GlobalEnv)
        sys.source(path, envir = env)
        env
    }
    read_patient_tpm <- function(path, uid) {
        if (!is.character(path) || length(path) != 1L || !file.exists(path)) {
            stop("Invalid RNA input file.")
        }
        dat <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
        if (ncol(dat) != 2L) {
            stop("Patient RNA input must contain ", "exactly 2 columns: Gene and TPM.")
        }
        genes <- trimws(as.character(dat[[1]]))
        values <- suppressWarnings(as.numeric(dat[[2]]))
        if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) || any(!is.finite(values)) || 
            any(values < 0)) {
            stop("Invalid gene symbols or TPM values.")
        }
        matrix(values, ncol = 1L, dimnames = list(genes, uid))
    }
    validate_mutation <- function(x) {
        if (!is.atomic(x) || is.null(names(x)) || anyNA(names(x)) || any(!nzchar(names(x))) || anyDuplicated(names(x))) {
            stop("Mutation input must be a named ", "0/1 vector.")
        }
        if (anyNA(x) || !all(as.character(x) %in% c("0", "1"))) {
            stop("Mutation status must contain only 0/1.")
        }
        setNames(as.integer(x), names(x))
    }
    functional_env <- load_engine("HCC_FunctionalScorer_v1.R")
    functional_ref <- readRDS(file.path(project_dir, "00_config", "Stage3_3C_module_reference_69_v1.rds"))
    score_function <- get("score_HCC_modules_fixed_v1", envir = functional_env, inherits = FALSE)
    provenance <- list(agent_version = "HCC_OmicsAgent_v2_core", mode = mode, functional_reference = functional_ref$version, 
        timestamp = as.character(Sys.time()), mutation_definition = "1=recorded mutation; 0=no recorded mutation")
    if (mode %in% c("patient_rna", "patient_multiomics")) {
        if (mode == "patient_rna") {
            rna_file <- input
        }
        else {
            if (!is.list(input) || !all(c("rna_file", "mutation_status", "mutation_patient_id") %in% 
                names(input))) {
                stop("patient_multiomics requires ", "rna_file, mutation_status, ", "mutation_patient_id.")
            }
            rna_file <- input$rna_file
        }
        uid <- sample_uid
        if (is.null(uid) && mode == "patient_multiomics") {
            uid <- input$mutation_patient_id
        }
        if (is.null(uid)) {
            uid <- sub("\\.[^.]*$", "", basename(rna_file))
        }
        stopifnot(is.character(uid), length(uid) == 1L, !is.na(uid), nzchar(uid))
        mutation_annotation <- NULL
        if (mode == "patient_multiomics") {
            if (!identical(as.character(input$mutation_patient_id), uid)) {
                stop("RNA and mutation patient IDs ", "do not match.")
            }
            mutation_annotation <- validate_mutation(input$mutation_status)
        }
        tpm <- read_patient_tpm(rna_file, uid)
        functional_try <- tryCatch(score_function(tpm = tpm, reference = functional_ref), error = function(e) e)
        functional_ok <- !inherits(functional_try, "error")
        stage2_try <- tryCatch({
            stage2_env <- load_engine("HCC_OmicsAgent_v1.R")
            original_runner <- get("run_HCC_OmicsAgent", envir = stage2_env, inherits = FALSE)
            original_runner(input_file = rna_file, sample_uid = uid, project_dir = project_dir)
        }, error = function(e) e)
        stage2_ok <- !inherits(stage2_try, "error")
        status <- if (functional_ok && stage2_ok) {
            "success"
        }
        else if (functional_ok || stage2_ok) {
            "partial"
        }
        else {
            "failed"
        }
        qc <- list(n_patients = 1L, n_input_genes = nrow(tpm), functional_status = if (functional_ok) {
            "success"
        } else {
            "unavailable"
        }, functional_reason = if (functional_ok) {
            "NA"
        } else {
            conditionMessage(functional_try)
        }, stage2_status = if (stage2_ok) {
            "success"
        } else {
            "unavailable"
        }, stage2_reason = if (stage2_ok) {
            "NA"
        } else {
            conditionMessage(stage2_try)
        })
        result <- list(status = status, mode = mode, sample_uid = uid, results = list(stage2 = if (stage2_ok) {
            stage2_try
        } else {
            NULL
        }, functional = if (functional_ok) {
            functional_try
        } else {
            NULL
        }, mutation_annotation = mutation_annotation), qc = qc, provenance = provenance)
        class(result) <- c("HCC_OmicsAgent_v2_Result", "list")
        return(result)
    }
    joint <- if (is.character(input) && length(input) == 1L) {
        readRDS(input)
    }
    else {
        input
    }
    required_fields <- c("cohort", "patient_ids", "clinical", "mechanism_TPM", "wnt_mutation")
    if (!is.list(joint) || !all(required_fields %in% names(joint))) {
        stop("Cohort input must be a ", "Stage3-compatible joint dataset.")
    }
    ids <- as.character(joint$patient_ids)
    tpm <- joint$mechanism_TPM
    stopifnot(is.matrix(tpm), is.numeric(tpm), identical(colnames(tpm), ids), !anyDuplicated(ids), is.matrix(joint$wnt_mutation), 
        all(ids %in% colnames(joint$wnt_mutation)), all(is.finite(tpm)), all(tpm >= 0))
    if (!mutation_gene %in% rownames(joint$wnt_mutation)) {
        stop("Mutation gene not found: ", mutation_gene)
    }
    frozen_try <- tryCatch(score_function(tpm, functional_ref), error = function(e) e)
    frozen_ok <- !inherits(frozen_try, "error")
    association_try <- NULL
    association_ok <- FALSE
    if (frozen_ok) {
        association_try <- tryCatch({
            assoc_env <- load_engine("HCC_AssociationEngine_v1.R")
            make_cov <- get("make_HCC_covariates_v1", envir = assoc_env, inherits = FALSE)
            fit_assoc <- get("hcc_fit_association_v1", envir = assoc_env, inherits = FALSE)
            module_names <- names(functional_ref$module_definitions)
            feature_matrix <- t(as.matrix(frozen_try$scores[, module_names, drop = FALSE]))
            colnames(feature_matrix) <- ids
            covariates <- make_cov(joint)
            fit_assoc(feature_matrix = feature_matrix, mutation_matrix = joint$wnt_mutation, covariates = covariates, 
                mutation_gene = mutation_gene)
        }, error = function(e) e)
        association_ok <- !inherits(association_try, "error")
    }
    status <- if (frozen_ok && association_ok) {
        "success"
    }
    else if (frozen_ok) {
        "partial"
    }
    else {
        "failed"
    }
    mutation_status <- as.integer(joint$wnt_mutation[mutation_gene, ids])
    result <- list(status = status, mode = mode, results = list(functional = if (frozen_ok) {
        frozen_try
    } else {
        NULL
    }, association = if (association_ok) {
        association_try
    } else {
        NULL
    }), qc = list(cohort = joint$cohort, n_patients = length(ids), n_input_genes = nrow(tpm), mutation_gene = mutation_gene, 
        n_mutation = sum(mutation_status == 1L, na.rm = TRUE), functional_status = if (frozen_ok) {
            "success"
        } else {
            "unavailable"
        }, functional_reason = if (frozen_ok) {
            "NA"
        } else {
            conditionMessage(frozen_try)
        }, association_status = if (association_ok) {
            "success"
        } else {
            "unavailable"
        }, association_reason = if (association_ok) {
            "NA"
        } else if (frozen_ok) {
            conditionMessage(association_try)
        } else {
            "Functional scoring unavailable."
        }), provenance = provenance)
    class(result) <- c("HCC_OmicsAgent_v2_Result", "list")
    result
}

v2_file <- file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v2.R")

dump("run_HCC_OmicsAgent_v2", file = v2_file, envir = .GlobalEnv)

stopifnot(file.exists(v2_file), length(parse(file = v2_file)) > 0L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v2.R"))

patient_result <- run_HCC_OmicsAgent_v2(mode = "patient_rna", input = file.path(PROJECT_DIR, "03_input", 
    "Patient_TPM.csv"), sample_uid = "PRJEB34724__S02-T-RNA", project_dir = PROJECT_DIR)

if (!is.null(patient_result$results$functional)) {
    print(patient_result$results$functional$scores, row.names = FALSE)
    print(patient_result$results$functional$percentiles, row.names = FALSE)
}

cohort_result <- run_HCC_OmicsAgent_v2(mode = "cohort_analysis", input = file.path(PROJECT_DIR, "04_processed", 
    "Stage3_2A_TCGA_LIHC_joint_v1.rds"), mutation_gene = "CTNNB1", project_dir = PROJECT_DIR)

if (!is.null(cohort_result$results$association)) {
    print(cohort_result$results$association[, c("Feature", "Mutation_Gene", "Beta_Adjusted", "P_Value", 
        "FDR_BH")], digits = 5, row.names = FALSE)
}

saveRDS(list(patient = patient_result, cohort = cohort_result), file.path(PROJECT_DIR, "08_reports", 
    "Stage4_1B_initial_integration_results.rds"))

if (identical(patient_result$status, "success") && identical(cohort_result$status, "success")) {
    cat("\nStage 4.1B: PASS")
} else {
    cat("\nStage 4.1B: NEEDS REVIEW", "(inspect QC reasons)")
}

PROJECT_DIR <- "HCC_OmicsAgent"

test_file <- file.path(PROJECT_DIR, "07_agent", "Stage4_1C_v2_acceptance_test.R")

log_file <- file.path(PROJECT_DIR, "08_reports", "Stage4_1C_v2_acceptance_log.txt")

test_code <- "\nPROJECT_DIR <- \"HCC_OmicsAgent\"\n\nsource(file.path(\n  PROJECT_DIR, \"07_agent\", \"HCC_OmicsAgent_v2.R\"\n))\n\nref <- readRDS(file.path(\n  PROJECT_DIR, \"00_config\",\n  \"Stage3_3C_module_reference_69_v1.rds\"\n))\n\nrna_file <- file.path(\n  PROJECT_DIR, \"03_input\", \"Patient_TPM.csv\"\n)\n\nuid <- \"PRJEB34724__S02-T-RNA\"\n\nstopifnot(file.exists(rna_file))\n\n\n\n\n\nrna_result <- run_HCC_OmicsAgent_v2(\n  mode = \"patient_rna\",\n  input = rna_file,\n  sample_uid = uid,\n  project_dir = PROJECT_DIR\n)\n\nstopifnot(\n  identical(rna_result$status, \"success\"),\n  identical(rna_result$qc$stage2_status, \"success\"),\n  identical(rna_result$qc$functional_status, \"success\"),\n  nrow(rna_result$results$functional$scores) == 1L\n)\n\nwriteLines(\"TEST 1 - patient_rna: PASS\")\n\n\n\n\n\n\nsynthetic_mutations <- c(\n  SYNTHETIC_GENE_A = 1L,\n  SYNTHETIC_GENE_B = 0L\n)\n\nmulti_result <- run_HCC_OmicsAgent_v2(\n  mode = \"patient_multiomics\",\n  input = list(\n    rna_file = rna_file,\n    mutation_status = synthetic_mutations,\n    mutation_patient_id = uid\n  ),\n  sample_uid = uid,\n  project_dir = PROJECT_DIR\n)\n\nstopifnot(\n  identical(multi_result$status, \"success\"),\n  identical(\n    multi_result$results$mutation_annotation,\n    synthetic_mutations\n  ),\n  identical(\n    multi_result$qc$functional_status,\n    \"success\"\n  )\n)\n\nmods <- names(ref$module_definitions)\n\nmulti_scores <- as.numeric(\n  multi_result$results$functional$scores[1, mods]\n)\n\nrna_scores <- as.numeric(\n  rna_result$results$functional$scores[1, mods]\n)\n\nstopifnot(\n  max(abs(multi_scores - rna_scores)) < 1e-12\n)\n\nwriteLines(\"TEST 2 - patient_multiomics: PASS\")\n\n\n\n\n\nbad_id <- try(\n  run_HCC_OmicsAgent_v2(\n    mode = \"patient_multiomics\",\n    input = list(\n      rna_file = rna_file,\n      mutation_status = synthetic_mutations,\n      mutation_patient_id = \"DIFFERENT_PATIENT\"\n    ),\n    sample_uid = uid,\n    project_dir = PROJECT_DIR\n  ),\n  silent = TRUE\n)\n\nstopifnot(inherits(bad_id, \"try-error\"))\n\nwriteLines(\"TEST 3 - patient ID mismatch rejected: PASS\")\n\n\n\n\n\nbad_mutation <- try(\n  run_HCC_OmicsAgent_v2(\n    mode = \"patient_multiomics\",\n    input = list(\n      rna_file = rna_file,\n      mutation_status = c(SYNTHETIC_GENE_A = 2L),\n      mutation_patient_id = uid\n    ),\n    sample_uid = uid,\n    project_dir = PROJECT_DIR\n  ),\n  silent = TRUE\n)\n\nstopifnot(inherits(bad_mutation, \"try-error\"))\n\nwriteLines(\"TEST 4 - invalid mutation rejected: PASS\")\n\n\n\n\n\n\n\nfull_df <- read.csv(\n  rna_file,\n  check.names = FALSE,\n  stringsAsFactors = FALSE\n)\n\nrequired <- ref$required_genes\n\nidx <- match(required, as.character(full_df[[1]]))\n\nstopifnot(!anyNA(idx))\n\nshort_df <- full_df[\n  idx, ,\n  drop = FALSE\n]\n\nshort_file <- file.path(\n  PROJECT_DIR,\n  \"08_reports\",\n  \"Stage4_1C_TEST_ONLY_16gene_TPM.csv\"\n)\n\nwrite.csv(\n  short_df,\n  short_file,\n  row.names = FALSE\n)\n\nshort_result <- run_HCC_OmicsAgent_v2(\n  mode = \"patient_rna\",\n  input = short_file,\n  sample_uid = \"TEST_ONLY_16GENES\",\n  project_dir = PROJECT_DIR\n)\n\nstopifnot(\n  identical(short_result$status, \"partial\"),\n  identical(\n    short_result$qc$functional_status,\n    \"success\"\n  ),\n  identical(\n    short_result$qc$stage2_status,\n    \"unavailable\"\n  ),\n  is.null(short_result$results$stage2),\n  !is.null(short_result$qc$stage2_reason)\n)\n\nwriteLines(\"TEST 5 - partial input safe degradation: PASS\")\n\n\n\n\n\nexpected_mut_counts <- c(\n  TCGA_LIHC = 27L,\n  ICGC_LICA_FR = 14L\n)\n\nfor (cohort in names(expected_mut_counts)) {\n\n  joint_file <- file.path(\n    PROJECT_DIR,\n    \"04_processed\",\n    paste0(\n      \"Stage3_2A_\",\n      cohort,\n      \"_joint_v1.rds\"\n    )\n  )\n\n  result <- run_HCC_OmicsAgent_v2(\n    mode = \"cohort_analysis\",\n    input = joint_file,\n    mutation_gene = \"AXIN1\",\n    project_dir = PROJECT_DIR\n  )\n\n  stopifnot(\n    identical(result$status, \"success\"),\n    result$qc$n_mutation ==\n      unname(expected_mut_counts[[cohort]]),\n    nrow(result$results$association) == 4L,\n    all(is.finite(\n      result$results$association$Beta_Adjusted\n    ))\n  )\n\n  writeLines(paste(\n    \"TEST 6 -\", cohort,\n    \"AXIN1 cohort analysis: PASS\"\n  ))\n}\n\n\n\n\n\nbad_gene <- try(\n  run_HCC_OmicsAgent_v2(\n    mode = \"cohort_analysis\",\n    input = file.path(\n      PROJECT_DIR,\n      \"04_processed\",\n      \"Stage3_2A_TCGA_LIHC_joint_v1.rds\"\n    ),\n    mutation_gene = \"UNAVAILABLE_GENE\",\n    project_dir = PROJECT_DIR\n  ),\n  silent = TRUE\n)\n\nstopifnot(inherits(bad_gene, \"try-error\"))\n\nwriteLines(\"TEST 7 - missing driver rejected: PASS\")\n\nwriteLines(\"Stage 4.1C COLD_START: PASS\")"

writeLines(test_code, test_file, useBytes = TRUE)

stopifnot(file.exists(test_file), length(parse(file = test_file)) > 0L)

rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") {
    "Rscript.exe"
} else {
    "Rscript"
})

stopifnot(file.exists(rscript_bin))

output <- suppressWarnings(system2(command = rscript_bin, args = c("--vanilla", shQuote(test_file)), 
    stdout = TRUE, stderr = TRUE))

writeLines(output, log_file, useBytes = TRUE)

status_code <- attr(output, "status")

if (is.null(status_code)) {
    status_code <- 0L
}

passed <- any(grepl("Stage 4.1C COLD_START: PASS", output, fixed = TRUE))

if (status_code != 0L || !passed) {
    stop("Stage 4.1C failed. Inspect log: ", log_file)
}

PROJECT_DIR <- "HCC_OmicsAgent"


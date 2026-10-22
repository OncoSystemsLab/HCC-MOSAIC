if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

run_HCC_OmicsAgent_auto_v1 <- function(input, project_dir, mode = "auto", sample_uid = NULL, mutation_gene = "CTNNB1", 
    dry_run = FALSE) {
    mode <- match.arg(mode, c("auto", "patient_rna", "patient_multiomics", "cohort_analysis"))
    stopifnot(is.character(project_dir), length(project_dir) == 1L, dir.exists(project_dir), is.logical(dry_run), 
        length(dry_run) == 1L, !is.na(dry_run))
    load_function <- function(filename, function_name) {
        path <- file.path(project_dir, "07_agent", filename)
        if (!file.exists(path)) {
            stop("Engine not found: ", path)
        }
        env <- new.env(parent = globalenv())
        sys.source(path, envir = env)
        get(function_name, envir = env, inherits = FALSE)
    }
    router <- load_function("HCC_TaskRouter_v1.R", "plan_HCC_task_v1")
    plan <- router(input = input, project_dir = project_dir, mode = mode, sample_uid = sample_uid, mutation_gene = mutation_gene)
    new_result <- function(status, results = list(), qc = list(), provenance = list()) {
        out <- list(status = status, mode = plan$detected_mode, plan = plan, results = results, qc = qc, 
            provenance = c(list(agent_version = "HCC_AutoExecutor_v1", router_version = plan$provenance$router_version, 
                execution_time = as.character(Sys.time())), provenance))
        class(out) <- c("HCC_OmicsAgent_Auto_Result", "list")
        out
    }
    if (identical(plan$status, "blocked")) {
        return(new_result(status = "blocked", qc = list(reason = plan$qc$reason), provenance = list(execution_path = "not_executed")))
    }
    if (dry_run) {
        return(new_result(status = "planned", qc = plan$qc, provenance = list(execution_path = "dry_run")))
    }
    cap_status <- function(name) {
        i <- match(name, plan$capabilities$Capability)
        if (is.na(i)) 
            return("not_applicable")
        as.character(plan$capabilities$Status[i])
    }
    cap_reason <- function(name) {
        i <- match(name, plan$capabilities$Capability)
        if (is.na(i)) 
            return("")
        as.character(plan$capabilities$Reason[i])
    }
    resolved_mode <- plan$detected_mode
    functional_ready <- identical(cap_status("functional_modules"), "ready")
    ba4_candidate <- identical(cap_status("BA4_stage2"), "candidate")
    association_ready <- identical(cap_status("mutation_association"), "ready")
    is_patient <- resolved_mode %in% c("patient_rna", "patient_multiomics")
    patient_module_only <- (is_patient && functional_ready && !ba4_candidate)
    cohort_module_only <- (identical(resolved_mode, "cohort_analysis") && functional_ready && !association_ready)
    execution <- tryCatch({
        if (patient_module_only || cohort_module_only) {
            scorer <- load_function("HCC_FunctionalScorer_v1.R", "score_HCC_modules_fixed_v1")
            ref <- readRDS(file.path(project_dir, "00_config", "Stage3_3C_module_reference_69_v1.rds"))
            if (patient_module_only) {
                rna_file <- if (identical(resolved_mode, "patient_rna")) 
                  input
                else input$rna_file
                uid <- sample_uid
                if (is.null(uid) && identical(resolved_mode, "patient_multiomics")) {
                  uid <- input$mutation_patient_id
                }
                if (is.null(uid)) {
                  uid <- sub("\\.[^.]*$", "", basename(rna_file))
                }
                df <- read.csv(rna_file, check.names = FALSE, stringsAsFactors = FALSE)
                tpm <- matrix(as.numeric(df[[2]]), ncol = 1L, dimnames = list(trimws(as.character(df[[1]])), 
                  uid))
                functional <- scorer(tpm = tpm, reference = ref)
                mutation_annotation <- NULL
                if (identical(resolved_mode, "patient_multiomics")) {
                  x <- input$mutation_status
                  mutation_annotation <- setNames(as.integer(x), names(x))
                }
                return(new_result(status = "partial", results = list(stage2 = NULL, functional = functional, 
                  mutation_annotation = mutation_annotation), qc = list(n_patients = 1L, n_input_genes = nrow(tpm), 
                  functional_status = "success", stage2_status = "skipped", stage2_reason = cap_reason("BA4_stage2"), 
                  biological_sample_match_verified = FALSE), provenance = list(execution_path = "functional_only", 
                  functional_reference = ref$version)))
            }
            joint <- if (is.character(input) && length(input) == 1L) 
                readRDS(input)
            else input
            functional <- scorer(tpm = joint$mechanism_TPM, reference = ref)
            return(new_result(status = "partial", results = list(functional = functional, association = NULL), 
                qc = list(cohort = joint$cohort, n_patients = length(joint$patient_ids), mutation_gene = mutation_gene, 
                  functional_status = "success", association_status = "skipped", association_reason = cap_reason("mutation_association")), 
                provenance = list(execution_path = "functional_only", functional_reference = ref$version)))
        }
        v2 <- load_function("HCC_OmicsAgent_v2.R", "run_HCC_OmicsAgent_v2")
        actual <- v2(mode = resolved_mode, input = input, project_dir = project_dir, sample_uid = sample_uid, 
            mutation_gene = mutation_gene)
        actual$plan <- plan
        actual$provenance$auto_executor_version <- "HCC_AutoExecutor_v1"
        actual$provenance$execution_path <- "full_v2"
        class(actual) <- unique(c("HCC_OmicsAgent_Auto_Result", class(actual)))
        actual
    }, error = function(e) {
        new_result(status = "failed", qc = list(execution_error = conditionMessage(e)), provenance = list(execution_path = "execution_failed"))
    })
    execution
}

executor_file <- file.path(PROJECT_DIR, "07_agent", "HCC_AutoExecutor_v1.R")

dump("run_HCC_OmicsAgent_auto_v1", file = executor_file, envir = .GlobalEnv)

stopifnot(file.exists(executor_file), length(parse(file = executor_file)) > 0L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_AutoExecutor_v1.R"))

full_file <- file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv")

short_file <- file.path(PROJECT_DIR, "08_reports", "Stage4_1C_TEST_ONLY_16gene_TPM.csv")

tcga_file <- file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")

stopifnot(all(file.exists(c(full_file, short_file, tcga_file))))

r1 <- run_HCC_OmicsAgent_auto_v1(input = full_file, sample_uid = "PRJEB34724__S02-T-RNA", project_dir = PROJECT_DIR)

stopifnot(identical(r1$status, "success"), identical(r1$mode, "patient_rna"), !is.null(r1$results$functional), 
    !is.null(r1$results$stage2))

r2 <- run_HCC_OmicsAgent_auto_v1(input = short_file, project_dir = PROJECT_DIR)

stopifnot(identical(r2$status, "partial"), identical(r2$mode, "patient_rna"), identical(r2$qc$stage2_status, 
    "skipped"), !is.null(r2$results$functional), is.null(r2$results$stage2))

r3 <- run_HCC_OmicsAgent_auto_v1(input = list(rna_file = full_file, mutation_status = c(SYNTHETIC_A = 1L, 
    SYNTHETIC_B = 0L), mutation_patient_id = "PRJEB34724__S02-T-RNA"), sample_uid = "PRJEB34724__S02-T-RNA", 
    project_dir = PROJECT_DIR)

stopifnot(identical(r3$status, "success"), identical(r3$mode, "patient_multiomics"), length(r3$results$mutation_annotation) == 
    2L)

r4 <- run_HCC_OmicsAgent_auto_v1(input = tcga_file, project_dir = PROJECT_DIR, mutation_gene = "CTNNB1")

stopifnot(identical(r4$status, "success"), identical(r4$mode, "cohort_analysis"), r4$qc$n_patients == 
    361L, nrow(r4$results$association) == 4L)

r5 <- run_HCC_OmicsAgent_auto_v1(input = tcga_file, project_dir = PROJECT_DIR, mutation_gene = "UNAVAILABLE_GENE")

stopifnot(identical(r5$status, "partial"), !is.null(r5$results$functional), is.null(r5$results$association), 
    identical(r5$qc$association_status, "skipped"))

r6 <- run_HCC_OmicsAgent_auto_v1(input = file.path(PROJECT_DIR, "03_input", "NONEXISTENT_FILE.csv"), 
    project_dir = PROJECT_DIR)

stopifnot(identical(r6$status, "blocked"), is.null(r6$results$functional))

r7 <- run_HCC_OmicsAgent_auto_v1(input = full_file, project_dir = PROJECT_DIR, dry_run = TRUE)

stopifnot(identical(r7$status, "planned"), identical(r7$mode, "patient_rna"), length(r7$results) == 0L)

all_results <- list(Full_RNA = r1, Partial_RNA = r2, Multiomics = r3, Cohort = r4, Missing_Driver = r5, 
    Missing_File = r6, Dry_Run = r7)

summary_table <- do.call(rbind, lapply(names(all_results), function(nm) {
    x <- all_results[[nm]]
    data.frame(Test = nm, Mode = if (length(x$mode) == 0L) {
        "NA"
    }
    else {
        x$mode
    }, Status = x$status, Execution_Path = if (is.null(x$provenance$execution_path)) {
        "NA"
    }
    else {
        x$provenance$execution_path
    }, stringsAsFactors = FALSE)
}))

write.csv(summary_table, file.path(PROJECT_DIR, "08_reports", "Stage4_2B_auto_execution_summary.csv"), 
    row.names = FALSE)

saveRDS(all_results, file.path(PROJECT_DIR, "08_reports", "Stage4_2B_auto_execution_results.rds"))

PROJECT_DIR <- "HCC_OmicsAgent"

write_HCC_report_v1 <- function(result, project_dir, run_label = "AgentRun") {
    stopifnot(is.list(result), length(result$status) == 1L, result$status %in% c("success", "partial", 
        "failed", "blocked", "planned"), is.character(project_dir), length(project_dir) == 1L, dir.exists(project_dir))
    safe_text <- function(x) {
        if (is.null(x) || length(x) == 0L) {
            return("NA")
        }
        x <- as.character(x)
        x[is.na(x)] <- "NA"
        paste(gsub("[\r\n|]", " ", x), collapse = "; ")
    }
    root <- file.path(project_dir, "08_reports", "AgentRuns")
    dir.create(root, recursive = TRUE, showWarnings = FALSE)
    label <- gsub("[^A-Za-z0-9_-]", "_", as.character(run_label)[1])
    if (!nzchar(label)) 
        label <- "AgentRun"
    base_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", label)
    run_id <- base_id
    run_dir <- file.path(root, run_id)
    counter <- 1L
    while (dir.exists(run_dir)) {
        counter <- counter + 1L
        run_id <- paste0(base_id, "_", sprintf("%02d", counter))
        run_dir <- file.path(root, run_id)
    }
    stopifnot(dir.create(run_dir))
    res <- result$results
    if (is.null(res) || !is.list(res) || result$status %in% c("blocked", "planned")) {
        res <- list()
    }
    has_stage2 <- !is.null(res$stage2)
    has_functional <- !is.null(res$functional)
    has_association <- is.data.frame(res$association)
    has_mutation <- !is.null(res$mutation_annotation)
    mode <- safe_text(result$mode)
    execution_path <- safe_text(result$provenance$execution_path)
    summary_df <- data.frame(Run_ID = run_id, Status = result$status, Mode = mode, Execution_Path = execution_path, 
        Stage2_Available = has_stage2, Functional_Available = has_functional, Association_Available = has_association, 
        Mutation_Annotation_Available = has_mutation, stringsAsFactors = FALSE)
    write.csv(summary_df, file.path(run_dir, "summary.csv"), row.names = FALSE)
    caps <- result$plan$capabilities
    if (is.data.frame(caps)) {
        write.csv(caps, file.path(run_dir, "capabilities.csv"), row.names = FALSE)
    }
    flatten_fields <- function(x, prefix = "") {
        if (is.null(x) || length(x) == 0L) {
            return(data.frame(Field = character(), Value = character()))
        }
        if (!is.list(x)) {
            x <- list(value = x)
        }
        nms <- names(x)
        if (is.null(nms)) {
            nms <- paste0("field_", seq_along(x))
        }
        data.frame(Field = paste0(prefix, nms), Value = vapply(x, safe_text, character(1)), stringsAsFactors = FALSE)
    }
    qc_df <- flatten_fields(result$qc)
    write.csv(qc_df, file.path(run_dir, "qc.csv"), row.names = FALSE)
    provenance_df <- rbind(flatten_fields(result$provenance, "execution."), flatten_fields(result$plan$provenance, 
        "preflight."))
    write.csv(provenance_df, file.path(run_dir, "provenance.csv"), row.names = FALSE)
    if (has_functional) {
        functional <- res$functional
        stopifnot(is.data.frame(functional$scores), is.data.frame(functional$percentiles), nrow(functional$scores) == 
            nrow(functional$percentiles))
        write.csv(functional$scores, file.path(run_dir, "functional_scores.csv"), row.names = FALSE)
        write.csv(functional$percentiles, file.path(run_dir, "functional_percentiles.csv"), row.names = FALSE)
    }
    if (has_mutation) {
        mutation <- res$mutation_annotation
        mutation_df <- data.frame(Gene = names(mutation), Recorded_Mutation_0_1 = as.integer(mutation), 
            stringsAsFactors = FALSE)
        write.csv(mutation_df, file.path(run_dir, "mutation_annotation.csv"), row.names = FALSE)
    }
    if (has_association) {
        write.csv(res$association, file.path(run_dir, "association_results.csv"), row.names = FALSE)
    }
    lines <- c(paste0("", run_id), "", "", "", paste0("- Status: ", result$status), paste0("- Mode: ", 
        mode), paste0("- Execution path: ", execution_path), paste0("- Stage2 available: ", has_stage2), 
        paste0("- Functional scores available: ", has_functional), paste0("- Association available: ", 
            has_association), "", "", "")
    if (nrow(qc_df) > 0L) {
        lines <- c(lines, sprintf("- %s: %s", qc_df$Field, qc_df$Value))
    }
    else {
        lines <- c(lines, "No additional QC fields.")
    }
    lines <- c(lines, "", "", "")
    if (has_functional) {
        scores <- res$functional$scores
        pcts <- res$functional$percentiles
        module_names <- setdiff(names(scores), "Patient_ID")
        if (nrow(scores) == 1L) {
            lines <- c(lines, "| Module | Fixed score | Reference percentile |", "|---|---:|---:|")
            for (m in module_names) {
                lines <- c(lines, sprintf("| %s | %.4f | %.2f |", m, as.numeric(scores[[m]][1]), as.numeric(pcts[[m]][1])))
            }
        }
        else {
            lines <- c(lines, paste0("Batch scoring completed for ", nrow(scores), " patients. See functional_scores.csv ", 
                "and functional_percentiles.csv."))
        }
    }
    else {
        lines <- c(lines, "Functional scores were not produced.")
    }
    lines <- c(lines, "", "", "")
    if (has_mutation) {
        lines <- c(lines, paste0("Mutation annotations: ", length(res$mutation_annotation), " genes. Source identity is not independently verified."))
    }
    if (has_association) {
        lines <- c(lines, paste0("Association analysis completed for ", nrow(res$association), " features. See association_results.csv."))
    }
    else {
        lines <- c(lines, "No association results were produced.")
    }
    lines <- c(lines, "", "", "")
    if (has_stage2) {
        lines <- c(lines, paste("Stage2 engine returned a result.", "The original object is preserved in result.rds.", 
            "Exact ICBScore extraction will be added after", "the legacy return schema is verified."))
    }
    else {
        lines <- c(lines, "Stage2 prediction was not produced.")
    }
    lines <- c(lines, "", "", "", paste("- Functional module scores describe relative", "RNA expression, not bile acid concentration,", 
        "metabolic flux, or pathway activity measured directly."), paste("- Reference percentiles use a 69-patient", 
        "development-cohort distribution and are not", "clinical normal ranges."), paste("- Mutation value 0 means no recorded mutation", 
        "in the supplied matrix; it is not necessarily", "sequence-verified wild type."), paste("- Mutation-expression associations do not", 
        "establish causal mechanisms or actual ICB response."), paste("- Synthetic mutation fixtures, when used,", 
        "must not be interpreted as patient WES findings."))
    writeLines(lines, file.path(run_dir, "report.md"), useBytes = TRUE)
    saveRDS(result, file.path(run_dir, "result.rds"))
    files <- list.files(run_dir, full.names = FALSE)
    stopifnot(file.exists(file.path(run_dir, "summary.csv")), file.exists(file.path(run_dir, "report.md")), 
        file.exists(file.path(run_dir, "result.rds")))
    list(status = "success", run_id = run_id, run_dir = run_dir, files = files, report_file = file.path(run_dir, 
        "report.md"))
}

writer_file <- file.path(PROJECT_DIR, "07_agent", "HCC_ReportWriter_v1.R")

dump("write_HCC_report_v1", file = writer_file, envir = .GlobalEnv)

stopifnot(file.exists(writer_file), length(parse(file = writer_file)) > 0L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_ReportWriter_v1.R"))

previous_results <- readRDS(file.path(PROJECT_DIR, "08_reports", "Stage4_2B_auto_execution_results.rds"))

stopifnot(length(previous_results) == 7L)

exported <- list()

for (nm in names(previous_results)) {
    exported[[nm]] <- write_HCC_report_v1(result = previous_results[[nm]], project_dir = PROJECT_DIR, 
        run_label = nm)
    cat("", nm, "| Export:", exported[[nm]]$status, "| Files:", length(exported[[nm]]$files), "")
}

for (nm in c("Full_RNA", "Partial_RNA", "Multiomics", "Cohort", "Missing_Driver")) {
    stopifnot(file.exists(file.path(exported[[nm]]$run_dir, "functional_scores.csv")))
}

for (nm in c("Missing_File", "Dry_Run")) {
    stopifnot(!file.exists(file.path(exported[[nm]]$run_dir, "functional_scores.csv")), !file.exists(file.path(exported[[nm]]$run_dir, 
        "association_results.csv")))
}

stopifnot(file.exists(file.path(exported$Cohort$run_dir, "association_results.csv")), !file.exists(file.path(exported$Missing_Driver$run_dir, 
    "association_results.csv")), file.exists(file.path(exported$Multiomics$run_dir, "mutation_annotation.csv")))

stage2_result <- previous_results$Full_RNA$results$stage2

if (is.list(stage2_result)) {
    cat("Top-level fields:", paste(names(stage2_result), collapse = ", "), "")
}


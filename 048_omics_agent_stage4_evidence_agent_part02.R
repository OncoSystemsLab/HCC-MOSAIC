if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

plan_HCC_task_v1 <- function(input, project_dir, mode = "auto", sample_uid = NULL, mutation_gene = "CTNNB1") {
    valid_modes <- c("auto", "patient_rna", "patient_multiomics", "cohort_analysis")
    stopifnot(mode %in% valid_modes, dir.exists(project_dir))
    new_capabilities <- function() {
        data.frame(Capability = c("functional_modules", "BA4_stage2", "mutation_annotation", "mutation_association"), 
            Status = rep("not_applicable", 4L), Reason = rep("", 4L), stringsAsFactors = FALSE)
    }
    caps <- new_capabilities()
    set_cap <- function(name, status, reason) {
        i <- match(name, caps$Capability)
        caps$Status[i] <<- status
        caps$Reason[i] <<- reason
    }
    result <- tryCatch({
        reference <- readRDS(file.path(project_dir, "00_config", "Stage3_3C_module_reference_69_v1.rds"))
        required_genes <- reference$required_genes
        joint_cached <- NULL
        detected <- NULL
        if (is.list(input) && all(c("rna_file", "mutation_status", "mutation_patient_id") %in% names(input))) {
            detected <- "patient_multiomics"
        }
        else if (is.list(input) && all(c("cohort", "patient_ids", "clinical", "mechanism_TPM", "wnt_mutation") %in% 
            names(input))) {
            detected <- "cohort_analysis"
            joint_cached <- input
        }
        else if (is.character(input) && length(input) == 1L && file.exists(input)) {
            if (grepl("\\.csv$", input, ignore.case = TRUE)) {
                detected <- "patient_rna"
            }
            else if (grepl("\\.rds$", input, ignore.case = TRUE)) {
                obj <- readRDS(input)
                if (is.list(obj) && all(c("cohort", "patient_ids", "clinical", "mechanism_TPM", "wnt_mutation") %in% 
                  names(obj))) {
                  detected <- "cohort_analysis"
                  joint_cached <- obj
                }
                else {
                  stop("RDS is not a supported Stage3 joint dataset.")
                }
            }
            else {
                stop("Unsupported input file extension.")
            }
        }
        else {
            stop("Unsupported input structure or missing file.")
        }
        if (!identical(mode, "auto") && !identical(mode, detected)) {
            stop("Requested mode does not match input type. ", "Requested: ", mode, "; detected: ", detected)
        }
        resolved_mode <- detected
        qc <- list()
        if (resolved_mode %in% c("patient_rna", "patient_multiomics")) {
            rna_file <- if (resolved_mode == "patient_rna") 
                input
            else input$rna_file
            if (!is.character(rna_file) || length(rna_file) != 1L || !file.exists(rna_file)) {
                stop("RNA file missing.")
            }
            dat <- read.csv(rna_file, check.names = FALSE, stringsAsFactors = FALSE)
            if (ncol(dat) != 2L) {
                stop("Patient TPM must have exactly two columns.")
            }
            genes <- trimws(as.character(dat[[1]]))
            values <- suppressWarnings(as.numeric(dat[[2]]))
            if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes) || any(!is.finite(values)) || 
                any(values < 0)) {
                stop("Invalid gene identifiers or TPM values.")
            }
            missing_modules <- setdiff(required_genes, genes)
            if (length(missing_modules) == 0L) {
                set_cap("functional_modules", "ready", "All 16 required module genes present.")
            }
            else {
                set_cap("functional_modules", "unavailable", paste("Missing module genes:", paste(missing_modules, 
                  collapse = ", ")))
            }
            if (length(genes) == 21308L) {
                set_cap("BA4_stage2", "candidate", paste("Gene count matches development background.", 
                  "Exact gene identity and model compatibility", "must be checked by Stage2."))
            }
            else {
                set_cap("BA4_stage2", "unavailable", paste("Input has", length(genes), "genes; exact frozen BA4 input", 
                  "requirements are not met by this preflight."))
            }
            qc$n_input_genes <- length(genes)
            qc$n_required_module_genes <- length(required_genes)
            qc$n_missing_module_genes <- length(missing_modules)
            qc$rna_identity_verified <- FALSE
            if (resolved_mode == "patient_multiomics") {
                mutation <- input$mutation_status
                if (!(is.integer(mutation) || (is.numeric(mutation) && !is.factor(mutation)))) {
                  stop("Mutation status must be numeric 0/1, ", "not character or factor.")
                }
                if (is.null(names(mutation)) || !length(mutation) || anyNA(names(mutation)) || any(!nzchar(names(mutation))) || 
                  anyDuplicated(names(mutation)) || anyNA(mutation) || !all(mutation %in% c(0, 1))) {
                  stop("Invalid named mutation vector.")
                }
                mut_id <- input$mutation_patient_id
                if (!is.character(mut_id) || length(mut_id) != 1L || is.na(mut_id) || !nzchar(mut_id)) {
                  stop("Invalid mutation patient ID.")
                }
                if (!is.null(sample_uid) && !identical(sample_uid, mut_id)) {
                  stop("Supplied RNA and mutation IDs mismatch.")
                }
                set_cap("mutation_annotation", "ready", paste("Mutation annotations are structurally valid.", 
                  "Biological sample identity remains unverified."))
                qc$mutation_patient_id <- mut_id
                qc$mutation_genes <- names(mutation)
            }
        }
        else {
            joint <- joint_cached
            ids <- as.character(joint$patient_ids)
            tpm <- joint$mechanism_TPM
            mut <- joint$wnt_mutation
            if (!length(ids) || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids) || !is.matrix(tpm) || 
                !is.numeric(tpm) || !identical(colnames(tpm), ids) || is.null(rownames(tpm)) || anyDuplicated(rownames(tpm)) || 
                any(!is.finite(tpm)) || any(tpm < 0) || !is.matrix(mut) || !all(ids %in% colnames(mut)) || 
                is.null(rownames(mut))) {
                stop("Invalid Stage3 joint dataset structure.")
            }
            missing_modules <- setdiff(required_genes, rownames(tpm))
            modules_ready <- length(missing_modules) == 0L
            if (modules_ready) {
                set_cap("functional_modules", "ready", "Frozen 16-gene scoring is available.")
            }
            else {
                set_cap("functional_modules", "unavailable", paste("Missing genes:", paste(missing_modules, 
                  collapse = ", ")))
            }
            qc$cohort <- as.character(joint$cohort)
            qc$n_patients <- length(ids)
            qc$n_input_genes <- nrow(tpm)
            qc$mutation_gene <- mutation_gene
            if (!mutation_gene %in% rownames(mut)) {
                set_cap("mutation_association", "unavailable", "Requested mutation gene is not in the matrix.")
            }
            else {
                status <- mut[mutation_gene, ids]
                if (!is.numeric(status) || anyNA(status) || !all(status %in% c(0, 1))) {
                  stop("Invalid cohort mutation encoding.")
                }
                n_mut <- sum(status == 1)
                n_ref <- sum(status == 0)
                qc$n_mutation <- n_mut
                qc$n_no_recorded_mutation <- n_ref
                cl <- joint$clinical
                clinical_ok <- is.data.frame(cl) && "Patient_ID" %in% names(cl) && identical(as.character(cl$Patient_ID), 
                  ids)
                if (clinical_ok) {
                  if (identical(joint$cohort, "TCGA_LIHC")) {
                    age <- suppressWarnings(as.numeric(cl$age_at_index))
                    sex <- cl$sex_at_birth
                  }
                  else if (identical(joint$cohort, "ICGC_LICA_FR")) {
                    age <- suppressWarnings(as.numeric(cl$Age))
                    sex <- cl$Gender
                  }
                  else {
                    age <- numeric()
                    sex <- character()
                  }
                  sex <- toupper(trimws(as.character(sex)))
                  clinical_ok <- (length(age) == length(ids) && length(sex) == length(ids) && all(is.finite(age)) && 
                    all(age > 0 & age < 120) && all(sex %in% c("M", "MALE", "F", "FEMALE")) && any(sex %in% 
                    c("M", "MALE")) && any(sex %in% c("F", "FEMALE")))
                }
                if (modules_ready && clinical_ok && min(n_mut, n_ref) >= 10L) {
                  set_cap("mutation_association", "ready", "Modules, mutation groups and clinical fields valid.")
                }
                else {
                  set_cap("mutation_association", "unavailable", paste("Requires module scores, supported", 
                    "age/sex fields and at least", "10 patients in each mutation group."))
                }
            }
        }
        active <- caps$Status[caps$Status != "not_applicable"]
        plan_status <- if (!any(active %in% c("ready", "candidate"))) {
            "blocked"
        }
        else if (any(active == "unavailable")) {
            "partial"
        }
        else {
            "ready"
        }
        list(status = plan_status, detected_mode = resolved_mode, capabilities = caps, qc = qc, provenance = list(router_version = "Stage4_TaskRouter_v1", 
            functional_reference = reference$version, preflight_only = TRUE))
    }, error = function(e) {
        list(status = "blocked", detected_mode = "NA", capabilities = new_capabilities(), qc = list(reason = conditionMessage(e)), 
            provenance = list(router_version = "Stage4_TaskRouter_v1", preflight_only = TRUE))
    })
    class(result) <- c("HCC_OmicsAgent_TaskPlan", "list")
    result
}

router_file <- file.path(PROJECT_DIR, "07_agent", "HCC_TaskRouter_v1.R")

dump("plan_HCC_task_v1", file = router_file, envir = .GlobalEnv)

stopifnot(file.exists(router_file), length(parse(file = router_file)) > 0L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_TaskRouter_v1.R"))

plan1 <- plan_HCC_task_v1(input = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"), project_dir = PROJECT_DIR, 
    sample_uid = "PRJEB34724__S02-T-RNA")

plan2 <- plan_HCC_task_v1(input = file.path(PROJECT_DIR, "08_reports", "Stage4_1C_TEST_ONLY_16gene_TPM.csv"), 
    project_dir = PROJECT_DIR)

plan3 <- plan_HCC_task_v1(input = list(rna_file = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"), 
    mutation_status = c(SYNTHETIC_A = 1L, SYNTHETIC_B = 0L), mutation_patient_id = "PRJEB34724__S02-T-RNA"), 
    project_dir = PROJECT_DIR, sample_uid = "PRJEB34724__S02-T-RNA")

plan4 <- plan_HCC_task_v1(input = file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds"), 
    project_dir = PROJECT_DIR, mutation_gene = "CTNNB1")

plan5 <- plan_HCC_task_v1(input = file.path(PROJECT_DIR, "03_input", "NONEXISTENT_FILE.csv"), project_dir = PROJECT_DIR)

plans <- list(Full_RNA = plan1, Partial_RNA = plan2, Synthetic_Multiomics = plan3, TCGA_Cohort = plan4, 
    Missing_File = plan5)

summary_table <- do.call(rbind, lapply(names(plans), function(nm) {
    p <- plans[[nm]]
    data.frame(Test = nm, Detected_Mode = p$detected_mode, Plan_Status = p$status, stringsAsFactors = FALSE)
}))

for (nm in names(plans)) {
    cat("\nTEST:", nm, "")
    print(plans[[nm]]$capabilities, row.names = FALSE, right = FALSE)
    if (identical(plans[[nm]]$status, "blocked")) {
        print(plans[[nm]]$qc)
    }
}

stopifnot(identical(plan1$detected_mode, "patient_rna"), identical(plan1$status, "ready"), identical(plan2$detected_mode, 
    "patient_rna"), identical(plan2$status, "partial"), identical(plan3$detected_mode, "patient_multiomics"), 
    identical(plan3$status, "ready"), identical(plan4$detected_mode, "cohort_analysis"), identical(plan4$status, 
        "ready"), identical(plan5$status, "blocked"))

saveRDS(plans, file.path(PROJECT_DIR, "08_reports", "Stage4_2A_task_plans.rds"))

write.csv(summary_table, file.path(PROJECT_DIR, "08_reports", "Stage4_2A_routing_summary.csv"), row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"


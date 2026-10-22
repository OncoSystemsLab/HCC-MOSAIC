if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

score_HCC_modules_fixed_v1 <- function(tpm, reference) {
    stopifnot(is.matrix(tpm), is.numeric(tpm), !is.null(rownames(tpm)), !is.null(colnames(tpm)), !anyDuplicated(rownames(tpm)), 
        !anyDuplicated(colnames(tpm)), all(is.finite(tpm)), all(tpm >= 0))
    genes <- reference$required_genes
    missing_genes <- setdiff(genes, rownames(tpm))
    if (length(missing_genes) > 0L) {
        stop("Required module genes are missing: ", paste(missing_genes, collapse = ", "))
    }
    parameters <- reference$gene_parameters
    mu <- setNames(parameters$Mean_logTPM, parameters$Gene)
    sigma <- setNames(parameters$SD_logTPM, parameters$Gene)
    expr <- log2(tpm[genes, , drop = FALSE] + 1)
    z <- sweep(expr, 1, mu[genes], "-")
    z <- sweep(z, 1, sigma[genes], "/")
    scores <- do.call(cbind, lapply(reference$module_definitions, function(g) {
        colMeans(z[g, , drop = FALSE])
    }))
    colnames(scores) <- names(reference$module_definitions)
    rownames(scores) <- colnames(tpm)
    ref_scores <- reference$reference_module_scores
    percentiles <- matrix(NA_real_, nrow = nrow(scores), ncol = ncol(scores), dimnames = dimnames(scores))
    for (m in colnames(scores)) {
        r <- ref_scores[, m]
        for (i in seq_len(nrow(scores))) {
            v <- scores[i, m]
            percentiles[i, m] <- 100 * (sum(r < v) + 0.5 * sum(r == v))/length(r)
        }
    }
    stopifnot(all(is.finite(scores)), all(is.finite(percentiles)), all(percentiles >= 0), all(percentiles <= 
        100))
    list(scores = data.frame(Patient_ID = rownames(scores), scores, check.names = FALSE, row.names = NULL), 
        percentiles = data.frame(Patient_ID = rownames(percentiles), percentiles, check.names = FALSE, 
            row.names = NULL), reference_version = reference$version)
}

test_file <- file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv")

stopifnot(file.exists(test_file))

test_df <- read.csv(test_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(ncol(test_df) == 2L)

test_genes <- as.character(test_df[[1]])

test_values <- suppressWarnings(as.numeric(test_df[[2]]))

stopifnot(!anyDuplicated(test_genes), all(is.finite(test_values)), all(test_values >= 0))

test_tpm <- matrix(test_values, ncol = 1L, dimnames = list(test_genes, "PRJEB34724__S02-T-RNA"))

test_result <- score_HCC_modules_fixed_v1(tpm = test_tpm, reference = module_reference_v1)

uid <- "PRJEB34724__S02-T-RNA"

if (uid %in% rownames(ref_module_scores)) {
    expected <- ref_module_scores[uid, , drop = TRUE]
    observed <- as.numeric(test_result$scores[1, names(module_sets), drop = TRUE])
    error <- max(abs(observed - as.numeric(expected)))
    cat("\nS02 reference reproduction error:", error, "")
    stopifnot(error < 1e-08)
} else {
    cat("\nReference row ID does not match S02 UID;", "direct score reproduction check skipped.")
}

function_file <- file.path(PROJECT_DIR, "07_agent", "HCC_FunctionalScorer_v1.R")

dump(list = "score_HCC_modules_fixed_v1", file = function_file, envir = .GlobalEnv)

invisible(parse(file = function_file))

stopifnot(file.exists(reference_file), file.exists(function_file), nrow(test_result$scores) == 1L, ncol(test_result$scores) == 
    5L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_FunctionalScorer_v1.R"))

ref <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage3_3C_module_reference_69_v1.rds"))

modules <- names(ref$module_definitions)

stopifnot(length(modules) == 4L, length(ref$required_genes) == 16L, nrow(ref$reference_module_scores) == 
    69L, is.function(score_HCC_modules_fixed_v1))

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

all_results <- list()

qc_list <- list()

for (cohort in names(joint_data)) {
    cat("\n===== EXTERNAL FROZEN SCORING:", cohort, "=====")
    joint <- joint_data[[cohort]]
    tpm <- joint$mechanism_TPM
    ids <- joint$patient_ids
    stopifnot(identical(colnames(tpm), ids), all(ref$required_genes %in% rownames(tpm)), all(is.finite(tpm)), 
        all(tpm >= 0))
    status <- as.integer(joint$wnt_mutation["CTNNB1", ids])
    result <- score_HCC_modules_fixed_v1(tpm = tpm, reference = ref)
    scores <- result$scores
    percentiles <- result$percentiles
    stopifnot(identical(scores$Patient_ID, ids), identical(percentiles$Patient_ID, ids), all(modules %in% 
        names(scores)), all(modules %in% names(percentiles)))
    probe_idx <- unique(c(1L, ceiling(length(ids)/2), length(ids)))
    probe <- score_HCC_modules_fixed_v1(tpm = tpm[, probe_idx, drop = FALSE], reference = ref)
    batch_values <- as.matrix(scores[probe_idx, modules, drop = FALSE])
    single_values <- as.matrix(probe$scores[, modules, drop = FALSE])
    batch_p <- as.matrix(percentiles[probe_idx, modules, drop = FALSE])
    single_p <- as.matrix(probe$percentiles[, modules, drop = FALSE])
    max_score_error <- max(abs(batch_values - single_values))
    max_percentile_error <- max(abs(batch_p - single_p))
    stopifnot(max_score_error < 1e-12, max_percentile_error < 1e-12)
    old_scores <- read.csv(file.path(PROJECT_DIR, "04_processed", paste0("Stage3_3B_", cohort, "_module_scores.csv")), 
        check.names = FALSE, stringsAsFactors = FALSE)
    old_scores <- old_scores[match(ids, old_scores$Patient_ID), , drop = FALSE]
    stopifnot(identical(old_scores$Patient_ID, ids))
    one_qc <- do.call(rbind, lapply(modules, function(m) {
        s <- as.numeric(scores[[m]])
        p <- as.numeric(percentiles[[m]])
        reference_scores <- as.numeric(ref$reference_module_scores[, m])
        old <- as.numeric(old_scores[[m]])
        data.frame(Cohort = cohort, Module = m, N_Patients = length(s), Mean_FrozenScore = mean(s), SD_FrozenScore = sd(s), 
            Median_ReferencePercentile = median(p), N_Below_Reference_Min = sum(s < min(reference_scores)), 
            N_Above_Reference_Max = sum(s > max(reference_scores)), Fraction_Outside_Reference = mean(s < 
                min(reference_scores) | s > max(reference_scores)), Spearman_Old_vs_Frozen = unname(cor(s, 
                old, method = "spearman")), Mean_Difference_Mut_NoRecorded = mean(s[status == 1L]) - 
                mean(s[status == 0L]), Batch_Single_Max_Error = max_score_error, Percentile_Max_Error = max_percentile_error, 
            stringsAsFactors = FALSE)
    }))
    qc_list[[cohort]] <- one_qc
    pct_cols <- percentiles[, modules, drop = FALSE]
    names(pct_cols) <- paste0("ReferencePctl_", modules)
    patient_output <- data.frame(Cohort = cohort, Patient_ID = ids, CTNNB1_Recorded_Mutation = status, 
        scores[, modules, drop = FALSE], pct_cols, check.names = FALSE, stringsAsFactors = FALSE)
    all_results[[cohort]] <- patient_output
    write.csv(patient_output, file.path(PROJECT_DIR, "04_processed", paste0("Stage3_3C_", cohort, "_frozen_module_scores.csv")), 
        row.names = FALSE)
    cat("Patients:", length(ids), "")
    cat("Required genes:", length(ref$required_genes), "")
    cat("Single/batch error:", max_score_error, "")
    cat("Percentile error:", max_percentile_error, "")
    print(one_qc[, c("Module", "Mean_FrozenScore", "Median_ReferencePercentile", "Fraction_Outside_Reference", 
        "Spearman_Old_vs_Frozen", "Mean_Difference_Mut_NoRecorded")], digits = 4, row.names = FALSE)
}

qc_all <- do.call(rbind, qc_list)

rownames(qc_all) <- NULL

patient_all <- do.call(rbind, all_results)

rownames(patient_all) <- NULL

write.csv(qc_all, file.path(PROJECT_DIR, "08_reports", "Stage3_3C_external_frozen_score_QC.csv"), row.names = FALSE)

saveRDS(list(reference_version = ref$version, cohort_scores = all_results, qc = qc_all), file.path(PROJECT_DIR, 
    "04_processed", "Stage3_3C_external_frozen_scores_v1.rds"))

stopifnot(nrow(patient_all) == 457L, nrow(qc_all) == 8L, sum(patient_all$CTNNB1_Recorded_Mutation) == 
    133L, all(is.finite(qc_all$Mean_FrozenScore)), all(qc_all$Batch_Single_Max_Error < 1e-12), all(qc_all$Percentile_Max_Error < 
    1e-12))

PROJECT_DIR <- "HCC_OmicsAgent"

make_HCC_covariates_v1 <- function(joint) {
    ids <- as.character(joint$patient_ids)
    cl <- joint$clinical
    stopifnot(identical(as.character(cl$Patient_ID), ids))
    if (joint$cohort == "TCGA_LIHC") {
        age <- as.numeric(cl$age_at_index)
        raw_sex <- cl$sex_at_birth
    }
    else if (joint$cohort == "ICGC_LICA_FR") {
        age <- as.numeric(cl$Age)
        raw_sex <- cl$Gender
    }
    else {
        stop("Unsupported clinical schema: ", joint$cohort)
    }
    raw_sex <- toupper(trimws(as.character(raw_sex)))
    sex <- ifelse(raw_sex %in% c("M", "MALE"), "Male", ifelse(raw_sex %in% c("F", "FEMALE"), "Female", 
        "NA"))
    out <- data.frame(Patient_ID = ids, Age = age, Sex = factor(sex, levels = c("Female", "Male")), stringsAsFactors = FALSE)
    stopifnot(all(is.finite(out$Age)), all(out$Age > 0 & out$Age < 120), !anyNA(out$Sex), nlevels(droplevels(out$Sex)) == 
        2L)
    out
}

hcc_fit_association_v1 <- function(feature_matrix, mutation_matrix, covariates, mutation_gene = "CTNNB1", 
    min_group_n = 10L) {
    stopifnot(is.matrix(feature_matrix), is.numeric(feature_matrix), is.matrix(mutation_matrix), !is.null(rownames(feature_matrix)), 
        !is.null(colnames(feature_matrix)), !anyDuplicated(rownames(feature_matrix)), !anyDuplicated(colnames(feature_matrix)), 
        all(is.finite(feature_matrix)), all(c("Patient_ID", "Age", "Sex") %in% names(covariates)))
    ids <- colnames(feature_matrix)
    stopifnot(identical(as.character(covariates$Patient_ID), ids), all(ids %in% colnames(mutation_matrix)))
    if (!mutation_gene %in% rownames(mutation_matrix)) {
        stop("Mutation gene unavailable in input matrix: ", mutation_gene)
    }
    mutation <- as.integer(mutation_matrix[mutation_gene, ids])
    stopifnot(!anyNA(mutation), all(mutation %in% c(0L, 1L)))
    n_mut <- sum(mutation == 1L)
    n_ref <- sum(mutation == 0L)
    if (min(n_mut, n_ref) < min_group_n) {
        stop("Mutation group too small for analysis: ", mutation_gene, " (", n_mut, " vs ", n_ref, ")")
    }
    dat <- data.frame(Mutation = mutation, Age = as.numeric(covariates$Age), Sex = factor(as.character(covariates$Sex), 
        levels = c("Female", "Male")))
    stopifnot(all(is.finite(dat$Age)), !anyNA(dat$Sex))
    results <- lapply(rownames(feature_matrix), function(feature) {
        df <- dat
        df$Y <- as.numeric(feature_matrix[feature, ids])
        fit <- lm(Y ~ Mutation + Age + Sex, data = df)
        X <- model.matrix(fit)
        if (fit$rank != ncol(X)) {
            stop("Rank deficient: ", feature)
        }
        A <- solve(crossprod(X))
        h <- hatvalues(fit)
        r <- residuals(fit)/pmax(1 - h, 1e-10)
        V <- A %*% crossprod(X * r) %*% A
        beta <- unname(coef(fit)["Mutation"])
        se <- sqrt(V["Mutation", "Mutation"])
        if (!is.finite(se) || se <= 0) {
            stop("Invalid robust SE: ", feature)
        }
        df_res <- df.residual(fit)
        p_value <- 2 * pt(-abs(beta/se), df = df_res)
        crit <- qt(0.975, df_res)
        data.frame(Feature = feature, Mutation_Gene = mutation_gene, N = nrow(df), N_Mutation = n_mut, 
            N_NoRecorded = n_ref, Beta_Adjusted = beta, SE_HC3 = se, CI_Lower = beta - crit * se, CI_Upper = beta + 
                crit * se, P_Value = p_value, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, results)
    out$FDR_BH <- p.adjust(out$P_Value, method = "BH")
    rownames(out) <- NULL
    out
}

engine_file <- file.path(PROJECT_DIR, "07_agent", "HCC_AssociationEngine_v1.R")

dump(c("make_HCC_covariates_v1", "hcc_fit_association_v1"), file = engine_file, envir = .GlobalEnv)

stopifnot(file.exists(engine_file))

test_file <- file.path(PROJECT_DIR, "07_agent", "Stage3_3C3_cold_start_test.R")

test_lines <- c("PROJECT_DIR <- \"HCC_OmicsAgent\"", "source(file.path(PROJECT_DIR, \"07_agent\",", "  \"HCC_FunctionalScorer_v1.R\"))", 
    "source(file.path(PROJECT_DIR, \"07_agent\",", "  \"HCC_AssociationEngine_v1.R\"))", "ref <- readRDS(file.path(PROJECT_DIR, \"00_config\",", 
    "  \"Stage3_3C_module_reference_69_v1.rds\"))", "previous <- read.csv(file.path(PROJECT_DIR, \"05_results\",", 
    "  \"Stage3_3B2_module_associations.csv\"))", "modules <- names(ref$module_definitions)", "for (cohort in c(\"TCGA_LIHC\", \"ICGC_LICA_FR\")) {", 
    "  cat(\"\\nTEST COHORT:\", cohort, \"\\n\")", "  joint <- readRDS(file.path(PROJECT_DIR,", "    \"04_processed\",", 
    "    paste0(\"Stage3_2A_\", cohort, \"_joint_v1.rds\")))", "  old <- read.csv(file.path(PROJECT_DIR,", 
    "    \"04_processed\",", "    paste0(\"Stage3_3B_\", cohort, \"_module_scores.csv\")),", "    check.names = FALSE)", 
    "  stopifnot(identical(as.character(old$Patient_ID),", "    joint$patient_ids))", "  x <- t(as.matrix(old[, modules, drop = FALSE]))", 
    "  colnames(x) <- old$Patient_ID", "  cov <- make_HCC_covariates_v1(joint)", "  fit <- hcc_fit_association_v1(", 
    "    feature_matrix = x,", "    mutation_matrix = joint$wnt_mutation,", "    covariates = cov,", 
    "    mutation_gene = \"CTNNB1\")", "  expected <- subset(previous, Cohort == cohort)", "  expected <- expected[", 
    "    match(fit$Feature, expected$Module), ]", "  stopifnot(", "    !anyNA(expected$Module),", "    identical(as.character(fit$Feature),", 
    "      as.character(expected$Module)))", "  association_error <- max(", "    abs(fit$Beta_Adjusted - expected$Beta_Adjusted),", 
    "    abs(fit$P_Value - expected$P_Value),", "    abs(fit$FDR_BH - expected$FDR_4))", "  stopifnot(association_error < 1e-9)", 
    "  axin1 <- hcc_fit_association_v1(", "    feature_matrix = x,", "    mutation_matrix = joint$wnt_mutation,", 
    "    covariates = cov,", "    mutation_gene = \"AXIN1\")", "  stopifnot(nrow(axin1) == length(modules))", 
    "", "  batch <- score_HCC_modules_fixed_v1(", "    joint$mechanism_TPM, ref)", "  saved <- read.csv(file.path(PROJECT_DIR,", 
    "    \"04_processed\",", "    paste0(\"Stage3_3C_\", cohort,", "      \"_frozen_module_scores.csv\")),", 
    "    check.names = FALSE)", "  stopifnot(identical(batch$scores$Patient_ID,", "    as.character(saved$Patient_ID)))", 
    "  frozen_error <- max(abs(", "    as.matrix(batch$scores[, modules]) -", "    as.matrix(saved[, modules])))", 
    "  stopifnot(frozen_error < 1e-12)", "", "  probe <- unique(c(1L, ceiling(length(joint$patient_ids)/2),", 
    "    length(joint$patient_ids)))", "  single_errors <- vapply(probe, function(i) {", "    single <- score_HCC_modules_fixed_v1(", 
    "      joint$mechanism_TPM[, i, drop = FALSE], ref)", "    max(abs(", "      as.numeric(single$scores[1, modules]) -", 
    "      as.numeric(batch$scores[i, modules])))", "  }, numeric(1))", "  stopifnot(all(single_errors < 1e-12))", 
    "", "  bad_gene <- try(", "    hcc_fit_association_v1(", "      x, joint$wnt_mutation, cov,", "      mutation_gene = \"UNAVAILABLE_GENE\"),", 
    "    silent = TRUE)", "  stopifnot(inherits(bad_gene, \"try-error\"))", "  cat(\"Association reproduction error:\",", 
    "    association_error, \"\\n\")", "  cat(\"Frozen score reproduction error:\",", "    frozen_error, \"\\n\")", 
    "  cat(\"Single patient max error:\",", "    max(single_errors), \"\\n\")", "  cat(\"AXIN1 mutated patients:\",", 
    "    unique(axin1$N_Mutation), \"\\n\")", "  cat(cohort, \": PASS\\n\")", "}", "cat(\"\\nStage 3.3C.3 COLD_START: PASS\\n\")")

writeLines(test_lines, con = test_file, useBytes = TRUE)

stopifnot(file.exists(test_file), length(parse(file = test_file)) > 0L)

rscript_bin <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") {
    "Rscript.exe"
} else {
    "Rscript"
})

stopifnot(file.exists(rscript_bin))

cold_output <- suppressWarnings(system2(command = rscript_bin, args = c("--vanilla", shQuote(test_file)), 
    stdout = TRUE, stderr = TRUE))

exit_status <- attr(cold_output, "status")

if (is.null(exit_status)) {
    exit_status <- 0L
}

log_file <- file.path(PROJECT_DIR, "08_reports", "Stage3_3C3_cold_start_log.txt")

writeLines(cold_output, con = log_file, useBytes = TRUE)

if (exit_status != 0L || !any(grepl("Stage 3.3C.3 COLD_START: PASS", cold_output, fixed = TRUE))) {
    stop("Cold-start acceptance failed. Check log: ", log_file)
}


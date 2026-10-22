if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

stopifnot(nrow(plot_df) == 22L, all(is.finite(plot_df$Delta_Beta)), all(file.exists(file.path(PROJECT_DIR, 
    "06_figures", "Stage3_2D_CTNNB1_model_stability.png"))))

PROJECT_DIR <- "HCC_OmicsAgent"

joint_data <- list(TCGA = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

genes_tcga <- rownames(joint_data$TCGA$mechanism_TPM)

genes_icgc <- rownames(joint_data$ICGC$mechanism_TPM)

stopifnot(length(genes_tcga) == 63L, length(genes_icgc) == 63L, setequal(genes_tcga, genes_icgc))

available_genes <- intersect(genes_tcga, genes_icgc)

module_candidates <- list(WNT_Response = c("NOTUM", "AXIN2", "NKD1", "TBX3", "LGR5", "GLUL"), BA_Classical_CA_Branch = c("CYP7A1", 
    "CYP8B1", "HSD3B7", "AKR1D1"), BA_Alternative_Synthesis = c("CYP27A1", "CYP7B1", "CYP39A1"), BA_Conjugation_Activation = c("BAAT", 
    "SLC27A5"), BA_Hepatic_Uptake = c("SLC10A1", "SLCO1B1", "SLCO1B3"), BA_Canalicular_Export = c("ABCB11", 
    "ABCC2"), BA_Basolateral_Export = c("ABCC3", "ABCC4"), BA_FXR_Feedback = c("NR1H4", "NR0B2", "FGF19", 
    "FGFR4", "KLB"), Immune_Cytotoxic = c("CD8A", "GZMB", "PRF1", "NKG7", "IFNG", "GNLY"))

stopifnot(!anyDuplicated(names(module_candidates)), all(lengths(module_candidates) >= 2L))

module_membership <- do.call(rbind, lapply(names(module_candidates), function(nm) {
    genes <- unique(module_candidates[[nm]])
    data.frame(Module = nm, Gene = genes, Present_in_63_Panel = genes %in% available_genes, stringsAsFactors = FALSE)
}))

rownames(module_membership) <- NULL

module_summary <- do.call(rbind, lapply(names(module_candidates), function(nm) {
    genes <- unique(module_candidates[[nm]])
    matched <- intersect(genes, available_genes)
    missing <- setdiff(genes, available_genes)
    data.frame(Module = nm, N_Candidate = length(genes), N_Matched = length(matched), N_Missing = length(missing), 
        Coverage_Pct = round(100 * length(matched)/length(genes), 2), Ready_for_Exploratory_Score = length(matched) >= 
            2L, Matched_Genes = paste(matched, collapse = ";"), Missing_Genes = paste(missing, collapse = ";"), 
        stringsAsFactors = FALSE)
}))

rownames(module_summary) <- NULL

original_panel_summary <- data.frame(Panel = names(config$gene_panels), N_Genes = lengths(config$gene_panels), 
    stringsAsFactors = FALSE)

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(module_summary, file.path(report_dir, "Stage3_3A_BA_module_coverage.csv"), row.names = FALSE)

write.csv(module_membership, file.path(report_dir, "Stage3_3A_BA_module_membership.csv"), row.names = FALSE)

saveRDS(module_candidates, file.path(report_dir, "Stage3_3A_BA_module_candidates_v0.rds"))

stopifnot(nrow(module_summary) == 9L, length(available_genes) == 63L, file.exists(file.path(report_dir, 
    "Stage3_3A_BA_module_coverage.csv")))

PROJECT_DIR <- "HCC_OmicsAgent"

RAW_DIR <- "raw_data"

modules <- readRDS(file.path(PROJECT_DIR, "08_reports", "Stage3_3A_BA_module_candidates_v0.rds"))

joint_data <- list(TCGA = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

stopifnot(length(modules) == 9L)

gene_qc_list <- list()

module_qc_list <- list()

for (cohort in names(joint_data)) {
    tpm <- joint_data[[cohort]]$mechanism_TPM
    stopifnot(nrow(tpm) == 63L, all(is.finite(tpm)), all(tpm >= 0))
    for (module_name in names(modules)) {
        candidates <- unique(modules[[module_name]])
        present <- intersect(candidates, rownames(tpm))
        missing <- setdiff(candidates, rownames(tpm))
        expr_tpm <- tpm[present, , drop = FALSE]
        expr_log <- log2(expr_tpm + 1)
        gene_qc <- data.frame(Cohort = cohort, Module = module_name, Gene = present, Median_TPM = apply(expr_tpm, 
            1, median), Zero_Fraction = rowMeans(expr_tpm == 0), SD_logTPM = apply(expr_log, 1, sd), 
            IQR_logTPM = apply(expr_log, 1, IQR), stringsAsFactors = FALSE)
        rownames(gene_qc) <- NULL
        gene_qc_list[[paste(cohort, module_name, sep = "__")]] <- gene_qc
        usable <- gene_qc$Gene[is.finite(gene_qc$SD_logTPM) & gene_qc$SD_logTPM > 1e-08]
        pair_rho <- numeric(0)
        if (length(usable) >= 2L) {
            mat <- expr_log[usable, , drop = FALSE]
            cor_mat <- suppressWarnings(cor(t(mat), method = "spearman"))
            pair_rho <- cor_mat[upper.tri(cor_mat)]
            pair_rho <- pair_rho[is.finite(pair_rho)]
        }
        median_rho <- if (length(pair_rho) > 0L) {
            median(pair_rho)
        }
        else {
            NA_real_
        }
        negative_fraction <- if (length(pair_rho) > 0L) {
            mean(pair_rho < 0)
        }
        else {
            NA_real_
        }
        module_qc <- data.frame(Cohort = cohort, Module = module_name, N_Candidates = length(candidates), 
            N_Present = length(present), N_Missing = length(missing), N_Variable_Genes = length(usable), 
            N_Genes_Zero_Over90 = sum(gene_qc$Zero_Fraction > 0.9), N_Pairs = length(pair_rho), Median_Spearman_Rho = median_rho, 
            Negative_Pair_Fraction = negative_fraction, Missing_Genes = paste(missing, collapse = ";"), 
            stringsAsFactors = FALSE)
        module_qc_list[[paste(cohort, module_name, sep = "__")]] <- module_qc
    }
}

gene_qc_table <- do.call(rbind, gene_qc_list)

module_qc_table <- do.call(rbind, module_qc_list)

rownames(gene_qc_table) <- NULL

rownames(module_qc_table) <- NULL

low_genes <- subset(gene_qc_table, Zero_Fraction > 0.9 | SD_logTPM <= 1e-08)

find_gene_in_gz <- function(path, target) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    readLines(con, n = 1L, warn = FALSE)
    n_found <- 0L
    repeat {
        lines <- readLines(con, n = 3000L, warn = FALSE)
        if (length(lines) == 0L) 
            break
        genes <- sub("\t.*$", "", lines)
        genes <- trimws(gsub("^\"|\"$", "", genes))
        n_found <- n_found + sum(genes == target)
    }
    n_found
}

raw_paths <- c(TCGA = file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_expression_TPM.tsv.gz"), ICGC = file.path(RAW_DIR, 
    "ICGC_LICA_FR_FINAL", "ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz"))

stopifnot(all(file.exists(raw_paths)))

cyp39a1_audit <- data.frame(Cohort = names(raw_paths), Gene = "CYP39A1", N_Rows_Found = vapply(raw_paths, 
    find_gene_in_gz, integer(1), target = "CYP39A1"), stringsAsFactors = FALSE)

cyp39a1_audit$Present_In_Raw_TPM <- cyp39a1_audit$N_Rows_Found > 0L

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(module_qc_table, file.path(report_dir, "Stage3_3A2_module_coexpression_QC.csv"), row.names = FALSE)

write.csv(gene_qc_table, file.path(report_dir, "Stage3_3A2_module_gene_expression_QC.csv"), row.names = FALSE)

write.csv(cyp39a1_audit, file.path(report_dir, "Stage3_3A2_CYP39A1_raw_TPM_audit.csv"), row.names = FALSE)

stopifnot(nrow(module_qc_table) == 18L, all(module_qc_table$N_Present >= 2L), nrow(cyp39a1_audit) == 
    2L)

PROJECT_DIR <- "HCC_OmicsAgent"

modules_v0 <- readRDS(file.path(PROJECT_DIR, "08_reports", "Stage3_3A_BA_module_candidates_v0.rds"))

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

selected_modules <- c("WNT_Response", "BA_Conjugation_Activation", "Immune_Cytotoxic", "BA_Canalicular_Export")

stopifnot(all(selected_modules %in% names(modules_v0)))

module_sets <- modules_v0[selected_modules]

score_module <- function(expr_log, genes) {
    stopifnot(all(genes %in% rownames(expr_log)), length(genes) >= 2L)
    x <- expr_log[genes, , drop = FALSE]
    gene_mean <- rowMeans(x)
    gene_sd <- apply(x, 1, sd)
    if (any(!is.finite(gene_sd)) || any(gene_sd <= 1e-08)) {
        stop("Module genes have zero variance: ", paste(genes[gene_sd <= 1e-08], collapse = ", "))
    }
    z <- sweep(x, 1, gene_mean, "-")
    z <- sweep(z, 1, gene_sd, "/")
    score <- colMeans(z)
    stopifnot(all(is.finite(score)), length(score) == ncol(x))
    list(score = score, gene_mean = gene_mean, gene_sd = gene_sd)
}

module_score_data <- list()

module_score_qc <- list()

gene_reference <- list()

for (cohort in names(joint_data)) {
    joint <- joint_data[[cohort]]
    ids <- joint$patient_ids
    expr <- joint$mechanism_logTPM
    mutation <- as.integer(joint$wnt_mutation["CTNNB1", ids])
    stopifnot(identical(colnames(expr), ids), all(mutation %in% c(0L, 1L)), nrow(expr) == 63L)
    score_mat <- matrix(NA_real_, nrow = length(ids), ncol = length(module_sets), dimnames = list(ids, 
        names(module_sets)))
    for (module_name in names(module_sets)) {
        genes <- module_sets[[module_name]]
        result <- score_module(expr_log = expr, genes = genes)
        score_mat[, module_name] <- result$score[ids]
        gene_reference[[paste(cohort, module_name, sep = "__")]] <- data.frame(Cohort = cohort, Module = module_name, 
            Gene = genes, Mean_logTPM = as.numeric(result$gene_mean[genes]), SD_logTPM = as.numeric(result$gene_sd[genes]), 
            stringsAsFactors = FALSE)
    }
    stopifnot(all(is.finite(score_mat)), max(abs(colMeans(score_mat))) < 1e-08, all(apply(score_mat, 
        2, sd) > 0))
    patient_scores <- data.frame(Cohort = cohort, Patient_ID = ids, CTNNB1_Recorded_Mutation = mutation, 
        score_mat, check.names = FALSE, row.names = NULL)
    module_score_data[[cohort]] <- patient_scores
    one_qc <- do.call(rbind, lapply(names(module_sets), function(nm) {
        s <- score_mat[, nm]
        data.frame(Cohort = cohort, Module = nm, N_Patients = length(s), N_Genes = length(module_sets[[nm]]), 
            Score_Mean = mean(s), Score_SD = sd(s), Score_Min = min(s), Score_Max = max(s), Mean_Mutation = mean(s[mutation == 
                1L]), Mean_NoRecorded = mean(s[mutation == 0L]), Mean_Difference = mean(s[mutation == 
                1L]) - mean(s[mutation == 0L]), stringsAsFactors = FALSE)
    }))
    module_score_qc[[cohort]] <- one_qc
    cat("\n===== MODULE SCORES:", cohort, "=====")
    cat("Patients:", length(ids), "| Mut:", sum(mutation == 1L), "| No recorded:", sum(mutation == 0L), 
        "")
    print(one_qc[, c("Module", "N_Genes", "Score_SD", "Mean_Mutation", "Mean_NoRecorded", "Mean_Difference")], 
        digits = 4, row.names = FALSE)
}

qc_all <- do.call(rbind, module_score_qc)

reference_all <- do.call(rbind, gene_reference)

rownames(qc_all) <- NULL

rownames(reference_all) <- NULL

outdir <- file.path(PROJECT_DIR, "04_processed")

report_dir <- file.path(PROJECT_DIR, "08_reports")

for (cohort in names(module_score_data)) {
    write.csv(module_score_data[[cohort]], file.path(outdir, paste0("Stage3_3B_", cohort, "_module_scores.csv")), 
        row.names = FALSE)
}

saveRDS(list(module_definitions = module_sets, cohort_scores = module_score_data, gene_standardization = reference_all, 
    note = paste("Exploratory cohort-standardized expression scores;", "not ssGSEA or measured pathway activity.")), 
    file.path(outdir, "Stage3_3B_module_scores_v1.rds"))

write.csv(qc_all, file.path(report_dir, "Stage3_3B_module_score_QC.csv"), row.names = FALSE)

write.csv(reference_all, file.path(report_dir, "Stage3_3B_gene_standardization.csv"), row.names = FALSE)

stopifnot(nrow(module_score_data$TCGA_LIHC) == 361L, nrow(module_score_data$ICGC_LICA_FR) == 96L, nrow(qc_all) == 
    8L, length(module_sets) == 4L, all(is.finite(qc_all$Mean_Difference)))

PROJECT_DIR <- "HCC_OmicsAgent"

score_obj <- readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_3B_module_scores_v1.rds"))

module_names <- names(score_obj$module_definitions)

joint_paths <- c(TCGA_LIHC = "Stage3_2A_TCGA_LIHC_joint_v1.rds", ICGC_LICA_FR = "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")

association_results <- list()

for (cohort in names(joint_paths)) {
    joint <- readRDS(file.path(PROJECT_DIR, "04_processed", joint_paths[[cohort]]))
    scores <- score_obj$cohort_scores[[cohort]]
    clinical <- joint$clinical
    ids <- joint$patient_ids
    stopifnot(identical(as.character(scores$Patient_ID), ids), identical(as.character(clinical$Patient_ID), 
        ids), all(module_names %in% names(scores)))
    age <- if (cohort == "TCGA_LIHC") {
        as.numeric(clinical$age_at_index)
    }
    else {
        as.numeric(clinical$Age)
    }
    sex_raw <- if (cohort == "TCGA_LIHC") {
        clinical$sex_at_birth
    }
    else {
        clinical$Gender
    }
    sex_raw <- toupper(trimws(as.character(sex_raw)))
    sex <- ifelse(sex_raw %in% c("M", "MALE"), "Male", ifelse(sex_raw %in% c("F", "FEMALE"), "Female", 
        "NA"))
    df <- data.frame(Patient_ID = ids, Mutation = as.integer(joint$wnt_mutation["CTNNB1", ids]), Age = age, 
        Sex = factor(sex, levels = c("Female", "Male")))
    stopifnot(all(is.finite(df$Age)), !anyNA(df$Sex), all(df$Mutation %in% c(0L, 1L)), identical(df$Mutation, 
        as.integer(scores$CTNNB1_Recorded_Mutation)))
    one <- do.call(rbind, lapply(module_names, function(m) {
        dat <- df
        dat$Y <- as.numeric(scores[[m]])
        fit <- lm(Y ~ Mutation + Age + Sex, data = dat)
        X <- model.matrix(fit)
        stopifnot(fit$rank == ncol(X))
        A <- solve(crossprod(X))
        h <- hatvalues(fit)
        r_hc3 <- residuals(fit)/pmax(1 - h, 1e-10)
        V <- A %*% crossprod(X * r_hc3) %*% A
        beta <- unname(coef(fit)["Mutation"])
        se <- sqrt(V["Mutation", "Mutation"])
        stopifnot(is.finite(se), se > 0)
        df_res <- df.residual(fit)
        tcrit <- qt(0.975, df_res)
        pval <- 2 * pt(-abs(beta/se), df = df_res)
        data.frame(Cohort = cohort, Module = m, N = nrow(dat), N_Mutation = sum(dat$Mutation == 1L), 
            Beta_Adjusted = beta, SE_HC3 = se, CI_Lower = beta - tcrit * se, CI_Upper = beta + tcrit * 
                se, P_Value = pval, stringsAsFactors = FALSE)
    }))
    one$FDR_4 <- p.adjust(one$P_Value, method = "BH")
    association_results[[cohort]] <- one
}

association_table <- do.call(rbind, association_results)

rownames(association_table) <- NULL

keep <- c("Module", "Beta_Adjusted", "CI_Lower", "CI_Upper", "P_Value", "FDR_4")

cross_module <- merge(association_results$TCGA_LIHC[, keep], association_results$ICGC_LICA_FR[, keep], 
    by = "Module", suffixes = c("_TCGA", "_ICGC"))

cross_module$Same_Direction <- with(cross_module, sign(Beta_Adjusted_TCGA) == sign(Beta_Adjusted_ICGC) & 
    Beta_Adjusted_TCGA != 0 & Beta_Adjusted_ICGC != 0)

cross_module$Strict_Replication <- with(cross_module, Same_Direction & FDR_4_TCGA < 0.05 & FDR_4_ICGC < 
    0.05)

write.csv(association_table, file.path(PROJECT_DIR, "05_results", "Stage3_3B2_module_associations.csv"), 
    row.names = FALSE)

write.csv(cross_module, file.path(PROJECT_DIR, "05_results", "Stage3_3B2_module_replication.csv"), row.names = FALSE)

stopifnot(nrow(association_table) == 8L, nrow(cross_module) == 4L, all(is.finite(association_table$Beta_Adjusted)), 
    all(is.finite(association_table$P_Value)))

PROJECT_DIR <- "HCC_OmicsAgent"

module_obj <- readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_3B_module_scores_v1.rds"))

ref_obj <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))

module_sets <- module_obj$module_definitions

stopifnot(length(module_sets) == 4L, !anyDuplicated(names(module_sets)))

required_genes <- unique(unlist(module_sets, use.names = FALSE))

ref_values <- ref_obj$reference_values

stopifnot(nrow(ref_values) == 69L, all(required_genes %in% colnames(ref_values)))

ref_expr <- as.matrix(ref_values[, required_genes, drop = FALSE])

storage.mode(ref_expr) <- "double"

stopifnot(all(is.finite(ref_expr)))

ref_mean <- colMeans(ref_expr)

ref_sd <- apply(ref_expr, 2, sd)

stopifnot(all(is.finite(ref_mean)), all(is.finite(ref_sd)), all(ref_sd > 1e-08))

gene_parameters <- data.frame(Gene = required_genes, Mean_logTPM = unname(ref_mean[required_genes]), 
    SD_logTPM = unname(ref_sd[required_genes]), stringsAsFactors = FALSE)

ref_z <- sweep(ref_expr, 2, ref_mean, "-")

ref_z <- sweep(ref_z, 2, ref_sd, "/")

ref_module_scores <- do.call(cbind, lapply(module_sets, function(genes) {
    rowMeans(ref_z[, genes, drop = FALSE])
}))

colnames(ref_module_scores) <- names(module_sets)

rownames(ref_module_scores) <- rownames(ref_expr)

stopifnot(nrow(ref_module_scores) == 69L, ncol(ref_module_scores) == 4L, all(is.finite(ref_module_scores)), 
    max(abs(colMeans(ref_module_scores))) < 1e-10)

module_reference_v1 <- list(version = "Stage3_ModuleReference_v1", reference_name = "HCC_Development_69", 
    reference_type = "Exploratory_Development_Cohorts", preprocessing = "log2(TPM+1)", scoring = "Mean_Gene_Zscore", 
    module_definitions = module_sets, required_genes = required_genes, gene_parameters = gene_parameters, 
    reference_module_scores = ref_module_scores, n_reference = 69L, note = paste("Fixed development-cohort expression reference.", 
        "Not clinical norms or measured pathway activity.", "Do not recalibrate on a single incoming patient."))

reference_file <- file.path(PROJECT_DIR, "00_config", "Stage3_3C_module_reference_69_v1.rds")

saveRDS(module_reference_v1, reference_file)


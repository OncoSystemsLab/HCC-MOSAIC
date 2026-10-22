if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

wnt_genes <- c("CTNNB1", "AXIN1", "APC", "AXIN2", "RNF43", "ZNRF3")

stopifnot(length(mechanism_genes) == 63L)

clean_id <- function(x) {
    x <- trimws(x)
    x <- sub("^﻿", "", x)
    sub("^\"(.*)\"$", "\\1", x)
}

extract_gene_matrix <- function(path, target_genes, data_type = c("TPM", "Mutation")) {
    data_type <- match.arg(data_type)
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    header <- clean_id(strsplit(readLines(con, n = 1L, warn = FALSE), "\t", fixed = TRUE)[[1]])
    stopifnot(header[1] == "Gene", length(header) > 1L, !anyDuplicated(header[-1]))
    patient_ids <- header[-1]
    out <- matrix(NA_real_, nrow = length(target_genes), ncol = length(patient_ids), dimnames = list(target_genes, 
        patient_ids))
    found <- rep(FALSE, length(target_genes))
    names(found) <- target_genes
    column_totals <- setNames(numeric(length(patient_ids)), patient_ids)
    total_gene_rows <- 0L
    repeat {
        lines <- readLines(con, n = 1000L, warn = FALSE)
        if (length(lines) == 0L) 
            break
        for (line in lines) {
            fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
            if (length(fields) != length(header)) {
                stop("Inconsistent field count in data rows: ", basename(path), ", line ", total_gene_rows + 
                  2L)
            }
            gene <- clean_id(fields[1])
            values <- suppressWarnings(as.numeric(clean_id(fields[-1])))
            if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
                stop("Invalid numeric value: ", basename(path), ", Gene=", gene)
            }
            if (data_type == "Mutation" && any(!(values %in% c(0, 1)))) {
                stop("Mutation data are not strictly binary (0/1): ", gene)
            }
            if (data_type == "TPM") {
                column_totals <- column_totals + values
            }
            j <- match(gene, target_genes, nomatch = 0L)
            if (j > 0L) {
                if (found[j]) {
                  stop("Duplicate target gene: ", gene)
                }
                out[j, ] <- values
                found[j] <- TRUE
            }
            total_gene_rows <- total_gene_rows + 1L
        }
    }
    if (!all(found)) {
        stop("Missing target genes: ", paste(target_genes[!found], collapse = ", "))
    }
    list(matrix = out, column_totals = column_totals, total_gene_rows = total_gene_rows, n_patients = length(patient_ids))
}

cohort_files <- list(TCGA_LIHC = list(rna = file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_expression_TPM.tsv.gz"), 
    mutation = file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_mutation_matrix.tsv.gz"), clinical = file.path(RAW_DIR, 
        "TCGA_LIHC_final", "TCGA_LIHC_clinical.tsv")), ICGC_LICA_FR = list(rna = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz"), mutation = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"), clinical = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_clinical_harmonized.tsv")))

joint_data <- list()

for (cohort_name in names(cohort_files)) {
    cat("\n===== BUILD:", cohort_name, "=====")
    f <- cohort_files[[cohort_name]]
    rna <- extract_gene_matrix(path = f$rna, target_genes = mechanism_genes, data_type = "TPM")
    mutation <- extract_gene_matrix(path = f$mutation, target_genes = wnt_genes, data_type = "Mutation")
    clinical <- read.delim(f$clinical, check.names = FALSE, stringsAsFactors = FALSE)
    clinical$Patient_ID <- clean_id(as.character(clinical$Patient_ID))
    patient_ids <- colnames(rna$matrix)
    stopifnot(identical(patient_ids, colnames(mutation$matrix)), !anyDuplicated(clinical$Patient_ID), 
        all(patient_ids %in% clinical$Patient_ID))
    clinical <- clinical[match(patient_ids, clinical$Patient_ID), , drop = FALSE]
    stopifnot(identical(patient_ids, clinical$Patient_ID))
    joint <- list(cohort = cohort_name, patient_ids = patient_ids, mechanism_TPM = rna$matrix, mechanism_logTPM = log2(rna$matrix + 
        1), wnt_mutation = mutation$matrix, clinical = clinical, TPM_column_sums = rna$column_totals, 
        RNA_gene_rows = rna$total_gene_rows, Mutation_gene_rows = mutation$total_gene_rows, source_files = f)
    joint_data[[cohort_name]] <- joint
    saveRDS(joint, file.path(PROJECT_DIR, "04_processed", paste0("Stage3_2A_", cohort_name, "_joint_v1.rds")))
    cat("Patients:", length(patient_ids), "", "Mechanism genes:", nrow(joint$mechanism_TPM), "", "Wnt mutation genes:", 
        nrow(joint$wnt_mutation), "", "CTNNB1 recorded mutations:", sum(joint$wnt_mutation["CTNNB1", 
            ] == 1), "")
}

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("joint_data"))

qc_tables <- list()

patient_status <- list()

for (cohort_name in names(joint_data)) {
    x <- joint_data[[cohort_name]]
    tpm <- x$mechanism_TPM
    mut <- x$wnt_mutation
    ids <- x$patient_ids
    stopifnot(identical(colnames(tpm), ids), identical(colnames(mut), ids), identical(x$clinical$Patient_ID, 
        ids), nrow(tpm) == 63L, nrow(mut) == 6L, all(is.finite(tpm)), all(tpm >= 0), all(mut %in% c(0, 
        1)))
    sums <- x$TPM_column_sums
    outside_nominal <- (sums < 9e+05 | sums > 1100000)
    qc_tables[[cohort_name]] <- data.frame(Cohort = cohort_name, Patients = length(ids), Mechanism_Genes = nrow(tpm), 
        Mutation_Genes = nrow(mut), TPM_Sum_Min = min(sums), TPM_Sum_Median = median(sums), TPM_Sum_Max = max(sums), 
        N_TPM_Sum_Outside_0_9_1_1M = sum(outside_nominal), CTNNB1_Recorded_Mutation = sum(mut["CTNNB1", 
            ] == 1), CTNNB1_No_Recorded_Mutation = sum(mut["CTNNB1", ] == 0))
    patient_status[[cohort_name]] <- data.frame(Cohort = cohort_name, Patient_ID = ids, CTNNB1_Recorded_Mutation = as.integer(mut["CTNNB1", 
        ]), Other_WNT_Recorded_Mutation = as.integer(colSums(mut[setdiff(rownames(mut), "CTNNB1"), , 
        drop = FALSE]) > 0), TPM_Column_Sum = as.numeric(sums[ids]), stringsAsFactors = FALSE)
}

qc_summary <- do.call(rbind, qc_tables)

rownames(qc_summary) <- NULL

status_table <- do.call(rbind, patient_status)

rownames(status_table) <- NULL

write.csv(qc_summary, file.path(PROJECT_DIR, "08_reports", "Stage3_2A_joint_dataset_QC.csv"), row.names = FALSE)

write.csv(status_table, file.path(PROJECT_DIR, "04_processed", "Stage3_2A_patient_mutation_status.csv"), 
    row.names = FALSE)

stopifnot(nrow(status_table) == 457L, sum(status_table$CTNNB1_Recorded_Mutation) == 133L)

PROJECT_DIR <- "HCC_OmicsAgent"

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

compare_gene <- function(gene, expression, tpm, status) {
    x_mut <- as.numeric(expression[gene, status == 1])
    x_ref <- as.numeric(expression[gene, status == 0])
    stopifnot(length(x_mut) >= 2, length(x_ref) >= 2, all(is.finite(c(x_mut, x_ref))))
    p_value <- wilcox.test(x_mut, x_ref, exact = FALSE)$p.value
    n1 <- length(x_mut)
    n0 <- length(x_ref)
    combined <- c(x_mut, x_ref)
    ranks <- rank(combined, ties.method = "average")
    U <- sum(ranks[seq_len(n1)]) - n1 * (n1 + 1)/2
    cliff_delta <- 2 * U/(n1 * n0) - 1
    data.frame(Gene = gene, N_Mutation = n1, N_No_Recorded_Mutation = n0, Median_logTPM_Mutation = median(x_mut), 
        Median_logTPM_NoMutation = median(x_ref), Median_Difference_logTPM = median(x_mut) - median(x_ref), 
        Median_TPM_Mutation = median(tpm[gene, status == 1]), Median_TPM_NoMutation = median(tpm[gene, 
            status == 0]), Cliff_Delta = cliff_delta, P_Value = p_value, stringsAsFactors = FALSE)
}

association_list <- list()

for (cohort_name in names(joint_data)) {
    x <- joint_data[[cohort_name]]
    expr <- x$mechanism_logTPM
    tpm <- x$mechanism_TPM
    ids <- x$patient_ids
    status <- as.integer(x$wnt_mutation["CTNNB1", ids])
    stopifnot(identical(colnames(expr), ids), identical(colnames(tpm), ids), length(status) == length(ids), 
        all(status %in% c(0L, 1L)), nrow(expr) == 63L, all(is.finite(expr)))
    res <- do.call(rbind, lapply(rownames(expr), compare_gene, expression = expr, tpm = tpm, status = status))
    res$FDR_63 <- p.adjust(res$P_Value, method = "BH")
    res$Cohort <- cohort_name
    res <- res[order(res$FDR_63, res$P_Value), , drop = FALSE]
    rownames(res) <- NULL
    association_list[[cohort_name]] <- res
    cat("\n===== GENE ASSOCIATION:", cohort_name, "=====")
    cat("Patients:", length(ids), "")
    cat("CTNNB1 mutation:", sum(status == 1), "")
    cat("FDR < 0.05:", sum(res$FDR_63 < 0.05), "")
    cat("\nTOP 15:")
    print(head(res[, c("Gene", "Median_Difference_logTPM", "Cliff_Delta", "P_Value", "FDR_63")], 15), 
        digits = 5, row.names = FALSE)
}

tcga <- association_list$TCGA_LIHC

icgc <- association_list$ICGC_LICA_FR

cross_cohort <- merge(tcga[, c("Gene", "Median_Difference_logTPM", "Cliff_Delta", "P_Value", "FDR_63")], 
    icgc[, c("Gene", "Median_Difference_logTPM", "Cliff_Delta", "P_Value", "FDR_63")], by = "Gene", suffixes = c("_TCGA", 
        "_ICGC"))

cross_cohort$Same_Direction <- with(cross_cohort, sign(Cliff_Delta_TCGA) == sign(Cliff_Delta_ICGC) & 
    Cliff_Delta_TCGA != 0 & Cliff_Delta_ICGC != 0)

cross_cohort$Directional_Replication <- with(cross_cohort, Same_Direction & FDR_63_TCGA < 0.05 & P_Value_ICGC < 
    0.05)

cross_cohort$Strict_Replication <- with(cross_cohort, Same_Direction & FDR_63_TCGA < 0.05 & FDR_63_ICGC < 
    0.05)

cross_cohort <- cross_cohort[order(cross_cohort$FDR_63_TCGA, cross_cohort$FDR_63_ICGC), , drop = FALSE]

rownames(cross_cohort) <- NULL

focus_genes <- c("CTNNB1", "AXIN2", "GLUL", "TBX3", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", "SLC10A1", 
    "ABCB11", "NR1H4", "NR0B2", "CD8A", "GZMB", "IFNG", "CD274", "B2M", "NKG7")

focus_results <- cross_cohort[match(focus_genes, cross_cohort$Gene), , drop = FALSE]

outdir <- file.path(PROJECT_DIR, "05_results")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(tcga, file.path(outdir, "Stage3_2B_TCGA_CTNNB1_gene_association.csv"), row.names = FALSE)

write.csv(icgc, file.path(outdir, "Stage3_2B_ICGC_CTNNB1_gene_association.csv"), row.names = FALSE)

write.csv(cross_cohort, file.path(outdir, "Stage3_2B_cross_cohort_replication.csv"), row.names = FALSE)

write.csv(focus_results, file.path(outdir, "Stage3_2B_CTNNB1_BA_Immune_focus.csv"), row.names = FALSE)

stopifnot(nrow(tcga) == 63L, nrow(icgc) == 63L, nrow(cross_cohort) == 63L, all(is.finite(cross_cohort$Cliff_Delta_TCGA)), 
    all(is.finite(cross_cohort$Cliff_Delta_ICGC)))

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("association_list"), exists("cross_cohort"), exists("focus_genes"), nrow(cross_cohort) == 
    63L)

available_genes <- cross_cohort$Gene

focus_present <- focus_genes[focus_genes %in% available_genes]

focus_missing <- focus_genes[!focus_genes %in% available_genes]

focus_results <- cross_cohort[match(focus_present, cross_cohort$Gene), , drop = FALSE]

rownames(focus_results) <- NULL

stopifnot(nrow(focus_results) == length(focus_present), identical(focus_results$Gene, focus_present), 
    !anyNA(focus_results$Gene))

focus_availability <- data.frame(Gene = focus_genes, In_63_Gene_Panel = focus_genes %in% available_genes, 
    stringsAsFactors = FALSE)

outdir <- file.path(PROJECT_DIR, "05_results")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

write.csv(association_list$TCGA_LIHC, file.path(outdir, "Stage3_2B_TCGA_CTNNB1_gene_association.csv"), 
    row.names = FALSE)

write.csv(association_list$ICGC_LICA_FR, file.path(outdir, "Stage3_2B_ICGC_CTNNB1_gene_association.csv"), 
    row.names = FALSE)

write.csv(cross_cohort, file.path(outdir, "Stage3_2B_cross_cohort_replication.csv"), row.names = FALSE)

write.csv(focus_results, file.path(outdir, "Stage3_2B_CTNNB1_BA_Immune_focus.csv"), row.names = FALSE)

write.csv(focus_availability, file.path(outdir, "Stage3_2B_focus_gene_availability.csv"), row.names = FALSE)

stopifnot(nrow(association_list$TCGA_LIHC) == 63L, nrow(association_list$ICGC_LICA_FR) == 63L, nrow(cross_cohort) == 
    63L, all(is.finite(cross_cohort$Cliff_Delta_TCGA)), all(is.finite(cross_cohort$Cliff_Delta_ICGC)))

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("joint_data"))

clinical_audit <- list()

pattern <- paste(c("age", "sex", "gender", "stage", "grade", "hepat", "hbv", "hcv", "alcohol", "fibrosis", 
    "etiolog", "purity", "tumor", "tumour", "vascular", "cirrho"), collapse = "|")

for (nm in names(joint_data)) {
    clinical <- joint_data[[nm]]$clinical
    selected <- grep(pattern, colnames(clinical), ignore.case = TRUE, value = TRUE)
    if (length(selected) == 0L) 
        next
    clinical_audit[[nm]] <- do.call(rbind, lapply(selected, function(v) {
        values <- trimws(as.character(clinical[[v]]))
        missing <- is.na(values) | toupper(values) %in% c("", "NA", "N/A", "UNKNOWN", "NOT REPORTED", 
            "NOT AVAILABLE", "[NOT AVAILABLE]", "[NOT REPORTED]")
        present_values <- values[!missing]
        data.frame(Cohort = nm, Variable = v, N_Patients = nrow(clinical), N_Available = sum(!missing), 
            Missing_Pct = round(100 * mean(missing), 2), N_Unique = length(unique(present_values)), Example_Values = paste(head(unique(present_values), 
                3), collapse = " | "), stringsAsFactors = FALSE)
    }))
}

clinical_audit_table <- do.call(rbind, clinical_audit)

rownames(clinical_audit_table) <- NULL

write.csv(clinical_audit_table, file.path(PROJECT_DIR, "08_reports", "Stage3_2B_clinical_covariate_audit.csv"), 
    row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

make_covariates <- function(joint) {
    cohort <- joint$cohort
    clinical <- joint$clinical
    ids <- joint$patient_ids
    if (cohort == "TCGA_LIHC") {
        age_raw <- clinical$age_at_index
        sex_raw <- clinical$sex_at_birth
    }
    else if (cohort == "ICGC_LICA_FR") {
        age_raw <- clinical$Age
        sex_raw <- clinical$Gender
    }
    else {
        stop("Unknown cohort: ", cohort)
    }
    age <- suppressWarnings(as.numeric(as.character(age_raw)))
    sex_raw <- toupper(trimws(as.character(sex_raw)))
    sex <- ifelse(sex_raw %in% c("M", "MALE"), "Male", ifelse(sex_raw %in% c("F", "FEMALE"), "Female", 
        "NA"))
    status <- as.integer(joint$wnt_mutation["CTNNB1", ids])
    dat <- data.frame(Patient_ID = ids, Mutation = status, Age = age, Sex = factor(sex, levels = c("Female", 
        "Male")))
    stopifnot(identical(clinical$Patient_ID, ids), length(age) == length(ids), all(is.finite(age)), all(age > 
        0 & age < 120), !anyNA(dat$Sex), length(unique(dat$Sex)) == 2L, all(status %in% c(0L, 1L)))
    dat
}

fit_adjusted_gene <- function(gene, expr, covariates) {
    df <- covariates
    df$Y <- as.numeric(expr[gene, df$Patient_ID])
    fit <- lm(Y ~ Mutation + Age + Sex, data = df)
    X <- model.matrix(fit)
    if (fit$rank != ncol(X)) {
        stop("Design matrix rank deficient: ", gene)
    }
    inv_xtx <- solve(crossprod(X))
    h <- hatvalues(fit)
    adjusted_residual <- residuals(fit)/pmax(1 - h, 1e-10)
    meat <- crossprod(X * adjusted_residual)
    vcov_hc3 <- inv_xtx %*% meat %*% inv_xtx
    beta <- unname(coef(fit)["Mutation"])
    se <- sqrt(vcov_hc3["Mutation", "Mutation"])
    if (!is.finite(se) || se <= 0) {
        stop("Invalid robust SE for gene: ", gene)
    }
    df_res <- df.residual(fit)
    t_value <- beta/se
    p_value <- 2 * pt(-abs(t_value), df = df_res)
    tcrit <- qt(0.975, df = df_res)
    data.frame(Gene = gene, N = nrow(df), N_Mutation = sum(df$Mutation == 1), N_Reference = sum(df$Mutation == 
        0), Beta_Adjusted = beta, SE_HC3 = se, CI_Lower = beta - tcrit * se, CI_Upper = beta + tcrit * 
        se, P_Value = p_value, stringsAsFactors = FALSE)
}

adjusted_results <- list()

for (cohort in names(joint_data)) {
    joint <- joint_data[[cohort]]
    covariates <- make_covariates(joint)
    expr <- joint$mechanism_logTPM
    stopifnot(identical(colnames(expr), covariates$Patient_ID), nrow(expr) == 63L)
    cat("\n===== ADJUSTED MODEL:", cohort, "=====")
    cat("N:", nrow(covariates), "| Mut:", sum(covariates$Mutation == 1), "| No recorded mutation:", sum(covariates$Mutation == 
        0), "")
    cat("Age median:", median(covariates$Age), "")
    print(table(covariates$Mutation, covariates$Sex))
    res <- do.call(rbind, lapply(rownames(expr), fit_adjusted_gene, expr = expr, covariates = covariates))
    res$FDR_63 <- p.adjust(res$P_Value, method = "BH")
    res$Cohort <- cohort
    res <- res[order(res$FDR_63, res$P_Value), , drop = FALSE]
    rownames(res) <- NULL
    adjusted_results[[cohort]] <- res
    cat("\nFDR < 0.05:", sum(res$FDR_63 < 0.05), "")
    print(head(res[, c("Gene", "Beta_Adjusted", "CI_Lower", "CI_Upper", "FDR_63")], 12), digits = 5, 
        row.names = FALSE)
}

keep <- c("Gene", "Beta_Adjusted", "CI_Lower", "CI_Upper", "P_Value", "FDR_63")

adjusted_cross <- merge(adjusted_results$TCGA_LIHC[, keep], adjusted_results$ICGC_LICA_FR[, keep], by = "Gene", 
    suffixes = c("_TCGA", "_ICGC"))

adjusted_cross$Same_Direction <- with(adjusted_cross, sign(Beta_Adjusted_TCGA) == sign(Beta_Adjusted_ICGC) & 
    Beta_Adjusted_TCGA != 0 & Beta_Adjusted_ICGC != 0)

adjusted_cross$Strict_Replication <- with(adjusted_cross, Same_Direction & FDR_63_TCGA < 0.05 & FDR_63_ICGC < 
    0.05)

focus_genes <- c("NOTUM", "AXIN2", "GLUL", "TBX3", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", "SLC10A1", 
    "ABCB11", "ABCC2", "ABCC4", "SLCO1B3", "CYP27A1", "NR1H4", "NR0B2", "CD8A", "GZMB", "IFNG", "NKG7", 
    "HAVCR2")

focus_present <- intersect(focus_genes, adjusted_cross$Gene)

focus_adjusted <- adjusted_cross[match(focus_present, adjusted_cross$Gene), , drop = FALSE]

rownames(focus_adjusted) <- NULL

outdir <- file.path(PROJECT_DIR, "05_results")

write.csv(adjusted_results$TCGA_LIHC, file.path(outdir, "Stage3_2C_TCGA_age_sex_adjusted.csv"), row.names = FALSE)

write.csv(adjusted_results$ICGC_LICA_FR, file.path(outdir, "Stage3_2C_ICGC_age_sex_adjusted.csv"), row.names = FALSE)


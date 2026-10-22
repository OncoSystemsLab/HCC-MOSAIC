if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

BG <- BA4_BUNDLE$background_genes

signature_source <- file.path(PROJECT_DIR, "02_reference", "HCC_ICB_57_gene_sets.rds")

if (!file.exists(signature_source)) {
    cat("\nSelect the required source file: ", "HCC_ICB_57_gene_sets.rds")
    chosen_file <- hcc_input_file("HCC_INPUT_FILE_004")
    original_sets <- readRDS(chosen_file)
    if (!is.list(original_sets) || length(original_sets) != 57L) {
        stop("The selected file is not the expected 57-gene-set RDS")
    }
    if (!file.copy(chosen_file, signature_source)) {
        stop("Could not copy the 57-gene-set file")
    }
}

GENESETS_57 <- readRDS(signature_source)

stopifnot(is.list(GENESETS_57), length(GENESETS_57) == 57L, !is.null(names(GENESETS_57)), !anyDuplicated(names(GENESETS_57)))

ba4_match <- vapply(BA4_BUNDLE$features, function(feature) {
    if (!feature %in% names(GENESETS_57)) {
        return(FALSE)
    }
    setequal(GENESETS_57[[feature]], BA4_BUNDLE$gene_sets[[feature]])
}, logical(1))

if (!all(ba4_match)) {
    stop("The 57 gene sets differ from the original BA4 model; stopped")
}

signature_modules <- list(Driver_Signaling = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_P53_PATHWAY", 
    "HALLMARK_MYC_TARGETS_V1", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING"), Metabolic_Reprogramming = c("HALLMARK_BILE_ACID_METABOLISM", 
    "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", 
    "HALLMARK_GLYCOLYSIS", "HALLMARK_PEROXISOME", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_XENOBIOTIC_METABOLISM", 
    "HALLMARK_HYPOXIA"), Immune_Ecosystem = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", 
    "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_ALLOGRAFT_REJECTION", 
    "HALLMARK_COMPLEMENT", "ICI__GEP18_T_CELL_INFLAMED", "ICI__IFNG6", "ICI__CYT2", "ICI__HCC_INFLAMMATORY4", 
    "ICI__HCC_TEFF3", "ICI__HCC_ISS10"))

selected_names <- unlist(signature_modules, use.names = FALSE)

stopifnot(!anyDuplicated(selected_names))

missing_signatures <- setdiff(selected_names, names(GENESETS_57))

if (length(missing_signatures) > 0) {
    stop("Gene sets not found: ", paste(missing_signatures, collapse = ", "))
}

selected_gene_sets <- GENESETS_57[selected_names]

signature_registry <- data.frame(Module = rep(names(signature_modules), lengths(signature_modules)), 
    Signature = selected_names, N_Genes = vapply(selected_gene_sets, length, integer(1)), N_BG_Matched = vapply(selected_gene_sets, 
        function(gs) sum(gs %in% BG), integer(1)), stringsAsFactors = FALSE)

signature_registry$BG_Coverage <- (signature_registry$N_BG_Matched/signature_registry$N_Genes)

signature_registry$Source <- ifelse(grepl("^HALLMARK_", signature_registry$Signature), "MSigDB_Hallmark_frozen", 
    "Original_ICI_signature_frozen")

stopifnot(all(signature_registry$N_Genes >= 2), all(signature_registry$BG_Coverage == 1), all(BA4_BUNDLE$features %in% 
    selected_names))

saveRDS(selected_gene_sets, file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_signatures_v1.rds"))

write.csv(signature_registry, file.path(PROJECT_DIR, "08_reports", "Stage2_1_signature_registry.csv"), 
    row.names = FALSE)

molecular_panels <- list(BA_Classical_Synthesis = c("CYP7A1", "CYP8B1", "HSD3B7", "AKR1D1", "CYP27A1", 
    "ACOX2"), BA_Alternative_Synthesis = c("CYP27A1", "CYP7B1", "HSD3B7", "AKR1D1"), BA_Conjugation = c("BAAT", 
    "SLC27A5"), BA_Sinusoidal_Uptake = c("SLC10A1", "SLCO1B1", "SLCO1B3"), BA_Canalicular_Export = c("ABCB11", 
    "ABCC2"), BA_Basolateral_Export = c("ABCC3", "ABCC4", "SLC51A", "SLC51B"), BA_FXR_FGF19_Axis = c("NR1H4", 
    "NR0B2", "FGF19", "FGFR4", "KLB"), BA_Detoxification = c("SULT2A1", "UGT2B4", "UGT2B7"), WNT_CTNNB1_Response = c("AXIN2", 
    "GLUL", "TBX3", "NOTUM", "RNF43", "APCDD1", "NKD1", "LGR5"), T_Cell_Cytotoxic = c("CD8A", "CD8B", 
    "GZMB", "PRF1", "NKG7", "IFNG"), NK_Cell_Markers = c("NCR1", "KLRD1", "FCGR3A", "GNLY", "NKG7", "PRF1"), 
    T_Cell_Checkpoints = c("PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX"), Myeloid_Regulatory = c("CD163", 
        "MRC1", "MSR1", "CSF1R", "IL10", "TGFB1"), Antigen_Presentation = c("B2M", "HLA-A", "HLA-B", 
        "HLA-C", "TAP1", "TAP2"))

panel_registry <- do.call(rbind, lapply(names(molecular_panels), function(name) {
    gs <- unique(molecular_panels[[name]])
    present <- intersect(gs, BG)
    missing <- setdiff(gs, BG)
    data.frame(Panel = name, N_Genes = length(gs), N_Matched = length(present), N_Missing = length(missing), 
        Coverage = length(present)/length(gs), Missing_Genes = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
}))

rownames(panel_registry) <- NULL

panel_gene_registry <- do.call(rbind, lapply(names(molecular_panels), function(name) {
    gs <- unique(molecular_panels[[name]])
    data.frame(Panel = name, Gene = gs, In_Training_Background = gs %in% BG, stringsAsFactors = FALSE)
}))

rownames(panel_gene_registry) <- NULL

priority_genes <- c("BAAT", "SLC27A5", "CYP7A1", "CYP8B1", "SLC10A1", "ABCB11", "NR1H4", "NR0B2", "AXIN2", 
    "GLUL")

priority_qc <- data.frame(Gene = priority_genes, Available = priority_genes %in% BG)

if (!all(priority_qc$Available)) {
    warning("Some key genes are absent from the training background; ", "handle these separately in downstream molecular profiling.")
}

saveRDS(molecular_panels, file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_gene_panels_v1.rds"))

write.csv(panel_registry, file.path(PROJECT_DIR, "08_reports", "Stage2_1_gene_panel_coverage.csv"), row.names = FALSE)

write.csv(panel_gene_registry, file.path(PROJECT_DIR, "08_reports", "Stage2_1_gene_panel_mapping.csv"), 
    row.names = FALSE)

stage2_config <- list(version = "HCC_Molecular_Fingerprint_v1", input_unit = "TPM", transformation = "log2(TPM + 1)", 
    background_genes = BG, signature_modules = signature_modules, signature_sets = selected_gene_sets, 
    gene_panels = molecular_panels, source_57_gene_sets_md5 = unname(tools::md5sum(signature_source)), 
    preprocessing = BA4_BUNDLE$preprocessing)

saveRDS(stage2_config, file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

library(GSVA)

library(BiocParallel)

library(glmnet)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

signature_sets <- config$signature_sets

signature_names <- names(signature_sets)

bg <- config$background_genes

pre <- config$preprocessing

stopifnot(length(signature_sets) == 26, length(bg) == 21308, !anyDuplicated(signature_names), identical(bg, 
    BA4_BUNDLE$background_genes), all(vapply(signature_sets, function(gs) all(gs %in% bg), logical(1))))

tpm_file <- file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds")

if (!file.exists(tpm_file)) {
    stop("Three-cohort TPM cache is missing; review Stage 1.5B")
}

tpm_reference <- readRDS(tpm_file)

tpm_reference <- lapply(tpm_reference, function(x) {
    mat <- as.matrix(x)
    storage.mode(mat) <- "numeric"
    stopifnot(!is.null(rownames(mat)), !is.null(colnames(mat)), !anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)), 
        all(bg %in% rownames(mat)))
    mat <- mat[bg, , drop = FALSE]
    stopifnot(all(is.finite(mat)), all(mat >= 0))
    return(mat)
})

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

stopifnot(nrow(reference_map) == 69, !anyDuplicated(reference_map$patient_uid), all(reference_map$Cohort %in% 
    names(tpm_reference)))

signature_matrix <- matrix(NA_real_, nrow = nrow(reference_map), ncol = length(signature_names), dimnames = list(reference_map$patient_uid, 
    signature_names))

for (i in seq_len(nrow(reference_map))) {
    uid <- reference_map$patient_uid[i]
    cohort <- reference_map$Cohort[i]
    sample_id <- reference_map$Raw_Sample_ID[i]
    cat("Processing", i, "/69 |", uid, "")
    mat <- tpm_reference[[cohort]]
    if (!sample_id %in% colnames(mat)) {
        stop("Patient not found: ", uid)
    }
    one_tpm <- mat[, sample_id, drop = FALSE]
    expr <- log2(one_tpm + 1)
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = signature_sets, minSize = pre$ssgsea_minSize, 
        maxSize = pre$ssgsea_maxSize, alpha = pre$ssgsea_alpha, normalize = pre$ssgsea_normalize)
    es <- suppressMessages(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()))
    if (!all(signature_names %in% rownames(es))) {
        stop("Missing ssGSEA features: ", uid)
    }
    signature_matrix[i, ] <- as.numeric(es[signature_names, 1])
}

stopifnot(identical(dim(signature_matrix), c(69L, 26L)), all(is.finite(signature_matrix)))

matrix_file <- file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_26signatures.rds")

saveRDS(signature_matrix, matrix_file)

signature_df <- data.frame(patient_uid = reference_map$patient_uid, Cohort = reference_map$Cohort, signature_matrix, 
    check.names = FALSE, row.names = NULL)

write.csv(signature_df, file.path(PROJECT_DIR, "04_processed", "Stage2_69patients_26signatures.csv"), 
    row.names = FALSE)

feature_ref <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage1_training_57features.rds"))

signature_matrix <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_26signatures.rds"))

signature_names <- colnames(signature_matrix)

stopifnot(all(signature_names %in% colnames(feature_ref)), "patient_uid" %in% colnames(feature_ref), 
    !anyDuplicated(feature_ref$patient_uid))

idx <- match(rownames(signature_matrix), feature_ref$patient_uid)

stopifnot(!anyNA(idx))

original_matrix <- as.matrix(feature_ref[idx, signature_names, drop = FALSE])

storage.mode(original_matrix) <- "numeric"

rownames(original_matrix) <- rownames(signature_matrix)

stopifnot(identical(dim(original_matrix), dim(signature_matrix)), all(is.finite(original_matrix)))

error_matrix <- abs(signature_matrix - original_matrix)

max_signature_error <- max(error_matrix)

feature_error_summary <- data.frame(Signature = signature_names, Max_Absolute_Error = apply(error_matrix, 
    2, max), Mean_Absolute_Error = colMeans(error_matrix), row.names = NULL)

ba4_features <- BA4_BUNDLE$features

ba4_feature_error <- max(error_matrix[, ba4_features, drop = FALSE])

X_BA4 <- signature_matrix[, ba4_features, drop = FALSE]

recalculated_icb <- as.numeric(predict(BA4_BUNDLE$model, newx = X_BA4, s = BA4_BUNDLE$lambda, type = "response"))

score_idx <- match(rownames(signature_matrix), reference_map$patient_uid)

stopifnot(!anyNA(score_idx))

original_icb <- reference_map$BA4_ReproducedScore[score_idx]

icb_error <- abs(recalculated_icb - original_icb)

max_icb_error <- max(icb_error)

patient_qc <- data.frame(patient_uid = rownames(signature_matrix), Cohort = reference_map$Cohort[score_idx], 
    Original_ICBScore = original_icb, Recalculated_ICBScore = recalculated_icb, Max_Signature_Error = apply(error_matrix, 
        1, max), ICBScore_Error = icb_error, row.names = NULL)

write.csv(feature_error_summary, file.path(PROJECT_DIR, "08_reports", "Stage2_2A_signature_validation.csv"), 
    row.names = FALSE)

write.csv(patient_qc, file.path(PROJECT_DIR, "08_reports", "Stage2_2A_patient_validation.csv"), row.names = FALSE)

passed <- (nrow(signature_matrix) == 69 && ncol(signature_matrix) == 26 && max_signature_error < 1e-08 && 
    max_icb_error < 1e-08)

if (passed) {
    cat("\nStage 2.2A: PASS")
} else {
    cat("\nStage 2.2A: CHECK REQUIRED")
    cat("Identify the features responsible for these differences; ", "keep the original model and gene sets unchanged during review.")
}

PROJECT_DIR <- "HCC_OmicsAgent"

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

panels <- config$gene_panels

bg <- config$background_genes

stopifnot(length(panels) == 14, length(bg) == 21308)

panel_genes <- unique(unlist(panels, use.names = FALSE))

extra_genes <- c("CTNNB1")

mechanism_genes <- unique(c(panel_genes, extra_genes))

missing_genes <- setdiff(mechanism_genes, bg)

if (length(missing_genes) > 0) {
    stop("Mechanistic genes absent from the fixed background: ", paste(missing_genes, collapse = ", "))
}

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

tpm_reference <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds"))

stopifnot(nrow(reference_map) == 69, !anyDuplicated(reference_map$patient_uid), all(reference_map$Cohort %in% 
    names(tpm_reference)))

mechanism_tpm <- matrix(NA_real_, nrow = nrow(reference_map), ncol = length(mechanism_genes), dimnames = list(reference_map$patient_uid, 
    mechanism_genes))

for (cohort in unique(reference_map$Cohort)) {
    cat("\nProcessing cohort:", cohort, "")
    idx <- which(reference_map$Cohort == cohort)
    sample_ids <- as.character(reference_map$Raw_Sample_ID[idx])
    mat <- as.matrix(tpm_reference[[cohort]])
    storage.mode(mat) <- "numeric"
    stopifnot(!anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)), all(mechanism_genes %in% 
        rownames(mat)), all(sample_ids %in% colnames(mat)))
    cohort_expr <- t(mat[mechanism_genes, sample_ids, drop = FALSE])
    mechanism_tpm[idx, ] <- cohort_expr
    cat("Patients:", length(idx), "")
}

stopifnot(nrow(mechanism_tpm) == 69, ncol(mechanism_tpm) == length(mechanism_genes), all(is.finite(mechanism_tpm)), 
    all(mechanism_tpm >= 0))

mechanism_logtpm <- log2(mechanism_tpm + 1)

saveRDS(mechanism_tpm, file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_mechanism_TPM.rds"))

saveRDS(mechanism_logtpm, file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_mechanism_logTPM.rds"))

make_expression_df <- function(mat) {
    data.frame(patient_uid = rownames(mat), Cohort = reference_map$Cohort, mat, check.names = FALSE, 
        row.names = NULL)
}

write.csv(make_expression_df(mechanism_tpm), file.path(PROJECT_DIR, "04_processed", "Stage2_69patients_mechanism_TPM.csv"), 
    row.names = FALSE)

write.csv(make_expression_df(mechanism_logtpm), file.path(PROJECT_DIR, "04_processed", "Stage2_69patients_mechanism_logTPM.csv"), 
    row.names = FALSE)

stopifnot(identical(rownames(mechanism_tpm), reference_map$patient_uid), identical(colnames(mechanism_tpm), 
    mechanism_genes))

transformation_error <- max(abs(mechanism_logtpm - log2(mechanism_tpm + 1)))

source_error <- numeric(nrow(reference_map))

for (i in seq_len(nrow(reference_map))) {
    cohort <- reference_map$Cohort[i]
    sample_id <- reference_map$Raw_Sample_ID[i]
    original <- as.numeric(tpm_reference[[cohort]][mechanism_genes, sample_id])
    extracted <- as.numeric(mechanism_tpm[i, ])
    source_error[i] <- max(abs(original - extracted))
}

max_source_error <- max(source_error)

panel_expression_qc <- do.call(rbind, lapply(names(panels), function(panel) {
    genes <- unique(panels[[panel]])
    x <- mechanism_tpm[, genes, drop = FALSE]
    detected_per_patient <- rowSums(x > 0)
    data.frame(Panel = panel, N_Genes = length(genes), Median_Detected_Genes = median(detected_per_patient), 
        Min_Detected_Genes = min(detected_per_patient), Max_Detected_Genes = max(detected_per_patient), 
        Zero_Fraction = mean(x == 0), stringsAsFactors = FALSE)
}))

rownames(panel_expression_qc) <- NULL

priority_genes <- c("CTNNB1", "AXIN2", "GLUL", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", "SLC10A1", "ABCB11", 
    "NR1H4")

stopifnot(all(priority_genes %in% colnames(mechanism_tpm)))

priority_gene_qc <- data.frame(Gene = priority_genes, Median_TPM = apply(mechanism_tpm[, priority_genes, 
    drop = FALSE], 2, median), Zero_Fraction = colMeans(mechanism_tpm[, priority_genes, drop = FALSE] == 
    0), row.names = NULL)

write.csv(panel_expression_qc, file.path(PROJECT_DIR, "08_reports", "Stage2_2B_panel_expression_QC.csv"), 
    row.names = FALSE)

write.csv(priority_gene_qc, file.path(PROJECT_DIR, "08_reports", "Stage2_2B_priority_gene_QC.csv"), row.names = FALSE)

passed <- (nrow(mechanism_tpm) == 69 && length(unique(reference_map$Cohort)) == 3 && transformation_error < 
    1e-12 && max_source_error < 1e-12 && nrow(panel_expression_qc) == 14)

if (passed) {
    cat("\nStage 2.2B: PASS")
} else {
    cat("\nStage 2.2B: CHECK REQUIRED")
}

PROJECT_DIR <- "HCC_OmicsAgent"

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

signature_mat <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_26signatures.rds"))

gene_mat <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_mechanism_logTPM.rds"))

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

patient_ids <- reference_map$patient_uid

stopifnot(length(patient_ids) == 69L, !anyDuplicated(patient_ids), all(patient_ids %in% rownames(signature_mat)), 
    all(patient_ids %in% rownames(gene_mat)))

signature_mat <- signature_mat[patient_ids, , drop = FALSE]

gene_mat <- gene_mat[patient_ids, , drop = FALSE]

stopifnot(ncol(signature_mat) == 26L, ncol(gene_mat) == 63L, all(is.finite(signature_mat)), all(is.finite(gene_mat)), 
    !anyDuplicated(c(colnames(signature_mat), colnames(gene_mat))))

reference_matrix <- cbind(signature_mat, gene_mat)

stopifnot(identical(dim(reference_matrix), c(69L, 89L)))

module_names <- unlist(config$signature_modules, use.names = FALSE)

module_labels <- rep(names(config$signature_modules), lengths(config$signature_modules))

module_lookup <- setNames(module_labels, module_names)

feature_types <- c(rep("ssGSEA", ncol(signature_mat)), rep("log2TPM", ncol(gene_mat)))

feature_modules <- c(unname(module_lookup[colnames(signature_mat)]), rep("Mechanism_Gene", ncol(gene_mat)))

stopifnot(!anyNA(feature_modules))

summarize_reference <- function(x) {
    c(Mean = mean(x), SD = sd(x), Median = median(x), MAD = mad(x), P05 = unname(quantile(x, 0.05)), 
        P25 = unname(quantile(x, 0.25)), P75 = unname(quantile(x, 0.75)), P95 = unname(quantile(x, 0.95)), 
        Zero_Fraction = mean(x == 0))
}

reference_stats <- t(vapply(seq_len(ncol(reference_matrix)), function(j) {
    summarize_reference(reference_matrix[, j])
}, numeric(9)))

rownames(reference_stats) <- colnames(reference_matrix)

reference_registry <- data.frame(Feature = colnames(reference_matrix), Type = feature_types, Module = feature_modules, 
    N_Reference = nrow(reference_matrix), reference_stats, check.names = FALSE, row.names = NULL)

reference_registry$Scale_Valid <- (is.finite(reference_registry$SD) & reference_registry$SD > 1e-12)

reference_bundle <- list(version = "HCC_Molecular_Reference_v1", reference_type = "Exploratory_Development_Cohorts", 
    n_patients = nrow(reference_matrix), cohort_counts = table(reference_map$Cohort), feature_names = colnames(reference_matrix), 
    feature_types = feature_types, reference_values = reference_matrix, reference_registry = reference_registry, 
    percentile_method = "Empirical_midrank", standardization = "Reference_mean_SD")

saveRDS(reference_bundle, file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))

write.csv(reference_registry, file.path(PROJECT_DIR, "08_reports", "Stage2_3_reference_registry.csv"), 
    row.names = FALSE)

reference_bundle <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))

X <- reference_bundle$reference_values

registry <- reference_bundle$reference_registry

stopifnot(identical(dim(X), c(69L, 89L)), identical(colnames(X), registry$Feature), identical(rownames(X), 
    reference_map$patient_uid))

percentile_matrix <- vapply(seq_len(ncol(X)), function(j) {
    100 * (rank(X[, j], ties.method = "average") - 0.5)/nrow(X)
}, numeric(nrow(X)))

dimnames(percentile_matrix) <- dimnames(X)

z_matrix <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X), dimnames = dimnames(X))


if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(GSVA)

library(BiocParallel)

library(glmnet)

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

predict_BA4_from_TPM <- function(tpm_matrix, bundle = BA4_BUNDLE) {
    tpm <- as.matrix(tpm_matrix)
    storage.mode(tpm) <- "numeric"
    if (is.null(rownames(tpm)) || is.null(colnames(tpm))) {
        stop("The TPM matrix must have gene and sample names")
    }
    if (anyDuplicated(rownames(tpm))) {
        stop("Duplicate gene symbols found; merge them first")
    }
    if (anyDuplicated(colnames(tpm))) {
        stop("Duplicate sample names found")
    }
    if (any(!is.finite(tpm)) || any(tpm < 0)) {
        stop("TPM contains NA, Inf or negative values")
    }
    bg <- bundle$background_genes
    missing_genes <- setdiff(bg, rownames(tpm))
    cat("Training background:", length(bg), "")
    cat("Missing background genes:", length(missing_genes), "")
    if (length(missing_genes) > 0) {
        stop("The gene background is incomplete; missing ", length(missing_genes), " genes. Prediction stopped.")
    }
    tpm_fixed <- tpm[bg, , drop = FALSE]
    expr <- log2(tpm_fixed + 1)
    pre <- bundle$preprocessing
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = bundle$gene_sets, minSize = pre$ssgsea_minSize, 
        maxSize = pre$ssgsea_maxSize, alpha = pre$ssgsea_alpha, normalize = pre$ssgsea_normalize)
    es <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
    stopifnot(all(bundle$features %in% rownames(es)))
    X <- t(es[bundle$features, , drop = FALSE])
    storage.mode(X) <- "numeric"
    if (any(!is.finite(X))) {
        stop("BA4 features contain nonfinite values")
    }
    score <- as.numeric(predict(bundle$model, newx = X, type = "response", s = bundle$lambda))
    result <- data.frame(Sample_ID = rownames(X), ICBScore = score, X, check.names = FALSE, row.names = NULL)
    return(list(prediction = result, features = X, ssgsea = es, n_background_genes = length(bg)))
}

test_sample <- "S02-T-RNA"

if (exists("tpm_background") && is.list(tpm_background) && "PRJEB34724" %in% names(tpm_background)) {
    tpm_test_source <- as.matrix(tpm_background$PRJEB34724)
    cat("Using existing tpm_background")
} else if (exists("PRJEB34724.tpm")) {
    df <- as.data.frame(PRJEB34724.tpm, check.names = FALSE)
    stopifnot("Symbol" %in% colnames(df))
    df$Symbol <- trimws(as.character(df$Symbol))
    df <- df[!is.na(df$Symbol) & df$Symbol != "", , drop = FALSE]
    sample_cols <- setdiff(colnames(df), "Symbol")
    mat <- as.matrix(df[, sample_cols, drop = FALSE])
    storage.mode(mat) <- "numeric"
    if (any(!is.finite(mat))) {
        stop("The original TPM matrix contains NA or Inf")
    }
    tpm_test_source <- rowsum(mat, group = df$Symbol, reorder = TRUE)
    cat("Using PRJEB34724.tpm")
} else {
    stop("Neither tpm_background nor PRJEB34724.tpm was found.", "Load the TPM objects from the preceding data-preparation stage first.")
}

stopifnot(test_sample %in% colnames(tpm_test_source))

single_tpm <- tpm_test_source[, test_sample, drop = FALSE]

single_result <- predict_BA4_from_TPM(single_tpm)

write.csv(single_result$prediction, file.path(PROJECT_DIR, "05_results", "Stage1_S02_RNA_to_ICBScore.csv"), 
    row.names = FALSE)

original_features <- readRDS(features_file)

patient_uid <- paste0("PRJEB34724__", test_sample)

idx <- match(patient_uid, original_features$patient_uid)

stopifnot(!is.na(idx))

old_features <- as.numeric(original_features[idx, BA4_BUNDLE$features, drop = TRUE])

new_features <- as.numeric(single_result$features[1, BA4_BUNDLE$features])

feature_comparison <- data.frame(Feature = BA4_BUNDLE$features, Original = old_features, Recalculated = new_features, 
    Absolute_Difference = abs(old_features - new_features))

max_feature_difference <- max(feature_comparison$Absolute_Difference)

old_score <- reproduced_df$BA4_ReproducedScore[match(patient_uid, reproduced_df$patient_uid)]

new_score <- single_result$prediction$ICBScore[1]

score_difference <- abs(old_score - new_score)

if (max_feature_difference < 1e-08 && score_difference < 1e-08) {
    cat("\nStage 1.1: PASS")
} else {
    cat("\nStage 1.1: CHECK REQUIRED")
    cat("Review differences between single-sample and original batch ssGSEA calculations.")
}

write.csv(feature_comparison, file.path(PROJECT_DIR, "08_reports", "Stage1_BA4_feature_reproducibility.csv"), 
    row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

read_HCC_TPM <- function(file_path) {
    if (!file.exists(file_path)) {
        stop("Input file not found: ", file_path)
    }
    ext <- tolower(tools::file_ext(file_path))
    if (ext == "csv") {
        df <- read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
    else if (ext %in% c("tsv", "txt")) {
        df <- read.delim(file_path, check.names = FALSE, stringsAsFactors = FALSE)
    }
    else {
        stop("Only CSV, TSV or tab-delimited TXT files are supported")
    }
    if (ncol(df) < 2) {
        stop("A gene column and at least one sample column are required")
    }
    genes <- trimws(as.character(df[[1]]))
    if (anyNA(genes) || any(genes == "")) {
        stop("Gene symbols contain blank values")
    }
    if (anyDuplicated(genes)) {
        stop("Duplicate gene symbols found; review them before continuing")
    }
    if (mean(grepl("^ENSG[0-9]+", genes)) > 0.5) {
        stop("Ensembl IDs detected; this model requires gene symbols")
    }
    sample_names <- colnames(df)[-1]
    if (anyNA(sample_names) || any(trimws(sample_names) == "") || anyDuplicated(sample_names)) {
        stop("Sample names contain blank or duplicate values")
    }
    expression_df <- df[, -1, drop = FALSE]
    for (j in seq_along(expression_df)) {
        values <- suppressWarnings(as.numeric(expression_df[[j]]))
        if (any(!is.finite(values))) {
            stop("Sample ", sample_names[j], " contains nonnumeric, NA or Inf values")
        }
        if (any(values < 0)) {
            stop("Sample ", sample_names[j], " contains negative expression values")
        }
        expression_df[[j]] <- values
    }
    tpm <- as.matrix(expression_df)
    rownames(tpm) <- genes
    colnames(tpm) <- sample_names
    cat("\nRNA input loaded successfully")
    cat("Genes:", nrow(tpm), "")
    cat("Samples:", ncol(tpm), "")
    return(tpm)
}

QC_HCC_TPM <- function(tpm, bundle = BA4_BUNDLE) {
    bg <- bundle$background_genes
    matched <- intersect(bg, rownames(tpm))
    missing <- setdiff(bg, rownames(tpm))
    qc <- data.frame(Sample_ID = colnames(tpm), N_Genes = nrow(tpm), N_BG_Matched = length(matched), 
        N_BG_Missing = length(missing), BG_Coverage = length(matched)/length(bg), TPM_Sum = colSums(tpm), 
        Zero_Fraction = colMeans(tpm == 0), row.names = NULL)
    cat("\n===== RNA INPUT QC =====")
    print(qc)
    cat("\nTraining background:", length(bg), "")
    cat("Missing background genes:", length(missing), "")
    return(list(summary = qc, missing_genes = missing, matched_genes = matched))
}

if (!exists("predict_BA4_from_TPM")) {
    stop("Run Stage 1.1.1 to define predict_BA4_from_TPM first")
}

predict_HCC_RNA_batch <- function(tpm, bundle = BA4_BUNDLE) {
    qc_result <- QC_HCC_TPM(tpm, bundle)
    if (length(qc_result$missing_genes) > 0) {
        stop("Strict mode does not allow missing background genes.", "Missing count: ", length(qc_result$missing_genes))
    }
    result_list <- vector("list", ncol(tpm))
    for (i in seq_len(ncol(tpm))) {
        sample_name <- colnames(tpm)[i]
        cat("\nProcessing:", sample_name, "(", i, "/", ncol(tpm), ")")
        one_sample <- tpm[, i, drop = FALSE]
        result_list[[i]] <- predict_BA4_from_TPM(one_sample, bundle = bundle)$prediction
    }
    prediction_df <- do.call(rbind, result_list)
    rownames(prediction_df) <- NULL
    return(list(prediction = prediction_df, qc = qc_result$summary))
}

stopifnot(exists("tpm_test_source"))

test_samples <- c("S02-T-RNA", "S03-T-RNA")

stopifnot(all(test_samples %in% colnames(tpm_test_source)))

test_tpm <- tpm_test_source[, test_samples, drop = FALSE]

input_df <- data.frame(GeneSymbol = rownames(test_tpm), test_tpm, check.names = FALSE, row.names = NULL)

test_file <- file.path(PROJECT_DIR, "03_input", "Stage1_demo_2patients_TPM.csv")

write.csv(input_df, test_file, row.names = FALSE)

new_tpm <- read_HCC_TPM(test_file)

stopifnot(isTRUE(all.equal(test_tpm, new_tpm, tolerance = 1e-12)))

batch_result <- predict_HCC_RNA_batch(new_tpm)

stopifnot(exists("reproduced_df"))

patient_ids <- paste0("PRJEB34724__", batch_result$prediction$Sample_ID)

idx <- match(patient_ids, reproduced_df$patient_uid)

stopifnot(!anyNA(idx))

validation <- data.frame(Sample_ID = batch_result$prediction$Sample_ID, Original_ICBScore = reproduced_df$BA4_ReproducedScore[idx], 
    New_ICBScore = batch_result$prediction$ICBScore)

validation$Absolute_Difference <- abs(validation$Original_ICBScore - validation$New_ICBScore)

max_diff <- max(validation$Absolute_Difference)

if (max_diff < 1e-08) {
    cat("Stage 1.2: PASS")
} else {
    cat("Stage 1.2: CHECK REQUIRED")
}

write.csv(batch_result$prediction, file.path(PROJECT_DIR, "05_results", "Stage1_2_batch_ICBScore.csv"), 
    row.names = FALSE)

write.csv(batch_result$qc, file.path(PROJECT_DIR, "08_reports", "Stage1_2_input_QC.csv"), row.names = FALSE)

write.csv(validation, file.path(PROJECT_DIR, "08_reports", "Stage1_2_batch_validation.csv"), row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

audit_BA4_input <- function(tpm, bundle = BA4_BUNDLE) {
    tpm <- as.matrix(tpm)
    storage.mode(tpm) <- "numeric"
    stopifnot(!is.null(rownames(tpm)), !is.null(colnames(tpm)), !anyDuplicated(rownames(tpm)), all(is.finite(tpm)), 
        all(tpm >= 0))
    bg <- bundle$background_genes
    features <- bundle$features
    gene_sets <- bundle$gene_sets
    matched <- intersect(bg, rownames(tpm))
    missing <- setdiff(bg, rownames(tpm))
    extra <- setdiff(rownames(tpm), bg)
    background_qc <- data.frame(Input_Genes = nrow(tpm), Training_BG_Genes = length(bg), Matched_BG_Genes = length(matched), 
        Missing_BG_Genes = length(missing), Extra_Genes = length(extra), BG_Coverage = length(matched)/length(bg))
    signature_qc <- do.call(rbind, lapply(features, function(feature) {
        genes <- unique(gene_sets[[feature]])
        present <- intersect(genes, rownames(tpm))
        data.frame(Feature = feature, Total_Genes = length(genes), Matched_Genes = length(present), Missing_Genes = length(genes) - 
            length(present), Coverage = length(present)/length(genes))
    }))
    rownames(signature_qc) <- NULL
    sample_qc <- data.frame(Sample_ID = colnames(tpm), TPM_Sum = colSums(tpm), Zero_Fraction = colMeans(tpm == 
        0), Detected_Genes = colSums(tpm > 0))
    exact_compatible <- (length(missing) == 0 && all(signature_qc$Missing_Genes == 0))
    cat("\n===== BA4 INPUT COMPATIBILITY =====")
    print(background_qc)
    cat("\n===== SIGNATURE COVERAGE =====")
    print(signature_qc)
    cat("\n===== SAMPLE QC =====")
    print(sample_qc)
    cat("\nExact background compatible:", exact_compatible, "")
    return(list(background = background_qc, signatures = signature_qc, samples = sample_qc, missing_genes = missing, 
        exact_compatible = exact_compatible))
}

demo_file <- file.path(PROJECT_DIR, "03_input", "Stage1_demo_2patients_TPM.csv")

new_tpm <- read_HCC_TPM(demo_file)

full_audit <- audit_BA4_input(new_tpm)

stopifnot(isTRUE(full_audit$exact_compatible))

write.csv(full_audit$background, file.path(PROJECT_DIR, "08_reports", "Stage1_3_full_background_QC.csv"), 
    row.names = FALSE)

write.csv(full_audit$signatures, file.path(PROJECT_DIR, "08_reports", "Stage1_3_full_signature_QC.csv"), 
    row.names = FALSE)

set.seed(2026)

BA4_all_genes <- unique(unlist(BA4_BUNDLE$gene_sets, use.names = FALSE))

non_signature_genes <- setdiff(BA4_BUNDLE$background_genes, BA4_all_genes)

stopifnot(length(non_signature_genes) >= 100)

remove_genes <- sample(non_signature_genes, size = 100)

masked_tpm <- new_tpm[!rownames(new_tpm) %in% remove_genes, , drop = FALSE]

masked_audit <- audit_BA4_input(masked_tpm)

stopifnot(length(masked_audit$missing_genes) == 100, all(masked_audit$signatures$Missing_Genes == 0), 
    !masked_audit$exact_compatible)

write.csv(masked_audit$background, file.path(PROJECT_DIR, "08_reports", "Stage1_3_masked_background_QC.csv"), 
    row.names = FALSE)

write.csv(masked_audit$signatures, file.path(PROJECT_DIR, "08_reports", "Stage1_3_masked_signature_QC.csv"), 
    row.names = FALSE)

library(GSVA)

library(BiocParallel)

library(glmnet)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

predict_BA4_missing_bg_experiment <- function(tpm, bundle = BA4_BUNDLE) {
    tpm <- as.matrix(tpm)
    storage.mode(tpm) <- "numeric"
    stopifnot(!is.null(rownames(tpm)), !is.null(colnames(tpm)), !anyDuplicated(rownames(tpm)), !anyDuplicated(colnames(tpm)), 
        all(is.finite(tpm)), all(tpm >= 0))
    bg <- bundle$background_genes
    available_bg <- bg[bg %in% rownames(tpm)]
    missing_bg <- setdiff(bg, rownames(tpm))
    gs <- bundle$gene_sets[bundle$features]
    missing_signature <- lapply(gs, function(g) {
        setdiff(g, available_bg)
    })
    if (any(lengths(missing_signature) > 0)) {
        stop("This stage permits removal of noncore BA4 genes only; ", "review the removed gene set.")
    }
    cat("\nBackground genes available:", length(available_bg), "")
    cat("Background genes missing:", length(missing_bg), "")
    results <- vector("list", ncol(tpm))
    pre <- bundle$preprocessing
    for (i in seq_len(ncol(tpm))) {
        sample_id <- colnames(tpm)[i]
        expr <- log2(tpm[available_bg, i, drop = FALSE] + 1)
        param <- GSVA::ssgseaParam(exprData = expr, geneSets = gs, minSize = pre$ssgsea_minSize, maxSize = pre$ssgsea_maxSize, 
            alpha = pre$ssgsea_alpha, normalize = pre$ssgsea_normalize)
        es <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
        X <- t(es[bundle$features, , drop = FALSE])
        storage.mode(X) <- "numeric"
        stopifnot(all(is.finite(X)))
        score <- as.numeric(predict(bundle$model, newx = X, s = bundle$lambda, type = "response"))
        one <- data.frame(Sample_ID = sample_id, ICBScore = score, check.names = FALSE)
        for (feature in bundle$features) {
            one[[feature]] <- as.numeric(X[1, feature])
        }
        results[[i]] <- one
    }
    prediction <- do.call(rbind, results)
    rownames(prediction) <- NULL
    return(list(prediction = prediction, missing_background = missing_bg, n_background_used = length(available_bg)))
}

demo_file <- file.path(PROJECT_DIR, "03_input", "Stage1_demo_2patients_TPM.csv")

new_tpm <- read_HCC_TPM(demo_file)

baseline_file <- file.path(PROJECT_DIR, "05_results", "Stage1_2_batch_ICBScore.csv")

baseline <- read.csv(baseline_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(all(BA4_BUNDLE$background_genes %in% rownames(new_tpm)), all(colnames(new_tpm) %in% baseline$Sample_ID))

BA4_all_genes <- unique(unlist(BA4_BUNDLE$gene_sets, use.names = FALSE))

non_signature_genes <- setdiff(BA4_BUNDLE$background_genes, BA4_all_genes)

set.seed(2026)

remove_genes <- sample(non_signature_genes, size = 100)

masked_tpm <- new_tpm[!rownames(new_tpm) %in% remove_genes, , drop = FALSE]

masked_result <- predict_BA4_missing_bg_experiment(masked_tpm)

write.csv(masked_result$prediction, file.path(PROJECT_DIR, "05_results", "Stage1_4A_masked_ICBScore.csv"), 
    row.names = FALSE)

features <- BA4_BUNDLE$features

masked <- masked_result$prediction

idx <- match(masked$Sample_ID, baseline$Sample_ID)

stopifnot(!anyNA(idx))

baseline_aligned <- baseline[idx, , drop = FALSE]

score_comparison <- data.frame(Sample_ID = masked$Sample_ID, Original_ICBScore = baseline_aligned$ICBScore, 
    Masked_ICBScore = masked$ICBScore)

score_comparison$Delta <- (score_comparison$Masked_ICBScore - score_comparison$Original_ICBScore)

score_comparison$Absolute_Difference <- abs(score_comparison$Delta)


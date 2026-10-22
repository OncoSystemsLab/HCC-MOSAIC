if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

folders <- c("00_config", "01_models", "02_reference", "03_input", "04_processed", "05_results", "06_figures", 
    "07_agent", "08_reports")

for (folder in folders) {
    dir.create(file.path(PROJECT_DIR, folder), recursive = TRUE, showWarnings = FALSE)
}

data.frame(Folder = folders, Exists = dir.exists(file.path(PROJECT_DIR, folders)))

writeLines(capture.output(sessionInfo()), file.path(PROJECT_DIR, "00_config", "Stage0_R_sessionInfo.txt"))

if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("First run install.packages('glmnet')")
}

library(glmnet)

bundle_path <- normalizePath(hcc_input_file("HCC_INPUT_FILE_001"), winslash = "/", mustWork = TRUE)

BA4_BUNDLE <- readRDS(bundle_path)

required_fields <- c("model_name", "model", "lambda", "features", "gene_sets", "background_genes", "feature_reference", 
    "preprocessing", "n_training")

missing_fields <- setdiff(required_fields, names(BA4_BUNDLE))

if (length(missing_fields) > 0) {
    stop("The bundle is missing required fields: ", paste(missing_fields, collapse = ", "))
}

EXPECTED_FEATURES <- c("ICI__HCC_INFLAMMATORY4", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_BILE_ACID_METABOLISM", 
    "HALLMARK_FATTY_ACID_METABOLISM")

stopifnot(identical(BA4_BUNDLE$features, EXPECTED_FEATURES))

stopifnot(inherits(BA4_BUNDLE$model, "glmnet"), length(BA4_BUNDLE$lambda) == 1, is.finite(BA4_BUNDLE$lambda), 
    BA4_BUNDLE$lambda > 0)

background <- BA4_BUNDLE$background_genes

gene_sets <- BA4_BUNDLE$gene_sets

stopifnot(is.character(background), !anyNA(background), !anyDuplicated(background), all(EXPECTED_FEATURES %in% 
    names(gene_sets)))

gene_set_QC <- data.frame(Feature = EXPECTED_FEATURES, N_genes = vapply(gene_sets[EXPECTED_FEATURES], 
    length, integer(1)), N_outside_background = vapply(gene_sets[EXPECTED_FEATURES], function(gs) {
    length(setdiff(gs, background))
}, integer(1)))

stopifnot(all(gene_set_QC$N_genes >= 2), all(gene_set_QC$N_outside_background == 0))

pre <- BA4_BUNDLE$preprocessing

stopifnot(identical(pre$expression_input, "TPM"), identical(pre$transformation, "log2(TPM + 1)"), identical(pre$enrichment, 
    "ssGSEA"), isTRUE(all.equal(pre$ssgsea_alpha, 0.25)), identical(pre$ssgsea_normalize, FALSE))

model_coef <- as.matrix(coef(BA4_BUNDLE$model, s = BA4_BUNDLE$lambda))

stopifnot(BA4_BUNDLE$n_training == 69)

write.csv(gene_set_QC, file.path(PROJECT_DIR, "08_reports", "Stage0_BA4_gene_set_QC.csv"), row.names = FALSE)

write.csv(data.frame(Term = rownames(model_coef), Coefficient = as.numeric(model_coef)), file.path(PROJECT_DIR, 
    "08_reports", "Stage0_BA4_model_coefficients.csv"), row.names = FALSE)

model_copy_path <- file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds")

if (!file.exists(model_copy_path)) {
    ok <- file.copy(bundle_path, model_copy_path)
    if (!ok) 
        stop("Model backup failed")
}

if (!identical(unname(tools::md5sum(bundle_path)), unname(tools::md5sum(model_copy_path)))) {
    stop("The model backup differs from the original file; review before continuing")
}

source_dir <- dirname(bundle_path)

features_file <- file.path(source_dir, "HCC_3cohort_ssGSEA_57features.rds")

original_score_file <- file.path(source_dir, "FINAL_BA4_training_scores_69patients.csv")

if (file.exists(features_file)) {
    feature_data <- readRDS(features_file)
    stopifnot(all(EXPECTED_FEATURES %in% colnames(feature_data)), "patient_uid" %in% colnames(feature_data), 
        !anyDuplicated(feature_data$patient_uid))
    X_test <- as.matrix(feature_data[, EXPECTED_FEATURES, drop = FALSE])
    storage.mode(X_test) <- "numeric"
    stopifnot(all(is.finite(X_test)))
    scores_reproduced <- as.numeric(predict(BA4_BUNDLE$model, newx = X_test, type = "response", s = BA4_BUNDLE$lambda))
    reproduced_df <- data.frame(patient_uid = feature_data$patient_uid, BA4_ReproducedScore = scores_reproduced)
    print(head(reproduced_df))
    if (file.exists(original_score_file)) {
        original_df <- read.csv(original_score_file, check.names = FALSE)
        stopifnot(all(c("patient_uid", "BA4_FinalScore") %in% colnames(original_df)), !anyDuplicated(original_df$patient_uid))
        idx <- match(reproduced_df$patient_uid, original_df$patient_uid)
        stopifnot(!anyNA(idx))
        max_difference <- max(abs(reproduced_df$BA4_ReproducedScore - original_df$BA4_FinalScore[idx]))
        cat("\nMaximum score difference:", max_difference, "")
        stopifnot(max_difference < 1e-08)
        cat("BA4 prediction reproduction: PASS")
    }
    else {
        cat("\nOriginal score CSV not found.")
        cat("Prediction generated, comparison pending.")
    }
    write.csv(reproduced_df, file.path(PROJECT_DIR, "05_results", "Stage0_BA4_reproduced_scores.csv"), 
        row.names = FALSE)
} else {
    cat("\nThe original 57-feature RDS is not in the bundle directory.")
    cat("Stage 0.3 comparison skipped.")
    cat("Use the Stage 0.2 results to assess model integrity.")
}


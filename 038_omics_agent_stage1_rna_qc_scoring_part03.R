if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

for (i in seq_len(nrow(selected_patients))) {
    patient <- selected_patients[i, ]
    uid <- as.character(patient$patient_uid)
    cohort <- as.character(patient$Cohort)
    sample_id <- as.character(patient$Raw_Sample_ID)
    cat("\nProcessing", i, "/", nrow(selected_patients), "|", cohort, "|", sample_id, "")
    if (!cohort %in% names(tpm_background)) {
        stop("Cohort not found: ", cohort)
    }
    mat <- as.matrix(tpm_background[[cohort]])
    storage.mode(mat) <- "numeric"
    if (!sample_id %in% colnames(mat)) {
        stop("Patient not found: ", uid)
    }
    if (!all(bg %in% rownames(mat))) {
        stop("Training gene background is incomplete: ", uid)
    }
    if (any(!is.finite(mat[, sample_id])) || any(mat[, sample_id] < 0)) {
        stop("TPM contains nonfinite or negative values: ", uid)
    }
    one_tpm <- mat[bg, sample_id, drop = FALSE]
    prediction <- suppressMessages(predict_BA4_from_TPM(one_tpm, bundle = BA4_BUNDLE))
    new_score <- prediction$prediction$ICBScore[1]
    score_idx <- match(uid, reference_map$patient_uid)
    if (is.na(score_idx)) {
        stop("Patient not found in the original model scores: ", uid)
    }
    old_score <- reference_map$BA4_ReproducedScore[score_idx]
    feature_idx <- match(uid, feature_reference$patient_uid)
    if (is.na(feature_idx)) {
        stop("Patient not found in the original features: ", uid)
    }
    old_features <- as.numeric(feature_reference[feature_idx, features, drop = TRUE])
    new_features <- as.numeric(prediction$features[1, features])
    feature_diff <- abs(new_features - old_features)
    score_results[[i]] <- data.frame(patient_uid = uid, Cohort = cohort, Raw_Sample_ID = sample_id, Original_ICBScore = old_score, 
        Recalculated_ICBScore = new_score, Score_Abs_Error = abs(new_score - old_score), Max_Feature_Abs_Error = max(feature_diff))
    feature_results[[i]] <- data.frame(patient_uid = uid, Cohort = cohort, Feature = features, Original = old_features, 
        Recalculated = new_features, Absolute_Difference = feature_diff)
}

score_validation <- do.call(rbind, score_results)

feature_validation <- do.call(rbind, feature_results)

rownames(score_validation) <- NULL

rownames(feature_validation) <- NULL

cohort_summary <- aggregate(cbind(Score_Abs_Error, Max_Feature_Abs_Error) ~ Cohort, data = score_validation, 
    FUN = max)

cohort_summary$N_Patients <- as.integer(table(score_validation$Cohort)[cohort_summary$Cohort])

max_score_error <- max(score_validation$Score_Abs_Error)

max_feature_error <- max(feature_validation$Absolute_Difference)

write.csv(score_validation, file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cross_cohort_score_validation.csv"), 
    row.names = FALSE)

write.csv(feature_validation, file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cross_cohort_feature_validation.csv"), 
    row.names = FALSE)

write.csv(cohort_summary, file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cross_cohort_summary.csv"), 
    row.names = FALSE)

if (nrow(score_validation) == 11 && length(unique(score_validation$Cohort)) == 3 && max_score_error < 
    1e-08 && max_feature_error < 1e-08) {
    cat("\nStage 1.5A.3: PASS")
} else {
    cat("\nStage 1.5A.3: CHECK REQUIRED")
}

library(GSVA)

library(BiocParallel)

library(glmnet)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

features <- BA4_BUNDLE$features

bg <- BA4_BUNDLE$background_genes

if (!exists("predict_BA4_from_TPM")) {
    stop("The predict_BA4_from_TPM function is missing.", "Run the function definitions in Stage 1.1.1 first.")
}

tpm_cache <- file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds")

if (file.exists(tpm_cache)) {
    tpm_reference <- readRDS(tpm_cache)
} else if (exists("tpm_background")) {
    tpm_reference <- tpm_background
    saveRDS(tpm_reference, tpm_cache, compress = "gzip")
} else {
    stop("TPM reference data for the three cohorts were not found.", "Run Stage 1.5A.3 first.")
}

stopifnot(is.list(tpm_reference), all(c("PRJEB34724", "GSE181946", "GSE202069") %in% names(tpm_reference)))

feature_cache <- file.path(PROJECT_DIR, "02_reference", "Stage1_training_57features.rds")

if (file.exists(feature_cache)) {
    feature_ref <- readRDS(feature_cache)
} else if (exists("features_file") && file.exists(features_file)) {
    feature_ref <- readRDS(features_file)
    saveRDS(feature_ref, feature_cache)
} else if (exists("original_features")) {
    feature_ref <- original_features
    saveRDS(feature_ref, feature_cache)
} else {
    cat("\nSelect the required source file: ", "HCC_3cohort_ssGSEA_57features.rds")
    feature_ref <- readRDS(hcc_input_file("HCC_INPUT_FILE_003"))
    saveRDS(feature_ref, feature_cache)
}

stopifnot("patient_uid" %in% colnames(feature_ref), all(features %in% colnames(feature_ref)), !anyDuplicated(feature_ref$patient_uid))

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

stopifnot(nrow(reference_map) == 69, !anyDuplicated(reference_map$patient_uid), all(reference_map$patient_uid %in% 
    feature_ref$patient_uid))

score_results <- vector("list", nrow(reference_map))

feature_results <- vector("list", nrow(reference_map))

for (i in seq_len(nrow(reference_map))) {
    uid <- reference_map$patient_uid[i]
    cohort <- reference_map$Cohort[i]
    sample_id <- reference_map$Raw_Sample_ID[i]
    cat("\nPatient", i, "/69", "|", uid, "")
    mat <- as.matrix(tpm_reference[[cohort]])
    storage.mode(mat) <- "numeric"
    stopifnot(sample_id %in% colnames(mat), !anyDuplicated(rownames(mat)), all(bg %in% rownames(mat)))
    one_tpm <- mat[bg, sample_id, drop = FALSE]
    stopifnot(all(is.finite(one_tpm)), all(one_tpm >= 0))
    prediction <- suppressWarnings(suppressMessages(predict_BA4_from_TPM(one_tpm, bundle = BA4_BUNDLE)))
    new_score <- prediction$prediction$ICBScore[1]
    old_score <- reference_map$BA4_ReproducedScore[i]
    fi <- match(uid, feature_ref$patient_uid)
    old_features <- as.numeric(unlist(feature_ref[fi, features, drop = FALSE], use.names = FALSE))
    new_features <- as.numeric(prediction$features[1, features])
    feature_error <- abs(new_features - old_features)
    score_results[[i]] <- data.frame(patient_uid = uid, Cohort = cohort, Original_ICBScore = old_score, 
        Recalculated_ICBScore = new_score, Score_Abs_Error = abs(new_score - old_score), Max_Feature_Abs_Error = max(feature_error))
    feature_results[[i]] <- data.frame(patient_uid = uid, Cohort = cohort, Feature = features, Original = old_features, 
        Recalculated = new_features, Absolute_Difference = feature_error)
}

score_validation_69 <- do.call(rbind, score_results)

feature_validation_69 <- do.call(rbind, feature_results)

rownames(score_validation_69) <- NULL

rownames(feature_validation_69) <- NULL

cohort_summary_69 <- aggregate(cbind(Score_Abs_Error, Max_Feature_Abs_Error) ~ Cohort, data = score_validation_69, 
    FUN = max)

n_by_cohort <- table(score_validation_69$Cohort)

cohort_summary_69$N_Patients <- as.integer(n_by_cohort[cohort_summary_69$Cohort])

max_score_error <- max(score_validation_69$Score_Abs_Error)

max_feature_error <- max(feature_validation_69$Absolute_Difference)

write.csv(score_validation_69, file.path(PROJECT_DIR, "08_reports", "Stage1_5B_69_score_validation.csv"), 
    row.names = FALSE)

write.csv(feature_validation_69, file.path(PROJECT_DIR, "08_reports", "Stage1_5B_69_feature_validation.csv"), 
    row.names = FALSE)

write.csv(cohort_summary_69, file.path(PROJECT_DIR, "08_reports", "Stage1_5B_cohort_summary.csv"), row.names = FALSE)

passed <- (nrow(score_validation_69) == 69 && nrow(feature_validation_69) == 276 && length(unique(score_validation_69$Cohort)) == 
    3 && max_score_error < 1e-08 && max_feature_error < 1e-08)

if (passed) {
    cat("\nStage 1.5B: PASS")
} else {
    cat("\nStage 1.5B: CHECK REQUIRED")
}

writeLines(capture.output(sessionInfo()), file.path(PROJECT_DIR, "00_config", "Stage1_final_sessionInfo.txt"))


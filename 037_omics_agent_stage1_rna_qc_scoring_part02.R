if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

score_comparison$Difference_Percentage_Points <- (score_comparison$Delta * 100)

feature_comparison <- do.call(rbind, lapply(features, function(feature) {
    old <- baseline_aligned[[feature]]
    new <- masked[[feature]]
    data.frame(Sample_ID = masked$Sample_ID, Feature = feature, Original = old, Masked = new, Delta = new - 
        old, Absolute_Difference = abs(new - old), Relative_Difference_Percent = 100 * abs(new - old)/pmax(abs(old), 
        1e-08))
}))

rownames(feature_comparison) <- NULL

write.csv(score_comparison, file.path(PROJECT_DIR, "08_reports", "Stage1_4A_ICBScore_sensitivity.csv"), 
    row.names = FALSE)

write.csv(feature_comparison, file.path(PROJECT_DIR, "08_reports", "Stage1_4A_feature_sensitivity.csv"), 
    row.names = FALSE)

features <- BA4_BUNDLE$features

coef_mat <- as.matrix(coef(BA4_BUNDLE$model, s = BA4_BUNDLE$lambda))

beta <- setNames(as.numeric(coef_mat[features, 1]), features)

masked <- masked_result$prediction

idx <- match(masked$Sample_ID, baseline$Sample_ID)

stopifnot(!anyNA(idx))

original <- baseline[idx, , drop = FALSE]

delta_features <- (as.matrix(masked[, features, drop = FALSE]) - as.matrix(original[, features, drop = FALSE]))

contributions <- sweep(delta_features, 2, beta, FUN = "*")

result <- data.frame(Sample_ID = masked$Sample_ID, Expected_Delta_Logit = rowSums(contributions), Actual_Delta_Logit = qlogis(masked$ICBScore) - 
    qlogis(original$ICBScore), row.names = NULL)

result$Difference <- (result$Expected_Delta_Logit - result$Actual_Delta_Logit)

stopifnot(max(abs(result$Difference)) < 1e-08)

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("read_HCC_TPM"), exists("predict_BA4_missing_bg_experiment"))

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

new_tpm <- read_HCC_TPM(file.path(PROJECT_DIR, "03_input", "Stage1_demo_2patients_TPM.csv"))

baseline <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage1_2_batch_ICBScore.csv"), check.names = FALSE)

features <- BA4_BUNDLE$features

bg <- BA4_BUNDLE$background_genes

stopifnot(all(bg %in% rownames(new_tpm)), all(colnames(new_tpm) %in% baseline$Sample_ID))

missing_rates <- c(0.005, 0.01, 0.02, 0.05)

N_REPEATS <- 10L

set.seed(20261011)

BA4_genes <- unique(unlist(BA4_BUNDLE$gene_sets, use.names = FALSE))

non_signature_genes <- setdiff(bg, BA4_genes)

stopifnot(length(non_signature_genes) >= round(max(missing_rates) * length(bg)))

all_results <- vector("list", length(missing_rates) * N_REPEATS)

counter <- 0L

for (rate in missing_rates) {
    n_remove <- round(length(bg) * rate)
    for (rep in seq_len(N_REPEATS)) {
        counter <- counter + 1L
        cat("\nExperiment:", counter, "/", length(all_results), "| Missing:", rate * 100, "% | Repeat:", 
            rep, "")
        removed <- sample(non_signature_genes, size = n_remove)
        masked_tpm <- new_tpm[!rownames(new_tpm) %in% removed, , drop = FALSE]
        invisible(capture.output({
            experiment <- suppressWarnings(suppressMessages(predict_BA4_missing_bg_experiment(masked_tpm, 
                bundle = BA4_BUNDLE)))
        }))
        pred <- experiment$prediction
        idx <- match(pred$Sample_ID, baseline$Sample_ID)
        stopifnot(!anyNA(idx))
        original <- baseline[idx, , drop = FALSE]
        score_delta <- (pred$ICBScore - original$ICBScore)
        old_X <- as.matrix(original[, features, drop = FALSE])
        new_X <- as.matrix(pred[, features, drop = FALSE])
        relative_feature_error <- (100 * abs(new_X - old_X)/pmax(abs(old_X), 1e-08))
        all_results[[counter]] <- data.frame(Missing_Rate = rate, Missing_Percent = rate * 100, Missing_Genes = n_remove, 
            Repeat = rep, Sample_ID = pred$Sample_ID, Original_ICBScore = original$ICBScore, New_ICBScore = pred$ICBScore, 
            Score_Delta = score_delta, Abs_Score_Delta = abs(score_delta), Max_Feature_Relative_Error = apply(relative_feature_error, 
                1, max), Flip_at_0_5 = (pred$ICBScore >= 0.5) != (original$ICBScore >= 0.5))
    }
}

simulation_results <- do.call(rbind, all_results)

rownames(simulation_results) <- NULL

write.csv(simulation_results, file.path(PROJECT_DIR, "08_reports", "Stage1_4B_simulation_results.csv"), 
    row.names = FALSE)

groups <- split(simulation_results, simulation_results$Missing_Rate)

summary_df <- do.call(rbind, lapply(groups, function(df) {
    data.frame(Missing_Percent = unique(df$Missing_Percent), N_Records = nrow(df), Median_Abs_Error = median(df$Abs_Score_Delta), 
        Mean_Abs_Error = mean(df$Abs_Score_Delta), P95_Abs_Error = as.numeric(quantile(df$Abs_Score_Delta, 
            probs = 0.95)), Max_Abs_Error = max(df$Abs_Score_Delta), Median_Max_Feature_Error = median(df$Max_Feature_Relative_Error), 
        N_Flip_0_5 = sum(df$Flip_at_0_5))
}))

rownames(summary_df) <- NULL

summary_df <- summary_df[order(summary_df$Missing_Percent), , drop = FALSE]

write.csv(summary_df, file.path(PROJECT_DIR, "08_reports", "Stage1_4B_sensitivity_summary.csv"), row.names = FALSE)

simulation_results$Rate_Label <- factor(paste0(simulation_results$Missing_Percent, "%"), levels = paste0(missing_rates * 
    100, "%"))

figure_file <- file.path(PROJECT_DIR, "06_figures", "Stage1_4B_ICBScore_missingness.png")

png(filename = figure_file, width = 1200, height = 800, res = 150)

boxplot(Abs_Score_Delta ~ Rate_Label, data = simulation_results, xlab = "Missing background genes (%)", 
    ylab = "Absolute ICBScore difference", main = "BA4 ICBScore sensitivity to missing genes", outline = TRUE)

stripchart(Abs_Score_Delta ~ Rate_Label, data = simulation_results, vertical = TRUE, method = "jitter", 
    add = TRUE, pch = 16, cex = 0.55)

dev.off()

library(GSVA)

library(BiocParallel)

library(glmnet)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

new_tpm <- read_HCC_TPM(file.path(PROJECT_DIR, "03_input", "Stage1_demo_2patients_TPM.csv"))

baseline <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage1_2_batch_ICBScore.csv"), check.names = FALSE)

features <- BA4_BUNDLE$features

gene_sets <- BA4_BUNDLE$gene_sets[features]

inflammatory_genes <- gene_sets[[features[1]]]

stopifnot(length(inflammatory_genes) == 4)

scenarios <- list()

for (gene in inflammatory_genes) {
    scenarios[[length(scenarios) + 1L]] <- list(Target_Feature = features[1], Loss_Type = "Inflammatory_single", 
        Removed_Genes = gene)
}

gene_pairs <- combn(inflammatory_genes, 2, simplify = FALSE)

for (pair in gene_pairs) {
    scenarios[[length(scenarios) + 1L]] <- list(Target_Feature = features[1], Loss_Type = "Inflammatory_pair", 
        Removed_Genes = pair)
}

for (feature in features[-1]) {
    genes <- intersect(gene_sets[[feature]], rownames(new_tpm))
    mean_expression <- rowMeans(log2(new_tpm[genes, , drop = FALSE] + 1))
    top_genes <- names(sort(mean_expression, decreasing = TRUE))[seq_len(min(3L, length(genes)))]
    for (gene in top_genes) {
        scenarios[[length(scenarios) + 1L]] <- list(Target_Feature = feature, Loss_Type = "Other_single", 
            Removed_Genes = gene)
    }
}

scenario_info <- do.call(rbind, lapply(seq_along(scenarios), function(i) {
    s <- scenarios[[i]]
    data.frame(Scenario_ID = sprintf("S%02d", i), Target_Feature = s$Target_Feature, Loss_Type = s$Loss_Type, 
        N_Removed = length(s$Removed_Genes), Removed_Genes = paste(s$Removed_Genes, collapse = ";"))
}))

rownames(scenario_info) <- NULL

stopifnot(nrow(scenario_info) == 19)

write.csv(scenario_info, file.path(PROJECT_DIR, "08_reports", "Stage1_4C_experiment_scenarios.csv"), 
    row.names = FALSE)

predict_BA4_core_missing <- function(tpm, bundle = BA4_BUNDLE) {
    bg <- bundle$background_genes
    gs <- bundle$gene_sets[bundle$features]
    pre <- bundle$preprocessing
    available <- bg[bg %in% rownames(tpm)]
    remaining <- vapply(gs, function(x) sum(x %in% available), integer(1))
    if (any(remaining < pre$ssgsea_minSize)) {
        stop("Some gene sets contain fewer genes than the ssGSEA minimum")
    }
    outputs <- vector("list", ncol(tpm))
    for (i in seq_len(ncol(tpm))) {
        expr <- log2(tpm[available, i, drop = FALSE] + 1)
        param <- GSVA::ssgseaParam(exprData = expr, geneSets = gs, minSize = pre$ssgsea_minSize, maxSize = pre$ssgsea_maxSize, 
            alpha = pre$ssgsea_alpha, normalize = pre$ssgsea_normalize)
        es <- GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())
        if (!all(bundle$features %in% rownames(es))) {
            stop("The ssGSEA results are missing required model features")
        }
        X <- t(es[bundle$features, , drop = FALSE])
        score <- as.numeric(predict(bundle$model, newx = X, s = bundle$lambda, type = "response"))
        outputs[[i]] <- data.frame(Sample_ID = colnames(tpm)[i], ICBScore = score, X, check.names = FALSE, 
            row.names = NULL)
    }
    do.call(rbind, outputs)
}

results <- list()

errors <- list()

for (i in seq_along(scenarios)) {
    s <- scenarios[[i]]
    cat("\nRunning scenario", i, "/", length(scenarios), ":", paste(s$Removed_Genes, collapse = ","), 
        "")
    masked_tpm <- new_tpm[!rownames(new_tpm) %in% s$Removed_Genes, , drop = FALSE]
    pred <- tryCatch(suppressMessages(suppressWarnings(predict_BA4_core_missing(masked_tpm))), error = function(e) e)
    if (inherits(pred, "error")) {
        errors[[length(errors) + 1L]] <- data.frame(Scenario_ID = sprintf("S%02d", i), Error = conditionMessage(pred))
        next
    }
    idx <- match(pred$Sample_ID, baseline$Sample_ID)
    stopifnot(!anyNA(idx))
    original <- baseline[idx, , drop = FALSE]
    target <- s$Target_Feature
    results[[length(results) + 1L]] <- data.frame(Scenario_ID = sprintf("S%02d", i), Target_Feature = target, 
        Loss_Type = s$Loss_Type, Removed_Genes = paste(s$Removed_Genes, collapse = ";"), Sample_ID = pred$Sample_ID, 
        Original_ICBScore = original$ICBScore, New_ICBScore = pred$ICBScore, Score_Delta = pred$ICBScore - 
            original$ICBScore, Abs_Score_Delta = abs(pred$ICBScore - original$ICBScore), Target_Feature_Relative_Error = 100 * 
            abs(pred[[target]] - original[[target]])/pmax(abs(original[[target]]), 1e-08))
}

if (length(results) == 0L) {
    stop("No experiments completed successfully; review errors")
}

core_results <- do.call(rbind, results)

rownames(core_results) <- NULL

write.csv(core_results, file.path(PROJECT_DIR, "08_reports", "Stage1_4C_core_gene_results.csv"), row.names = FALSE)

if (length(errors) > 0) {
    print(do.call(rbind, errors))
}

groups <- split(core_results, interaction(core_results$Target_Feature, core_results$Loss_Type, drop = TRUE))

core_summary <- do.call(rbind, lapply(groups, function(df) {
    data.frame(Target_Feature = df$Target_Feature[1], Loss_Type = df$Loss_Type[1], N_Records = nrow(df), 
        Mean_Abs_ICB_Error = mean(df$Abs_Score_Delta), Max_Abs_ICB_Error = max(df$Abs_Score_Delta), Mean_Target_Feature_Error = mean(df$Target_Feature_Relative_Error), 
        Max_Target_Feature_Error = max(df$Target_Feature_Relative_Error))
}))

rownames(core_summary) <- NULL

top_results <- core_results[order(core_results$Abs_Score_Delta, decreasing = TRUE), , drop = FALSE]

write.csv(core_summary, file.path(PROJECT_DIR, "08_reports", "Stage1_4C_core_gene_summary.csv"), row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

reference_scores <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage0_BA4_reproduced_scores.csv"), 
    check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(nrow(reference_scores) == 69, !anyDuplicated(reference_scores$patient_uid), all(is.finite(reference_scores$BA4_ReproducedScore)))

if (!exists("tpm_background")) {
    stop("The current R session is missing tpm_background; ", "load the three-cohort TPM background matrices from the preceding preparation stage.")
}

stopifnot(is.list(tpm_background), !is.null(names(tpm_background)))

bg <- BA4_BUNDLE$background_genes

cohort_inventory <- do.call(rbind, lapply(names(tpm_background), function(cohort) {
    mat <- as.matrix(tpm_background[[cohort]])
    if (is.null(rownames(mat)) || is.null(colnames(mat))) {
        stop("Cohort is missing gene or sample names: ", cohort)
    }
    if (anyDuplicated(colnames(mat))) {
        stop("Duplicate sample names within cohort: ", cohort)
    }
    data.frame(Cohort = cohort, Raw_Sample_ID = colnames(mat), patient_uid = paste0(cohort, "__", colnames(mat)), 
        Input_Genes = nrow(mat), Matched_BG_Genes = sum(bg %in% rownames(mat)), stringsAsFactors = FALSE)
}))

rownames(cohort_inventory) <- NULL

stopifnot(!anyDuplicated(cohort_inventory$patient_uid))

idx <- match(reference_scores$patient_uid, cohort_inventory$patient_uid)

reference_map <- cbind(reference_scores, cohort_inventory[idx, c("Cohort", "Raw_Sample_ID", "Input_Genes", 
    "Matched_BG_Genes"), drop = FALSE])

if (anyNA(idx)) {
    cat("\nUnmatched patient IDs:")
    print(reference_scores$patient_uid[is.na(idx)])
}

write.csv(reference_map, file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), row.names = FALSE)

if (anyNA(reference_map$Cohort)) {
    stop("Unmatched patients found; review Stage 1.5A.1 first")
}

scores <- reference_map$BA4_ReproducedScore

targets <- as.numeric(quantile(scores, probs = c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)))

selected_index <- vapply(targets, function(target) {
    which.min(abs(scores - target))
}, integer(1))

selected_index <- c(selected_index, which.min(abs(scores - 0.5)))

for (cohort in unique(reference_map$Cohort)) {
    cohort_idx <- which(reference_map$Cohort == cohort)
    cohort_median <- median(scores[cohort_idx])
    nearest <- cohort_idx[which.min(abs(scores[cohort_idx] - cohort_median))]
    selected_index <- c(selected_index, nearest)
}

selected_index <- unique(selected_index)

selected_patients <- reference_map[selected_index, , drop = FALSE]

selected_patients <- selected_patients[order(selected_patients$BA4_ReproducedScore), , drop = FALSE]

rownames(selected_patients) <- NULL

write.csv(selected_patients, file.path(PROJECT_DIR, "08_reports", "Stage1_5A_selected_patients.csv"), 
    row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

BA4_BUNDLE <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

if (!exists("predict_BA4_from_TPM")) {
    stop("Run Stage 1.1.1 first to ", "define predict_BA4_from_TPM()")
}

tpm_cache <- file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds")

if (!exists("tpm_background")) {
    if (file.exists(tpm_cache)) {
        tpm_background <- readRDS(tpm_cache)
    }
    else {
        stop("Missing tpm_background; ", "load TPM data for all three cohorts from the preceding preparation stage")
    }
}

if (!file.exists(tpm_cache)) {
    saveRDS(tpm_background, file = tpm_cache, compress = "gzip")
}

selected_patients <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_selected_patients.csv"), 
    check.names = FALSE, stringsAsFactors = FALSE)

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

stopifnot(nrow(selected_patients) == 11, nrow(reference_map) == 69, !anyDuplicated(selected_patients$patient_uid))

if (exists("original_features")) {
    feature_reference <- original_features
} else if (exists("features_file") && file.exists(features_file)) {
    feature_reference <- readRDS(features_file)
} else {
    cat("\nSelect the required source file: ", "HCC_3cohort_ssGSEA_57features.rds")
    feature_reference <- readRDS(hcc_input_file("HCC_INPUT_FILE_002"))
}

features <- BA4_BUNDLE$features

bg <- BA4_BUNDLE$background_genes

stopifnot("patient_uid" %in% colnames(feature_reference), all(features %in% colnames(feature_reference)), 
    !anyDuplicated(feature_reference$patient_uid))

score_results <- vector("list", nrow(selected_patients))

feature_results <- vector("list", nrow(selected_patients))


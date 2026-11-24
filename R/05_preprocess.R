# Export rank features, encoder inputs, and observation masks.
PROJECT_DIR <- CFG$root
INTERSECTION_FILE <- artifact("stage03_intersection_candidate_bundle.rds")
ALIGNED_FILE <- artifact("stage03_ici_aligned_with_mask_bundle.rds")
if (!all(file.exists(c(INTERSECTION_FILE, ALIGNED_FILE)))) {
    stop("\u4EA4\u96C6\u77E9\u9635\u6216\u5E94\u7B54\u961F\u5217\u5BF9\u9F50\u7ED3\u679C\u4E0D\u5B58\u5728\u3002")
}
intersection_bundle <- readRDS(INTERSECTION_FILE)
aligned_bundle <- readRDS(ALIGNED_FILE)
current_intersection_md5 <- unname(tools::md5sum(INTERSECTION_FILE))
if (current_intersection_md5 != aligned_bundle$provenance$intersection_md5) {
    stop("\u4EA4\u96C6\u6587\u4EF6\u4E0E\u5E94\u7B54\u961F\u5217\u5BF9\u9F50\u65F6\u4F7F\u7528\u7684\u7248\u672C\u4E0D\u4E00\u81F4\u3002")
}
fixed_genes <- intersection_bundle$intersection_genes
ssl_log <- intersection_bundle$ssl_expression_log2tpm
ssl_meta <- intersection_bundle$ssl_patient_metadata
stopifnot(identical(rownames(ssl_log), fixed_genes), identical(colnames(ssl_log), ssl_meta$patient_uid), identical(aligned_bundle$fixed_genes, 
    fixed_genes), ncol(ssl_log) == 616L, all(is.finite(ssl_log)))
PREP_DIR <- out_dir("05_preprocess")
if (dir.exists(PREP_DIR) || !dir.create(PREP_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_prep_tsv <- function(x, filename) {
    write.table(x, file.path(PREP_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
SPLIT_SEED <- 20261116L
VALIDATION_FRACTION <- 0.1
SD_FLOOR <- 0.1
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(SPLIT_SEED)
ssl_meta$split <- "train"
for (cohort in c("TCGA", "ICGC", "CPTAC")) {
    cohort_rows <- which(ssl_meta$cohort == cohort)
    cohort_rows <- cohort_rows[order(ssl_meta$patient_uid[cohort_rows], method = "radix")]
    n_validation <- max(1L, as.integer(round(length(cohort_rows) * VALIDATION_FRACTION)))
    validation_rows <- sample(cohort_rows, size = n_validation, replace = FALSE)
    ssl_meta$split[validation_rows] <- "validation"
}
train_index <- which(ssl_meta$split == "train")
validation_index <- which(ssl_meta$split == "validation")
stopifnot(length(intersect(train_index, validation_index)) == 0L, length(train_index) + length(validation_index) == 
    616L, !anyDuplicated(ssl_meta$patient_uid))
train_log <- ssl_log[, train_index, drop = FALSE]
gene_mean <- rowMeans(train_log)
gene_sd_raw <- apply(train_log, 1L, sd)
gene_median <- apply(train_log, 1L, median)
gene_sd <- pmax(gene_sd_raw, SD_FLOOR)
if (any(!is.finite(gene_mean)) || any(!is.finite(gene_sd)) || any(!is.finite(gene_median))) {
    stop("\u8BAD\u7EC3\u96C6\u9884\u5904\u7406\u53C2\u6570\u5305\u542B\u65E0\u6548\u503C\u3002")
}
preprocessing_parameters <- data.frame(ensembl_gene_id = fixed_genes, gene_symbol = intersection_bundle$gene_manifest$gene_symbol, 
    train_mean_log2tpm = gene_mean, train_sd_log2tpm = gene_sd_raw, scale_used = gene_sd, train_median_log2tpm = gene_median, 
    sd_floor_applied = gene_sd_raw < SD_FLOOR, zero_in_all_SSL_training_samples = rowSums(train_log > 0) == 0L, 
    stringsAsFactors = FALSE)
write_prep_tsv(ssl_meta, "SSL_candidate_train_validation_split.tsv")
write_prep_tsv(preprocessing_parameters, "SSL_training_preprocessing_parameters.tsv")
cat("\n===== SSL \u5019\u9009\u5212\u5206 =====\n")
print(with(ssl_meta, table(cohort, split)))
cat("\n\u6807\u51C6\u5DEE\u4E0B\u9650\u4F5C\u7528\u7684\u57FA\u56E0\u6570\uFF1A", sum(preprocessing_parameters$sd_floor_applied), "\n")
cat("SSL \u8BAD\u7EC3\u5B50\u96C6\u4E2D\u5168\u90E8\u4E3A\u96F6\u7684\u57FA\u56E0\u6570\uFF1A", sum(preprocessing_parameters$zero_in_all_SSL_training_samples), 
    "\n")
transform_expression <- function(log_expression, observed_mask) {
    stopifnot(identical(rownames(log_expression), fixed_genes), identical(dimnames(log_expression), dimnames(observed_mask)))
    if (anyNA(observed_mask)) {
        stop("\u89C2\u6D4B\u63A9\u7801\u5305\u542B NA\u3002")
    }
    if (any(!is.finite(log_expression[observed_mask]))) {
        stop("\u5DF2\u89C2\u6D4B\u4F4D\u7F6E\u5B58\u5728\u65E0\u6548\u8868\u8FBE\u503C\u3002")
    }
    if (any(!is.na(log_expression[!observed_mask]))) {
        stop("\u7F3A\u6D4B\u4F4D\u7F6E\u5E94\u4FDD\u7559\u4E3A NA\uFF0C\u8F93\u5165\u72B6\u6001\u4E0D\u7B26\u5408\u9884\u671F\u3002")
    }
    filled <- log_expression
    for (j in seq_len(ncol(filled))) {
        missing <- !observed_mask[, j]
        filled[missing, j] <- gene_median[missing]
    }
    z <- sweep(filled, 1L, gene_mean, "-")
    z <- sweep(z, 1L, gene_sd, "/")
    if (any(!is.finite(z))) {
        stop("\u9884\u5904\u7406\u540E\u4ECD\u5B58\u5728\u65E0\u6548\u503C\u3002")
    }
    z
}
make_sample_qc <- function(cohort, log_expression, observed_mask, standardized, sample_ids) {
    rows <- lapply(seq_len(ncol(log_expression)), function(j) {
        observed <- observed_mask[, j]
        values <- log_expression[observed, j]
        z_observed <- standardized[observed, j]
        if (length(values) == 0L) {
            stop("\u5B58\u5728\u5B8C\u5168\u6CA1\u6709\u89C2\u6D4B\u57FA\u56E0\u7684\u6837\u672C\u3002")
        }
        data.frame(cohort = if (length(cohort) == 1L) 
            cohort
        else cohort[j], sample_id = sample_ids[j], n_observed_genes = sum(observed), n_missing_genes = sum(!observed), 
            n_detected_observed_genes = sum(values > 0), zero_fraction_among_observed = mean(values == 0), median_observed_log2tpm = median(values), 
            p99_observed_log2tpm = unname(quantile(values, probs = 0.99, names = FALSE)), mean_abs_z_among_observed = mean(abs(z_observed)), 
            fraction_observed_abs_z_gt5 = mean(abs(z_observed) > 5), stringsAsFactors = FALSE)
    })
    do.call(rbind, rows)
}
ssl_mask <- matrix(TRUE, nrow = nrow(ssl_log), ncol = ncol(ssl_log), dimnames = dimnames(ssl_log))
ssl_z <- transform_expression(ssl_log, ssl_mask)
if (max(abs(rowMeans(ssl_z[, train_index, drop = FALSE]))) > 1e-08) {
    stop("\u8BAD\u7EC3\u96C6\u6807\u51C6\u5316\u540E\u7684\u57FA\u56E0\u5747\u503C\u672A\u63A5\u8FD1\u96F6\uFF0C\u8BF7\u68C0\u67E5\u3002")
}
ssl_sample_qc <- make_sample_qc(cohort = ssl_meta$cohort, log_expression = ssl_log, observed_mask = ssl_mask, standardized = ssl_z, 
    sample_ids = ssl_meta$patient_uid)
ssl_sample_qc$split <- ssl_meta$split
ici_prepared <- list()
ici_qc_rows <- list()
for (cohort in names(aligned_bundle$ici_aligned)) {
    obj <- aligned_bundle$ici_aligned[[cohort]]
    z <- transform_expression(obj$expression_log2tpm, obj$observed_mask)
    ici_prepared[[cohort]] <- list(expression_z = z, observed_mask = obj$observed_mask, metadata = obj$metadata)
    ici_qc_rows[[cohort]] <- make_sample_qc(cohort = cohort, log_expression = obj$expression_log2tpm, observed_mask = obj$observed_mask, 
        standardized = z, sample_ids = obj$metadata$sample_id)
}
ici_sample_qc <- do.call(rbind, ici_qc_rows)
rownames(ici_sample_qc) <- NULL
write_prep_tsv(ssl_sample_qc, "SSL_sample_preprocessing_QC.tsv")
write_prep_tsv(ici_sample_qc, "ICI_sample_preprocessing_QC.tsv")
saveRDS(list(fixed_genes = fixed_genes, gene_manifest = intersection_bundle$gene_manifest, ssl_expression_z = ssl_z, 
    ssl_patient_metadata = ssl_meta, train_index = train_index, validation_index = validation_index, preprocessing_parameters = preprocessing_parameters, 
    ici_prepared = ici_prepared, token_gene_indices = aligned_bundle$token_gene_indices, token_gene_ids = aligned_bundle$token_gene_ids, 
    ssl_sample_qc = ssl_sample_qc, ici_sample_qc = ici_sample_qc, settings = list(split_seed = SPLIT_SEED, validation_fraction = VALIDATION_FRACTION, 
        transform = "log2(TPM+1)", scaler_fit_scope = "SSL_training_subset_only", missing_fill = "SSL_training_gene_median_log2tpm", 
        sd_floor = SD_FLOOR, clipping_applied = FALSE, response_labels_used_for_preprocessing = FALSE, patient_identity_audit_complete = FALSE, 
        final_training_input_locked = FALSE), provenance = list(intersection_file = INTERSECTION_FILE, intersection_md5 = current_intersection_md5, 
        aligned_file = ALIGNED_FILE, aligned_md5 = unname(tools::md5sum(ALIGNED_FILE)))), artifact("stage04_SSL_prepared_candidate_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(PREP_DIR, "sessionInfo.txt"))
combined_qc <- rbind(ssl_sample_qc[setdiff(names(ssl_sample_qc), "split")], ici_sample_qc)
qc_summary_rows <- lapply(unique(combined_qc$cohort), function(cohort) {
    q <- combined_qc[combined_qc$cohort == cohort, , drop = FALSE]
    data.frame(cohort = cohort, n_samples = nrow(q), min_detected = min(q$n_detected_observed_genes), max_detected = max(q$n_detected_observed_genes), 
        min_mean_abs_z = min(q$mean_abs_z_among_observed), max_mean_abs_z = max(q$mean_abs_z_among_observed), max_fraction_abs_z_gt5 = max(q$fraction_observed_abs_z_gt5), 
        stringsAsFactors = FALSE)
})
preprocessing_qc_summary <- do.call(rbind, qc_summary_rows)
write_prep_tsv(preprocessing_qc_summary, "preprocessing_QC_summary.tsv")
cat("\n===== \u516D\u4E2A\u961F\u5217\u9884\u5904\u7406\u8D28\u63A7\u8303\u56F4 =====\n")
print(preprocessing_qc_summary, row.names = FALSE)
cat("\n===== \u5E94\u7B54\u961F\u5217\u504F\u79BB\u53C2\u8003\u8F83\u5927\u7684\u6837\u672C\uFF1A\u524D 10 \u4F8B =====\n")
print(head(ici_sample_qc[order(-ici_sample_qc$mean_abs_z_among_observed), , drop = FALSE], 10L), row.names = FALSE)
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", PREP_DIR, "\n")

# Export rank features, encoder inputs, and observation masks.
PROJECT_DIR <- CFG$root
PREP_FILE <- artifact("stage04_SSL_prepared_candidate_bundle.rds")
if (!file.exists(PREP_FILE)) {
    stop("\u627E\u4E0D\u5230\u4FEE\u590D\u540E\u7684\u5019\u9009\u8F93\u5165\u5305\u3002")
}
prepared <- readRDS(PREP_FILE)
PCA_DIR <- out_dir("05_PCA_QC")
if (dir.exists(PCA_DIR) || !dir.create(PCA_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_pca_tsv <- function(x, filename) {
    write.table(x, file.path(PCA_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
ssl_z <- prepared$ssl_expression_z
ssl_meta <- prepared$ssl_patient_metadata
parameters <- prepared$preprocessing_parameters
train_index <- prepared$train_index
validation_index <- prepared$validation_index
stopifnot(identical(rownames(ssl_z), prepared$fixed_genes), identical(colnames(ssl_z), ssl_meta$patient_uid), identical(parameters$ensembl_gene_id, 
    prepared$fixed_genes), all(ssl_meta$split[train_index] == "train"), all(ssl_meta$split[validation_index] == 
    "validation"))
eligible <- which(parameters$train_sd_log2tpm > 0)
ranked <- eligible[order(-parameters$train_sd_log2tpm[eligible], parameters$ensembl_gene_id[eligible], method = "radix")]
pca_gene_index <- head(ranked, 2000L)
if (length(pca_gene_index) < 10L) {
    stop("PCA \u53EF\u7528\u53D8\u5F02\u57FA\u56E0\u8FC7\u5C11\u3002")
}
pca_genes <- prepared$fixed_genes[pca_gene_index]
train_matrix <- t(ssl_z[pca_gene_index, train_index, drop = FALSE])
pca_fit <- prcomp(train_matrix, center = TRUE, scale. = FALSE, rank. = 10L)
explained_fraction <- pca_fit$sdev^2/sum(pca_fit$sdev^2)
ssl_scores <- predict(pca_fit, newdata = t(ssl_z[pca_gene_index, , drop = FALSE]))
ssl_score_table <- data.frame(cohort = ssl_meta$cohort, sample_id = ssl_meta$patient_id, patient_uid = ssl_meta$patient_uid, 
    data_role = ifelse(ssl_meta$split == "train", "SSL_train", "SSL_validation"), ssl_scores, check.names = FALSE, 
    stringsAsFactors = FALSE)
ici_score_rows <- list()
for (cohort in names(prepared$ici_prepared)) {
    obj <- prepared$ici_prepared[[cohort]]
    stopifnot(identical(rownames(obj$expression_z), prepared$fixed_genes), identical(colnames(obj$expression_z), 
        obj$metadata$sample_id))
    scores <- predict(pca_fit, newdata = t(obj$expression_z[pca_gene_index, , drop = FALSE]))
    ici_score_rows[[cohort]] <- data.frame(cohort = cohort, sample_id = obj$metadata$sample_id, patient_uid = obj$metadata$patient_uid, 
        data_role = "ICI_projection", scores, check.names = FALSE, stringsAsFactors = FALSE)
}
ici_score_table <- do.call(rbind, ici_score_rows)
all_scores <- rbind(ssl_score_table, ici_score_table)
rownames(all_scores) <- NULL
stopifnot(nrow(all_scores) == 683L, !anyDuplicated(all_scores$patient_uid))
n_distance_pc <- min(10L, ncol(pca_fit$rotation))
pc_names <- paste0("PC", seq_len(n_distance_pc))
pc_scale <- pca_fit$sdev[seq_len(n_distance_pc)]
if (any(!is.finite(pc_scale)) || any(pc_scale <= 0)) {
    stop("PCA \u4E3B\u6210\u5206\u5C3A\u5EA6\u65E0\u6548\u3002")
}
scaled_scores <- sweep(as.matrix(all_scores[pc_names]), 2L, pc_scale, "/")
all_scores$PCA_reference_distance <- sqrt(rowSums(scaled_scores^2))
write_pca_tsv(all_scores, "all_samples_PCA_scores.tsv")
write_pca_tsv(data.frame(ensembl_gene_id = pca_genes, stringsAsFactors = FALSE), "PCA_training_selected_genes.tsv")
variance_report <- data.frame(PC = paste0("PC", seq_len(ncol(pca_fit$rotation))), explained_fraction = explained_fraction[seq_len(ncol(pca_fit$rotation))], 
    stringsAsFactors = FALSE)
write_pca_tsv(variance_report, "PCA_explained_variance.tsv")
saveRDS(list(pca_fit = pca_fit, pca_genes = pca_genes, all_scores = all_scores, explained_variance = variance_report, 
    prepared_file = PREP_FILE, prepared_md5 = unname(tools::md5sum(PREP_FILE)), purpose = "QC_only", sample_exclusion_applied = FALSE), 
    file.path(PCA_DIR, "stage04_PCA_QC_bundle.rds"))
cat("\n===== \u524D 10 \u4E2A\u4E3B\u6210\u5206\u89E3\u91CA\u7387 =====\n")
print(variance_report, row.names = FALSE)
cat("\n===== \u5E94\u7B54\u961F\u5217\u504F\u79BB\u8BAD\u7EC3\u53C2\u8003\u8F83\u5927\u7684\u6837\u672C =====\n")
ici_rows <- all_scores[all_scores$data_role == "ICI_projection", , drop = FALSE]
print(head(ici_rows[order(-ici_rows$PCA_reference_distance), c("cohort", "sample_id", "PC1", "PC2", "PCA_reference_distance")], 
    10L), row.names = FALSE)
cohort_order <- c("TCGA", "ICGC", "CPTAC", "GSE181946", "GSE202069", "PRJEB34724")
cohort_colors <- c(TCGA = "#4477AA", ICGC = "#228833", CPTAC = "#AA3377", GSE181946 = "#EE7733", GSE202069 = "#66CCEE", 
    PRJEB34724 = "#CCBB44")
axis_label <- function(pc) {
    paste0("PC", pc, " (", round(100 * explained_fraction[pc], 1), "% training variance)")
}
draw_pca <- function(scores, title) {
    point_shape <- ifelse(scores$data_role == "SSL_train", 16L, ifelse(scores$data_role == "SSL_validation", 1L, 
        17L))
    plot(scores$PC1, scores$PC2, col = unname(cohort_colors[scores$cohort]), pch = point_shape, cex = 0.8, xlab = axis_label(1L), 
        ylab = axis_label(2L), main = title)
    present_cohorts <- cohort_order[cohort_order %in% scores$cohort]
    legend("topright", legend = present_cohorts, col = cohort_colors[present_cohorts], pch = 16, bty = "n", cex = 0.8)
    legend("bottomright", legend = c("SSL train", "SSL validation", "ICI projection"), pch = c(16, 1, 17), bty = "n", 
        cex = 0.8)
}
save_png <- function(filename, draw_function) {
    png(file.path(PCA_DIR, filename), width = 1800, height = 1400, res = 180)
    on.exit(dev.off(), add = TRUE)
    draw_function()
}
save_png("PCA_all_cohorts.png", function() draw_pca(all_scores, "PCA fitted on SSL training samples"))
save_png("PCA_SSL_only.png", function() draw_pca(all_scores[all_scores$data_role != "ICI_projection", , drop = FALSE], 
    "SSL training and validation samples"))
qc <- rbind(prepared$ssl_sample_qc[setdiff(names(prepared$ssl_sample_qc), "split")], prepared$ici_sample_qc)
qc$cohort <- factor(qc$cohort, levels = cohort_order)
save_png("sample_QC_overview.png", function() {
    par(mfrow = c(2, 1), mar = c(5, 5, 3, 1))
    boxplot(n_detected_observed_genes ~ cohort, data = qc, col = cohort_colors[cohort_order], las = 2, xlab = "", 
        ylab = "Detected observed genes", main = "Observed protein-coding gene detection", cex.axis = 0.8)
    boxplot(mean_abs_z_among_observed ~ cohort, data = qc, col = cohort_colors[cohort_order], las = 2, xlab = "", 
        ylab = "Mean absolute standardized expression", main = "Deviation from SSL training reference", cex.axis = 0.8)
})
writeLines(capture.output(sessionInfo()), file.path(PCA_DIR, "sessionInfo.txt"))
cat("\n\u56FE\u7247\u5DF2\u4FDD\u5B58\uFF1A\n")
cat("PCA_all_cohorts.png\n")
cat("PCA_SSL_only.png\n")
cat("sample_QC_overview.png\n")
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", PCA_DIR, "\n")

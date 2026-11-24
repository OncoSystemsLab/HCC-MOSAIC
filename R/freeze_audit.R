INPUT <- out_dir("06_inputs")
SSL_ROOT <- out_dir("07_SSL")
CONTEXT <- out_dir("07_context")
LOCO <- out_dir("08_LOCO")
BA37_FINAL_DIR <- out_dir("09_final")
pred <- read.delim(file.path(LOCO, "patient_predictions.tsv"), stringsAsFactors = FALSE)
pred <- pred[pred$panel_id == "P009_BA37" & pred$mode == "Early", ]
stopifnot(nrow(pred) == 67L * 4L, !anyDuplicated(pred[c("patient_uid", "seed")]))
ensemble <- pred[as.character(pred$seed) == "ensemble", ]
stopifnot(nrow(ensemble) == 67L, !anyDuplicated(ensemble$patient_uid))
seed_matrix <- sapply(c("24", "42", "64"), function(seed) {
    p <- pred[as.character(pred$seed) == seed, ]
    stopifnot(nrow(p) == 67L, setequal(p$patient_uid, ensemble$patient_uid))
    p <- p[match(ensemble$patient_uid, p$patient_uid), ]
    stopifnot(identical(p$response, ensemble$response), identical(p$cohort, ensemble$cohort))
    p$prediction
})
stopifnot(max(abs(rowMeans(seed_matrix) - ensemble$prediction)) < 1e-06)
auc <- function(y, p) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    stopifnot(n1 > 0, n0 > 0, all(is.finite(p)))
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}
audit <- do.call(rbind, lapply(c("24", "42", "64", "ensemble"), function(seed) {
    p <- pred[as.character(pred$seed) == seed, ]
    do.call(rbind, lapply(c("GSE181946", "GSE202069", "PRJEB34724"), function(cc) {
        x <- p[p$cohort == cc, ]
        data.frame(seed = seed, cohort = cc, AUC = auc(x$response, x$prediction))
    }))
}))
summary <- do.call(rbind, lapply(unique(audit$seed), function(seed) {
    a <- audit$AUC[audit$seed == seed]
    data.frame(seed = seed, macro_AUC = mean(a), minimum_AUC = min(a), target_met = mean(a) > 0.8 && min(a) > 0.7)
}))
outer <- read.delim(file.path(LOCO, "outer_cohort_metrics.tsv"))
outer <- outer[outer$panel_id == "P009_BA37" & outer$mode == "Early", ]
ens_audit <- audit[audit$seed == "ensemble", ]
stopifnot(nrow(outer) == 3L, max(abs(outer$AUC[match(ens_audit$cohort, outer$cohort)] - ens_audit$AUC)) < 1e-08, 
    all(is.finite(outer$mean_abs_SSL_ablation_logit_change)), all(outer$mean_abs_SSL_ablation_logit_change > 1e-08))
checkpoint_names <- unlist(lapply(c("GSE181946", "GSE202069", "PRJEB34724"), function(cc) paste0("P009_BA37_Early_", 
    cc, "_seed", c(24, 42, 64), ".pt")), use.names = FALSE)
stopifnot(all(file.exists(file.path(LOCO, checkpoint_names))))
cap_counts <- table(outer$selected_cap)
cap <- as.numeric(names(cap_counts)[which.max(cap_counts)])
stopifnot(cap %in% c(0.25, 0.75), max(cap_counts) >= 2L)
dir.create(BA37_FINAL_DIR, recursive = TRUE)
freeze_dir <- file.path(BA37_FINAL_DIR, "00_frozen_development")
dir.create(freeze_dir)
sources <- c(inputs = INPUT, SSL_training = SSL_ROOT, context = CONTEXT, LOCO = LOCO)
copy_rows <- list()
for (name in names(sources)) {
    src <- normalizePath(sources[[name]], winslash = "/", mustWork = TRUE)
    files <- list.files(src, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    files <- files[!file.info(files)$isdir]
    rel <- substring(normalizePath(files, winslash = "/", mustWork = TRUE), nchar(src) + 2L)
    targets <- file.path(freeze_dir, name, rel)
    for (d in unique(dirname(targets))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
    hashes <- unname(tools::md5sum(files))
    if (anyNA(hashes) || !all(file.copy(files, targets, overwrite = FALSE)) || any(hashes != unname(tools::md5sum(targets)))) 
        stop("\u51BB\u7ED3\u6587\u4EF6\u590D\u5236\u5931\u8D25\uFF1A", name)
    copy_rows[[name]] <- data.frame(filename = paste0(name, "/", rel), md5 = hashes)
}
write.table(do.call(rbind, copy_rows), file.path(freeze_dir, "copied_files_manifest.tsv"), sep = "\t", quote = FALSE, 
    row.names = FALSE)
write.table(audit, file.path(freeze_dir, "LOCO_auc_audit.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summary, file.path(freeze_dir, "seed_stability.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(outer, file.path(freeze_dir, "SSL_dependence_and_outer_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
FULL <- file.path(BA37_FINAL_DIR, "01_all67_fixed_model")
SCORES <- file.path(BA37_FINAL_DIR, "02_scores616")
ATTR <- file.path(BA37_FINAL_DIR, "03_attribution")
for (d in c(FULL, SCORES, ATTR)) dir.create(d)

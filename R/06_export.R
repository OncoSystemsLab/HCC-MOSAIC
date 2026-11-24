PROJECT_DIR <- CFG$root
AUDIT37_FILE <- artifact("stage09_token37_audit_bundle.rds")
PREP_FILE <- artifact("stage04_SSL_prepared_candidate_bundle.rds")
stopifnot(file.exists(AUDIT37_FILE), file.exists(PREP_FILE))
audit37 <- readRDS(AUDIT37_FILE)
prepared <- readRDS(PREP_FILE)
sources <- audit37$source_manifest
current_md5 <- unname(tools::md5sum(sources$path))
if (anyNA(current_md5) || any(current_md5 != sources$md5)) {
    stop("37-token\u5BA1\u8BA1\u7684\u6E90\u6587\u4EF6\u5DF2\u6539\u53D8\uFF0C\u8BF7\u5148\u91CD\u65B0\u5BA1\u8BA1\u3002")
}
get_source_path <- function(name) {
    hit <- which(sources$source == name)
    if (length(hit) != 1L) 
        stop("\u6E90\u6587\u4EF6\u8BB0\u5F55\u4E0D\u552F\u4E00\uFF1A", name)
    sources$path[hit]
}
intersection_bundle <- readRDS(get_source_path("pretraining_intersection"))
aligned_bundle <- readRDS(get_source_path("ici_aligned"))
background_genes <- audit37$fixed_background_genes
union_indices <- audit37$member_union_indices
union_genes <- background_genes[union_indices]
tokens37 <- audit37$tokens
ssl_log <- intersection_bundle$ssl_expression_log2tpm
ssl_meta <- prepared$ssl_patient_metadata
stopifnot(length(background_genes) == 17101L, length(union_genes) == 1753L, nrow(tokens37) == 37L, identical(rownames(ssl_log), 
    background_genes), identical(background_genes, prepared$fixed_genes), identical(colnames(ssl_log), ssl_meta$patient_uid))
train_indices <- which(ssl_meta$split == "train")
validation_indices <- which(ssl_meta$split == "validation")
stopifnot(length(train_indices) == 554L, length(validation_indices) == 62L, identical(train_indices, prepared$train_index))
DUAL_INPUT_DIR <- out_dir("06_inputs")
if (dir.exists(DUAL_INPUT_DIR) || !dir.create(DUAL_INPUT_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_dual_tsv <- function(x, filename) {
    write.table(x, file.path(DUAL_INPUT_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", 
        fileEncoding = "UTF-8")
}
write_float32 <- function(x, filename) {
    if (any(!is.finite(x))) {
        stop("\u4E8C\u8FDB\u5236\u5BFC\u51FA\u5305\u542B\u65E0\u6548\u503C\uFF1A", filename)
    }
    con <- file(file.path(DUAL_INPUT_DIR, filename), "wb")
    on.exit(close(con), add = TRUE)
    writeBin(as.vector(x), con, size = 4L, endian = "little")
}
write_mask <- function(x, filename) {
    stopifnot(is.logical(x), !anyNA(x))
    con <- file(file.path(DUAL_INPUT_DIR, filename), "wb")
    on.exit(close(con), add = TRUE)
    writeBin(as.raw(as.vector(x)), con)
}
train_background_log <- ssl_log[, train_indices, drop = FALSE]
background_median <- apply(train_background_log, 1L, median)
train_union_log <- train_background_log[union_indices, , drop = FALSE]
union_mean <- rowMeans(train_union_log)
union_sd_raw <- apply(train_union_log, 1L, sd)
SSL_SD_FLOOR <- 0.25
SSL_CLIP <- 6
union_scale <- pmax(union_sd_raw, SSL_SD_FLOOR)
scaler_table <- data.frame(ensembl_gene_id = union_genes, train_mean_log2tpm = union_mean, train_sd_log2tpm = union_sd_raw, 
    scale_used = union_scale, stringsAsFactors = FALSE)
write_dual_tsv(scaler_table, "SSL_union_training_scaler.tsv")
write_dual_tsv(data.frame(ensembl_gene_id = background_genes, train_median_log2tpm = background_median), "background_training_imputation.tsv")
union_gene_table <- data.frame(gene_index0 = seq_along(union_genes) - 1L, ensembl_gene_id = union_genes, gene_symbol = intersection_bundle$gene_manifest$gene_symbol[union_indices], 
    stringsAsFactors = FALSE)
write_dual_tsv(union_gene_table, "union_genes.tsv")
write_dual_tsv(intersection_bundle$gene_manifest, "fixed_background_genes.tsv")
hierarchy37 <- tokens37[c("token_name", "concept_id", "concept_name", "module_id", "module_name")]
hierarchy37$token_index0 <- seq_len(nrow(hierarchy37)) - 1L
write_dual_tsv(hierarchy37, "token_hierarchy.tsv")
membership_rows <- lapply(seq_len(nrow(tokens37)), function(i) {
    token <- tokens37$token_name[i]
    background_index <- audit37$token_gene_indices[[token]]
    ids <- background_genes[background_index]
    union_index <- match(ids, union_genes)
    if (length(ids) == 0L || anyNA(union_index)) {
        stop(token, "\uFF1A\u6210\u5458\u7D22\u5F15\u65E0\u6548\u3002")
    }
    data.frame(token_index0 = i - 1L, token_name = token, gene_index0 = union_index - 1L, ensembl_gene_id = ids, 
        stringsAsFactors = FALSE)
})
membership37 <- do.call(rbind, membership_rows)
write_dual_tsv(membership37, "token_gene_membership.tsv")
dataset_order <- c("ssl", "GSE181946", "GSE202069", "PRJEB34724")
shape_rows <- list()
rank_qc_rows <- list()
imputation_rows <- list()
value_check_rows <- list()
for (dataset in dataset_order) {
    cat("\n\u6B63\u5728\u751F\u6210\uFF1A", dataset, "\n")
    if (dataset == "ssl") {
        log_expression <- ssl_log
        observed_mask <- matrix(TRUE, nrow = nrow(log_expression), ncol = ncol(log_expression), dimnames = dimnames(log_expression))
        metadata <- ssl_meta
        sample_uids <- metadata$patient_uid
    }
    else {
        obj <- aligned_bundle$ici_aligned[[dataset]]
        log_expression <- obj$expression_log2tpm
        observed_mask <- obj$observed_mask
        metadata <- obj$metadata
        sample_uids <- metadata$patient_uid
    }
    stopifnot(identical(rownames(log_expression), background_genes), identical(dimnames(log_expression), dimnames(observed_mask)), 
        length(sample_uids) == ncol(log_expression), !anyDuplicated(sample_uids), all(is.finite(log_expression[observed_mask])), 
        all(is.na(log_expression[!observed_mask])))
    filled_log <- log_expression
    for (j in seq_len(ncol(filled_log))) {
        missing <- !observed_mask[, j]
        filled_log[missing, j] <- background_median[missing]
    }
    percentile_rank <- vapply(seq_len(ncol(filled_log)), function(j) {
        rank(filled_log[, j], ties.method = "average")/length(background_genes)
    }, numeric(length(background_genes)))
    dim(percentile_rank) <- c(length(background_genes), ncol(filled_log))
    token_rank <- matrix(NA_real_, nrow = nrow(tokens37), ncol = ncol(filled_log), dimnames = list(tokens37$token_name, 
        sample_uids))
    for (i in seq_len(nrow(tokens37))) {
        token <- tokens37$token_name[i]
        indices <- audit37$token_gene_indices[[token]]
        token_rank[i, ] <- colMeans(percentile_rank[indices, , drop = FALSE])
        observed_fraction <- colMeans(observed_mask[indices, , drop = FALSE])
        observed_only_score <- vapply(seq_len(ncol(filled_log)), function(j) {
            keep <- observed_mask[indices, j]
            if (!any(keep)) 
                return(NA_real_)
            mean(percentile_rank[indices[keep], j])
        }, numeric(1))
        key <- paste(dataset, token, sep = "::")
        imputation_rows[[key]] <- data.frame(dataset = dataset, token_name = token, n_fixed_members = length(indices), 
            min_observed_fraction = min(observed_fraction), max_missing_members = max(colSums(!observed_mask[indices, 
                , drop = FALSE])), max_abs_rank_difference_observed_only = if (all(is.finite(observed_only_score))) {
                max(abs(token_rank[i, ] - observed_only_score))
            }
            else {
                NA_real_
            }, stringsAsFactors = FALSE)
    }
    if (any(!is.finite(token_rank)) || any(token_rank < 0 | token_rank > 1)) {
        stop(dataset, "\uFF1Atoken\u767E\u5206\u4F4D\u79E9\u65E0\u6548\u3002")
    }
    union_log <- filled_log[union_indices, , drop = FALSE]
    union_z <- sweep(union_log, 1L, union_mean, "-")
    union_z <- sweep(union_z, 1L, union_scale, "/")
    clipped_fraction <- mean(abs(union_z) > SSL_CLIP)
    union_z <- pmax(union_z, -SSL_CLIP)
    union_z <- pmin(union_z, SSL_CLIP)
    stopifnot(is.matrix(union_z), identical(dim(union_z), dim(union_log)), identical(dimnames(union_z), dimnames(union_log)), 
        all(is.finite(union_z)))
    union_mask <- observed_mask[union_indices, , drop = FALSE]
    write_float32(union_z, paste0(dataset, "_gene_z.float32.bin"))
    write_float32(token_rank, paste0(dataset, "_token_rank.float32.bin"))
    write_mask(union_mask, paste0(dataset, "_gene_observed.uint8.bin"))
    write_dual_tsv(metadata, paste0(dataset, "_metadata.tsv"))
    for (kind in c("gene_z", "token_rank")) {
        mat <- if (kind == "gene_z") 
            union_z
        else token_rank
        checks <- expand.grid(feature_index1 = unique(c(1L, nrow(mat))), sample_index1 = unique(c(1L, ncol(mat))))
        checks$value <- mat[cbind(checks$feature_index1, checks$sample_index1)]
        checks$dataset <- dataset
        checks$kind <- kind
        value_check_rows[[paste(dataset, kind)]] <- checks
    }
    shape_rows[[dataset]] <- data.frame(dataset = dataset, n_samples = ncol(filled_log), n_background_genes = length(background_genes), 
        n_SSL_genes = length(union_genes), n_tokens = nrow(tokens37), stringsAsFactors = FALSE)
    rank_qc_rows[[dataset]] <- data.frame(dataset = dataset, min_token_rank = min(token_rank), max_token_rank = max(token_rank), 
        n_tokens_constant_within_dataset = sum(apply(token_rank, 1L, sd) < 1e-10), SSL_expression_fraction_clipped = clipped_fraction, 
        max_missing_union_genes_per_sample = max(colSums(!union_mask)), stringsAsFactors = FALSE)
}
shapes37 <- do.call(rbind, shape_rows)
rank_qc37 <- do.call(rbind, rank_qc_rows)
imputation_qc37 <- do.call(rbind, imputation_rows)
value_checks37 <- do.call(rbind, value_check_rows)
write_dual_tsv(shapes37, "shapes.tsv")
write_dual_tsv(rank_qc37, "rank_and_SSL_numeric_QC.tsv")
write_dual_tsv(imputation_qc37, "token_imputation_sensitivity.tsv")
write_dual_tsv(value_checks37, "R_Python_value_checks.tsv")
settings37 <- list(background_gene_count = length(background_genes), SSL_gene_count = length(union_genes), n_tokens = 37L, 
    n_concepts = 22L, n_modules = 8L, rank_definition = "mean member percentile rank within fixed background", 
    ties_method = "average", missing_fill = "554-case SSL training gene median log2TPM", SSL_scaler_fit = "554-case SSL training subset only", 
    SSL_sd_floor = SSL_SD_FLOOR, SSL_clip = SSL_CLIP, response_labels_used_for_transform = FALSE, panel_selection_performed = FALSE, 
    audit37_file = AUDIT37_FILE, audit37_md5 = unname(tools::md5sum(AUDIT37_FILE)), prepared_file = PREP_FILE, 
    prepared_md5 = unname(tools::md5sum(PREP_FILE)))
saveRDS(list(settings = settings37, shapes = shapes37, hierarchy = hierarchy37, membership = membership37, scaler = scaler_table, 
    background_gene_medians = background_median), file.path(DUAL_INPUT_DIR, "stage09_token37_input_bundle.rds"))
exported_paths <- list.files(DUAL_INPUT_DIR, full.names = TRUE)
write_dual_tsv(data.frame(filename = basename(exported_paths), size_bytes = file.info(exported_paths)$size, md5 = unname(tools::md5sum(exported_paths))), 
    "export_manifest.tsv")
writeLines(capture.output(sessionInfo()), file.path(DUAL_INPUT_DIR, "R_sessionInfo.txt"))
cat("\n===== \u53CC\u8DEF\u8F93\u5165\u7EF4\u5EA6 =====\n")
print(shapes37, row.names = FALSE)
cat("\n===== Rank\u53CASSL\u6570\u503C\u8D28\u63A7 =====\n")
print(rank_qc37, row.names = FALSE)
cat("\n===== \u586B\u8865\u654F\u611F\u6027\u8F83\u5927\u7684token\uFF1A\u524D15\u6761 =====\n")
affected <- imputation_qc37[imputation_qc37$max_missing_members > 0L, , drop = FALSE]
affected <- affected[order(-affected$max_abs_rank_difference_observed_only, affected$dataset, affected$token_name, 
    na.last = TRUE), , drop = FALSE]
print(head(affected, 15L), row.names = FALSE)
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", DUAL_INPUT_DIR, "\n")

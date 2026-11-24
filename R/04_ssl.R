# Create the SSL split, preprocess expression, and audit samples.
PROJECT_DIR <- CFG$root
RAW_DIR <- file.path(CFG$root, "data/raw_data")
SSL_INPUT_DIR <- out_dir("04_ssl_raw")
if (dir.exists(SSL_INPUT_DIR) || !dir.create(SSL_INPUT_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
source_files <- data.frame(cohort = rep(c("TCGA", "ICGC", "CPTAC"), each = 2L), role = rep(c("expression_TPM", 
    "mutation_matrix"), 3L), path = c(file.path(RAW_DIR, "TCGA_LIHC_final/TCGA_LIHC_expression_TPM.tsv.gz"), file.path(RAW_DIR, 
    "TCGA_LIHC_final/TCGA_LIHC_mutation_matrix.tsv.gz"), file.path(RAW_DIR, "ICGC_LICA_FR_FINAL/ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz"), 
    file.path(RAW_DIR, "ICGC_LICA_FR_FINAL/ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"), paste0("data/CPTAC/", 
        "CPTAC_validation_TPM_v19_159/", "CPTAC_159_TPM_Ensembl_v19.rds"), paste0("data/CPTAC/", 
        "CPTAC_validation_TPM_v19_159/", "CPTAC_159_mutation_matrix.rds")), stringsAsFactors = FALSE)
if (any(!file.exists(source_files$path))) {
    stop("\u4EE5\u4E0B\u8F93\u5165\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A\n", paste(source_files$path[!file.exists(source_files$path)], collapse = "\n"))
}
write_ssl_tsv <- function(x, filename) {
    write.table(x, file.path(SSL_INPUT_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", 
        fileEncoding = "UTF-8")
}
validate_matrix <- function(x, label) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop(label, "\uFF1A\u4E0D\u662F\u6570\u503C\u77E9\u9635\uFF0C\u9700\u8981\u68C0\u67E5\u5BF9\u8C61\u7ED3\u6784\u3002")
    }
    if (nrow(x) == 0L || ncol(x) == 0L) {
        stop(label, "\uFF1A\u77E9\u9635\u4E3A\u7A7A\u3002")
    }
    for (ids in list(rownames(x), colnames(x))) {
        if (is.null(ids) || anyNA(ids) || any(!nzchar(ids)) || any(ids != trimws(ids)) || anyDuplicated(ids)) {
            stop(label, "\uFF1A\u884C\u540D\u6216\u5217\u540D\u7F3A\u5931\u3001\u91CD\u590D\u6216\u5305\u542B\u4E24\u7AEF\u7A7A\u683C\u3002")
        }
    }
    if (any(!is.finite(x)) || any(x < 0)) {
        stop(label, "\uFF1A\u5B58\u5728\u7F3A\u5931\u3001\u65E0\u7A77\u6216\u8D1F\u503C\u3002")
    }
    x
}
read_source_matrix <- function(path, label) {
    if (grepl("\\.rds$", path, ignore.case = TRUE)) {
        obj <- readRDS(path)
        if (is.data.frame(obj)) {
            if (!all(vapply(obj, is.numeric, logical(1)))) {
                stop(label, "\uFF1ARDS \u6570\u636E\u6846\u542B\u975E\u6570\u503C\u5217\uFF0C\u9700\u8981\u5355\u72EC\u68C0\u67E5\u3002")
            }
            obj <- as.matrix(obj)
        }
        return(validate_matrix(obj, label))
    }
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    dat <- read.delim(con, header = TRUE, sep = "\t", quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, 
        colClasses = "character", fill = FALSE)
    if (sum(names(dat) == "Gene") != 1L) {
        stop(label, "\uFF1A\u672A\u627E\u5230\u552F\u4E00 Gene \u5217\u3002")
    }
    gene_ids <- dat$Gene
    values <- as.matrix(dat[names(dat) != "Gene"])
    numeric_values <- suppressWarnings(as.numeric(values))
    x <- matrix(numeric_values, nrow = nrow(dat), ncol = ncol(values), dimnames = list(gene_ids, colnames(values)))
    validate_matrix(x, label)
}
unlabeled_raw <- list()
matrix_summaries <- list()
pairing_reports <- list()
rna_qc_reports <- list()
for (cohort in c("TCGA", "ICGC", "CPTAC")) {
    cat("\n\u6B63\u5728\u8BFB\u53D6\uFF1A", cohort, "\n")
    expr_path <- source_files$path[source_files$cohort == cohort & source_files$role == "expression_TPM"]
    mut_path <- source_files$path[source_files$cohort == cohort & source_files$role == "mutation_matrix"]
    expr <- read_source_matrix(expr_path, paste(cohort, "TPM"))
    mut <- read_source_matrix(mut_path, paste(cohort, "mutation"))
    if (any(abs(mut - round(mut)) > 1e-08)) {
        stop(cohort, "\uFF1A\u7A81\u53D8\u77E9\u9635\u5B58\u5728\u975E\u6574\u6570\u503C\uFF0C\u9700\u8981\u6838\u5B9E\u5B9A\u4E49\u3002")
    }
    mutation_is_binary <- all(mut %in% c(0, 1))
    patient_ids <- union(colnames(expr), colnames(mut))
    pairing <- data.frame(cohort = cohort, patient_id = patient_ids, has_RNA_column = patient_ids %in% colnames(expr), 
        has_WES_column = patient_ids %in% colnames(mut), stringsAsFactors = FALSE)
    rna_qc <- data.frame(cohort = cohort, patient_id = colnames(expr), expression_sum = colSums(expr), n_detected_genes = colSums(expr > 
        0), zero_fraction = colMeans(expr == 0), all_zero = colSums(expr > 0) == 0L, stringsAsFactors = FALSE)
    summary_row <- data.frame(cohort = cohort, n_RNA_genes = nrow(expr), n_RNA_samples = ncol(expr), n_mutation_genes = nrow(mut), 
        n_WES_columns = ncol(mut), n_RNA_WES_matched = length(intersect(colnames(expr), colnames(mut))), n_RNA_without_WES_column = length(setdiff(colnames(expr), 
            colnames(mut))), n_WES_without_RNA_column = length(setdiff(colnames(mut), colnames(expr))), mutation_is_binary = mutation_is_binary, 
        min_RNA_sum = min(colSums(expr)), max_RNA_sum = max(colSums(expr)), stringsAsFactors = FALSE)
    unlabeled_raw[[cohort]] <- list(expression_tpm = expr, mutation_matrix = mut, pairing = pairing)
    matrix_summaries[[cohort]] <- summary_row
    pairing_reports[[cohort]] <- pairing
    rna_qc_reports[[cohort]] <- rna_qc
    write_ssl_tsv(pairing, paste0(cohort, "_RNA_WES_column_matching.tsv"))
    write_ssl_tsv(rna_qc, paste0(cohort, "_RNA_numeric_qc.tsv"))
    cat("\u5B8C\u6210\uFF1A", cohort, " | RNA\uFF1A", ncol(expr), " | WES \u5217\uFF1A", ncol(mut), " | \u5339\u914D\uFF1A", summary_row$n_RNA_WES_matched, 
        "\n", sep = "")
}
CPTAC_MAPPING_FILE <- paste0("data/CPTAC/", "CPTAC_source_review/sample_identity_mapping.tsv")
if (!file.exists(CPTAC_MAPPING_FILE)) {
    stop("CPTAC \u6837\u672C\u6620\u5C04\u6587\u4EF6\u4E0D\u5B58\u5728\u3002")
}
cptac_mapping <- read.delim(CPTAC_MAPPING_FILE, sep = "\t", quote = "", comment.char = "", colClasses = "character", 
    check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("Count_column", "Tissue_type", "Patient_tumor_ID")
if (!all(required_columns %in% names(cptac_mapping))) {
    stop("CPTAC \u6620\u5C04\u8868\u7F3A\u5C11\u5FC5\u8981\u5B57\u6BB5\u3002")
}
tumor_mapping <- cptac_mapping[!is.na(cptac_mapping$Tissue_type) & cptac_mapping$Tissue_type == "Tumor", , drop = FALSE]
if (nrow(tumor_mapping) == 0L) {
    stop("CPTAC \u6620\u5C04\u8868\u6CA1\u6709\u8BC6\u522B\u5230 Tumor \u8BB0\u5F55\u3002")
}
if (anyNA(tumor_mapping$Patient_tumor_ID) || any(!nzchar(tumor_mapping$Patient_tumor_ID)) || anyDuplicated(tumor_mapping$Patient_tumor_ID)) {
    stop("CPTAC \u80BF\u7624\u6620\u5C04 ID \u7F3A\u5931\u6216\u91CD\u590D\uFF0C\u9700\u8981\u68C0\u67E5\u3002")
}
cptac_rna_ids <- colnames(unlabeled_raw$CPTAC$expression_tpm)
mapping_index <- match(cptac_rna_ids, tumor_mapping$Patient_tumor_ID)
cptac_identity_check <- data.frame(RNA_patient_id = cptac_rna_ids, matched_tumor_mapping = !is.na(mapping_index), 
    source_count_column = tumor_mapping$Count_column[mapping_index], stringsAsFactors = FALSE)
write_ssl_tsv(cptac_identity_check, "CPTAC_RNA_tumor_identity_check.tsv")
unlabeled_summary <- do.call(rbind, matrix_summaries)
unlabeled_rna_qc <- do.call(rbind, rna_qc_reports)
rownames(unlabeled_summary) <- NULL
rownames(unlabeled_rna_qc) <- NULL
source_files$md5 <- unname(tools::md5sum(source_files$path))
if (anyNA(source_files$md5)) {
    stop("\u90E8\u5206\u8F93\u5165\u6587\u4EF6\u65E0\u6CD5\u8BA1\u7B97 MD5\u3002")
}
write_ssl_tsv(source_files, "unlabeled_input_manifest.tsv")
write_ssl_tsv(unlabeled_summary, "unlabeled_cohort_summary.tsv")
saveRDS(list(unlabeled_raw = unlabeled_raw, unlabeled_summary = unlabeled_summary, unlabeled_rna_qc = unlabeled_rna_qc, 
    cptac_identity_check = cptac_identity_check, input_files = source_files, cptac_mapping_file = CPTAC_MAPPING_FILE, 
    cptac_mapping_md5 = unname(tools::md5sum(CPTAC_MAPPING_FILE)), response_labels_available = FALSE, WES_quality_verified = FALSE, 
    final_pretraining_samples_locked = FALSE), artifact("stage03_unlabeled_raw_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(SSL_INPUT_DIR, "sessionInfo.txt"))
cat("\n===== \u65E0\u5E94\u7B54\u6807\u7B7E\u961F\u5217\u6C47\u603B =====\n")
print(unlabeled_summary, row.names = FALSE)
cat("\n===== CPTAC \u80BF\u7624\u6620\u5C04\u5339\u914D\u60C5\u51B5 =====\n")
print(table(cptac_identity_check$matched_tumor_mapping))
cat("\n===== \u57FA\u56E0 ID \u793A\u4F8B =====\n")
for (cohort in names(unlabeled_raw)) {
    cat(cohort, "\uFF1A", paste(head(rownames(unlabeled_raw[[cohort]]$expression_tpm)), collapse = ", "), "\n")
}
cat("\nRNA \u6837\u672C\u5217\u5408\u8BA1\uFF1A", sum(unlabeled_summary$n_RNA_samples), "\n")
cat("\u672C\u6B21\u8F93\u51FA\u76EE\u5F55\uFF1A", SSL_INPUT_DIR, "\n")
if (any(!cptac_identity_check$matched_tumor_mapping)) {
    warning("\u90E8\u5206 CPTAC RNA ID \u672A\u5339\u914D\u80BF\u7624\u6620\u5C04\uFF0C\u9700\u5148\u6838\u67E5\u3002")
}
if (any(unlabeled_rna_qc$all_zero)) {
    warning("\u53D1\u73B0\u5168\u96F6 RNA \u6837\u672C\uFF0C\u9700\u5148\u6838\u67E5\u3002")
}
if (any(!unlabeled_summary$mutation_is_binary)) {
    warning("\u90E8\u5206\u7A81\u53D8\u77E9\u9635\u4E0D\u662F\u4E8C\u503C\u77E9\u9635\uFF0C\u9700\u6838\u5B9E\u5176\u6570\u503C\u5B9A\u4E49\u3002")
}
PROJECT_DIR <- CFG$root
UNLABELED_FILE <- artifact("stage03_unlabeled_raw_bundle.rds")
REFERENCE_FILE <- artifact("stage01_reference_bundle.rds")
if (!all(file.exists(c(UNLABELED_FILE, REFERENCE_FILE)))) {
    stop("\u65E0\u6807\u7B7E\u961F\u5217\u7ED3\u679C\u6216\u53C2\u8003 RDS \u4E0D\u5B58\u5728\u3002")
}
unlabeled_bundle <- readRDS(UNLABELED_FILE)
reference_bundle <- readRDS(REFERENCE_FILE)
unlabeled_raw <- unlabeled_bundle$unlabeled_raw
gene_ref <- reference_bundle$gene_reference
stopifnot(sum(vapply(unlabeled_raw, function(obj) ncol(obj$expression_tpm), integer(1))) == 616L, !anyDuplicated(gene_ref$ensembl_gene_id))
source_files <- unlabeled_bundle$input_files
current_md5 <- unname(tools::md5sum(source_files$path))
if (anyNA(current_md5) || any(current_md5 != source_files$md5)) {
    stop("\u6E90\u6587\u4EF6\u7F3A\u5931\u6216\u5DF2\u6539\u53D8\uFF0C\u8BF7\u91CD\u65B0\u5B8C\u6210\u7B2C 25\u201326 \u6BB5\u3002")
}
SSL_MAP_DIR <- out_dir("04_ssl_mapping")
if (dir.exists(SSL_MAP_DIR) || !dir.create(SSL_MAP_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_ssl_map <- function(x, filename) {
    write.table(x, file.path(SSL_MAP_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
rows_by_symbol <- split(seq_len(nrow(gene_ref)), gene_ref$gene_symbol)
ssl_gene_mappings <- list()
for (cohort in names(unlabeled_raw)) {
    x <- unlabeled_raw[[cohort]]$expression_tpm
    source_ids <- rownames(x)
    if (cohort %in% c("TCGA", "ICGC")) {
        matches <- lapply(source_ids, function(s) rows_by_symbol[[s]])
        n_matches <- lengths(matches)
        selected_index <- vapply(matches, function(rows) {
            if (length(rows) == 1L) {
                as.integer(rows)
            }
            else {
                NA_integer_
            }
        }, integer(1))
        mapping_status <- ifelse(n_matches == 0L, "symbol_unmatched", ifelse(n_matches == 1L, "symbol_exact_unique_candidate", 
            "symbol_ambiguous"))
        source_id_type <- "gene_symbol"
    }
    else {
        ensembl_ids <- sub("\\.[0-9]+$", "", source_ids)
        selected_index <- match(ensembl_ids, gene_ref$ensembl_gene_id)
        n_matches <- as.integer(!is.na(selected_index))
        mapping_status <- ifelse(is.na(selected_index), "id_not_found_review", "id_exact")
        source_id_type <- "Ensembl_gene_ID"
    }
    mapping <- data.frame(original_gene_id = source_ids, source_id_type = source_id_type, n_reference_matches = n_matches, 
        mapping_status = mapping_status, target_ensembl_id = gene_ref$ensembl_gene_id[selected_index], target_symbol = gene_ref$gene_symbol[selected_index], 
        target_gene_type = gene_ref$gene_type[selected_index], stringsAsFactors = FALSE)
    mapping$protein_coding_candidate <- mapping$mapping_status %in% c("id_exact", "symbol_exact_unique_candidate") & 
        !is.na(mapping$target_gene_type) & mapping$target_gene_type == "protein_coding"
    candidate_index <- which(mapping$protein_coding_candidate)
    candidate_targets <- mapping$target_ensembl_id[candidate_index]
    mapping$target_has_multiple_source_rows <- FALSE
    mapping$target_has_multiple_source_rows[candidate_index] <- duplicated(candidate_targets) | duplicated(candidate_targets, 
        fromLast = TRUE)
    mapping$total_TPM_across_samples <- rowSums(x)
    ssl_gene_mappings[[cohort]] <- mapping
    write_ssl_map(mapping, paste0(cohort, "_GENCODE_v50_mapping.tsv"))
    needs_review <- mapping$mapping_status %in% c("symbol_unmatched", "symbol_ambiguous", "id_not_found_review") | 
        mapping$target_has_multiple_source_rows
    review <- mapping[needs_review, , drop = FALSE]
    review <- review[order(-review$total_TPM_across_samples, review$original_gene_id), , drop = FALSE]
    write_ssl_map(review, paste0(cohort, "_mapping_review.tsv"))
    cat("\n===== ", cohort, "\uFF1A\u5F85\u6838\u67E5\u8BB0\u5F55\uFF0C\u524D 15 \u884C =====\n", sep = "")
    print(head(review, 15L), row.names = FALSE)
}
status_rows <- list()
summary_rows <- list()
sample_rows <- list()
for (cohort in names(ssl_gene_mappings)) {
    mapping <- ssl_gene_mappings[[cohort]]
    x <- unlabeled_raw[[cohort]]$expression_tpm
    stopifnot(identical(rownames(x), mapping$original_gene_id))
    candidate <- mapping$protein_coding_candidate
    conservative <- candidate & !mapping$target_has_multiple_source_rows
    unresolved <- mapping$mapping_status %in% c("symbol_unmatched", "symbol_ambiguous", "id_not_found_review")
    status_table <- as.data.frame(table(mapping$mapping_status), stringsAsFactors = FALSE)
    names(status_table) <- c("mapping_status", "n_source_rows")
    status_table$cohort <- cohort
    status_rows[[cohort]] <- status_table[c("cohort", "mapping_status", "n_source_rows")]
    summary_rows[[cohort]] <- data.frame(cohort = cohort, n_samples = ncol(x), n_source_rows = nrow(x), n_pc_candidate_rows = sum(candidate), 
        n_unique_pc_targets = length(unique(mapping$target_ensembl_id[candidate])), n_pc_rows_with_shared_target = sum(candidate & 
            mapping$target_has_multiple_source_rows), n_conservative_pc_rows = sum(conservative), n_unresolved_rows = sum(unresolved), 
        stringsAsFactors = FALSE)
    sample_rows[[cohort]] <- data.frame(cohort = cohort, patient_id = colnames(x), conservative_pc_TPM_fraction = colSums(x[conservative, 
        , drop = FALSE])/colSums(x), unresolved_TPM_fraction = colSums(x[unresolved, , drop = FALSE])/colSums(x), 
        stringsAsFactors = FALSE)
}
ssl_mapping_status <- do.call(rbind, status_rows)
ssl_mapping_summary <- do.call(rbind, summary_rows)
ssl_mapping_sample_qc <- do.call(rbind, sample_rows)
rownames(ssl_mapping_status) <- NULL
rownames(ssl_mapping_summary) <- NULL
rownames(ssl_mapping_sample_qc) <- NULL
write_ssl_map(ssl_mapping_status, "unlabeled_mapping_status.tsv")
write_ssl_map(ssl_mapping_summary, "unlabeled_mapping_summary.tsv")
write_ssl_map(ssl_mapping_sample_qc, "unlabeled_mapping_sample_coverage.tsv")
saveRDS(list(unlabeled_file = UNLABELED_FILE, unlabeled_md5 = unname(tools::md5sum(UNLABELED_FILE)), reference_file = REFERENCE_FILE, 
    reference_md5 = unname(tools::md5sum(REFERENCE_FILE)), ssl_gene_mappings = ssl_gene_mappings, mapping_status = ssl_mapping_status, 
    mapping_summary = ssl_mapping_summary, sample_coverage = ssl_mapping_sample_qc, final_model_gene_list_locked = FALSE), 
    artifact("stage03_unlabeled_GENCODE_mapping_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(SSL_MAP_DIR, "sessionInfo.txt"))
cat("\n===== \u9884\u8BAD\u7EC3\u961F\u5217\u6620\u5C04\u72B6\u6001 =====\n")
print(ssl_mapping_status, row.names = FALSE)
cat("\n===== \u86CB\u767D\u7F16\u7801\u5019\u9009\u6C47\u603B =====\n")
print(ssl_mapping_summary, row.names = FALSE)
cat("\n===== \u6837\u672C TPM \u8986\u76D6\u8303\u56F4 =====\n")
for (cohort in names(ssl_gene_mappings)) {
    q <- ssl_mapping_sample_qc[ssl_mapping_sample_qc$cohort == cohort, , drop = FALSE]
    cat(cohort, " | \u4FDD\u5B88\u86CB\u767D\u7F16\u7801 TPM \u6BD4\u4F8B\uFF1A", paste(round(range(q$conservative_pc_TPM_fraction), 4), collapse = " ~ "), 
        " | \u5F85\u6838\u67E5 TPM \u6BD4\u4F8B\uFF1A", paste(round(range(q$unresolved_TPM_fraction), 4), collapse = " ~ "), "\n", 
        sep = "")
}
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", SSL_MAP_DIR, "\n")
PROJECT_DIR <- CFG$root
SSL_MAPPING_FILE <- artifact("stage03_unlabeled_GENCODE_mapping_bundle.rds")
ICI_PC_FILE <- artifact("stage02_ici_pc_candidate_bundle.rds")
TOKEN_AUDIT_FILE <- artifact("stage01_token_audit_bundle.rds")
required_files <- c(SSL_MAPPING_FILE, ICI_PC_FILE, TOKEN_AUDIT_FILE)
if (!all(file.exists(required_files))) {
    stop("\u4EE5\u4E0B\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}
ssl_mapping_bundle <- readRDS(SSL_MAPPING_FILE)
verify_md5 <- function(path, expected) {
    actual <- unname(tools::md5sum(path))
    if (length(actual) != 1L || is.na(actual) || actual != expected) {
        stop("\u4E0A\u6E38\u6587\u4EF6\u7F3A\u5931\u6216\u6307\u7EB9\u4E0D\u4E00\u81F4\uFF1A", path)
    }
}
verify_md5(ssl_mapping_bundle$unlabeled_file, ssl_mapping_bundle$unlabeled_md5)
verify_md5(ssl_mapping_bundle$reference_file, ssl_mapping_bundle$reference_md5)
unlabeled_bundle <- readRDS(ssl_mapping_bundle$unlabeled_file)
reference_bundle <- readRDS(ssl_mapping_bundle$reference_file)
unlabeled_raw <- unlabeled_bundle$unlabeled_raw
ssl_mappings <- ssl_mapping_bundle$ssl_gene_mappings
gene_ref <- reference_bundle$gene_reference
INTERSECTION_DIR <- out_dir("04_intersection")
if (dir.exists(INTERSECTION_DIR) || !dir.create(INTERSECTION_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_intersection_tsv <- function(x, filename) {
    write.table(x, file.path(INTERSECTION_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", 
        fileEncoding = "UTF-8")
}
cohort_order <- c("TCGA", "ICGC", "CPTAC")
eligible_mappings <- list()
cohort_gene_sets <- list()
for (cohort in cohort_order) {
    mapping <- ssl_mappings[[cohort]]
    keep <- mapping$protein_coding_candidate & !mapping$target_has_multiple_source_rows
    eligible <- mapping[keep, , drop = FALSE]
    if (anyNA(eligible$target_ensembl_id) || anyDuplicated(eligible$target_ensembl_id)) {
        stop(cohort, "\uFF1A\u5019\u9009\u76EE\u6807 ID \u7F3A\u5931\u6216\u91CD\u590D\u3002")
    }
    eligible_mappings[[cohort]] <- eligible
    cohort_gene_sets[[cohort]] <- eligible$target_ensembl_id
}
intersection_genes <- sort(Reduce(intersect, cohort_gene_sets), method = "radix")
if (length(intersection_genes) == 0L) {
    stop("\u4E09\u4E2A\u9884\u8BAD\u7EC3\u961F\u5217\u7684\u86CB\u767D\u7F16\u7801\u57FA\u56E0\u4EA4\u96C6\u4E3A\u7A7A\u3002")
}
reference_index <- match(intersection_genes, gene_ref$ensembl_gene_id)
if (anyNA(reference_index)) {
    stop("\u4EA4\u96C6\u57FA\u56E0\u672A\u80FD\u5728\u56FA\u5B9A\u53C2\u8003\u8868\u4E2D\u627E\u5230\u3002")
}
intersection_annotation <- gene_ref[reference_index, , drop = FALSE]
rownames(intersection_annotation) <- NULL
stopifnot(identical(intersection_genes, intersection_annotation$ensembl_gene_id), all(intersection_annotation$gene_type == 
    "protein_coding"))
gene_manifest <- data.frame(ensembl_gene_id = intersection_genes, gene_symbol = intersection_annotation$gene_symbol, 
    stringsAsFactors = FALSE)
ssl_tpm_list <- list()
patient_metadata_list <- list()
for (cohort in cohort_order) {
    source_matrix <- unlabeled_raw[[cohort]]$expression_tpm
    mapping <- eligible_mappings[[cohort]]
    map_index <- match(intersection_genes, mapping$target_ensembl_id)
    original_ids <- mapping$original_gene_id[map_index]
    expression_index <- match(original_ids, rownames(source_matrix))
    if (anyNA(expression_index)) {
        stop(cohort, "\uFF1A\u4EA4\u96C6\u57FA\u56E0\u65E0\u6CD5\u5B9A\u4F4D\u5230\u539F\u59CB\u8868\u8FBE\u884C\u3002")
    }
    x <- source_matrix[expression_index, , drop = FALSE]
    rownames(x) <- intersection_genes
    original_patient_ids <- colnames(x)
    patient_uid <- paste(cohort, original_patient_ids, sep = "::")
    colnames(x) <- patient_uid
    patient_metadata_list[[cohort]] <- data.frame(cohort = cohort, patient_id = original_patient_ids, patient_uid = patient_uid, 
        stringsAsFactors = FALSE)
    ssl_tpm_list[[cohort]] <- x
    gene_manifest[[paste0(cohort, "_source_gene_id")]] <- original_ids
}
ssl_tpm <- do.call(cbind, ssl_tpm_list)
ssl_patient_metadata <- do.call(rbind, patient_metadata_list)
rownames(ssl_patient_metadata) <- NULL
stopifnot(ncol(ssl_tpm) == 616L, !anyDuplicated(colnames(ssl_tpm)), identical(colnames(ssl_tpm), ssl_patient_metadata$patient_uid), 
    identical(rownames(ssl_tpm), intersection_genes), all(is.finite(ssl_tpm)), all(ssl_tpm >= 0))
ssl_log2tpm <- log2(ssl_tpm + 1)
intersection_summary <- data.frame(cohort = cohort_order, n_pc_genes_before_intersection = vapply(cohort_gene_sets, 
    length, integer(1)), n_intersection_genes = length(intersection_genes), stringsAsFactors = FALSE)
write_intersection_tsv(gene_manifest, "pretraining_intersection_gene_manifest.tsv")
write_intersection_tsv(ssl_patient_metadata, "pretraining_patient_metadata.tsv")
write_intersection_tsv(intersection_summary, "pretraining_intersection_summary.tsv")
cat("\n===== \u9884\u8BAD\u7EC3\u57FA\u56E0\u4EA4\u96C6 =====\n")
print(intersection_summary, row.names = FALSE)
cat("\u5019\u9009\u77E9\u9635\u7EF4\u5EA6\uFF1A", nrow(ssl_tpm), "\u57FA\u56E0 \u00D7", ncol(ssl_tpm), "\u6837\u672C\n")
token_bundle <- readRDS(TOKEN_AUDIT_FILE)
ici_bundle <- readRDS(ICI_PC_FILE)
tokens <- token_bundle$tokens
token_members <- unique(token_bundle$members_long[nzchar(token_bundle$members_long$gene_symbol), c("token_name", 
    "concept_id", "module_id", "gene_symbol")])
pc_reference_symbols <- unique(gene_ref$gene_symbol[gene_ref$gene_type == "protein_coding"])
intersection_ids_by_symbol <- split(intersection_annotation$ensembl_gene_id, intersection_annotation$gene_symbol)
member_target_lists <- lapply(token_members$gene_symbol, function(s) intersection_ids_by_symbol[[s]])
token_members$n_intersection_target_ids <- lengths(member_target_lists)
token_members$intersection_target_ids <- vapply(member_target_lists, function(ids) paste(sort(unique(ids)), collapse = ";"), 
    character(1))
token_members$in_strict_pc_reference <- token_members$gene_symbol %in% pc_reference_symbols
token_members$present_in_intersection <- token_members$n_intersection_target_ids > 0L
token_coverage_rows <- list()
for (i in seq_len(nrow(tokens))) {
    m <- token_members[token_members$token_name == tokens$token_name[i], , drop = FALSE]
    n_original <- nrow(m)
    n_pc_reference <- sum(m$in_strict_pc_reference)
    n_present <- sum(m$present_in_intersection)
    token_coverage_rows[[i]] <- data.frame(token_name = tokens$token_name[i], concept_id = tokens$concept_id[i], 
        module_id = tokens$module_id[i], n_original_members = n_original, n_strict_pc_reference_members = n_pc_reference, 
        n_members_in_intersection = n_present, n_members_with_multiple_target_ids = sum(m$n_intersection_target_ids > 
            1L), original_member_coverage = n_present/n_original, pc_reference_member_coverage = if (n_pc_reference > 
            0L) {
            n_present/n_pc_reference
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
}
token_intersection_coverage <- do.call(rbind, token_coverage_rows)
write_intersection_tsv(token_members, "token_members_in_pretraining_intersection.tsv")
write_intersection_tsv(token_intersection_coverage, "token_pretraining_intersection_coverage.tsv")
ici_input_coverage <- list()
ici_missing_genes <- list()
for (cohort in names(ici_bundle$pc_candidates)) {
    obj <- ici_bundle$pc_candidates[[cohort]]
    available_ids <- rownames(obj$expression_tpm)
    present <- intersection_genes %in% available_ids
    missing <- gene_manifest[!present, , drop = FALSE]
    ici_missing_genes[[cohort]] <- missing
    write_intersection_tsv(missing, paste0(cohort, "_missing_pretraining_intersection_genes.tsv"))
    ici_input_coverage[[cohort]] <- data.frame(cohort = cohort, n_samples = ncol(obj$expression_tpm), n_required_genes = length(intersection_genes), 
        n_available_genes = sum(present), n_missing_genes = sum(!present), gene_availability_fraction = mean(present), 
        stringsAsFactors = FALSE)
}
ici_input_coverage <- do.call(rbind, ici_input_coverage)
rownames(ici_input_coverage) <- NULL
write_intersection_tsv(ici_input_coverage, "ICI_fixed_input_gene_availability.tsv")
input_artifacts <- c(ssl_mapping = SSL_MAPPING_FILE, unlabeled_data = ssl_mapping_bundle$unlabeled_file, gencode_reference = ssl_mapping_bundle$reference_file, 
    ici_pc_candidates = ICI_PC_FILE, token_audit = TOKEN_AUDIT_FILE)
input_manifest <- data.frame(input_name = names(input_artifacts), path = unname(input_artifacts), md5 = unname(tools::md5sum(unname(input_artifacts))), 
    stringsAsFactors = FALSE)
if (anyNA(input_manifest$md5)) {
    stop("\u90E8\u5206\u8F93\u5165\u6587\u4EF6\u65E0\u6CD5\u8BA1\u7B97 MD5\u3002")
}
write_intersection_tsv(input_manifest, "intersection_input_manifest.tsv")
saveRDS(list(intersection_genes = intersection_genes, gene_manifest = gene_manifest, ssl_expression_tpm = ssl_tpm, 
    ssl_expression_log2tpm = ssl_log2tpm, ssl_patient_metadata = ssl_patient_metadata, intersection_summary = intersection_summary, 
    token_members = token_members, token_coverage = token_intersection_coverage, ici_input_coverage = ici_input_coverage, 
    ici_missing_genes = ici_missing_genes, input_manifest = input_manifest, policy = list(gene_selection_source = "TCGA_ICGC_CPTAC_annotation_intersection", 
        response_labels_used_for_selection = FALSE, expression_prevalence_filter_applied = FALSE, TPM_renormalized = FALSE, 
        ICI_missing_values_filled = FALSE, final_model_gene_list_locked = FALSE)), artifact("stage03_intersection_candidate_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(INTERSECTION_DIR, "sessionInfo.txt"))
cat("\n===== Token \u8986\u76D6\u8F83\u4F4E\u6216\u5B58\u5728\u591A\u76EE\u6807 ID \u7684\u8BB0\u5F55 =====\n")
flagged_tokens <- token_intersection_coverage[is.na(token_intersection_coverage$pc_reference_member_coverage) | 
    token_intersection_coverage$pc_reference_member_coverage < 0.8 | token_intersection_coverage$n_members_in_intersection <= 
    1L | token_intersection_coverage$n_members_with_multiple_target_ids > 0L, , drop = FALSE]
print(flagged_tokens, row.names = FALSE)
cat("\n\u5B8C\u5168\u65E0\u4EA4\u96C6\u6210\u5458\u7684 token \u6570\uFF1A", sum(token_intersection_coverage$n_members_in_intersection == 0L), 
    "\n")
cat("\n===== \u514D\u75AB\u6CBB\u7597\u961F\u5217\u5BF9\u56FA\u5B9A\u8F93\u5165\u7684\u8986\u76D6 =====\n")
print(ici_input_coverage, row.names = FALSE)
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", INTERSECTION_DIR, "\n")
PROJECT_DIR <- CFG$root
INTERSECTION_FILE <- artifact("stage03_intersection_candidate_bundle.rds")
if (!file.exists(INTERSECTION_FILE)) {
    stop("\u627E\u4E0D\u5230\u4EA4\u96C6\u5019\u9009\u7ED3\u679C\u3002")
}
intersection_bundle <- readRDS(INTERSECTION_FILE)
input_manifest <- intersection_bundle$input_manifest
ici_index <- which(input_manifest$input_name == "ici_pc_candidates")
if (length(ici_index) != 1L) {
    stop("\u4EA4\u96C6\u7ED3\u679C\u4E2D\u6CA1\u6709\u552F\u4E00\u7684\u514D\u75AB\u6CBB\u7597\u5019\u9009\u6587\u4EF6\u8BB0\u5F55\u3002")
}
ICI_PC_FILE <- input_manifest$path[ici_index]
current_md5 <- unname(tools::md5sum(ICI_PC_FILE))
if (is.na(current_md5) || current_md5 != input_manifest$md5[ici_index]) {
    stop("\u514D\u75AB\u6CBB\u7597\u5019\u9009 RDS \u7F3A\u5931\u6216\u5DF2\u53D1\u751F\u6539\u53D8\u3002")
}
ici_bundle <- readRDS(ICI_PC_FILE)
fixed_genes <- intersection_bundle$intersection_genes
gene_manifest <- intersection_bundle$gene_manifest
stopifnot(length(fixed_genes) == 17101L, !anyDuplicated(fixed_genes), identical(fixed_genes, gene_manifest$ensembl_gene_id))
ALIGNED_DIR <- out_dir("04_aligned")
if (dir.exists(ALIGNED_DIR) || !dir.create(ALIGNED_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_aligned_tsv <- function(x, filename) {
    write.table(x, file.path(ALIGNED_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
ici_aligned <- list()
availability_rows <- list()
for (cohort in names(ici_bundle$pc_candidates)) {
    obj <- ici_bundle$pc_candidates[[cohort]]
    source_tpm <- obj$expression_tpm
    source_index <- match(fixed_genes, rownames(source_tpm))
    observed_gene <- !is.na(source_index)
    aligned_tpm <- matrix(NA_real_, nrow = length(fixed_genes), ncol = ncol(source_tpm), dimnames = list(fixed_genes, 
        colnames(source_tpm)))
    aligned_tpm[observed_gene, ] <- source_tpm[source_index[observed_gene], , drop = FALSE]
    observed_mask <- matrix(rep(observed_gene, times = ncol(source_tpm)), nrow = length(fixed_genes), ncol = ncol(source_tpm), 
        dimnames = dimnames(aligned_tpm))
    if (any(!is.finite(aligned_tpm[observed_mask])) || any(aligned_tpm[observed_mask] < 0)) {
        stop(cohort, "\uFF1A\u5DF2\u89C2\u6D4B\u4F4D\u7F6E\u5B58\u5728\u65E0\u6548\u503C\u3002")
    }
    stopifnot(all(is.na(aligned_tpm[!observed_mask])), identical(colnames(aligned_tpm), obj$metadata$sample_id))
    aligned_log2tpm <- log2(aligned_tpm + 1)
    ici_aligned[[cohort]] <- list(expression_tpm = aligned_tpm, expression_log2tpm = aligned_log2tpm, observed_mask = observed_mask, 
        metadata = obj$metadata)
    availability <- data.frame(cohort = cohort, ensembl_gene_id = fixed_genes, gene_symbol = gene_manifest$gene_symbol, 
        observed = observed_gene, stringsAsFactors = FALSE)
    availability_rows[[cohort]] <- availability
    write_aligned_tsv(availability[!observed_gene, , drop = FALSE], paste0(cohort, "_missing_fixed_input_genes.tsv"))
}
gene_availability <- do.call(rbind, availability_rows)
rownames(gene_availability) <- NULL
cat("\n\u56FA\u5B9A\u8F93\u5165\u57FA\u56E0\u6570\uFF1A", length(fixed_genes), "\n")
for (cohort in names(ici_aligned)) {
    obj <- ici_aligned[[cohort]]
    cat(cohort, " | \u6837\u672C\u6570\uFF1A", ncol(obj$expression_tpm), " | \u7F3A\u6D4B\u57FA\u56E0\uFF1A", sum(!obj$observed_mask[, 1L]), 
        "\n", sep = "")
}
token_members <- intersection_bundle$token_members
if (any(token_members$n_intersection_target_ids > 1L)) {
    stop("\u5B58\u5728 token symbol \u5BF9\u5E94\u591A\u4E2A\u4EA4\u96C6\u76EE\u6807 ID\uFF0C\u9700\u8981\u5148\u5904\u7406\u3002")
}
token_names <- unique(token_members$token_name)
token_gene_indices <- list()
token_gene_ids <- list()
token_missing_rows <- list()
for (token_name in token_names) {
    members <- token_members[token_members$token_name == token_name & token_members$present_in_intersection, , 
        drop = FALSE]
    ids <- unique(members$intersection_target_ids)
    if (length(ids) == 0L || anyNA(ids) || any(!nzchar(ids))) {
        stop(token_name, "\uFF1A\u6CA1\u6709\u53EF\u7528\u7684\u4EA4\u96C6\u6210\u5458\u3002")
    }
    indices <- match(ids, fixed_genes)
    if (anyNA(indices)) {
        stop(token_name, "\uFF1A\u6210\u5458\u65E0\u6CD5\u5B9A\u4F4D\u5230\u56FA\u5B9A\u77E9\u9635\u3002")
    }
    token_gene_ids[[token_name]] <- ids
    token_gene_indices[[token_name]] <- indices
    for (cohort in names(ici_aligned)) {
        observed <- ici_aligned[[cohort]]$observed_mask[indices, 1L]
        missing_ids <- ids[!observed]
        missing_symbols <- gene_manifest$gene_symbol[match(missing_ids, gene_manifest$ensembl_gene_id)]
        key <- paste(cohort, token_name, sep = "::")
        token_missing_rows[[key]] <- data.frame(cohort = cohort, token_name = token_name, concept_id = members$concept_id[1L], 
            module_id = members$module_id[1L], n_fixed_token_members = length(ids), n_observed_members = sum(observed), 
            n_missing_members = sum(!observed), observed_member_fraction = mean(observed), missing_symbols = paste(missing_symbols, 
                collapse = ";"), missing_ensembl_ids = paste(missing_ids, collapse = ";"), stringsAsFactors = FALSE)
    }
}
token_missing_report <- do.call(rbind, token_missing_rows)
rownames(token_missing_report) <- NULL
write_aligned_tsv(gene_availability, "ICI_fixed_gene_availability.tsv")
write_aligned_tsv(token_missing_report, "ICI_token_missingness_report.tsv")
token_index_table <- do.call(rbind, lapply(names(token_gene_indices), function(nm) {
    data.frame(token_name = nm, gene_row_index = token_gene_indices[[nm]], ensembl_gene_id = token_gene_ids[[nm]], 
        stringsAsFactors = FALSE)
}))
write_aligned_tsv(token_index_table, "token_fixed_gene_indices.tsv")
saveRDS(list(fixed_genes = fixed_genes, gene_manifest = gene_manifest, ici_aligned = ici_aligned, gene_availability = gene_availability, 
    token_gene_indices = token_gene_indices, token_gene_ids = token_gene_ids, token_missing_report = token_missing_report, 
    provenance = list(intersection_file = INTERSECTION_FILE, intersection_md5 = unname(tools::md5sum(INTERSECTION_FILE)), 
        ici_pc_file = ICI_PC_FILE, ici_pc_md5 = current_md5, missing_values_filled = FALSE, missing_gene_policy = "NA_with_observed_mask", 
        final_model_gene_list_locked = FALSE)), artifact("stage03_ici_aligned_with_mask_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(ALIGNED_DIR, "sessionInfo.txt"))
affected_tokens <- token_missing_report[token_missing_report$n_missing_members > 0L, , drop = FALSE]
cat("\n===== \u53D7\u7F3A\u6D4B\u5F71\u54CD\u7684 token =====\n")
print(affected_tokens[c("cohort", "token_name", "n_fixed_token_members", "n_observed_members", "n_missing_members", 
    "observed_member_fraction", "missing_symbols")], row.names = FALSE)
cat("\n\u5B8C\u5168\u65E0\u89C2\u6D4B\u6210\u5458\u7684 cohort-token \u7EC4\u5408\u6570\uFF1A", sum(token_missing_report$n_observed_members == 0L), 
    "\n")
cat("\nToken \u6570\uFF1A", length(token_gene_indices), "\n")
cat("\u8F93\u51FA\u76EE\u5F55\uFF1A", ALIGNED_DIR, "\n")

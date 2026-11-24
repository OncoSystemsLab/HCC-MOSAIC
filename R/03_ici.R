# Align the HCC and ICI cohorts to a shared gene background.
PROJECT_DIR <- CFG$root
INVENTORY_FILE <- artifact("stage02_input_inventory.rds")
if (!file.exists(INVENTORY_FILE)) {
    stop("\u627E\u4E0D\u5230\u4E0A\u4E00\u90E8\u5206\u7684\u8F93\u5165\u6E05\u5355\u3002")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("\u8BF7\u5148\u5B89\u88C5 readxl\u3002")
}
inventory <- readRDS(INVENTORY_FILE)
manifest <- inventory$input_manifest
ICI_RUN_DIR <- out_dir("03_ici_raw")
if (dir.exists(ICI_RUN_DIR) || !dir.create(ICI_RUN_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\uFF0C\u8BF7\u91CD\u65B0\u8FD0\u884C\u672C\u6BB5\u3002")
}
get_input_path <- function(cohort, role) {
    rows <- manifest$cohort == cohort & manifest$role == role
    if (sum(rows) != 1L) {
        stop("\u6587\u4EF6\u8BB0\u5F55\u4E0D\u662F\u552F\u4E00\u5339\u914D\uFF1A", cohort, " / ", role)
    }
    path <- manifest$path[rows]
    if (!file.exists(path)) {
        stop("\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A", path)
    }
    path
}
read_tsv_character <- function(path) {
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(con), add = TRUE)
    read.delim(con, header = TRUE, sep = "\t", quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, 
        colClasses = "character", na.strings = c("", "NA"), fill = FALSE)
}
make_expression_matrix <- function(dat, id_column, cohort) {
    if (sum(names(dat) == id_column) != 1L) {
        stop(cohort, "\uFF1A\u57FA\u56E0 ID \u5217\u7F3A\u5931\u6216\u4E0D\u552F\u4E00\u3002")
    }
    if (anyNA(names(dat)) || any(!nzchar(trimws(names(dat)))) || anyDuplicated(names(dat))) {
        stop(cohort, "\uFF1A\u5217\u540D\u5B58\u5728\u7A7A\u503C\u6216\u91CD\u590D\u3002")
    }
    ids <- as.character(dat[[id_column]])
    if (anyNA(ids) || any(!nzchar(trimws(ids)))) {
        stop(cohort, "\uFF1A\u5B58\u5728\u7A7A\u57FA\u56E0 ID\u3002")
    }
    if (any(ids != trimws(ids)) || any(names(dat) != trimws(names(dat)))) {
        stop(cohort, "\uFF1A\u57FA\u56E0 ID \u6216\u5217\u540D\u542B\u4E24\u4FA7\u7A7A\u683C\uFF0C\u9700\u8981\u68C0\u67E5\u3002")
    }
    if (anyDuplicated(ids)) {
        stop(cohort, "\uFF1A\u539F\u59CB\u57FA\u56E0 ID \u91CD\u590D\uFF0C\u9700\u8981\u5148\u68C0\u67E5\uFF0C\u6682\u4E0D\u81EA\u52A8\u5408\u5E76\u3002")
    }
    sample_columns <- setdiff(names(dat), id_column)
    if (length(sample_columns) == 0L) {
        stop(cohort, "\uFF1A\u6CA1\u6709\u8868\u8FBE\u6837\u672C\u5217\u3002")
    }
    x_text <- as.matrix(dat[sample_columns])
    x_numeric <- suppressWarnings(as.numeric(x_text))
    x <- matrix(x_numeric, nrow = nrow(dat), ncol = length(sample_columns), dimnames = list(ids, sample_columns))
    if (any(!is.finite(x))) {
        stop(cohort, "\uFF1A\u8868\u8FBE\u77E9\u9635\u5305\u542B\u7F3A\u5931\u3001\u975E\u6570\u503C\u6216\u65E0\u7A77\u503C\u3002")
    }
    if (any(x < 0)) {
        stop(cohort, "\uFF1A\u8868\u8FBE\u77E9\u9635\u5B58\u5728\u8D1F\u6570\uFF0C\u9700\u8981\u6838\u5B9E\u6570\u636E\u5C3A\u5EA6\u3002")
    }
    x
}
save_ici_tsv <- function(x, filename) {
    write.table(x, file.path(ICI_RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
required_objects <- c("manifest", "INVENTORY_FILE", "ICI_RUN_DIR", "get_input_path", "read_tsv_character", "make_expression_matrix", 
    "save_ici_tsv")
missing_objects <- required_objects[!vapply(required_objects, function(name, ...) exists(name, inherits = TRUE), 
    logical(1), inherits = TRUE)]
if (length(missing_objects) > 0L) {
    stop("\u8BF7\u5148\u5B8C\u6574\u8FD0\u884C\u7B2C 13 \u6BB5\u3002\u7F3A\u5C11\uFF1A", paste(missing_objects, collapse = ", "))
}
expected_inventory_file <- artifact("stage02_input_inventory.rds")
if (normalizePath(INVENTORY_FILE, winslash = "/", mustWork = TRUE) != normalizePath(expected_inventory_file, winslash = "/", 
    mustWork = TRUE)) {
    stop("\u7B2C 13 \u6BB5\u4ECD\u5F15\u7528\u65E7\u8F93\u5165\u6E05\u5355\uFF0C\u8BF7\u4FEE\u6539\u8DEF\u5F84\u540E\u91CD\u8DD1\u7B2C 13 \u6BB5\u3002")
}
cohorts <- c("GSE181946", "GSE202069", "PRJEB34724")
ici_raw <- setNames(vector("list", length(cohorts)), cohorts)
sample_match_reports <- list()
label_reports <- list()
actual_input_records <- list()
expected_sample_counts <- c(GSE181946 = 15L, GSE202069 = 17L, PRJEB34724 = 35L)
for (cohort in cohorts) {
    cat("\n\u6B63\u5728\u8BFB\u53D6\uFF1A", cohort, "\n")
    if (cohort == "GSE181946") {
        source_path <- get_input_path(cohort, "expression_and_clinical")
        required_sheets <- c("Sheet1", "clinical")
        available_sheets <- readxl::excel_sheets(source_path)
        if (!all(required_sheets %in% available_sheets)) {
            stop(cohort, "\uFF1AExcel \u7F3A\u5C11 Sheet1 \u6216 clinical \u5DE5\u4F5C\u8868\u3002\u5B9E\u9645\u5DE5\u4F5C\u8868\uFF1A", paste(available_sheets, 
                collapse = ", "))
        }
        expr_df <- as.data.frame(readxl::read_excel(source_path, sheet = "Sheet1", col_types = "text", .name_repair = "check_unique"), 
            stringsAsFactors = FALSE)
        clin <- as.data.frame(readxl::read_excel(source_path, sheet = "clinical", col_types = "text", .name_repair = "check_unique"), 
            stringsAsFactors = FALSE)
        gene_id_column <- "GeneID"
        source_scale <- "counts_expected_not_yet_confirmed"
        if (!"GeneID" %in% names(expr_df)) {
            stop(cohort, "\uFF1ASheet1 \u672A\u627E\u5230 GeneID \u5217\u3002\u5B9E\u9645\u524D\u51E0\u5217\u4E3A\uFF1A", paste(head(names(expr_df)), 
                collapse = ", "))
        }
        gene_ids <- as.character(expr_df$GeneID)
        if (anyNA(gene_ids) || any(!grepl("^[0-9]+$", gene_ids))) {
            stop(cohort, "\uFF1AGeneID \u5217\u542B\u7A7A\u503C\u6216\u975E\u6570\u5B57 ID\uFF0C", "\u8BF7\u786E\u8BA4 Excel \u5DF2\u66FF\u6362\u4E3A NCBI GeneID \u8BA1\u6570\u77E9\u9635\u3002")
        }
        actual_input_records[[cohort]] <- data.frame(cohort = cohort, role = "expression_and_clinical", path = source_path, 
            expression_sheet = "Sheet1", clinical_sheet = "clinical", stringsAsFactors = FALSE)
    }
    else {
        expression_role <- if (cohort == "GSE202069") {
            "expression_counts"
        }
        else {
            "expression_TPM"
        }
        expression_path <- get_input_path(cohort, expression_role)
        clinical_path <- get_input_path(cohort, "clinical")
        expr_df <- read_tsv_character(expression_path)
        clin <- as.data.frame(readxl::read_excel(clinical_path, sheet = "Sheet1", col_types = "text", .name_repair = "check_unique"), 
            stringsAsFactors = FALSE)
        gene_id_column <- if (cohort == "GSE202069") {
            "GeneID"
        }
        else {
            "gene_id"
        }
        source_scale <- if (cohort == "GSE202069") {
            "counts_expected_not_yet_confirmed"
        }
        else {
            "TPM_as_provided"
        }
        actual_input_records[[cohort]] <- data.frame(cohort = cohort, role = c(expression_role, "clinical"), path = c(expression_path, 
            clinical_path), expression_sheet = c(NA_character_, NA_character_), clinical_sheet = c(NA_character_, 
            "Sheet1"), stringsAsFactors = FALSE)
    }
    x_all <- make_expression_matrix(expr_df, gene_id_column, cohort)
    if (cohort %in% c("GSE181946", "GSE202069")) {
        if (any(abs(x_all - round(x_all)) > 1e-08)) {
            stop(cohort, "\uFF1A\u8BA1\u6570\u77E9\u9635\u542B\u975E\u6574\u6570\u6570\u503C\uFF0C\u8BF7\u6838\u5B9E\u8868\u8FBE\u5C3A\u5EA6\u3002")
        }
    }
    if (!all(c("sample", "label") %in% names(clin))) {
        stop(cohort, "\uFF1A\u4E34\u5E8A\u8868\u7F3A\u5C11 sample \u6216 label \u5217\u3002")
    }
    sample_ids <- as.character(clin$sample)
    labels_original <- as.character(clin$label)
    if (anyNA(sample_ids) || any(!nzchar(trimws(sample_ids))) || anyDuplicated(sample_ids)) {
        stop(cohort, "\uFF1A\u4E34\u5E8A\u6837\u672C ID \u5B58\u5728\u7F3A\u5931\u3001\u7A7A\u503C\u6216\u91CD\u590D\u3002")
    }
    if (any(sample_ids != trimws(sample_ids))) {
        stop(cohort, "\uFF1A\u4E34\u5E8A\u6837\u672C ID \u4E24\u7AEF\u5B58\u5728\u7A7A\u683C\uFF0C\u8BF7\u68C0\u67E5\u3002")
    }
    labels_standard <- tolower(trimws(labels_original))
    label_report <- data.frame(cohort = cohort, sample_id = sample_ids, label_original = labels_original, label_standard = labels_standard, 
        stringsAsFactors = FALSE)
    label_reports[[cohort]] <- label_report
    save_ici_tsv(label_report, paste0(cohort, "_label_audit.tsv"))
    allowed_labels <- c("responder", "non-responder")
    if (anyNA(labels_standard) || any(!labels_standard %in% allowed_labels)) {
        stop(cohort, "\uFF1A\u5B58\u5728\u672A\u77E5\u6216\u7F3A\u5931\u6807\u7B7E\uFF0C\u8BF7\u67E5\u770B\u8BE5\u961F\u5217\u7684 label_audit.tsv\u3002")
    }
    response <- ifelse(labels_standard == "responder", 1L, 0L)
    all_sample_ids <- union(colnames(x_all), sample_ids)
    match_report <- data.frame(cohort = cohort, sample_id = all_sample_ids, in_expression = all_sample_ids %in% 
        colnames(x_all), in_clinical = all_sample_ids %in% sample_ids, stringsAsFactors = FALSE)
    sample_match_reports[[cohort]] <- match_report
    save_ici_tsv(match_report, paste0(cohort, "_sample_matching.tsv"))
    missing_expression <- setdiff(sample_ids, colnames(x_all))
    if (length(missing_expression) > 0L) {
        stop(cohort, "\uFF1A\u4EE5\u4E0B\u4E34\u5E8A\u6837\u672C\u7F3A\u5C11\u8868\u8FBE\uFF1A", paste(missing_expression, collapse = ", "))
    }
    extra_expression <- setdiff(colnames(x_all), sample_ids)
    if (cohort == "GSE181946" && length(extra_expression) > 0L) {
        stop(cohort, "\uFF1A\u8868\u8FBE\u4E2D\u51FA\u73B0\u4E0D\u5728\u4E34\u5E8A\u8868\u4E2D\u7684\u5217\uFF1A", paste(extra_expression, collapse = ", "))
    }
    x_selected <- x_all[, sample_ids, drop = FALSE]
    metadata <- data.frame(cohort = cohort, sample_id = sample_ids, patient_uid = paste(cohort, sample_ids, sep = "::"), 
        response = response, label_original = labels_original, source_scale = source_scale, stringsAsFactors = FALSE)
    stopifnot(identical(colnames(x_selected), metadata$sample_id))
    if (length(unique(metadata$response)) != 2L) {
        stop(cohort, "\uFF1A\u672A\u540C\u65F6\u5305\u542B\u5E94\u7B54\u548C\u4E0D\u5E94\u7B54\u6837\u672C\u3002")
    }
    if (ncol(x_selected) != expected_sample_counts[[cohort]]) {
        stop(cohort, "\uFF1A\u5B9E\u9645\u9009\u4E2D ", ncol(x_selected), " \u4F8B\uFF0C\u9884\u671F\u4E3A ", expected_sample_counts[[cohort]], 
            " \u4F8B\u3002\u8BF7\u5148\u6838\u67E5\uFF0C\u4E0D\u8981\u5F3A\u884C\u5220\u51CF\u3002")
    }
    ici_raw[[cohort]] <- list(expression = x_selected, metadata = metadata, n_expression_columns_before_selection = ncol(x_all), 
        gene_id_column = gene_id_column, source_scale = source_scale, input_sources = actual_input_records[[cohort]])
    saveRDS(ici_raw[[cohort]], file.path(ICI_RUN_DIR, paste0(cohort, "_raw_selected.rds")))
    cat("\u5B8C\u6210\uFF1A", cohort, " | \u57FA\u56E0\u884C\u6570\uFF1A", nrow(x_selected), " | \u6837\u672C\u6570\uFF1A", ncol(x_selected), " | \u5E94\u7B54\uFF1A", 
        sum(response == 1L), " | \u4E0D\u5E94\u7B54\uFF1A", sum(response == 0L), "\n", sep = "")
    rm(expr_df, x_all, x_selected)
}
actual_inputs <- do.call(rbind, actual_input_records)
rownames(actual_inputs) <- NULL
actual_inputs$md5 <- unname(tools::md5sum(actual_inputs$path))
if (anyNA(actual_inputs$md5)) {
    stop("\u90E8\u5206\u8F93\u5165\u6587\u4EF6\u65E0\u6CD5\u8BA1\u7B97 MD5\uFF0C\u8BF7\u68C0\u67E5\u6587\u4EF6\u662F\u5426\u53EF\u8BFB\u53D6\u3002")
}
save_ici_tsv(actual_inputs, "ici_actual_input_sources.tsv")
cat("\n\u7B2C 14 \u6BB5\u5B8C\u6210\uFF0C\u4E09\u4E2A\u961F\u5217\u5747\u5DF2\u91CD\u65B0\u8BFB\u53D6\u3002\n")
cat("\u8F93\u51FA\u76EE\u5F55\uFF1A", ICI_RUN_DIR, "\n")
sample_qc <- do.call(rbind, lapply(names(ici_raw), function(cohort) {
    x <- ici_raw[[cohort]]$expression
    data.frame(cohort = cohort, sample_id = colnames(x), source_scale = ici_raw[[cohort]]$source_scale, expression_sum = colSums(x), 
        n_detected_genes = colSums(x > 0), zero_fraction = colMeans(x == 0), all_zero = colSums(x > 0) == 0L, stringsAsFactors = FALSE)
}))
expected_n <- c(GSE181946 = 15L, GSE202069 = 17L, PRJEB34724 = 35L)
cohort_summary <- do.call(rbind, lapply(names(ici_raw), function(cohort) {
    obj <- ici_raw[[cohort]]
    x <- obj$expression
    meta <- obj$metadata
    data.frame(cohort = cohort, n_genes_raw = nrow(x), n_expression_before_selection = obj$n_expression_columns_before_selection, 
        n_selected = ncol(x), expected_n = unname(expected_n[cohort]), n_responders = sum(meta$response == 1L), 
        n_nonresponders = sum(meta$response == 0L), fraction_non_integer = mean(abs(x - round(x)) > 1e-08), min_sample_sum = min(colSums(x)), 
        max_sample_sum = max(colSums(x)), stringsAsFactors = FALSE)
}))
metadata_all <- do.call(rbind, lapply(ici_raw, function(x) x$metadata))
rownames(metadata_all) <- NULL
sample_matching <- do.call(rbind, sample_match_reports)
save_ici_tsv(sample_qc, "ici_sample_numeric_qc.tsv")
save_ici_tsv(cohort_summary, "ici_cohort_summary.tsv")
save_ici_tsv(metadata_all, "ici_metadata_audited.tsv")
save_ici_tsv(sample_matching, "ici_sample_matching_all.tsv")
used_files <- manifest[manifest$cohort %in% cohorts & manifest$role %in% c("expression_and_clinical", "expression_counts", 
    "expression_TPM", "clinical"), c("cohort", "role", "path")]
used_files$md5 <- unname(tools::md5sum(used_files$path))
if (anyNA(used_files$md5)) {
    stop("\u90E8\u5206\u8F93\u5165\u6587\u4EF6\u65E0\u6CD5\u8BA1\u7B97 MD5\u3002")
}
save_ici_tsv(used_files, "ici_inputs_md5.tsv")
saveRDS(list(ici_raw = ici_raw, metadata = metadata_all, cohort_summary = cohort_summary, sample_qc = sample_qc, 
    input_files = used_files, inventory_file = INVENTORY_FILE), artifact("stage02_ici_raw_audit_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(ICI_RUN_DIR, "sessionInfo.txt"))
cat("\n===== \u961F\u5217\u6C47\u603B =====\n")
print(cohort_summary, row.names = FALSE)
cat("\n===== \u6837\u672C\u6570\u503C\u8D28\u63A7 =====\n")
print(sample_qc, row.names = FALSE)
cat("\n===== \u539F\u59CB\u57FA\u56E0 ID \u793A\u4F8B =====\n")
for (cohort in names(ici_raw)) {
    cat(cohort, "\uFF1A", paste(head(rownames(ici_raw[[cohort]]$expression)), collapse = ", "), "\n")
}
cat("\n\u76EE\u6807\u6837\u672C\u603B\u6570\uFF1A", nrow(metadata_all), "\n")
cat("\u7ED3\u679C\u76EE\u5F55\uFF1A", ICI_RUN_DIR, "\n")
if (any(cohort_summary$n_selected != cohort_summary$expected_n)) {
    warning("\u5B9E\u9645\u6837\u672C\u6570\u4E0E\u9884\u671F\u4E0D\u4E00\u81F4\uFF0C\u8BF7\u5148\u53CD\u9988\u7ED3\u679C\uFF0C\u4E0D\u8981\u5F3A\u884C\u5220\u51CF\u81F3\u9884\u671F\u4EBA\u6570\u3002")
}
if (any(sample_qc$all_zero)) {
    warning("\u53D1\u73B0\u5168\u96F6\u6837\u672C\uFF0C\u9700\u8981\u5148\u68C0\u67E5\uFF0C\u6682\u4E0D\u8FDB\u5165\u540E\u7EED\u5904\u7406\u3002")
}
PROJECT_DIR <- CFG$root
RAW_BUNDLE_FILE <- artifact("stage02_ici_raw_audit_bundle.rds")
if (!file.exists(RAW_BUNDLE_FILE)) {
    stop("\u627E\u4E0D\u5230\u65B0\u7248\u539F\u59CB\u8868\u8FBE\u7ED3\u679C\uFF0C\u8BF7\u68C0\u67E5\u8DEF\u5F84\u3002")
}
raw_bundle <- readRDS(RAW_BUNDLE_FILE)
ici_raw <- raw_bundle$ici_raw
stopifnot(nrow(ici_raw$GSE181946$expression) == 39376L, ncol(ici_raw$GSE181946$expression) == 15L, ici_raw$GSE181946$gene_id_column == 
    "GeneID", nrow(raw_bundle$metadata) == 67L)
input_files <- raw_bundle$input_files
current_md5 <- unname(tools::md5sum(input_files$path))
if (anyNA(current_md5) || any(current_md5 != input_files$md5)) {
    stop("\u539F\u59CB\u8F93\u5165\u6587\u4EF6\u5DF2\u6539\u53D8\uFF0C\u8BF7\u5148\u91CD\u65B0\u5B8C\u6210\u7B2C 13\u201315 \u6BB5\u3002")
}
inventory <- readRDS(raw_bundle$inventory_file)
manifest <- inventory$input_manifest
MAP_RUN_DIR <- out_dir("03_ici_annotation")
if (dir.exists(MAP_RUN_DIR) || !dir.create(MAP_RUN_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_map_tsv <- function(x, filename) {
    write.table(x, file.path(MAP_RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
strip_version <- function(x) {
    sub("\\.[0-9]+$", "", x)
}
annotation_results <- list()
annotation_input_files <- list()
for (cohort in names(ici_raw)) {
    cat("\n\u6B63\u5728\u5339\u914D\u539F\u59CB\u6CE8\u91CA\uFF1A", cohort, "\n")
    rows <- manifest$cohort == cohort & manifest$role == "gene_annotation"
    if (sum(rows) != 1L) {
        stop(cohort, "\uFF1A\u6CE8\u91CA\u6587\u4EF6\u8BB0\u5F55\u4E0D\u552F\u4E00\u3002")
    }
    annotation_path <- manifest$path[rows]
    if (!file.exists(annotation_path)) {
        stop("\u6CE8\u91CA\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A", annotation_path)
    }
    ann <- read.delim(annotation_path, header = TRUE, sep = "\t", quote = "", comment.char = "", colClasses = "character", 
        check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"), fill = FALSE)
    x <- ici_raw[[cohort]]$expression
    is_count_cohort <- cohort %in% c("GSE181946", "GSE202069")
    if (is_count_cohort) {
        required_columns <- c("GeneID", "Symbol", "EnsemblGeneID", "Length")
        if (!all(required_columns %in% names(ann))) {
            stop(cohort, "\uFF1A\u6CE8\u91CA\u8868\u7F3A\u5C11\u5FC5\u8981\u5217\u3002")
        }
        if (ici_raw[[cohort]]$gene_id_column != "GeneID") {
            stop(cohort, "\uFF1A\u8868\u8FBE ID \u7C7B\u578B\u4E0D\u662F\u9884\u671F\u7684 GeneID\u3002")
        }
        standardized_ann <- data.frame(source_key = ann$GeneID, source_symbol = ann$Symbol, source_ensembl_id = strip_version(ann$EnsemblGeneID), 
            length_bp = suppressWarnings(as.numeric(ann$Length)), stringsAsFactors = FALSE)
        expression_key <- rownames(x)
    }
    else {
        required_columns <- c("gene_id", "gene_symbol", "gene_type")
        if (!all(required_columns %in% names(ann))) {
            stop(cohort, "\uFF1A\u6CE8\u91CA\u8868\u7F3A\u5C11\u5FC5\u8981\u5217\u3002")
        }
        standardized_ann <- data.frame(source_key = strip_version(ann$gene_id), source_symbol = ann$gene_symbol, 
            source_ensembl_id = strip_version(ann$gene_id), length_bp = NA_real_, stringsAsFactors = FALSE)
        expression_key <- strip_version(rownames(x))
    }
    if (anyNA(expression_key) || any(!nzchar(expression_key)) || anyDuplicated(expression_key)) {
        stop(cohort, "\uFF1A\u8868\u8FBE\u5339\u914D ID \u5B58\u5728\u7F3A\u5931\u3001\u7A7A\u503C\u6216\u91CD\u590D\u3002")
    }
    standardized_ann <- unique(standardized_ann)
    valid_key <- !is.na(standardized_ann$source_key) & nzchar(standardized_ann$source_key)
    ann_valid <- standardized_ann[valid_key, , drop = FALSE]
    key_counts <- table(ann_valid$source_key)
    n_records <- as.integer(key_counts[expression_key])
    n_records[is.na(n_records)] <- 0L
    annotation_index <- match(expression_key, ann_valid$source_key)
    annotation_index[n_records != 1L] <- NA_integer_
    mapping <- data.frame(original_gene_id = rownames(x), source_key = expression_key, n_annotation_records = n_records, 
        source_symbol = ann_valid$source_symbol[annotation_index], source_ensembl_id = ann_valid$source_ensembl_id[annotation_index], 
        length_bp = ann_valid$length_bp[annotation_index], stringsAsFactors = FALSE)
    mapping$mapping_status <- ifelse(n_records == 0L, "unmatched", ifelse(n_records == 1L, "unique", "multiple_annotations"))
    mapping$length_valid <- is.finite(mapping$length_bp) & mapping$length_bp > 0
    usable <- mapping$mapping_status == "unique"
    if (is_count_cohort) {
        usable <- usable & mapping$length_valid
    }
    conflict_keys <- names(key_counts)[key_counts > 1L]
    conflicts <- ann_valid[ann_valid$source_key %in% intersect(expression_key, conflict_keys), , drop = FALSE]
    annotation_results[[cohort]] <- list(mapping = mapping, conflicts = conflicts, usable_for_next_step = usable, 
        annotation_path = annotation_path, mapping_key_type = if (is_count_cohort) {
            "NCBI_GeneID"
        } else {
            "Ensembl_gene_ID_without_version"
        })
    annotation_input_files[[cohort]] <- data.frame(cohort = cohort, path = annotation_path, md5 = unname(tools::md5sum(annotation_path)), 
        stringsAsFactors = FALSE)
    write_map_tsv(mapping, paste0(cohort, "_source_annotation_mapping.tsv"))
    write_map_tsv(conflicts, paste0(cohort, "_conflicting_annotation_records.tsv"))
}
summary_rows <- list()
for (cohort in names(ici_raw)) {
    x <- ici_raw[[cohort]]$expression
    result <- annotation_results[[cohort]]
    mapping <- result$mapping
    usable <- result$usable_for_next_step
    stopifnot(identical(rownames(x), mapping$original_gene_id), length(usable) == nrow(x))
    total_expression <- colSums(x)
    if (any(total_expression <= 0)) {
        stop(cohort, "\uFF1A\u5B58\u5728\u603B\u8868\u8FBE\u91CF\u4E3A\u96F6\u7684\u6837\u672C\u3002")
    }
    excluded_fraction <- colSums(x[!usable, , drop = FALSE])/total_expression
    sample_report <- data.frame(cohort = cohort, sample_id = colnames(x), expression_fraction_not_usable = excluded_fraction, 
        stringsAsFactors = FALSE)
    write_map_tsv(sample_report, paste0(cohort, "_annotation_expression_coverage.tsv"))
    problem_index <- which(!usable)
    expression_total <- rowSums(x)
    problem_index <- problem_index[order(expression_total[problem_index], decreasing = TRUE)]
    top_index <- head(problem_index, 20L)
    top_problems <- mapping[top_index, , drop = FALSE]
    top_problems$total_expression <- expression_total[top_index]
    write_map_tsv(top_problems, paste0(cohort, "_top_unusable_expression_rows.tsv"))
    missing_symbol <- is.na(mapping$source_symbol) | !nzchar(mapping$source_symbol)
    missing_ensembl <- is.na(mapping$source_ensembl_id) | !nzchar(mapping$source_ensembl_id)
    is_count_cohort <- cohort %in% c("GSE181946", "GSE202069")
    summary_rows[[cohort]] <- data.frame(cohort = cohort, n_expression_rows = nrow(x), n_samples = ncol(x), n_unique_mapping = sum(mapping$mapping_status == 
        "unique"), n_unmatched = sum(mapping$mapping_status == "unmatched"), n_multiple_annotations = sum(mapping$mapping_status == 
        "multiple_annotations"), n_unique_with_invalid_length = if (is_count_cohort) {
        sum(mapping$mapping_status == "unique" & !mapping$length_valid)
    }
    else {
        NA_integer_
    }, n_usable_for_next_step = sum(usable), n_missing_symbol = sum(missing_symbol), n_missing_ensembl_id = sum(missing_ensembl), 
        max_expression_fraction_not_usable = max(excluded_fraction), stringsAsFactors = FALSE)
    cat("\n===== ", cohort, "\uFF1A\u9AD8\u8868\u8FBE\u5F85\u68C0\u67E5\u884C =====\n", sep = "")
    print(top_problems, row.names = FALSE)
}
mapping_summary <- do.call(rbind, summary_rows)
rownames(mapping_summary) <- NULL
annotation_files <- do.call(rbind, annotation_input_files)
rownames(annotation_files) <- NULL
if (anyNA(annotation_files$md5)) {
    stop("\u5B58\u5728\u65E0\u6CD5\u8BA1\u7B97 MD5 \u7684\u6CE8\u91CA\u6587\u4EF6\u3002")
}
write_map_tsv(mapping_summary, "source_mapping_summary.tsv")
write_map_tsv(annotation_files, "annotation_input_md5.tsv")
saveRDS(list(raw_bundle_file = RAW_BUNDLE_FILE, raw_bundle_md5 = unname(tools::md5sum(RAW_BUNDLE_FILE)), annotation_results = annotation_results, 
    mapping_summary = mapping_summary, annotation_input_files = annotation_files), artifact("stage02_ici_annotation_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(MAP_RUN_DIR, "sessionInfo.txt"))
cat("\n===== \u65B0\u7248\u539F\u59CB\u6CE8\u91CA\u6620\u5C04\u6C47\u603B =====\n")
print(mapping_summary, row.names = FALSE)
cat("\n\u539F\u59CB\u8868\u8FBE\u6765\u6E90\uFF1A", RAW_BUNDLE_FILE, "\n")
cat("\u672C\u6B21\u8F93\u51FA\u76EE\u5F55\uFF1A", MAP_RUN_DIR, "\n")
if (any(mapping_summary$n_unmatched > 0L | mapping_summary$n_multiple_annotations > 0L | mapping_summary$n_usable_for_next_step < 
    mapping_summary$n_expression_rows)) {
    warning("\u5B58\u5728\u6620\u5C04\u6216\u957F\u5EA6\u95EE\u9898\uFF0C\u8BF7\u5148\u53CD\u9988\u7ED3\u679C\uFF0C\u6682\u4E0D\u8F6C\u6362 TPM\u3002")
}
PROJECT_DIR <- CFG$root
ANNOTATION_BUNDLE_FILE <- artifact("stage02_ici_annotation_bundle.rds")
if (!file.exists(ANNOTATION_BUNDLE_FILE)) {
    stop("\u627E\u4E0D\u5230\u6700\u65B0\u7684\u6CE8\u91CA\u5BA1\u8BA1\u7ED3\u679C\u3002")
}
annotation_bundle <- readRDS(ANNOTATION_BUNDLE_FILE)
RAW_BUNDLE_FILE <- annotation_bundle$raw_bundle_file
raw_md5 <- unname(tools::md5sum(RAW_BUNDLE_FILE))
if (is.na(raw_md5) || raw_md5 != annotation_bundle$raw_bundle_md5) {
    stop("\u539F\u59CB\u8868\u8FBE RDS \u7F3A\u5931\u6216\u5185\u5BB9\u5DF2\u6539\u53D8\u3002")
}
raw_bundle <- readRDS(RAW_BUNDLE_FILE)
ici_raw <- raw_bundle$ici_raw
annotation_results <- annotation_bundle$annotation_results
annotation_files <- annotation_bundle$annotation_input_files
current_annotation_md5 <- unname(tools::md5sum(annotation_files$path))
if (anyNA(current_annotation_md5) || any(current_annotation_md5 != annotation_files$md5)) {
    stop("\u539F\u59CB\u6CE8\u91CA\u6587\u4EF6\u7F3A\u5931\u6216\u5DF2\u6539\u53D8\uFF0C\u8BF7\u91CD\u65B0\u5BA1\u8BA1\u3002")
}
expected_n <- c(GSE181946 = 15L, GSE202069 = 17L, PRJEB34724 = 35L)
stopifnot(setequal(names(ici_raw), names(expected_n)), nrow(raw_bundle$metadata) == 67L)
TPM_RUN_DIR <- out_dir("03_ici_TPM")
if (dir.exists(TPM_RUN_DIR) || !dir.create(TPM_RUN_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_tpm_tsv <- function(x, filename) {
    write.table(x, file.path(TPM_RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
ici_tpm <- list()
for (cohort in names(expected_n)) {
    obj <- ici_raw[[cohort]]
    x <- obj$expression
    ann <- annotation_results[[cohort]]$mapping
    stopifnot(ncol(x) == expected_n[[cohort]], identical(rownames(x), ann$original_gene_id), identical(colnames(x), 
        obj$metadata$sample_id))
    if (any(!is.finite(x)) || any(x < 0)) {
        stop(cohort, "\uFF1A\u8868\u8FBE\u77E9\u9635\u5B58\u5728\u65E0\u6548\u6216\u8D1F\u503C\u3002")
    }
    if (any(ann$mapping_status != "unique")) {
        stop(cohort, "\uFF1A\u4ECD\u5B58\u5728\u975E\u552F\u4E00\u6CE8\u91CA\uFF0C\u6682\u4E0D\u8F6C\u6362\u3002")
    }
    if (cohort %in% c("GSE181946", "GSE202069")) {
        if (any(!is.finite(ann$length_bp)) || any(ann$length_bp <= 0)) {
            stop(cohort, "\uFF1A\u5B58\u5728\u65E0\u6548\u957F\u5EA6\u3002")
        }
        rpk <- sweep(x, MARGIN = 1L, STATS = ann$length_bp/1000, FUN = "/")
        rpk_sum <- colSums(rpk)
        if (any(!is.finite(rpk_sum)) || any(rpk_sum <= 0)) {
            stop(cohort, "\uFF1ARPK \u5206\u6BCD\u65E0\u6548\u3002")
        }
        tpm <- sweep(rpk, MARGIN = 2L, STATS = rpk_sum, FUN = "/") * 1e+06
        method <- "derived_TPM_from_counts_and_source_Length"
    }
    else {
        tpm <- x
        rpk_sum <- rep(NA_real_, ncol(x))
        method <- "provided_TPM_unchanged"
    }
    if (any(!is.finite(tpm)) || any(tpm < 0)) {
        stop(cohort, "\uFF1ATPM \u4EA7\u751F\u65E0\u6548\u503C\u3002")
    }
    stopifnot(identical(dimnames(tpm), dimnames(x)))
    ici_tpm[[cohort]] <- list(expression_tpm = tpm, gene_annotation = ann, metadata = obj$metadata, expression_method = method, 
        rpk_denominator = rpk_sum)
    saveRDS(ici_tpm[[cohort]], file.path(TPM_RUN_DIR, paste0(cohort, "_TPM_sourceIDs.rds")))
    cat("\u5B8C\u6210\uFF1A", cohort, " | \u57FA\u56E0\u884C\u6570\uFF1A", nrow(tpm), " | \u6837\u672C\u6570\uFF1A", ncol(tpm), "\n", sep = "")
}
tpm_qc_rows <- list()
tpm_summary_rows <- list()
for (cohort in names(ici_tpm)) {
    result <- ici_tpm[[cohort]]
    tpm <- result$expression_tpm
    original <- ici_raw[[cohort]]$expression
    if (!identical(tpm == 0, original == 0)) {
        stop(cohort, "\uFF1A\u8F6C\u6362\u540E\u96F6\u503C\u7ED3\u6784\u53D1\u751F\u5F02\u5E38\u53D8\u5316\u3002")
    }
    if (cohort == "PRJEB34724" && !identical(tpm, original)) {
        stop("PRJEB34724 \u539F\u6709 TPM \u88AB\u610F\u5916\u6539\u53D8\u3002")
    }
    tpm_sum <- colSums(tpm)
    if (any(abs(tpm_sum - 1e+06) > 1)) {
        stop(cohort, "\uFF1A\u90E8\u5206\u6837\u672C TPM \u603B\u548C\u504F\u79BB\u4E00\u767E\u4E07\u8D85\u8FC7 1\uFF0C\u8BF7\u68C0\u67E5\u3002")
    }
    tpm_qc_rows[[cohort]] <- data.frame(cohort = cohort, sample_id = colnames(tpm), tpm_sum = tpm_sum, absolute_sum_error = abs(tpm_sum - 
        1e+06), n_detected_genes = colSums(tpm > 0), zero_fraction = colMeans(tpm == 0), max_gene_tpm = apply(tpm, 
        2L, max), stringsAsFactors = FALSE)
    tpm_summary_rows[[cohort]] <- data.frame(cohort = cohort, n_genes = nrow(tpm), n_samples = ncol(tpm), n_responders = sum(result$metadata$response == 
        1L), n_nonresponders = sum(result$metadata$response == 0L), min_tpm_sum = min(tpm_sum), max_tpm_sum = max(tpm_sum), 
        max_absolute_sum_error = max(abs(tpm_sum - 1e+06)), expression_method = result$expression_method, stringsAsFactors = FALSE)
}
tpm_sample_qc <- do.call(rbind, tpm_qc_rows)
tpm_summary <- do.call(rbind, tpm_summary_rows)
rownames(tpm_sample_qc) <- NULL
rownames(tpm_summary) <- NULL
write_tpm_tsv(tpm_sample_qc, "TPM_sample_numeric_qc.tsv")
write_tpm_tsv(tpm_summary, "TPM_cohort_summary.tsv")
conversion_settings <- list(count_cohorts = c("GSE181946", "GSE202069"), length_source = "original_annotation_Length", 
    denominator_gene_scope = "all_input_rows_before_protein_coding_filter", provided_TPM_cohort = "PRJEB34724", 
    provided_TPM_modified = FALSE, gene_filter_applied = FALSE, batch_correction_applied = FALSE, annotation_bundle_file = ANNOTATION_BUNDLE_FILE, 
    annotation_bundle_md5 = unname(tools::md5sum(ANNOTATION_BUNDLE_FILE)), raw_bundle_file = RAW_BUNDLE_FILE, raw_bundle_md5 = raw_md5)
saveRDS(list(ici_tpm = ici_tpm, metadata = raw_bundle$metadata, tpm_summary = tpm_summary, tpm_sample_qc = tpm_sample_qc, 
    conversion_settings = conversion_settings), artifact("stage02_ici_TPM_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(TPM_RUN_DIR, "sessionInfo.txt"))
cat("\n===== TPM \u8F6C\u6362\u9A8C\u8BC1\u6C47\u603B =====\n")
print(tpm_summary, row.names = FALSE, digits = 10)
cat("\n\u6837\u672C\u603B\u6570\uFF1A", nrow(raw_bundle$metadata), "\n")
cat("\u672C\u6B21\u8F93\u51FA\u76EE\u5F55\uFF1A", TPM_RUN_DIR, "\n")
PROJECT_DIR <- CFG$root
TPM_BUNDLE_FILE <- artifact("stage02_ici_TPM_bundle.rds")
REFERENCE_FILE <- artifact("stage01_reference_bundle.rds")
if (!all(file.exists(c(TPM_BUNDLE_FILE, REFERENCE_FILE)))) {
    stop("TPM \u7ED3\u679C\u6216 GENCODE \u53C2\u8003 RDS \u4E0D\u5B58\u5728\u3002")
}
tpm_bundle <- readRDS(TPM_BUNDLE_FILE)
reference_bundle <- readRDS(REFERENCE_FILE)
ici_tpm <- tpm_bundle$ici_tpm
gene_ref <- reference_bundle$gene_reference
if (anyDuplicated(gene_ref$ensembl_gene_id)) {
    stop("GENCODE \u53C2\u8003\u8868\u4E2D\u5B58\u5728\u91CD\u590D\u7684\u53BB\u7248\u672C\u53F7 Ensembl ID\u3002")
}
GENE_MAP_DIR <- out_dir("03_ici_mapping")
if (dir.exists(GENE_MAP_DIR) || !dir.create(GENE_MAP_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_gene_map <- function(x, filename) {
    write.table(x, file.path(GENE_MAP_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", 
        fileEncoding = "UTF-8")
}
reference_rows_by_symbol <- split(seq_len(nrow(gene_ref)), gene_ref$gene_symbol)
gencode_mappings <- list()
for (cohort in names(ici_tpm)) {
    obj <- ici_tpm[[cohort]]
    ann <- obj$gene_annotation
    x <- obj$expression_tpm
    stopifnot(identical(rownames(x), ann$original_gene_id))
    source_id <- ann$source_ensembl_id
    source_symbol <- ann$source_symbol
    id_present <- !is.na(source_id) & nzchar(source_id)
    id_index <- match(source_id, gene_ref$ensembl_gene_id)
    id_index[!id_present] <- NA_integer_
    id_found <- !is.na(id_index)
    symbol_matches <- lapply(source_symbol, function(s) {
        if (is.na(s) || !nzchar(s)) {
            integer()
        }
        else {
            reference_rows_by_symbol[[s]]
        }
    })
    n_symbol_matches <- lengths(symbol_matches)
    symbol_index <- vapply(symbol_matches, function(rows) {
        if (length(rows) == 1L) {
            as.integer(rows)
        }
        else {
            NA_integer_
        }
    }, integer(1))
    fallback <- !id_present & n_symbol_matches == 1L
    selected_index <- id_index
    selected_index[fallback] <- symbol_index[fallback]
    mapping_status <- rep("no_id_symbol_unmatched", nrow(ann))
    mapping_status[id_present & !id_found] <- "id_not_found_review"
    mapping_status[!id_present & n_symbol_matches > 1L] <- "no_id_symbol_ambiguous"
    mapping_status[fallback] <- "symbol_exact_unique_candidate"
    mapping_status[id_found] <- "id_exact"
    symbol_consistent_with_id <- vapply(seq_len(nrow(ann)), function(i) {
        id_found[i] && id_index[i] %in% symbol_matches[[i]]
    }, logical(1))
    id_symbol_conflict <- id_found & n_symbol_matches > 0L & !symbol_consistent_with_id
    mapping_status[id_found & n_symbol_matches == 0L] <- "id_exact_symbol_changed"
    mapping_status[id_symbol_conflict] <- "id_symbol_conflict_review"
    mapping <- data.frame(original_gene_id = ann$original_gene_id, source_symbol = source_symbol, source_ensembl_id = source_id, 
        n_symbol_reference_matches = n_symbol_matches, mapping_status = mapping_status, target_ensembl_id = gene_ref$ensembl_gene_id[selected_index], 
        target_symbol = gene_ref$gene_symbol[selected_index], target_gene_type = gene_ref$gene_type[selected_index], 
        stringsAsFactors = FALSE)
    accepted_status <- c("id_exact", "id_exact_symbol_changed", "symbol_exact_unique_candidate")
    mapping$protein_coding_candidate <- mapping$mapping_status %in% accepted_status & !is.na(mapping$target_gene_type) & 
        mapping$target_gene_type == "protein_coding"
    candidate_index <- which(mapping$protein_coding_candidate)
    candidate_targets <- mapping$target_ensembl_id[candidate_index]
    duplicated_target <- duplicated(candidate_targets) | duplicated(candidate_targets, fromLast = TRUE)
    mapping$target_has_multiple_source_rows <- FALSE
    mapping$target_has_multiple_source_rows[candidate_index] <- duplicated_target
    gencode_mappings[[cohort]] <- mapping
    write_gene_map(mapping, paste0(cohort, "_GENCODE_v50_mapping.tsv"))
    review_status <- c("id_not_found_review", "no_id_symbol_ambiguous", "id_symbol_conflict_review")
    review_rows <- mapping[mapping$mapping_status %in% review_status | mapping$target_has_multiple_source_rows, 
        , drop = FALSE]
    review_rows$total_TPM_across_samples <- rowSums(x)[match(review_rows$original_gene_id, rownames(x))]
    review_rows <- review_rows[order(-review_rows$total_TPM_across_samples), , drop = FALSE]
    write_gene_map(review_rows, paste0(cohort, "_GENCODE_mapping_review.tsv"))
}
status_tables <- list()
candidate_summaries <- list()
sample_coverage_tables <- list()
for (cohort in names(gencode_mappings)) {
    mapping <- gencode_mappings[[cohort]]
    x <- ici_tpm[[cohort]]$expression_tpm
    candidate <- mapping$protein_coding_candidate
    status_table <- as.data.frame(table(mapping$mapping_status), stringsAsFactors = FALSE)
    names(status_table) <- c("mapping_status", "n_source_rows")
    status_table$cohort <- cohort
    status_tables[[cohort]] <- status_table[c("cohort", "mapping_status", "n_source_rows")]
    duplicate_targets <- unique(mapping$target_ensembl_id[candidate & mapping$target_has_multiple_source_rows])
    candidate_summaries[[cohort]] <- data.frame(cohort = cohort, n_source_rows = nrow(mapping), n_pc_candidate_rows = sum(candidate), 
        n_unique_pc_target_ids = length(unique(mapping$target_ensembl_id[candidate])), n_pc_symbol_only_candidate_rows = sum(candidate & 
            mapping$mapping_status == "symbol_exact_unique_candidate"), n_pc_rows_with_shared_target = sum(candidate & 
            mapping$target_has_multiple_source_rows), n_shared_target_ids = length(duplicate_targets), stringsAsFactors = FALSE)
    sample_coverage_tables[[cohort]] <- data.frame(cohort = cohort, sample_id = colnames(x), pc_candidate_TPM_fraction = colSums(x[candidate, 
        , drop = FALSE])/colSums(x), stringsAsFactors = FALSE)
    cat("\n===== ", cohort, "\uFF1A\u4F18\u5148\u6838\u67E5\u8BB0\u5F55\uFF0C\u6700\u591A 10 \u884C =====\n", sep = "")
    review_file <- file.path(GENE_MAP_DIR, paste0(cohort, "_GENCODE_mapping_review.tsv"))
    review <- read.delim(review_file, check.names = FALSE, stringsAsFactors = FALSE)
    print(head(review, 10L), row.names = FALSE)
}
mapping_status_summary <- do.call(rbind, status_tables)
pc_candidate_summary <- do.call(rbind, candidate_summaries)
pc_sample_coverage <- do.call(rbind, sample_coverage_tables)
rownames(mapping_status_summary) <- NULL
rownames(pc_candidate_summary) <- NULL
rownames(pc_sample_coverage) <- NULL
write_gene_map(mapping_status_summary, "GENCODE_mapping_status_summary.tsv")
write_gene_map(pc_candidate_summary, "protein_coding_candidate_summary.tsv")
write_gene_map(pc_sample_coverage, "protein_coding_candidate_sample_coverage.tsv")
saveRDS(list(tpm_bundle_file = TPM_BUNDLE_FILE, tpm_bundle_md5 = unname(tools::md5sum(TPM_BUNDLE_FILE)), reference_file = REFERENCE_FILE, 
    reference_md5 = unname(tools::md5sum(REFERENCE_FILE)), gencode_mappings = gencode_mappings, mapping_status_summary = mapping_status_summary, 
    pc_candidate_summary = pc_candidate_summary, pc_sample_coverage = pc_sample_coverage, mapping_policy = list(primary_key = "Ensembl_gene_ID", 
        missing_ID_fallback = "unique_exact_symbol_candidate", nonempty_unmatched_ID_fallback = "none", id_symbol_conflict = "review_before_use", 
        protein_coding_rule = "gene_type == \"protein_coding\"", duplicate_target_action = "report_without_aggregation")), 
    artifact("stage02_ici_GENCODE_mapping_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(GENE_MAP_DIR, "sessionInfo.txt"))
cat("\n===== GENCODE \u6620\u5C04\u72B6\u6001 =====\n")
print(mapping_status_summary, row.names = FALSE)
cat("\n===== \u86CB\u767D\u7F16\u7801\u5019\u9009\u6C47\u603B =====\n")
print(pc_candidate_summary, row.names = FALSE)
cat("\n\u672C\u6B21\u8F93\u51FA\u76EE\u5F55\uFF1A", GENE_MAP_DIR, "\n")
PROJECT_DIR <- CFG$root
GENE_MAPPING_FILE <- artifact("stage02_ici_GENCODE_mapping_bundle.rds")
if (!file.exists(GENE_MAPPING_FILE)) {
    stop("\u627E\u4E0D\u5230\u6700\u65B0\u57FA\u56E0\u6620\u5C04\u7ED3\u679C\u3002")
}
mapping_bundle <- readRDS(GENE_MAPPING_FILE)
TPM_BUNDLE_FILE <- mapping_bundle$tpm_bundle_file
REFERENCE_FILE <- mapping_bundle$reference_file
check_md5 <- function(path, expected_md5) {
    actual_md5 <- unname(tools::md5sum(path))
    if (length(actual_md5) != 1L || is.na(actual_md5) || actual_md5 != expected_md5) {
        stop("\u4E0A\u6E38\u6587\u4EF6\u7F3A\u5931\u6216\u6307\u7EB9\u4E0D\u4E00\u81F4\uFF1A", path)
    }
}
check_md5(TPM_BUNDLE_FILE, mapping_bundle$tpm_bundle_md5)
check_md5(REFERENCE_FILE, mapping_bundle$reference_md5)
tpm_bundle <- readRDS(TPM_BUNDLE_FILE)
PC_RUN_DIR <- out_dir("03_ici_pc")
if (dir.exists(PC_RUN_DIR) || !dir.create(PC_RUN_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_pc_tsv <- function(x, filename) {
    write.table(x, file.path(PC_RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
pc_candidates <- list()
summary_rows <- list()
sample_qc_rows <- list()
for (cohort in names(tpm_bundle$ici_tpm)) {
    obj <- tpm_bundle$ici_tpm[[cohort]]
    x <- obj$expression_tpm
    mapping <- mapping_bundle$gencode_mappings[[cohort]]
    stopifnot(identical(rownames(x), mapping$original_gene_id), identical(colnames(x), obj$metadata$sample_id))
    candidate <- mapping$protein_coding_candidate
    duplicate_source <- candidate & mapping$target_has_multiple_source_rows
    keep <- candidate & !duplicate_source
    gene_metadata <- mapping[keep, , drop = FALSE]
    pc_tpm <- x[keep, , drop = FALSE]
    if (nrow(pc_tpm) == 0L) {
        stop(cohort, "\uFF1A\u6CA1\u6709\u4FDD\u7559\u7684\u5019\u9009\u86CB\u767D\u7F16\u7801\u57FA\u56E0\u3002")
    }
    if (anyNA(gene_metadata$target_ensembl_id) || anyDuplicated(gene_metadata$target_ensembl_id)) {
        stop(cohort, "\uFF1A\u5019\u9009\u77E9\u9635\u4ECD\u5B58\u5728\u7F3A\u5931\u6216\u91CD\u590D\u76EE\u6807 ID\u3002")
    }
    rownames(pc_tpm) <- gene_metadata$target_ensembl_id
    gene_order <- order(gene_metadata$target_ensembl_id, method = "radix")
    pc_tpm <- pc_tpm[gene_order, , drop = FALSE]
    gene_metadata <- gene_metadata[gene_order, , drop = FALSE]
    rownames(gene_metadata) <- NULL
    pc_log_tpm <- log2(pc_tpm + 1)
    stopifnot(identical(rownames(pc_tpm), gene_metadata$target_ensembl_id), identical(dimnames(pc_tpm), dimnames(pc_log_tpm)), 
        identical(colnames(pc_tpm), obj$metadata$sample_id), all(is.finite(pc_log_tpm)), all(gene_metadata$target_gene_type == 
            "protein_coding"))
    decision_table <- mapping
    decision_table$candidate_matrix_action <- ifelse(keep, "retain", ifelse(duplicate_source, "hold_multiple_source_rows", 
        "not_eligible_under_current_mapping"))
    write_pc_tsv(decision_table, paste0(cohort, "_source_row_decisions.tsv"))
    write_pc_tsv(decision_table[duplicate_source, , drop = FALSE], paste0(cohort, "_held_multiple_source_rows.tsv"))
    write_pc_tsv(gene_metadata, paste0(cohort, "_pc_candidate_gene_metadata.tsv"))
    pc_candidates[[cohort]] <- list(expression_tpm = pc_tpm, expression_log2tpm = pc_log_tpm, gene_metadata = gene_metadata, 
        metadata = obj$metadata, expression_method = obj$expression_method)
    saveRDS(pc_candidates[[cohort]], file.path(PC_RUN_DIR, paste0(cohort, "_pc_candidate.rds")))
    sample_qc_rows[[cohort]] <- data.frame(cohort = cohort, sample_id = colnames(pc_tpm), retained_TPM_sum = colSums(pc_tpm), 
        retained_TPM_fraction = colSums(pc_tpm)/colSums(x), n_detected_pc_genes = colSums(pc_tpm > 0), pc_zero_fraction = colMeans(pc_tpm == 
            0), all_zero = colSums(pc_tpm > 0) == 0L, stringsAsFactors = FALSE)
    summary_rows[[cohort]] <- data.frame(cohort = cohort, n_samples = ncol(pc_tpm), n_pc_candidate_rows_before_hold = sum(candidate), 
        n_held_source_rows = sum(duplicate_source), n_retained_pc_genes = nrow(pc_tpm), n_retained_symbol_only_candidates = sum(gene_metadata$mapping_status == 
            "symbol_exact_unique_candidate"), n_genes_zero_in_all_samples = sum(rowSums(pc_tpm > 0) == 0L), stringsAsFactors = FALSE)
}
pc_summary <- do.call(rbind, summary_rows)
pc_sample_qc <- do.call(rbind, sample_qc_rows)
rownames(pc_summary) <- NULL
rownames(pc_sample_qc) <- NULL
write_pc_tsv(pc_summary, "pc_candidate_summary.tsv")
write_pc_tsv(pc_sample_qc, "pc_candidate_sample_qc.tsv")
saveRDS(list(pc_candidates = pc_candidates, metadata = tpm_bundle$metadata, pc_summary = pc_summary, pc_sample_qc = pc_sample_qc, 
    provenance = list(gene_mapping_file = GENE_MAPPING_FILE, gene_mapping_md5 = unname(tools::md5sum(GENE_MAPPING_FILE)), 
        tpm_bundle_file = TPM_BUNDLE_FILE, tpm_bundle_md5 = mapping_bundle$tpm_bundle_md5, duplicate_target_policy = "hold_all_candidate_source_rows", 
        TPM_renormalized = FALSE, cohort_gene_intersection_applied = FALSE, expression_prevalence_filter_applied = FALSE, 
        final_model_gene_list_locked = FALSE)), artifact("stage02_ici_pc_candidate_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(PC_RUN_DIR, "sessionInfo.txt"))
cat("\n===== \u4FDD\u5B88\u7248\u86CB\u767D\u7F16\u7801\u5019\u9009\u77E9\u9635 =====\n")
print(pc_summary, row.names = FALSE)
cat("\n===== \u6837\u672C\u8D28\u63A7\u8303\u56F4 =====\n")
for (cohort in names(pc_candidates)) {
    q <- pc_sample_qc[pc_sample_qc$cohort == cohort, , drop = FALSE]
    cat(cohort, " | \u4FDD\u7559 TPM \u6BD4\u4F8B\uFF1A", paste(round(range(q$retained_TPM_fraction), 4), collapse = " ~ "), " | \u68C0\u51FA\u86CB\u767D\u7F16\u7801\u57FA\u56E0\u6570\uFF1A", 
        paste(range(q$n_detected_pc_genes), collapse = " ~ "), "\n", sep = "")
}
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", PC_RUN_DIR, "\n")
if (any(pc_sample_qc$all_zero)) {
    warning("\u5019\u9009\u86CB\u767D\u7F16\u7801\u77E9\u9635\u5B58\u5728\u5168\u96F6\u6837\u672C\uFF0C\u8BF7\u5148\u68C0\u67E5\u3002")
}

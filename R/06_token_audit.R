PROJECT_DIR <- CFG$root
TOKEN37_FILE <- file.path(CFG$root, "HCC_37token.xlsx")
REFERENCE_FILE <- artifact("stage01_reference_bundle.rds")
INTERSECTION_FILE <- artifact("stage03_intersection_candidate_bundle.rds")
ALIGNED_FILE <- artifact("stage03_ici_aligned_with_mask_bundle.rds")
required_files <- c(TOKEN37_FILE, REFERENCE_FILE, INTERSECTION_FILE, ALIGNED_FILE)
if (!all(file.exists(required_files))) {
    stop("\u4EE5\u4E0B\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A\n", paste(required_files[!file.exists(required_files)], collapse = "\n"))
}
if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("\u8BF7\u5148\u5B89\u88C5 readxl\u3002")
}
TOKEN37_DIR <- out_dir("06_token_audit")
if (dir.exists(TOKEN37_DIR) || !dir.create(TOKEN37_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u72EC\u7ACB\u8F93\u51FA\u76EE\u5F55\u3002")
}
write_37_tsv <- function(x, filename) {
    write.table(x, file.path(TOKEN37_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
raw37 <- as.data.frame(readxl::read_excel(TOKEN37_FILE, sheet = "Sheet1", col_types = "text", .name_repair = "check_unique"), 
    stringsAsFactors = FALSE)
required_columns <- c("\u6982\u5FF5ID", "\u4E00\u7EA7\u540D\u79F0", "\u4E8C\u7EA7\u540D\u79F0\uFF08\u4E2D\u6587\uFF09", "\u57FA\u56E0\u96C6\u51C6\u786E\u540D\u79F0", "\u6210\u5458\u6570\u91CF", 
    "\u57FA\u56E0\u6210\u5458\uFF08HGNC Symbol\uFF0C\u7528:\u5206\u9694\uFF09")
if (!all(required_columns %in% names(raw37))) {
    stop("37-token\u8868\u7F3A\u5C11\u5FC5\u8981\u5217\u3002")
}
module_lookup <- setNames(paste0("M", seq_len(8L)), c(
    "\u514D\u75AB\u4E0E\u708E\u75C7",
    "\u4EE3\u8C22\u91CD\u7F16\u7A0B",
    "\u81F4\u764C\u4FE1\u53F7\u4E0E\u589E\u6B96",
    "\u57FA\u56E0\u7EC4\u7A33\u5B9A\u6027\u4E0E\u7EC6\u80DE\u547D\u8FD0",
    "\u57FA\u8D28\u91CD\u5851\u4E0E\u4FB5\u88AD",
    "\u8840\u7BA1\u751F\u6210\u4E0E\u7F3A\u6C27",
    "\u7EC6\u80DE\u5206\u5316\u4E0E\u5E72\u6027",
    "\u7EC6\u80DE\u5E94\u6FC0\u4E0E\u8C03\u63A7"))
tokens37 <- data.frame(excel_row = seq_len(nrow(raw37)) + 1L, token_name = trimws(raw37[["\u57FA\u56E0\u96C6\u51C6\u786E\u540D\u79F0"]]), 
    concept_id = trimws(raw37[["\u6982\u5FF5ID"]]), concept_name = trimws(raw37[["\u4E8C\u7EA7\u540D\u79F0\uFF08\u4E2D\u6587\uFF09"]]), module_name = trimws(raw37[["\u4E00\u7EA7\u540D\u79F0"]]), 
    declared_n = suppressWarnings(as.numeric(raw37[["\u6210\u5458\u6570\u91CF"]])), members_raw = raw37[["\u57FA\u56E0\u6210\u5458\uFF08HGNC Symbol\uFF0C\u7528:\u5206\u9694\uFF09"]], 
    stringsAsFactors = FALSE)
tokens37$module_id <- unname(module_lookup[tokens37$module_name])
essential <- c("token_name", "concept_id", "concept_name", "module_id", "module_name", "members_raw")
if (any(vapply(tokens37[essential], function(x) any(is.na(x) | !nzchar(trimws(x))), logical(1)))) {
    stop("\u5FC5\u8981\u5B57\u6BB5\u5B58\u5728\u7A7A\u503C\u6216\u672A\u8BC6\u522B\u6A21\u5757\u3002")
}
if (nrow(tokens37) != 37L || anyDuplicated(tokens37$token_name)) {
    stop("Token\u6570\u4E0D\u662F37\uFF0C\u6216\u540D\u79F0\u91CD\u590D\u3002")
}
concept_mapping <- unique(tokens37[c("concept_id", "concept_name", "module_id")])
if (anyDuplicated(concept_mapping$concept_id)) {
    stop("\u540C\u4E00Concept ID\u5BF9\u5E94\u591A\u4E2A\u540D\u79F0\u6216\u6A21\u5757\uFF0C\u9700\u8981\u68C0\u67E5\u3002")
}
split_members <- function(x) {
    parts <- strsplit(paste0(x, ":__END__"), split = ":", fixed = TRUE)[[1L]]
    trimws(parts[-length(parts)])
}
member_lists37 <- lapply(tokens37$members_raw, split_members)
tokens37$n_empty <- vapply(member_lists37, function(x) sum(!nzchar(x)), integer(1))
tokens37$n_unique <- vapply(member_lists37, function(x) length(unique(x[nzchar(x)])), integer(1))
tokens37$n_duplicates <- lengths(member_lists37) - tokens37$n_empty - tokens37$n_unique
tokens37$count_ok <- !is.na(tokens37$declared_n) & tokens37$declared_n == tokens37$n_unique
write_37_tsv(tokens37, "token37_structure_audit.tsv")
if (any(tokens37$n_empty > 0L | tokens37$n_duplicates > 0L | !tokens37$count_ok)) {
    stop("\u5B58\u5728\u7A7A\u6210\u5458\u3001\u91CD\u590D\u6210\u5458\u6216\u6210\u5458\u6570\u91CF\u95EE\u9898\uFF0C\u8BF7\u68C0\u67E5\u5BA1\u8BA1\u8868\u3002")
}
members37 <- do.call(rbind, lapply(seq_len(nrow(tokens37)), function(i) {
    data.frame(token_name = tokens37$token_name[i], concept_id = tokens37$concept_id[i], module_id = tokens37$module_id[i], 
        gene_symbol = member_lists37[[i]], stringsAsFactors = FALSE)
}))
rownames(members37) <- NULL
module_summary37 <- do.call(rbind, lapply(paste0("M", 1:8), function(m) {
    x <- tokens37[tokens37$module_id == m, , drop = FALSE]
    data.frame(module_id = m, module_name = unique(x$module_name), n_tokens = nrow(x), n_concepts = length(unique(x$concept_id)), 
        stringsAsFactors = FALSE)
}))
write_37_tsv(module_summary37, "token37_module_summary.tsv")
write_37_tsv(concept_mapping, "token37_concept_mapping.tsv")
cat("\n===== 37-token\u5C42\u7EA7 =====\n")
print(module_summary37, row.names = FALSE)
cat("Token\u6570\uFF1A", nrow(tokens37), "\n")
cat("Concept\u6570\uFF1A", nrow(concept_mapping), "\n")
cat("\u552F\u4E00\u539F\u59CB\u57FA\u56E0\u540D\u6570\uFF1A", length(unique(members37$gene_symbol)), "\n")
reference_bundle <- readRDS(REFERENCE_FILE)
intersection_bundle <- readRDS(INTERSECTION_FILE)
aligned_bundle <- readRDS(ALIGNED_FILE)
gene_ref <- reference_bundle$gene_reference
pc_ref <- reference_bundle$protein_coding_reference
fixed_genes <- intersection_bundle$intersection_genes
stopifnot(identical(fixed_genes, aligned_bundle$fixed_genes), identical(fixed_genes, rownames(intersection_bundle$ssl_expression_tpm)))
current_intersection_md5 <- unname(tools::md5sum(INTERSECTION_FILE))
if (current_intersection_md5 != aligned_bundle$provenance$intersection_md5) {
    stop("ICI\u5BF9\u9F50\u5BF9\u8C61\u4F7F\u7528\u7684\u4EA4\u96C6\u7248\u672C\u4E0D\u4E00\u81F4\u3002")
}
pc_ids_by_symbol <- split(pc_ref$ensembl_gene_id, pc_ref$gene_symbol)
types_by_symbol <- split(gene_ref$gene_type, gene_ref$gene_symbol)
symbols37 <- sort(unique(members37$gene_symbol))
symbol_audit37 <- data.frame(gene_symbol = symbols37, n_pc_ids = vapply(symbols37, function(s) length(unique(pc_ids_by_symbol[[s]])), 
    integer(1)), matched_gene_types = vapply(symbols37, function(s) paste(sort(unique(types_by_symbol[[s]])), collapse = ";"), 
    character(1)), stringsAsFactors = FALSE)
symbol_audit37$target_ensembl_id <- vapply(symbols37, function(s) {
    ids <- unique(pc_ids_by_symbol[[s]])
    if (length(ids) == 1L) 
        ids
    else NA_character_
}, character(1))
symbol_audit37$status <- ifelse(symbol_audit37$n_pc_ids == 1L, "protein_coding_unique", ifelse(symbol_audit37$n_pc_ids > 
    1L, "protein_coding_ambiguous", ifelse(nzchar(symbol_audit37$matched_gene_types), "matched_other_GENCODE_type", 
    "not_found_exact_symbol")))
symbol_audit37$in_fixed_background <- !is.na(symbol_audit37$target_ensembl_id) & symbol_audit37$target_ensembl_id %in% 
    fixed_genes
symbol_index <- match(members37$gene_symbol, symbol_audit37$gene_symbol)
members37$status <- symbol_audit37$status[symbol_index]
members37$target_ensembl_id <- symbol_audit37$target_ensembl_id[symbol_index]
members37$in_fixed_background <- symbol_audit37$in_fixed_background[symbol_index]
token_gene_indices37 <- list()
coverage_rows37 <- list()
ici_coverage_rows37 <- list()
for (i in seq_len(nrow(tokens37))) {
    token_name <- tokens37$token_name[i]
    m <- members37[members37$token_name == token_name, , drop = FALSE]
    target_ids <- sort(unique(m$target_ensembl_id[m$in_fixed_background]))
    indices <- match(target_ids, fixed_genes)
    stopifnot(!anyNA(indices))
    token_gene_indices37[[token_name]] <- indices
    n_original <- nrow(m)
    n_unique_pc <- sum(m$status == "protein_coding_unique")
    n_in_background <- sum(m$in_fixed_background)
    coverage_rows37[[i]] <- data.frame(token_name = token_name, concept_id = tokens37$concept_id[i], module_id = tokens37$module_id[i], 
        n_original_members = n_original, n_unique_pc_members = n_unique_pc, n_ambiguous_pc_members = sum(m$status == 
            "protein_coding_ambiguous"), n_unmatched_symbols = sum(m$status == "not_found_exact_symbol"), n_members_in_background = n_in_background, 
        original_member_coverage = n_in_background/n_original, unique_pc_member_coverage = if (n_unique_pc > 0L) {
            n_in_background/n_unique_pc
        }
        else {
            NA_real_
        }, stringsAsFactors = FALSE)
    for (cohort in names(aligned_bundle$ici_aligned)) {
        obj <- aligned_bundle$ici_aligned[[cohort]]
        stopifnot(identical(rownames(obj$observed_mask), fixed_genes))
        if (!all(rowSums(obj$observed_mask) %in% c(0L, ncol(obj$observed_mask)))) {
            stop(cohort, "\uFF1A\u5B58\u5728\u60A3\u8005\u7279\u5F02\u7F3A\u6D4B\uFF0C\u9700\u5355\u72EC\u5BA1\u8BA1\u3002")
        }
        observed <- obj$observed_mask[indices, 1L]
        key <- paste(cohort, token_name, sep = "::")
        ici_coverage_rows37[[key]] <- data.frame(cohort = cohort, token_name = token_name, concept_id = tokens37$concept_id[i], 
            module_id = tokens37$module_id[i], n_background_members = length(indices), n_observed_members = sum(observed), 
            n_missing_members = sum(!observed), observed_member_fraction = if (length(indices) > 0L) {
                mean(observed)
            }
            else {
                NA_real_
            }, stringsAsFactors = FALSE)
    }
}
coverage37 <- do.call(rbind, coverage_rows37)
ici_coverage37 <- do.call(rbind, ici_coverage_rows37)
rownames(coverage37) <- NULL
rownames(ici_coverage37) <- NULL
coverage37$review_flag <- coverage37$n_members_in_background == 0L | is.na(coverage37$unique_pc_member_coverage) | 
    coverage37$unique_pc_member_coverage < 0.8 | coverage37$n_ambiguous_pc_members > 0L
write_37_tsv(symbol_audit37, "token37_symbol_annotation.tsv")
write_37_tsv(symbol_audit37[symbol_audit37$status != "protein_coding_unique", , drop = FALSE], "token37_annotation_issues.tsv")
write_37_tsv(members37, "token37_member_mapping.tsv")
write_37_tsv(coverage37, "token37_background_coverage.tsv")
write_37_tsv(ici_coverage37, "token37_ICI_observed_coverage.tsv")
member_union_indices <- sort(unique(unlist(token_gene_indices37, use.names = FALSE)))
source_paths <- c(token37 = TOKEN37_FILE, gencode_reference = REFERENCE_FILE, pretraining_intersection = INTERSECTION_FILE, 
    ici_aligned = ALIGNED_FILE)
source_manifest <- data.frame(source = names(source_paths), path = unname(source_paths), md5 = unname(tools::md5sum(unname(source_paths))), 
    stringsAsFactors = FALSE)
if (anyNA(source_manifest$md5)) {
    stop("\u90E8\u5206\u6E90\u6587\u4EF6\u65E0\u6CD5\u8BA1\u7B97MD5\u3002")
}
write_37_tsv(source_manifest, "token37_source_manifest.tsv")
saveRDS(list(tokens = tokens37, concept_mapping = concept_mapping, module_summary = module_summary37, members = members37, 
    symbol_audit = symbol_audit37, background_coverage = coverage37, ici_observed_coverage = ici_coverage37, fixed_background_genes = fixed_genes, 
    token_gene_indices = token_gene_indices37, member_union_indices = member_union_indices, member_union_gene_ids = fixed_genes[member_union_indices], 
    source_manifest = source_manifest, policy = list(response_labels_used = FALSE, automatic_token_exclusion = FALSE, 
        automatic_alias_replacement = FALSE, ambiguous_symbol_resolution = "pending", final_panel_locked = FALSE, 
        note = "Coverage audit in existing fixed pretraining gene universe")), artifact("stage09_token37_audit_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(TOKEN37_DIR, "sessionInfo.txt"))
cat("\n===== \u552F\u4E00\u6210\u5458\u7684\u6CE8\u91CA\u72B6\u6001 =====\n")
print(as.data.frame(table(symbol_audit37$status)), row.names = FALSE)
cat("\n===== \u9700\u91CD\u70B9\u6838\u67E5\u7684token =====\n")
print(coverage37[coverage37$review_flag, , drop = FALSE], row.names = FALSE)
cat("\n\u73B0\u6709\u80CC\u666F\u4E2D\u53EF\u5B9A\u4F4D\u7684token\u6210\u5458\u57FA\u56E0\u5E76\u96C6\uFF1A", length(member_union_indices), "\n")
cat("\u65E0\u80CC\u666F\u6210\u5458\u7684token\u6570\uFF1A", sum(coverage37$n_members_in_background == 0L), "\n")
cat("ICI\u4E2D\u65E0\u89C2\u6D4B\u6210\u5458\u7684cohort-token\u7EC4\u5408\u6570\uFF1A", sum(ici_coverage37$n_observed_members == 0L), "\n")
cat("\n\u8F93\u51FA\u76EE\u5F55\uFF1A", TOKEN37_DIR, "\n")

# Prepare ICI expression inputs and response labels.
PROJECT_DIR <- CFG$root
REFERENCE_FILE <- artifact("stage01_reference_bundle.rds")
TOKEN_FILE <- file.path(CFG$root, "HCC_37token.xlsx")
if (!file.exists(REFERENCE_FILE) || !file.exists(TOKEN_FILE)) {
    stop("\u53C2\u8003 RDS \u6216\u57FA\u56E0\u96C6 Excel \u4E0D\u5B58\u5728\uFF0C\u8BF7\u68C0\u67E5\u8DEF\u5F84\u3002")
}
if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("\u8BF7\u5148\u8FD0\u884C install.packages('readxl')\uFF0C\u518D\u91CD\u65B0\u8FD0\u884C\u672C\u6BB5\u3002")
}
ref_bundle <- readRDS(REFERENCE_FILE)
gene_ref <- ref_bundle$gene_reference
pc_ref <- ref_bundle$protein_coding_reference
old_token_md5 <- ref_bundle$input_manifest$md5[ref_bundle$input_manifest$input_name == "token_table"]
current_token_md5 <- unname(tools::md5sum(TOKEN_FILE))
if (length(old_token_md5) != 1L || is.na(current_token_md5) || current_token_md5 != old_token_md5) {
    stop("Excel \u4E0E\u4E0A\u4E00\u9636\u6BB5\u8BB0\u5F55\u4E0D\u4E00\u81F4\uFF0C\u8BF7\u5148\u786E\u8BA4\u6587\u4EF6\u662F\u5426\u53D1\u751F\u4FEE\u6539\u3002")
}
TOKEN_RUN_DIR <- out_dir("02_token_reference")
if (dir.exists(TOKEN_RUN_DIR) || !dir.create(TOKEN_RUN_DIR, recursive = TRUE)) {
    stop("\u8F93\u51FA\u76EE\u5F55\u5DF2\u5B58\u5728\u6216\u65E0\u6CD5\u521B\u5EFA\uFF0C\u8BF7\u91CD\u65B0\u8FD0\u884C\u672C\u6BB5\u3002")
}
token_raw <- as.data.frame(readxl::read_excel(TOKEN_FILE, sheet = "Sheet1", col_types = "text", .name_repair = "check_unique"), 
    stringsAsFactors = FALSE)
token_raw[["Module ID"]] <- paste0("M", match(token_raw[["\u4E00\u7EA7\u540D\u79F0"]], c("\u514D\u75AB\u4E0E\u708E\u75C7", "\u4EE3\u8C22\u91CD\u7F16\u7A0B", 
    "\u81F4\u764C\u4FE1\u53F7\u4E0E\u589E\u6B96", "\u57FA\u56E0\u7EC4\u7A33\u5B9A\u6027\u4E0E\u7EC6\u80DE\u547D\u8FD0", "\u57FA\u8D28\u91CD\u5851\u4E0E\u4FB5\u88AD", "\u8840\u7BA1\u751F\u6210\u4E0E\u7F3A\u6C27", 
    "\u7EC6\u80DE\u5206\u5316\u4E0E\u5E72\u6027", "\u7EC6\u80DE\u5E94\u6FC0\u4E0E\u8C03\u63A7")))
required_columns <- c("\u57FA\u56E0\u96C6\u51C6\u786E\u540D\u79F0", "Module ID", "\u4E00\u7EA7\u540D\u79F0", "\u6982\u5FF5ID", "\u6210\u5458\u6570\u91CF", "\u57FA\u56E0\u6210\u5458\uFF08HGNC Symbol\uFF0C\u7528:\u5206\u9694\uFF09")
if (!all(required_columns %in% names(token_raw))) {
    stop("Excel \u7F3A\u5C11\u5217\uFF1A", paste(setdiff(required_columns, names(token_raw)), collapse = ", "))
}
save_token_tsv <- function(x, filename) {
    write.table(x, file.path(TOKEN_RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", 
        fileEncoding = "UTF-8")
}
cat("Excel \u884C\u6570\uFF1A", nrow(token_raw), "\n")
cat("\u672C\u90E8\u5206\u8F93\u51FA\u76EE\u5F55\uFF1A", TOKEN_RUN_DIR, "\n")
token_raw[["Module ID"]] <- paste0("M", match(token_raw[["\u4E00\u7EA7\u540D\u79F0"]], c("\u514D\u75AB\u4E0E\u708E\u75C7", "\u4EE3\u8C22\u91CD\u7F16\u7A0B", 
    "\u81F4\u764C\u4FE1\u53F7\u4E0E\u589E\u6B96", "\u57FA\u56E0\u7EC4\u7A33\u5B9A\u6027\u4E0E\u7EC6\u80DE\u547D\u8FD0", "\u57FA\u8D28\u91CD\u5851\u4E0E\u4FB5\u88AD", "\u8840\u7BA1\u751F\u6210\u4E0E\u7F3A\u6C27", 
    "\u7EC6\u80DE\u5206\u5316\u4E0E\u5E72\u6027", "\u7EC6\u80DE\u5E94\u6FC0\u4E0E\u8C03\u63A7")))
tokens <- data.frame(excel_row = seq_len(nrow(token_raw)) + 1L, token_name = trimws(token_raw[["\u57FA\u56E0\u96C6\u51C6\u786E\u540D\u79F0"]]), 
    module_id = trimws(token_raw[["Module ID"]]), module_name = trimws(token_raw[["\u4E00\u7EA7\u540D\u79F0"]]), concept_id = trimws(token_raw[["\u6982\u5FF5ID"]]), 
    declared_n = suppressWarnings(as.numeric(token_raw[["\u6210\u5458\u6570\u91CF"]])), members_raw = token_raw[["\u57FA\u56E0\u6210\u5458\uFF08HGNC Symbol\uFF0C\u7528:\u5206\u9694\uFF09"]], 
    stringsAsFactors = FALSE)
essential <- c("token_name", "module_id", "module_name", "concept_id", "members_raw")
if (any(vapply(tokens[essential], function(x) any(is.na(x) | !nzchar(trimws(x))), logical(1)))) {
    stop("token \u540D\u79F0\u3001\u5C42\u7EA7\u5F52\u5C5E\u6216\u6210\u5458\u5B57\u6BB5\u5B58\u5728\u7A7A\u503C\uFF0C\u8BF7\u68C0\u67E5\u539F\u8868\u3002")
}
if (anyDuplicated(tokens$token_name)) {
    stop("Token Name \u5B58\u5728\u91CD\u590D\uFF0C\u9700\u8981\u5148\u786E\u8BA4\u662F\u5426\u4E3A\u91CD\u590D\u884C\u6216\u540C\u540D\u4E0D\u540C\u57FA\u56E0\u96C6\u3002")
}
concept_parent <- unique(tokens[c("concept_id", "module_id")])
if (anyDuplicated(concept_parent$concept_id)) {
    stop("\u540C\u4E00 Concept ID \u88AB\u5206\u914D\u5230\u591A\u4E2A\u6A21\u5757\uFF0C\u9700\u8981\u5148\u786E\u8BA4\u7ED3\u6784\u3002")
}
module_labels <- unique(tokens[c("module_id", "module_name")])
if (anyDuplicated(module_labels$module_id)) {
    stop("\u540C\u4E00 Module ID \u5BF9\u5E94\u591A\u4E2A\u6A21\u5757\u540D\u79F0\u3002")
}
split_members <- function(x) {
    parts <- strsplit(paste0(x, ":__END__"), split = ":", fixed = TRUE)[[1L]]
    trimws(parts[-length(parts)])
}
member_lists <- lapply(tokens$members_raw, split_members)
tokens$n_member_entries <- lengths(member_lists)
tokens$n_empty_entries <- vapply(member_lists, function(x) sum(!nzchar(x)), integer(1))
tokens$n_unique_members <- vapply(member_lists, function(x) length(unique(x[nzchar(x)])), integer(1))
tokens$n_duplicate_entries <- tokens$n_member_entries - tokens$n_empty_entries - tokens$n_unique_members
tokens$declared_count_ok <- !is.na(tokens$declared_n) & tokens$declared_n == tokens$n_unique_members
members_long <- do.call(rbind, lapply(seq_len(nrow(tokens)), function(i) {
    data.frame(token_name = tokens$token_name[i], concept_id = tokens$concept_id[i], module_id = tokens$module_id[i], 
        member_index = seq_along(member_lists[[i]]), gene_symbol = member_lists[[i]], stringsAsFactors = FALSE)
}))
save_token_tsv(tokens, "token_structure_audit.tsv")
structure_issues <- tokens[!tokens$declared_count_ok | tokens$n_empty_entries > 0L | tokens$n_duplicate_entries > 
    0L, , drop = FALSE]
save_token_tsv(structure_issues, "token_structure_issues.tsv")
cat("Token \u6570\uFF1A", nrow(tokens), "\n")
cat("Concept \u6570\uFF1A", length(unique(tokens$concept_id)), "\n")
cat("Module \u6570\uFF1A", length(unique(tokens$module_id)), "\n")
cat("\u5B58\u5728\u6210\u5458\u6570\u91CF\u6216\u6210\u5458\u683C\u5F0F\u95EE\u9898\u7684 token \u6570\uFF1A", nrow(structure_issues), "\n")
if (nrow(structure_issues) > 0L) {
    print(structure_issues[c("token_name", "declared_n", "n_unique_members", "n_empty_entries", "n_duplicate_entries")], 
        row.names = FALSE)
}
all_member_symbols <- sort(unique(members_long$gene_symbol[nzchar(members_long$gene_symbol)]))
pc_ids_by_symbol <- split(pc_ref$ensembl_gene_id, pc_ref$gene_symbol)
types_by_symbol <- split(gene_ref$gene_type, gene_ref$gene_symbol)
symbol_audit <- data.frame(gene_symbol = all_member_symbols, n_pc_ids = vapply(all_member_symbols, function(s) length(unique(pc_ids_by_symbol[[s]])), 
    integer(1)), matched_gene_types = vapply(all_member_symbols, function(s) paste(sort(unique(types_by_symbol[[s]])), 
    collapse = ";"), character(1)), pc_ensembl_ids = vapply(all_member_symbols, function(s) paste(sort(unique(pc_ids_by_symbol[[s]])), 
    collapse = ";"), character(1)), stringsAsFactors = FALSE)
symbol_audit$status <- ifelse(symbol_audit$n_pc_ids == 1L, "protein_coding_unique", ifelse(symbol_audit$n_pc_ids > 
    1L, "protein_coding_ambiguous", ifelse(nzchar(symbol_audit$matched_gene_types), "matched_non_protein_coding", 
    "not_found_exact_symbol")))
idx <- match(members_long$gene_symbol, symbol_audit$gene_symbol)
members_long$n_pc_ids <- symbol_audit$n_pc_ids[idx]
members_long$status <- symbol_audit$status[idx]
members_long$status[!nzchar(members_long$gene_symbol)] <- "empty_member"
members_unique <- unique(members_long[nzchar(members_long$gene_symbol), c("token_name", "concept_id", "module_id", 
    "gene_symbol", "n_pc_ids", "status")])
token_coverage <- do.call(rbind, lapply(seq_len(nrow(tokens)), function(i) {
    x <- members_unique[members_unique$token_name == tokens$token_name[i], , drop = FALSE]
    n_total <- nrow(x)
    n_pc <- sum(x$n_pc_ids > 0L)
    data.frame(token_name = tokens$token_name[i], concept_id = tokens$concept_id[i], module_id = tokens$module_id[i], 
        n_members = n_total, n_pc_matched = n_pc, n_pc_ambiguous = sum(x$n_pc_ids > 1L), n_non_pc = sum(x$status == 
            "matched_non_protein_coding"), n_not_found = sum(x$status == "not_found_exact_symbol"), pc_coverage = if (n_total > 
            0L) 
            n_pc/n_total
        else NA_real_, stringsAsFactors = FALSE)
}))
pc_candidate_sets <- setNames(lapply(tokens$token_name, function(nm) {
    sort(unique(members_unique$gene_symbol[members_unique$token_name == nm & members_unique$n_pc_ids > 0L]))
}), tokens$token_name)
save_token_tsv(symbol_audit, "member_symbol_annotation.tsv")
save_token_tsv(members_long, "token_gene_members_audited.tsv")
save_token_tsv(token_coverage, "token_protein_coding_coverage.tsv")
annotation_issues <- symbol_audit[symbol_audit$status != "protein_coding_unique", , drop = FALSE]
save_token_tsv(annotation_issues, "member_annotation_issues.tsv")
cat("\n===== \u552F\u4E00\u6210\u5458\u7684\u6CE8\u91CA\u72B6\u6001 =====\n")
print(as.data.frame(table(symbol_audit$status)), row.names = FALSE)
cat("\n===== \u86CB\u767D\u7F16\u7801\u8986\u76D6\u4E0D\u5B8C\u6574\u6216\u5B58\u5728\u6B67\u4E49\u7684 token =====\n")
print(token_coverage[token_coverage$n_pc_matched < token_coverage$n_members | token_coverage$n_pc_ambiguous > 0L, 
    , drop = FALSE], row.names = FALSE)
cat("\n===== \u9700\u8981\u6838\u67E5\u7684\u6210\u5458\uFF0C\u6700\u591A\u663E\u793A 50 \u4E2A =====\n")
print(head(annotation_issues, 50L), row.names = FALSE)
if (length(pc_candidate_sets) < 2L) {
    stop("\u4E0D\u8DB3\u4E24\u4E2A token\uFF0C\u65E0\u6CD5\u8FDB\u884C\u6210\u5BF9\u91CD\u53E0\u5BA1\u8BA1\u3002")
}
pair_indices <- combn(seq_along(pc_candidate_sets), 2L)
overlap_audit <- do.call(rbind, lapply(seq_len(ncol(pair_indices)), function(k) {
    i <- pair_indices[1L, k]
    j <- pair_indices[2L, k]
    a <- pc_candidate_sets[[i]]
    b <- pc_candidate_sets[[j]]
    n_intersection <- length(intersect(a, b))
    n_union <- length(union(a, b))
    data.frame(token_a = names(pc_candidate_sets)[i], token_b = names(pc_candidate_sets)[j], n_a = length(a), n_b = length(b), 
        n_intersection = n_intersection, jaccard = if (n_union > 0L) {
            n_intersection/n_union
        }
        else {
            NA_real_
        }, identical_nonempty_sets = length(a) > 0L && length(b) > 0L && setequal(a, b), stringsAsFactors = FALSE)
}))
high_overlap <- overlap_audit[!is.na(overlap_audit$jaccard) & overlap_audit$jaccard >= 0.8, , drop = FALSE]
module_summary <- do.call(rbind, lapply(sort(unique(tokens$module_id)), function(m) {
    x <- tokens[tokens$module_id == m, , drop = FALSE]
    data.frame(module_id = m, module_name = unique(x$module_name), n_tokens = nrow(x), n_concepts = length(unique(x$concept_id)), 
        stringsAsFactors = FALSE)
}))
save_token_tsv(overlap_audit, "token_pairwise_overlap.tsv")
save_token_tsv(high_overlap, "token_high_overlap.tsv")
save_token_tsv(module_summary, "module_structure_summary.tsv")
saveRDS(list(reference_file = REFERENCE_FILE, reference_md5 = unname(tools::md5sum(REFERENCE_FILE)), token_file = TOKEN_FILE, 
    token_md5 = current_token_md5, token_raw = token_raw, tokens = tokens, symbol_audit = symbol_audit, members_long = members_long, 
    token_coverage = token_coverage, pc_candidate_sets = pc_candidate_sets, overlap_audit = overlap_audit, module_summary = module_summary), 
    artifact("stage01_token_audit_bundle.rds"))
writeLines(capture.output(sessionInfo()), file.path(TOKEN_RUN_DIR, "sessionInfo.txt"))
cat("\n===== \u6A21\u5757\u7ED3\u6784 =====\n")
print(module_summary, row.names = FALSE)
cat("\n===== Jaccard >= 0.8 \u7684\u57FA\u56E0\u96C6\u5BF9 =====\n")
print(high_overlap, row.names = FALSE)
cat("\n\u552F\u4E00\u539F\u59CB\u6210\u5458\u6570\uFF1A", nrow(symbol_audit), "\n")
cat("\u81F3\u5C11\u5339\u914D\u4E00\u4E2A protein_coding ID \u7684\u552F\u4E00\u6210\u5458\u6570\uFF1A", sum(symbol_audit$n_pc_ids > 0L), "\n")
cat("\u672C\u90E8\u5206\u7ED3\u679C\u76EE\u5F55\uFF1A", TOKEN_RUN_DIR, "\n")

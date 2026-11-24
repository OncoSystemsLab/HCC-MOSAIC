# Build the GENCODE reference and audit the token panel.
PROJECT_DIR <- CFG$root
GTF_FILE <- file.path(PROJECT_DIR, "gencode.v50.annotation.gtf.gz")
TOKEN_FILE <- file.path(CFG$root, "HCC_37token.xlsx")
SEED <- 20261116L
set.seed(SEED)
options(warn = 1)
required_files <- c(gencode_gtf = GTF_FILE, token_table = TOKEN_FILE)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
    stop("\u4EE5\u4E0B\u8F93\u5165\u6587\u4EF6\u4E0D\u5B58\u5728\uFF1A\n", paste(missing_files, collapse = "\n"))
}
RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
RUN_DIR <- out_dir("01_reference")
if (dir.exists(RUN_DIR)) {
    stop("\u672C\u6B21\u8FD0\u884C\u76EE\u5F55\u5DF2\u5B58\u5728\uFF0C\u8BF7\u91CD\u65B0\u8FD0\u884C\u672C\u6BB5\u751F\u6210\u65B0\u7684 RUN_ID\u3002")
}
if (!dir.create(RUN_DIR, recursive = TRUE)) {
    stop("\u65E0\u6CD5\u521B\u5EFA\u8F93\u51FA\u76EE\u5F55\uFF1A", RUN_DIR)
}
writeLines(capture.output(sessionInfo()), file.path(RUN_DIR, "sessionInfo.txt"))
config <- list(project_dir = PROJECT_DIR, gtf_file = GTF_FILE, token_file = TOKEN_FILE, annotation_release = "GENCODE v50", 
    seed = SEED, run_id = RUN_ID, run_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), timezone = Sys.timezone())
saveRDS(config, file.path(RUN_DIR, "stage01_config.rds"))
cat("\u7B2C 1 \u6BB5\u5B8C\u6210\u3002\n")
cat("\u672C\u6B21\u8F93\u51FA\u76EE\u5F55\uFF1A", RUN_DIR, "\n")
input_info <- file.info(unname(required_files))
input_manifest <- data.frame(input_name = names(required_files), path = normalizePath(unname(required_files), winslash = "/", 
    mustWork = TRUE), size_bytes = input_info$size, modified_time = as.character(input_info$mtime), md5 = unname(tools::md5sum(unname(required_files))), 
    stringsAsFactors = FALSE)
if (anyNA(input_manifest$md5)) {
    stop("\u5B58\u5728\u65E0\u6CD5\u8BA1\u7B97 MD5 \u7684\u8F93\u5165\u6587\u4EF6\uFF0C\u8BF7\u68C0\u67E5\u6587\u4EF6\u8BFB\u53D6\u6743\u9650\u3002")
}
write.table(input_manifest, file.path(RUN_DIR, "input_manifest.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, 
    fileEncoding = "UTF-8")
extract_gtf_gene_lines <- function(path, chunk_size = 100000L) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    gene_chunks <- list()
    header_lines <- character()
    n_lines <- 0
    chunk_index <- 0L
    repeat {
        lines <- readLines(con, n = chunk_size, warn = FALSE)
        if (length(lines) == 0L) {
            break
        }
        n_lines <- n_lines + length(lines)
        is_comment <- startsWith(lines, "#")
        header_lines <- c(header_lines, lines[is_comment])
        is_gene <- !is_comment & grepl("^[^\t]+\t[^\t]+\tgene\t", lines)
        if (any(is_gene)) {
            chunk_index <- chunk_index + 1L
            gene_chunks[[chunk_index]] <- lines[is_gene]
        }
    }
    list(gene_lines = unlist(gene_chunks, use.names = FALSE), header = unique(header_lines), total_lines = n_lines)
}
gtf_raw <- extract_gtf_gene_lines(GTF_FILE)
if (length(gtf_raw$gene_lines) == 0L) {
    stop("\u672A\u63D0\u53D6\u5230 gene \u8BB0\u5F55\uFF0C\u8BF7\u68C0\u67E5 GTF \u6587\u4EF6\u683C\u5F0F\u3002")
}
writeLines(gtf_raw$header, file.path(RUN_DIR, "gencode_gtf_header.txt"), useBytes = TRUE)
cat("GTF \u603B\u884C\u6570\uFF1A", gtf_raw$total_lines, "\n")
cat("gene \u8BB0\u5F55\u6570\uFF1A", length(gtf_raw$gene_lines), "\n")
cat("GTF \u5934\u90E8\u4FE1\u606F\uFF1A\n")
cat(paste(head(gtf_raw$header, 10L), collapse = "\n"), "\n")
gtf_fields <- strsplit(gtf_raw$gene_lines, split = "\t", fixed = TRUE)
if (any(lengths(gtf_fields) != 9L)) {
    stop("\u5B58\u5728\u4E0D\u662F 9 \u5217\u7684 gene \u8BB0\u5F55\uFF0C\u8BF7\u68C0\u67E5 GTF\u3002")
}
gtf_matrix <- do.call(rbind, gtf_fields)
extract_gtf_attribute <- function(attributes, key) {
    pattern <- paste0("(^|;)[[:space:]]*", key, "[[:space:]]+\"([^\"]*)\"")
    matches <- regmatches(attributes, regexec(pattern, attributes, perl = TRUE))
    vapply(matches, function(x) {
        if (length(x) >= 3L) 
            x[[3L]]
        else NA_character_
    }, character(1))
}
gene_reference <- data.frame(ensembl_gene_id_version = extract_gtf_attribute(gtf_matrix[, 9L], "gene_id"), gene_symbol = extract_gtf_attribute(gtf_matrix[, 
    9L], "gene_name"), gene_type = extract_gtf_attribute(gtf_matrix[, 9L], "gene_type"), chromosome = gtf_matrix[, 
    1L], start = as.integer(gtf_matrix[, 4L]), end = as.integer(gtf_matrix[, 5L]), strand = gtf_matrix[, 7L], stringsAsFactors = FALSE)
gene_reference$ensembl_gene_id <- sub("\\.[0-9]+$", "", gene_reference$ensembl_gene_id_version)
required_annotation <- c("ensembl_gene_id_version", "ensembl_gene_id", "gene_symbol", "gene_type")
bad_annotation <- vapply(gene_reference[required_annotation], function(x) any(is.na(x) | !nzchar(x)), logical(1))
if (any(bad_annotation)) {
    stop("\u5173\u952E\u6CE8\u91CA\u5B57\u6BB5\u5B58\u5728\u7F3A\u5931\uFF1A", paste(names(bad_annotation)[bad_annotation], collapse = ", "))
}
if (anyDuplicated(gene_reference$ensembl_gene_id_version)) {
    stop("\u5B8C\u6574 Ensembl gene ID \u5B58\u5728\u91CD\u590D\uFF0C\u9700\u8981\u5148\u68C0\u67E5\u3002")
}
if (anyDuplicated(gene_reference$ensembl_gene_id)) {
    stop("\u53BB\u7248\u672C\u53F7 Ensembl gene ID \u5B58\u5728\u91CD\u590D\uFF0C\u9700\u8981\u5148\u68C0\u67E5\u6620\u5C04\u89C4\u5219\u3002")
}
if (anyNA(gene_reference$start) || anyNA(gene_reference$end) || any(gene_reference$end < gene_reference$start)) {
    stop("\u5B58\u5728\u65E0\u6548\u7684\u57FA\u56E0\u5750\u6807\u3002")
}
protein_coding_reference <- gene_reference[gene_reference$gene_type == "protein_coding", , drop = FALSE]
if (nrow(protein_coding_reference) == 0L) {
    stop("\u672A\u8BC6\u522B\u5230 protein_coding \u57FA\u56E0\u3002")
}
protein_coding_reference$symbol_ambiguous <- duplicated(protein_coding_reference$gene_symbol) | duplicated(protein_coding_reference$gene_symbol, 
    fromLast = TRUE)
gene_reference <- gene_reference[order(gene_reference$ensembl_gene_id), , drop = FALSE]
protein_coding_reference <- protein_coding_reference[order(protein_coding_reference$ensembl_gene_id), , drop = FALSE]
rownames(gene_reference) <- NULL
rownames(protein_coding_reference) <- NULL
cat("\u5168\u90E8 gene \u8BB0\u5F55\uFF1A", nrow(gene_reference), "\n")
cat("\u86CB\u767D\u7F16\u7801 gene ID \u6570\uFF1A", nrow(protein_coding_reference), "\n")
cat("\u86CB\u767D\u7F16\u7801\u552F\u4E00 symbol \u6570\uFF1A", length(unique(protein_coding_reference$gene_symbol)), "\n")
gene_type_counts <- as.data.frame(table(gene_reference$gene_type), stringsAsFactors = FALSE)
names(gene_type_counts) <- c("gene_type", "n_genes")
gene_type_counts <- gene_type_counts[order(-gene_type_counts$n_genes, gene_type_counts$gene_type), , drop = FALSE]
ambiguous_symbols <- protein_coding_reference[protein_coding_reference$symbol_ambiguous, , drop = FALSE]
qc_summary <- data.frame(metric = c("all_gene_records", "protein_coding_gene_ids", "protein_coding_unique_symbols", 
    "protein_coding_ids_with_ambiguous_symbol", "ambiguous_symbol_names", "protein_coding_chrM", "protein_coding_chrX", 
    "protein_coding_chrY"), value = c(nrow(gene_reference), nrow(protein_coding_reference), length(unique(protein_coding_reference$gene_symbol)), 
    nrow(ambiguous_symbols), length(unique(ambiguous_symbols$gene_symbol)), sum(protein_coding_reference$chromosome == 
        "chrM"), sum(protein_coding_reference$chromosome == "chrX"), sum(protein_coding_reference$chromosome == 
        "chrY")), stringsAsFactors = FALSE)
write_tsv <- function(x, filename) {
    write.table(x, file.path(RUN_DIR, filename), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA", fileEncoding = "UTF-8")
}
write_tsv(gene_reference, "gencode_v50_all_genes.tsv")
write_tsv(protein_coding_reference, "gencode_v50_protein_coding_genes.tsv")
write_tsv(ambiguous_symbols, "protein_coding_ambiguous_symbols.tsv")
write_tsv(gene_type_counts, "gene_type_counts.tsv")
write_tsv(qc_summary, "reference_qc_summary.tsv")
saveRDS(list(config = config, input_manifest = input_manifest, gtf_header = gtf_raw$header, gene_reference = gene_reference, 
    protein_coding_reference = protein_coding_reference, qc_summary = qc_summary), artifact("stage01_reference_bundle.rds"))
cat("\n===== \u672C\u9636\u6BB5\u8D28\u63A7\u6458\u8981 =====\n")
print(qc_summary, row.names = FALSE)
cat("\n===== \u6570\u91CF\u6700\u591A\u7684 15 \u7C7B\u57FA\u56E0 =====\n")
print(head(gene_type_counts, 15L), row.names = FALSE)
cat("\n===== \u86CB\u767D\u7F16\u7801\u53C2\u8003\u8868\u524D 6 \u884C =====\n")
print(head(protein_coding_reference), row.names = FALSE)
cat("\n===== \u540C\u540D\u86CB\u767D\u7F16\u7801\u57FA\u56E0\u793A\u4F8B\uFF0C\u6700\u591A\u663E\u793A 10 \u884C =====\n")
print(head(ambiguous_symbols, 10L), row.names = FALSE)
cat("\n\u7ED3\u679C\u4FDD\u5B58\u5230\uFF1A", RUN_DIR, "\n")
rm(gtf_raw, gtf_fields, gtf_matrix)
invisible(gc())

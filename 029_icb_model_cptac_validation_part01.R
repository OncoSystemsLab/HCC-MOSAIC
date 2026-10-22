if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

CPTAC.data_dir <- "resources/CPTAC_HCC/RNA-seq"

CPTAC.review_dir <- paste0("ICB.model/", "CPTAC_source_review")

CPTAC.out_dir <- paste0("ICB.model/", "CPTAC_validation_TPM_v19_159")

CPTAC.count_file <- file.path(CPTAC.data_dir, "CPTAC_159_readcount.tsv")

CPTAC.xlsx_file <- file.path(CPTAC.data_dir, "Table S1. Clinicopathoglogic and multi-omics data.xlsx")

CPTAC.mapping_file <- file.path(CPTAC.review_dir, "sample_identity_mapping.tsv")

CPTAC.gtf_file <- file.path(CPTAC.review_dir, "gencode.v19.annotation.gtf.gz")

CPTAC.n_expected <- 159L

CPTAC.wes_variants_file <- file.path(CPTAC.review_dir, "T777_WES_annotated_variants.tsv")

CPTAC.wes_qc_file <- file.path(CPTAC.review_dir, "T777_WES_annotation_QC.tsv")

CPTAC.OS_unit <- "months"

CPTAC.bootstrap_B <- 2000L

CPTAC.mutation_threshold <- 0.05

CPTAC.required_packages <- c("data.table", "readxl", "dplyr", "tidyr", "tibble", "rtracklayer", "GenomicRanges", 
    "IRanges", "ggplot2", "survival", "survminer", "glmnet", "GSVA", "BiocParallel")

CPTAC.missing_packages <- CPTAC.required_packages[!vapply(CPTAC.required_packages, requireNamespace, 
    quietly = TRUE, FUN.VALUE = logical(1))]

if (length(CPTAC.missing_packages)) {
    stop("Install the following dependencies before continuing: ", paste(CPTAC.missing_packages, collapse = ", "))
}

library(dplyr)

library(tidyr)

library(tibble)

library(ggplot2)

library(survival)

stopifnot(file.exists(CPTAC.count_file), file.exists(CPTAC.xlsx_file), file.exists(CPTAC.mapping_file), 
    file.exists(CPTAC.gtf_file))

dir.create(CPTAC.out_dir, recursive = TRUE, showWarnings = FALSE)

CPTAC.write <- function(x, filename) {
    data.table::fwrite(as.data.frame(x), file.path(CPTAC.out_dir, filename), sep = "\t", quote = FALSE, 
        na = "NA")
}

CPTAC.save_plot <- function(p, filename, width = 8, height = 5) {
    ggsave(file.path(CPTAC.out_dir, filename), plot = p, width = width, height = height, limitsize = FALSE)
}

CPTAC.expression_label <- "log2(TPM + 1)"

if (!exists("BA4_MODEL_BUNDLE")) {
    CPTAC.bundle_file <- paste0("ICB.model/", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds")
    if (!file.exists(CPTAC.bundle_file)) {
        stop("Load BA4_MODEL_BUNDLE or update the model bundle path first.")
    }
    BA4_MODEL_BUNDLE <- readRDS(CPTAC.bundle_file)
}

stopifnot(all(c("model", "lambda", "features", "background_genes", "gene_sets") %in% names(BA4_MODEL_BUNDLE)), 
    length(BA4_MODEL_BUNDLE$features) == 4L)

CPTAC.sample_mapping <- data.table::fread(CPTAC.mapping_file, data.table = FALSE, colClasses = "character", 
    check.names = FALSE)

stopifnot(all(c("Count_column", "Tissue_type", "Reference_sample_ID", "Clinical_sample_ID", "Patient_tumor_ID") %in% 
    names(CPTAC.sample_mapping)))

CPTAC.selected_mapping <- CPTAC.sample_mapping %>% filter(Tissue_type == "Tumor")

stopifnot(nrow(CPTAC.selected_mapping) == CPTAC.n_expected, !anyNA(CPTAC.selected_mapping$Count_column), 
    !anyNA(CPTAC.selected_mapping$Patient_tumor_ID), !anyDuplicated(CPTAC.selected_mapping$Count_column), 
    !anyDuplicated(CPTAC.selected_mapping$Patient_tumor_ID), all(CPTAC.selected_mapping$Reference_sample_ID == 
        CPTAC.selected_mapping$Patient_tumor_ID), all(CPTAC.selected_mapping$Clinical_sample_ID == CPTAC.selected_mapping$Patient_tumor_ID))

CPTAC.patient_ids <- CPTAC.selected_mapping$Patient_tumor_ID

CPTAC.write(CPTAC.selected_mapping, "01_selected_159_tumor_samples.tsv")

CPTAC.count_raw <- data.table::fread(CPTAC.count_file, header = TRUE, sep = " ", data.table = FALSE, 
    colClasses = "character", check.names = FALSE)

names(CPTAC.count_raw)[1] <- "Gene"

CPTAC.missing_count_columns <- setdiff(CPTAC.selected_mapping$Count_column, names(CPTAC.count_raw))

if (length(CPTAC.missing_count_columns)) {
    stop("The counts matrix is missing sample columns: ", paste(CPTAC.missing_count_columns, collapse = ", "))
}

stopifnot(!anyDuplicated(names(CPTAC.count_raw)))

CPTAC.raw_gene_ids <- trimws(CPTAC.count_raw$Gene)

CPTAC.raw_gene_ids <- trimws(CPTAC.count_raw$Gene)

CPTAC.is_gene <- grepl("^ENSGR?[0-9]+([.][0-9]+)?$", CPTAC.raw_gene_ids)

CPTAC.removed_rows <- data.frame(Data_row = which(!CPTAC.is_gene), First_column = CPTAC.raw_gene_ids[!CPTAC.is_gene])

CPTAC.write(CPTAC.removed_rows, "02_removed_internal_header_rows.tsv")

if (any(!CPTAC.is_gene & !(CPTAC.raw_gene_ids %in% c("Geneid", "GENES", "Gene")))) {
    stop("Unknown nongene rows remain; review the output file.")
}

stopifnot(sum(CPTAC.is_gene) == 57820L)

CPTAC.removed_rows <- data.frame(Data_row = which(!CPTAC.is_gene), First_column = CPTAC.raw_gene_ids[!CPTAC.is_gene])

CPTAC.write(CPTAC.removed_rows, "02_removed_internal_header_rows.tsv")

if (any(!CPTAC.is_gene & !CPTAC.raw_gene_ids %in% c("Geneid", "GENES", "Gene"))) {
    stop("Unknown non-Ensembl rows found; review 02_removed_internal_header_rows.tsv.")
}

CPTAC.counts <- as.matrix(CPTAC.count_raw[CPTAC.is_gene, CPTAC.selected_mapping$Count_column, drop = FALSE])

storage.mode(CPTAC.counts) <- "numeric"

rownames(CPTAC.counts) <- sub("\\.[0-9]+$", "", CPTAC.raw_gene_ids[CPTAC.is_gene])

colnames(CPTAC.counts) <- CPTAC.patient_ids

stopifnot(ncol(CPTAC.counts) == CPTAC.n_expected, nrow(CPTAC.counts) == 57820L, !anyDuplicated(rownames(CPTAC.counts)), 
    all(is.finite(CPTAC.counts)), all(CPTAC.counts >= 0), all(CPTAC.counts == floor(CPTAC.counts)), all(colSums(CPTAC.counts) > 
        0))

rm(CPTAC.count_raw)

CPTAC.gtf <- rtracklayer::import(CPTAC.gtf_file, format = "gtf")

stopifnot(all(c("type", "gene_id", "gene_name") %in% names(S4Vectors::mcols(CPTAC.gtf))))

CPTAC.gene_ranges <- CPTAC.gtf[as.character(CPTAC.gtf$type) == "gene"]

CPTAC.annotation <- data.frame(Ensembl_ID = sub("\\.[0-9]+$", "", as.character(CPTAC.gene_ranges$gene_id)), 
    Ensembl_ID_versioned = as.character(CPTAC.gene_ranges$gene_id), Gene_Symbol = as.character(CPTAC.gene_ranges$gene_name), 
    stringsAsFactors = FALSE)

stopifnot(!anyNA(CPTAC.annotation$Ensembl_ID), !anyDuplicated(CPTAC.annotation$Ensembl_ID), !anyNA(CPTAC.annotation$Gene_Symbol), 
    all(nzchar(CPTAC.annotation$Gene_Symbol)))

CPTAC.exons <- CPTAC.gtf[as.character(CPTAC.gtf$type) == "exon"]

CPTAC.exon_gene_ids <- sub("\\.[0-9]+$", "", as.character(CPTAC.exons$gene_id))

CPTAC.exons_by_gene <- split(CPTAC.exons, CPTAC.exon_gene_ids)

CPTAC.exons_by_gene <- GenomicRanges::GRangesList(CPTAC.exons_by_gene)

CPTAC.exon_union <- GenomicRanges::reduce(CPTAC.exons_by_gene, ignore.strand = FALSE)

CPTAC.exon_lengths <- vapply(CPTAC.exon_union, function(gr) sum(IRanges::width(gr)), FUN.VALUE = numeric(1))

stopifnot(!anyDuplicated(names(CPTAC.exon_lengths)), all(is.finite(CPTAC.exon_lengths)), all(CPTAC.exon_lengths > 
    0), all(rownames(CPTAC.counts) %in% names(CPTAC.exon_lengths)))

CPTAC.annotation$Exon_union_length_bp <- unname(CPTAC.exon_lengths[CPTAC.annotation$Ensembl_ID])

CPTAC.annotation_index <- match(rownames(CPTAC.counts), CPTAC.annotation$Ensembl_ID)

if (anyNA(CPTAC.annotation_index)) {
    stop("Some count-matrix genes do not match the v19 annotation; TPM calculation stopped.")
}

CPTAC.count_annotation <- CPTAC.annotation[CPTAC.annotation_index, , drop = FALSE]

stopifnot(identical(rownames(CPTAC.counts), CPTAC.count_annotation$Ensembl_ID), all(is.finite(CPTAC.count_annotation$Exon_union_length_bp)), 
    all(CPTAC.count_annotation$Exon_union_length_bp > 0))

CPTAC.write(CPTAC.count_annotation, "03_GENCODE_v19_gene_annotation_lengths.tsv")

CPTAC.annotation_QC <- data.frame(count_genes = nrow(CPTAC.counts), matched_genes = sum(!is.na(CPTAC.annotation_index)), 
    genes_with_positive_length = sum(CPTAC.count_annotation$Exon_union_length_bp > 0), annotation = "GENCODE v19", 
    length_definition = "Union of exonic intervals per gene")

CPTAC.write(CPTAC.annotation_QC, "03_annotation_QC.tsv")

CPTAC.RPK <- sweep(CPTAC.counts, MARGIN = 1, STATS = CPTAC.count_annotation$Exon_union_length_bp/1000, 
    FUN = "/")

CPTAC.RPK_sum <- colSums(CPTAC.RPK)

stopifnot(all(is.finite(CPTAC.RPK_sum)), all(CPTAC.RPK_sum > 0))

CPTAC.tpm_ensembl <- sweep(CPTAC.RPK, MARGIN = 2, STATS = CPTAC.RPK_sum, FUN = "/") * 1e+06

stopifnot(all(is.finite(CPTAC.tpm_ensembl)), all(CPTAC.tpm_ensembl >= 0), max(abs(colSums(CPTAC.tpm_ensembl) - 
    1e+06)) < 1e-05)

CPTAC.tpm <- rowsum(CPTAC.tpm_ensembl, group = CPTAC.count_annotation$Gene_Symbol, reorder = FALSE)

stopifnot(!anyDuplicated(rownames(CPTAC.tpm)), identical(colnames(CPTAC.tpm), CPTAC.patient_ids), max(abs(colSums(CPTAC.tpm) - 
    1e+06)) < 1e-05)

CPTAC.TPM_QC <- data.frame(Patient_ID = CPTAC.patient_ids, total_counts = colSums(CPTAC.counts), detected_Ensembl_genes = colSums(CPTAC.counts > 
    0), TPM_sum_Ensembl = colSums(CPTAC.tpm_ensembl), TPM_sum_Symbol = colSums(CPTAC.tpm))

CPTAC.write(CPTAC.TPM_QC, "04_TPM_sample_QC.tsv")

CPTAC.symbol_merge_audit <- CPTAC.count_annotation %>% group_by(Gene_Symbol) %>% summarise(n_Ensembl_ID = n(), 
    Ensembl_IDs = paste(Ensembl_ID, collapse = ";"), .groups = "drop") %>% filter(n_Ensembl_ID > 1)

CPTAC.write(CPTAC.symbol_merge_audit, "04_symbols_with_multiple_Ensembl_IDs.tsv")

saveRDS(CPTAC.counts, file.path(CPTAC.out_dir, "CPTAC_159_raw_counts_Ensembl.rds"))

saveRDS(CPTAC.tpm_ensembl, file.path(CPTAC.out_dir, "CPTAC_159_TPM_Ensembl_v19.rds"))

saveRDS(CPTAC.tpm, file.path(CPTAC.out_dir, "CPTAC_159_TPM_GeneSymbol_v19.rds"))

data.table::fwrite(data.table::as.data.table(CPTAC.tpm, keep.rownames = "Gene"), file.path(CPTAC.out_dir, 
    "CPTAC_159_TPM_GeneSymbol_v19.tsv.gz"), sep = "\t", quote = FALSE, compress = "gzip")

CPTAC.expr <- log2(CPTAC.tpm + 1)

CPTAC.log2_tpm <- CPTAC.expr

saveRDS(CPTAC.log2_tpm, file.path(CPTAC.out_dir, "CPTAC_159_log2TPM_GeneSymbol_v19.rds"))

rm(CPTAC.RPK, CPTAC.gtf, CPTAC.gene_ranges, CPTAC.exons, CPTAC.exons_by_gene, CPTAC.exon_union)

CPTAC.clinical_raw <- readxl::read_excel(CPTAC.xlsx_file, sheet = "1. Clinical and molecular data", .name_repair = "minimal")

CPTAC.clinical <- CPTAC.clinical_raw %>% transmute(Patient_ID = trimws(as.character(`Tumor (T) sample ID`)), 
    Normal_ID = trimws(as.character(`Adjacent liver tissue (N) sample ID`)), Age = as.numeric(Age), Gender = as.character(Gender), 
    OS_time_native = as.numeric(`Overall survial`), OS_status = as.numeric(`Survial status (1, dead; 0, alive)`), 
    RFS_time_native = as.numeric(`Recurrence-free survival`), RFS_status = as.numeric(`Recurrence status (1, yes; 0, no)`), 
    BCLC_stage = as.character(`BCLC stage`), TNM_stage = as.character(`TNM stage`), Tumor_size_cm = as.numeric(`Tumor size (cm)`), 
    Tumor_purity = as.numeric(`Tumor purity by HE staining`))

stopifnot(!anyDuplicated(CPTAC.clinical$Patient_ID), all(CPTAC.patient_ids %in% CPTAC.clinical$Patient_ID))

CPTAC.clinical <- CPTAC.clinical[match(CPTAC.patient_ids, CPTAC.clinical$Patient_ID), , drop = FALSE]

CPTAC.write(CPTAC.clinical, "05_clinical_159_patients.tsv")

CPTAC.variants_published <- readxl::read_excel(CPTAC.xlsx_file, sheet = "2. SNV & INDELs", .name_repair = "minimal")

if (!all(file.exists(c(CPTAC.wes_variants_file, CPTAC.wes_qc_file)))) {
    stop(paste0("159-patient RNA preprocessing saved. T777 WES annotation/QC is pending: ", CPTAC.wes_variants_file, 
        " ; ", CPTAC.wes_qc_file, ". Do not zero-fill T777 or reuse the 158-patient mutation matrix."))
}

CPTAC.wes_qc <- data.table::fread(CPTAC.wes_qc_file, data.table = FALSE, colClasses = "character")

stopifnot(all(c("Patient_ID", "WES_usable", "Annotation_comparable", "Variants_file_md5", "Reference_build", 
    "Annotation_method", "Concordance_report", "Concordance_report_md5", "Source_VCF", "Source_VCF_md5", 
    "PASS_records", "Retained_variants", "Record_audit", "Record_audit_md5") %in% names(CPTAC.wes_qc)))

CPTAC.wes_md5 <- unname(tools::md5sum(CPTAC.wes_variants_file))

stopifnot(nrow(CPTAC.wes_qc) == 1L, all(!is.na(CPTAC.wes_qc)), identical(CPTAC.wes_qc$Patient_ID, "T777"), 
    all(toupper(CPTAC.wes_qc$WES_usable) == "TRUE"), all(toupper(CPTAC.wes_qc$Annotation_comparable) == 
        "TRUE"), all(tolower(CPTAC.wes_qc$Variants_file_md5) == CPTAC.wes_md5), all(CPTAC.wes_qc$Reference_build %in% 
        c("hg19", "GRCh37")), all(nzchar(trimws(CPTAC.wes_qc$Annotation_method))), file.exists(CPTAC.wes_qc$Concordance_report), 
    all(tolower(CPTAC.wes_qc$Concordance_report_md5) == unname(tools::md5sum(CPTAC.wes_qc$Concordance_report))), 
    file.exists(CPTAC.wes_qc$Source_VCF), all(tolower(CPTAC.wes_qc$Source_VCF_md5) == unname(tools::md5sum(CPTAC.wes_qc$Source_VCF))), 
    file.exists(CPTAC.wes_qc$Record_audit), all(tolower(CPTAC.wes_qc$Record_audit_md5) == unname(tools::md5sum(CPTAC.wes_qc$Record_audit))))

CPTAC.T777_source_lines <- readLines(CPTAC.wes_qc$Source_VCF, warn = FALSE)

CPTAC.T777_record_audit <- data.table::fread(CPTAC.wes_qc$Record_audit, data.table = FALSE, colClasses = "character")

stopifnot("" %in% CPTAC.T777_source_lines, "" %in% CPTAC.T777_source_lines, sum(!startsWith(CPTAC.T777_source_lines, 
    "")) == as.integer(CPTAC.wes_qc$PASS_records), nrow(CPTAC.T777_record_audit) == as.integer(CPTAC.wes_qc$PASS_records), 
    nrow(CPTAC.T777_record_audit) > 0L, all(CPTAC.T777_record_audit$Patient_ID == "T777"), !anyNA(CPTAC.T777_record_audit$Final_retained))

CPTAC.variants_T777 <- data.table::fread(CPTAC.wes_variants_file, sep = "\t", data.table = FALSE, colClasses = "character", 
    check.names = FALSE)

CPTAC.variant_columns <- c("Gene symbol", "Sample ID", "Function in refGene", "Exonic Function in refGene", 
    "Chromosome", "Start", "End", "Reference (Ref.)", "Alteration (Alt.)", "VAF")

stopifnot(all(CPTAC.variant_columns %in% names(CPTAC.variants_T777)), all(CPTAC.variant_columns %in% 
    names(CPTAC.variants_published)), !anyNA(CPTAC.variants_T777[["Sample ID"]]), all(CPTAC.variants_T777[["Sample ID"]] == 
    "T777"), !("T777" %in% CPTAC.variants_published[["Sample ID"]]), setequal(unique(CPTAC.variants_published[["Sample ID"]]), 
    setdiff(CPTAC.patient_ids, "T777")), nrow(CPTAC.variants_T777) == as.integer(CPTAC.wes_qc$Retained_variants), 
    sum(toupper(CPTAC.T777_record_audit$Final_retained) == "TRUE") == nrow(CPTAC.variants_T777))

CPTAC.variants_raw <- dplyr::bind_rows(as.data.frame(lapply(CPTAC.variants_published[CPTAC.variant_columns], 
    as.character), check.names = FALSE), CPTAC.variants_T777[CPTAC.variant_columns])

CPTAC.variants <- CPTAC.variants_raw %>% transmute(Gene = trimws(as.character(`Gene symbol`)), Patient_ID = trimws(as.character(`Sample ID`)), 
    Function = trimws(as.character(`Function in refGene`)), Exonic_function = trimws(as.character(`Exonic Function in refGene`)), 
    Chromosome = as.character(Chromosome), Start = as.numeric(Start), End = as.numeric(End), Reference = as.character(`Reference (Ref.)`), 
    Alternate = as.character(`Alteration (Alt.)`), VAF = as.numeric(VAF))

stopifnot(!anyNA(CPTAC.variants$Patient_ID), all(CPTAC.variants$Patient_ID %in% CPTAC.patient_ids))

CPTAC.variant_functions <- c("nonsynonymous SNV", "stopgain", "stoploss", "frameshift deletion", "frameshift insertion", 
    "frameshift substitution", "nonframeshift deletion", "nonframeshift insertion", "nonframeshift substitution")

CPTAC.variant_filter_QC <- CPTAC.variants %>% filter(Patient_ID %in% CPTAC.patient_ids) %>% count(Function, 
    Exonic_function, name = "n_variants")

CPTAC.write(CPTAC.variant_filter_QC, "05_WES_function_counts_before_filter.tsv")

CPTAC.variants_filtered <- CPTAC.variants %>% filter(Patient_ID %in% CPTAC.patient_ids, !is.na(Gene), 
    nzchar(Gene), Gene != ".", Exonic_function %in% CPTAC.variant_functions | grepl("splicing", Function, 
        ignore.case = TRUE))

CPTAC.mut_pairs <- CPTAC.variants_filtered %>% distinct(Gene, Patient_ID)

CPTAC.mut <- matrix(0L, nrow = n_distinct(CPTAC.mut_pairs$Gene), ncol = length(CPTAC.patient_ids), dimnames = list(sort(unique(CPTAC.mut_pairs$Gene)), 
    CPTAC.patient_ids))

CPTAC.mut[cbind(match(CPTAC.mut_pairs$Gene, rownames(CPTAC.mut)), match(CPTAC.mut_pairs$Patient_ID, colnames(CPTAC.mut)))] <- 1L

stopifnot("CTNNB1" %in% rownames(CPTAC.mut), all(CPTAC.mut %in% c(0, 1)), identical(colnames(CPTAC.mut), 
    colnames(CPTAC.expr)), identical(colnames(CPTAC.expr), CPTAC.clinical$Patient_ID))

CPTAC.write(CPTAC.variants_filtered, "05_WES_variants_filtered.tsv")

saveRDS(CPTAC.mut, file.path(CPTAC.out_dir, "CPTAC_159_mutation_matrix.rds"))

CPTAC.WES_patient_records <- data.frame(Patient_ID = CPTAC.patient_ids, N_Retained_Variants = tabulate(match(CPTAC.variants_filtered$Patient_ID, 
    CPTAC.patient_ids), nbins = length(CPTAC.patient_ids)), N_Retained_Mutation_Genes = as.integer(colSums(CPTAC.mut)), 
    Source = ifelse(CPTAC.patient_ids == "T777", "Paired VCF reannotation; complete record-level audit", 
        "Published Table S1"))

CPTAC.write(CPTAC.WES_patient_records, "05_WES_patient_record_counts.tsv")

CPTAC.background_genes <- BA4_MODEL_BUNDLE$background_genes

CPTAC.BA4_feature_names <- BA4_MODEL_BUNDLE$features

CPTAC.BA4_gene_sets <- BA4_MODEL_BUNDLE$gene_sets[CPTAC.BA4_feature_names]

stopifnot(!anyDuplicated(CPTAC.background_genes), all(CPTAC.BA4_feature_names %in% names(CPTAC.BA4_gene_sets)))

CPTAC.symbol_mapping_used <- data.frame(source_gene = character(), target_gene = character(), stringsAsFactors = FALSE)

CPTAC.present_background <- CPTAC.background_genes[CPTAC.background_genes %in% rownames(CPTAC.expr)]

CPTAC.background_missing <- setdiff(CPTAC.background_genes, rownames(CPTAC.expr))

CPTAC.background_QC <- data.frame(training_background = length(CPTAC.background_genes), present = length(CPTAC.present_background), 
    missing = length(CPTAC.background_missing), coverage = length(CPTAC.present_background)/length(CPTAC.background_genes))

CPTAC.BA4_expr_input <- CPTAC.expr[CPTAC.present_background, CPTAC.patient_ids, drop = FALSE]

CPTAC.BA4_gene_sets_available <- lapply(CPTAC.BA4_gene_sets, function(z) {
    intersect(unique(z), rownames(CPTAC.BA4_expr_input))
})

CPTAC.BA4_gene_coverage <- bind_rows(lapply(CPTAC.BA4_feature_names, function(f) {
    genes <- unique(CPTAC.BA4_gene_sets[[f]])
    available <- CPTAC.BA4_gene_sets_available[[f]]
    data.frame(feature = f, total_genes = length(genes), present_genes = length(available), coverage = length(available)/length(genes), 
        missing_genes = paste(setdiff(genes, available), collapse = ";"))
}))

CPTAC.write(CPTAC.background_QC, "06_background_QC.tsv")

CPTAC.write(data.frame(Gene = CPTAC.background_missing), "06_missing_background_genes.tsv")

CPTAC.write(CPTAC.BA4_gene_coverage, "06_BA4_gene_coverage.tsv")

CPTAC.write(CPTAC.symbol_mapping_used, "06_symbol_mapping_used.tsv")

if (any(lengths(CPTAC.BA4_gene_sets_available) < 2L)) {
    stop("A BA4 gene set has fewer than two available genes; review coverage.")
}

CPTAC.constant_gene_QC <- bind_rows(lapply(CPTAC.BA4_feature_names, function(f) {
    genes <- CPTAC.BA4_gene_sets_available[[f]]
    x <- CPTAC.BA4_expr_input[genes, , drop = FALSE]
    constant <- genes[apply(x, 1, function(z) all(z == z[1]))]
    data.frame(feature = f, n_constant = length(constant), constant_genes = paste(constant, collapse = ";"))
}))

CPTAC.write(CPTAC.constant_gene_QC, "06_constant_gene_QC.tsv")

CPTAC.ssgsea_param <- GSVA::ssgseaParam(exprData = CPTAC.BA4_expr_input, geneSets = CPTAC.BA4_gene_sets_available, 
    minSize = 2, maxSize = 500, alpha = 0.25, normalize = FALSE)

CPTAC.ssgsea_BA4 <- GSVA::gsva(CPTAC.ssgsea_param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())

stopifnot(all(CPTAC.BA4_feature_names %in% rownames(CPTAC.ssgsea_BA4)), all(CPTAC.patient_ids %in% colnames(CPTAC.ssgsea_BA4)))

CPTAC.ssgsea_BA4 <- CPTAC.ssgsea_BA4[CPTAC.BA4_feature_names, CPTAC.patient_ids, drop = FALSE]

stopifnot(all(is.finite(CPTAC.ssgsea_BA4)))

CPTAC.BA4_features <- as.data.frame(t(CPTAC.ssgsea_BA4), check.names = FALSE) %>% rownames_to_column("Patient_ID")

CPTAC.X <- as.matrix(CPTAC.BA4_features[, CPTAC.BA4_feature_names, drop = FALSE])

CPTAC.BA4_features$BA4_ICBScore <- as.numeric(predict(BA4_MODEL_BUNDLE$model, newx = CPTAC.X, s = BA4_MODEL_BUNDLE$lambda, 
    type = "response"))

CPTAC.BA4_features$BA4_eta <- as.numeric(predict(BA4_MODEL_BUNDLE$model, newx = CPTAC.X, s = BA4_MODEL_BUNDLE$lambda, 
    type = "link"))

stopifnot(nrow(CPTAC.BA4_features) == CPTAC.n_expected, all(is.finite(CPTAC.BA4_features$BA4_ICBScore)), 
    all(CPTAC.BA4_features$BA4_ICBScore >= 0), all(CPTAC.BA4_features$BA4_ICBScore <= 1), all(is.finite(CPTAC.BA4_features$BA4_eta)))

CPTAC.CTNNB1_status <- data.frame(Patient_ID = CPTAC.patient_ids, CTNNB1_mut = as.integer(CPTAC.mut["CTNNB1", 
    CPTAC.patient_ids])) %>% mutate(CTNNB1_status = factor(ifelse(CTNNB1_mut == 1L, "MUT", "WT"), levels = c("WT", 
    "MUT")))

CPTAC.analysis <- CPTAC.BA4_features %>% left_join(CPTAC.CTNNB1_status, by = "Patient_ID") %>% left_join(CPTAC.clinical, 
    by = "Patient_ID")

stopifnot(nrow(CPTAC.analysis) == CPTAC.n_expected, identical(CPTAC.analysis$Patient_ID, CPTAC.patient_ids), 
    !anyNA(CPTAC.analysis$CTNNB1_mut))

CPTAC.write(CPTAC.analysis, "06_patient_BA4_clinical.tsv")

CPTAC.TPM_preprocessing_settings <- list(n_patients = CPTAC.n_expected, excluded_tumor = character(0), 
    annotation = "GENCODE v19", length_definition = "Per-gene exon union length", TPM_denominator_genes = nrow(CPTAC.counts), 
    gene_symbol_aggregation = "Sum TPM before log transformation", model_input_scale = "log2(TPM + 1)", 
    model_refitted = FALSE)

saveRDS(list(settings = CPTAC.TPM_preprocessing_settings, selected_mapping = CPTAC.selected_mapping, 
    annotation_QC = CPTAC.annotation_QC, TPM_QC = CPTAC.TPM_QC, background_QC = CPTAC.background_QC, 
    gene_coverage = CPTAC.BA4_gene_coverage, BA4_features = CPTAC.BA4_features, analysis = CPTAC.analysis), 
    file.path(CPTAC.out_dir, "CPTAC_159_TPM_preprocessing_and_BA4_results.rds"))

CPTAC.rank_biserial <- function(x, y) {
    nx <- length(x)
    ny <- length(y)
    if (!nx || !ny) 
        return(NA_real_)
    r <- rank(c(x, y), ties.method = "average")
    u <- sum(r[seq_len(nx)]) - nx * (nx + 1)/2
    2 * u/(nx * ny) - 1
}


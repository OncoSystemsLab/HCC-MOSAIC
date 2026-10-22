if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

search_dirs <- file.path(PROJECT_DIR, c("02_reference", "03_input"))

all_files <- unlist(lapply(search_dirs, function(d) {
    if (dir.exists(d)) {
        list.files(d, recursive = TRUE, full.names = TRUE)
    }
    else {
        character(0)
    }
}), use.names = FALSE)

pattern <- paste(c("TCGA", "LIHC", "ICGC", "MAF", "VCF", "WES", "mutation", "RNA", "TPM", "expression", 
    "clinical"), collapse = "|")

candidate_files <- all_files[grepl(pattern, basename(all_files), ignore.case = TRUE)]

file_inventory <- data.frame(File = basename(candidate_files), Path = candidate_files, Size_MB = round(file.info(candidate_files)$size/1024^2, 
    2), stringsAsFactors = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

RAW_DIR <- "raw_data"

external_files <- data.frame(Cohort = c(rep("TCGA_LIHC", 3), rep("ICGC_LICA_FR", 4)), Data_Type = c("Clinical", 
    "RNA_TPM", "Mutation", "Clinical", "RNA_TPM", "RNA_Normalized", "Mutation_WES"), Folder = c(rep("TCGA_LIHC_final", 
    3), rep("ICGC_LICA_FR_FINAL", 4)), Filename = c("TCGA_LIHC_clinical.tsv", "TCGA_LIHC_expression_TPM.tsv.gz", 
    "TCGA_LIHC_mutation_matrix.tsv.gz", "ICGC_LICA_FR_clinical_harmonized.tsv", "ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz", 
    "ICGC_LICA_FR_expression_normalized_GeneSymbol.tsv.gz", "ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"), 
    stringsAsFactors = FALSE)

external_files$Path <- file.path(RAW_DIR, external_files$Folder, external_files$Filename)

external_files$Exists <- file.exists(external_files$Path)

external_files$Size_MB <- round(file.info(external_files$Path)$size/1024^2, 2)

if (!all(external_files$Exists)) {
    print(external_files$Path[!external_files$Exists])
    stop("Some external data files are missing; review the paths")
}

inspect_external_file <- function(path) {
    is_gzip <- grepl("\\.gz$", path, ignore.case = TRUE)
    con <- if (is_gzip) {
        gzfile(path, open = "rt")
    }
    else {
        file(path, open = "rt")
    }
    on.exit(close(con), add = TRUE)
    lines <- readLines(con, n = 3L, warn = FALSE)
    if (length(lines) < 2L) {
        stop("The file contains fewer than two lines: ", path)
    }
    lines[1] <- sub("^﻿", "", lines[1])
    fields <- lapply(lines, function(x) {
        strsplit(x, "\t", fixed = TRUE)[[1]]
    })
    show_fields <- function(x, n = 8L) {
        paste(substr(head(x, n), 1, 45), collapse = " | ")
    }
    data.frame(N_Header_Fields = length(fields[[1]]), N_Row1_Fields = length(fields[[2]]), N_Row2_Fields = if (length(fields) >= 
        3L) {
        length(fields[[3]])
    }
    else {
        NA_integer_
    }, First_Header = fields[[1]][1], First_Row_ID = fields[[2]][1], Header_Preview = show_fields(fields[[1]]), 
        Row1_Preview = show_fields(fields[[2]]), stringsAsFactors = FALSE)
}

schema_list <- vector("list", nrow(external_files))

for (i in seq_len(nrow(external_files))) {
    f <- external_files[i, ]
    cat("\nInspecting:", f$Cohort, "|", f$Data_Type, "")
    schema_list[[i]] <- inspect_external_file(f$Path)
}

schema_audit <- cbind(external_files[, c("Cohort", "Data_Type", "Filename")], do.call(rbind, schema_list))

rownames(schema_audit) <- NULL

for (i in seq_len(nrow(schema_audit))) {
    x <- schema_audit[i, ]
    cat("\n-------------------------------------")
    cat("Dataset:", x$Cohort, "|", x$Data_Type, "")
    cat("Header columns:", x$N_Header_Fields, "")
    cat("First data row columns:", x$N_Row1_Fields, "")
    cat("First header field:", x$First_Header, "")
    cat("First row identifier:", x$First_Row_ID, "")
    cat("Header preview:", x$Header_Preview, "")
    cat("First data row:", x$Row1_Preview, "")
}

write.csv(external_files, file.path(PROJECT_DIR, "08_reports", "Stage3_1_external_file_inventory.csv"), 
    row.names = FALSE)

write.csv(schema_audit, file.path(PROJECT_DIR, "08_reports", "Stage3_1_external_schema_audit.csv"), row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

RAW_DIR <- "raw_data"

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

background <- config$background_genes

mechanism_genes <- unique(c(unlist(config$gene_panels, use.names = FALSE), "CTNNB1"))

stopifnot(length(background) == 21308L, length(mechanism_genes) == 63L)

files <- list(TCGA_RNA = file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_expression_TPM.tsv.gz"), TCGA_MUT = file.path(RAW_DIR, 
    "TCGA_LIHC_final", "TCGA_LIHC_mutation_matrix.tsv.gz"), ICGC_RNA = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz"), ICGC_NORMALIZED = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_expression_normalized_GeneSymbol.tsv.gz"), ICGC_MUT = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", 
    "ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"))

stopifnot(all(file.exists(unlist(files))))

clean_field <- function(x) {
    x <- sub("^﻿", "", x)
    x <- trimws(x)
    sub("^\"(.*)\"$", "\\1", x)
}

scan_gene_matrix <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    header_line <- readLines(con, n = 1L, warn = FALSE)
    header <- clean_field(strsplit(header_line, "\t", fixed = TRUE)[[1]])
    stopifnot(length(header) >= 2L, header[1] == "Gene")
    sample_ids <- header[-1]
    gene_chunks <- list()
    ctnnb1_rows <- list()
    repeat {
        lines <- readLines(con, n = 2000L, warn = FALSE)
        if (length(lines) == 0L) 
            break
        genes <- clean_field(sub("\t.*$", "", lines))
        gene_chunks[[length(gene_chunks) + 1L]] <- genes
        hits <- which(genes == "CTNNB1")
        if (length(hits) > 0L) {
            for (h in hits) {
                values <- clean_field(strsplit(lines[h], "\t", fixed = TRUE)[[1]])
                ctnnb1_rows[[length(ctnnb1_rows) + 1L]] <- values[-1]
            }
        }
    }
    all_genes <- unlist(gene_chunks, use.names = FALSE)
    list(genes = all_genes, samples = sample_ids, n_genes = length(all_genes), n_samples = length(sample_ids), 
        duplicated_genes = sum(duplicated(all_genes)), duplicated_samples = sum(duplicated(sample_ids)), 
        ctnnb1_rows = ctnnb1_rows)
}

matrix_info <- lapply(files, scan_gene_matrix)

dimension_qc <- do.call(rbind, lapply(names(matrix_info), function(nm) {
    x <- matrix_info[[nm]]
    data.frame(Dataset = nm, Gene_Rows = x$n_genes, Sample_Columns = x$n_samples, Duplicate_Genes = x$duplicated_genes, 
        Duplicate_Samples = x$duplicated_samples, CTNNB1_Present = "CTNNB1" %in% x$genes)
}))

rownames(dimension_qc) <- NULL

rna_names <- c("TCGA_RNA", "ICGC_RNA", "ICGC_NORMALIZED")

coverage_qc <- do.call(rbind, lapply(rna_names, function(nm) {
    genes <- unique(matrix_info[[nm]]$genes)
    data.frame(Dataset = nm, Background_Total = length(background), Background_Matched = sum(background %in% 
        genes), Background_Missing = sum(!background %in% genes), Background_Coverage = round(100 * mean(background %in% 
        genes), 3), Mechanism_Matched = sum(mechanism_genes %in% genes), Mechanism_Total = length(mechanism_genes))
}))

rownames(coverage_qc) <- NULL

missing_background <- do.call(rbind, lapply(rna_names, function(nm) {
    missing <- setdiff(background, matrix_info[[nm]]$genes)
    data.frame(Dataset = rep(nm, length(missing)), Missing_Gene = missing, stringsAsFactors = FALSE)
}))

write.csv(missing_background, file.path(PROJECT_DIR, "08_reports", "Stage3_1B_missing_background_genes.csv"), 
    row.names = FALSE)

tcga_clinical <- read.delim(file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_clinical.tsv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

icgc_clinical <- read.delim(file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", "ICGC_LICA_FR_clinical_harmonized.tsv"), 
    check.names = FALSE, stringsAsFactors = FALSE)

stopifnot("Patient_ID" %in% colnames(tcga_clinical), "Patient_ID" %in% colnames(icgc_clinical))

tcga_ids <- clean_field(as.character(tcga_clinical$Patient_ID))

icgc_ids <- clean_field(as.character(icgc_clinical$Patient_ID))

make_pairing <- function(cohort_name, rna_samples, mutation_samples, clinical_ids) {
    all_ids <- unique(c(rna_samples, mutation_samples, clinical_ids))
    pairing <- data.frame(Cohort = cohort_name, Patient_ID = all_ids, RNA_Available = all_ids %in% rna_samples, 
        Mutation_Matrix_Available = all_ids %in% mutation_samples, Clinical_Available = all_ids %in% 
            clinical_ids, stringsAsFactors = FALSE)
    pairing$All_Three_Available <- with(pairing, RNA_Available & Mutation_Matrix_Available & Clinical_Available)
    pairing
}

tcga_pairing <- make_pairing("TCGA_LIHC", matrix_info$TCGA_RNA$samples, matrix_info$TCGA_MUT$samples, 
    tcga_ids)

icgc_pairing <- make_pairing("ICGC_LICA_FR", matrix_info$ICGC_RNA$samples, matrix_info$ICGC_MUT$samples, 
    icgc_ids)

all_pairing <- rbind(tcga_pairing, icgc_pairing)

rownames(all_pairing) <- NULL

pairing_summary <- do.call(rbind, lapply(split(all_pairing, all_pairing$Cohort), function(df) {
    data.frame(Cohort = df$Cohort[1], RNA = sum(df$RNA_Available), Mutation = sum(df$Mutation_Matrix_Available), 
        Clinical = sum(df$Clinical_Available), RNA_Mutation = sum(df$RNA_Available & df$Mutation_Matrix_Available), 
        All_Three = sum(df$All_Three_Available))
}))

rownames(pairing_summary) <- NULL

for (nm in c("TCGA_MUT", "ICGC_MUT")) {
    x <- matrix_info[[nm]]
    cat("\nDataset:", nm, "")
    cat("CTNNB1 row count:", length(x$ctnnb1_rows), "")
    if (length(x$ctnnb1_rows) == 1L) {
        vals <- x$ctnnb1_rows[[1]]
        cat("Values:", paste(unique(vals), collapse = ", "), "")
        if (length(vals) == x$n_samples) {
            numeric_vals <- suppressWarnings(as.numeric(vals))
            is_binary <- all(!is.na(numeric_vals) & numeric_vals %in% c(0, 1))
            cat("Binary 0/1:", is_binary, "")
            if (is_binary) {
                cat("CTNNB1 recorded variant (1):", sum(numeric_vals == 1), "")
                cat("CTNNB1 no recorded variant (0):", sum(numeric_vals == 0), "")
            }
        }
        else {
            cat("WARNING: CTNNB1 row width mismatch")
        }
    }
    else {
        cat("WARNING: CTNNB1 row absent or duplicated")
    }
}

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(dimension_qc, file.path(report_dir, "Stage3_1B_matrix_dimensions.csv"), row.names = FALSE)

write.csv(coverage_qc, file.path(report_dir, "Stage3_1B_gene_coverage.csv"), row.names = FALSE)

write.csv(pairing_summary, file.path(report_dir, "Stage3_1B_pairing_summary.csv"), row.names = FALSE)

write.csv(all_pairing, file.path(report_dir, "Stage3_1B_patient_pairing.csv"), row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("matrix_info"), exists("background"), exists("config"))

ba4 <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

signature_sets <- config$signature_sets

model_features <- ba4$features

datasets <- c("TCGA_RNA", "ICGC_RNA", "ICGC_NORMALIZED")

signature_coverage <- do.call(rbind, lapply(datasets, function(ds) {
    available <- unique(matrix_info[[ds]]$genes)
    do.call(rbind, lapply(names(signature_sets), function(sig) {
        genes <- unique(signature_sets[[sig]])
        genes_bg <- intersect(genes, background)
        matched <- intersect(genes_bg, available)
        missing <- setdiff(genes_bg, available)
        data.frame(Dataset = ds, Feature = sig, BA4_Model_Feature = sig %in% model_features, N_Reference_Genes = length(genes_bg), 
            N_Matched = length(matched), N_Missing = length(missing), Coverage_Pct = round(100 * length(matched)/max(1, 
                length(genes_bg)), 2), Missing_Genes = paste(missing, collapse = ";"), stringsAsFactors = FALSE)
    }))
}))

rownames(signature_coverage) <- NULL

summary_table <- do.call(rbind, lapply(datasets, function(ds) {
    df <- subset(signature_coverage, Dataset == ds)
    data.frame(Dataset = ds, Signatures = nrow(df), Fully_Matched = sum(df$N_Missing == 0), Partially_Missing = sum(df$N_Missing > 
        0), Median_Coverage = median(df$Coverage_Pct), Min_Coverage = min(df$Coverage_Pct))
}))

model_coverage <- signature_coverage[signature_coverage$BA4_Model_Feature, c("Dataset", "Feature", "N_Reference_Genes", 
    "N_Matched", "N_Missing", "Coverage_Pct")]

for (ds in datasets) {
    cat("\nDataset:", ds, "")
    df <- signature_coverage[signature_coverage$Dataset == ds, ]
    df <- df[order(df$Coverage_Pct), ]
    print(head(df[, c("Feature", "N_Missing", "Coverage_Pct")], 10), row.names = FALSE)
}

priority_genes <- c("CTNNB1", "AXIN2", "GLUL", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", "SLC10A1", "NR1H4", 
    "NR0B2", "CD274", "CD8A", "LAG3", "STAT1", "GZMB", "IFNG", "NKG7")

priority_coverage <- data.frame(Gene = priority_genes, TCGA_TPM = priority_genes %in% matrix_info$TCGA_RNA$genes, 
    ICGC_TPM = priority_genes %in% matrix_info$ICGC_RNA$genes, ICGC_Normalized = priority_genes %in% 
        matrix_info$ICGC_NORMALIZED$genes)

icgc_tpm_genes <- matrix_info$ICGC_RNA$genes

icgc_norm_genes <- matrix_info$ICGC_NORMALIZED$genes

rescued_in_norm <- intersect(setdiff(background, icgc_tpm_genes), icgc_norm_genes)

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(signature_coverage, file.path(report_dir, "Stage3_1B_signature_coverage.csv"), row.names = FALSE)

write.csv(priority_coverage, file.path(report_dir, "Stage3_1B_priority_gene_coverage.csv"), row.names = FALSE)

write.csv(data.frame(Gene = rescued_in_norm), file.path(report_dir, "Stage3_1B_ICGC_norm_only_background.csv"), 
    row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

RAW_DIR <- "raw_data"

mutation_paths <- c(TCGA_LIHC = file.path(RAW_DIR, "TCGA_LIHC_final", "TCGA_LIHC_mutation_matrix.tsv.gz"), 
    ICGC_LICA_FR = file.path(RAW_DIR, "ICGC_LICA_FR_FINAL", "ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"))

read_mutation_matrix <- function(path) {
    con <- gzfile(path, open = "rt")
    on.exit(close(con), add = TRUE)
    df <- read.delim(con, check.names = FALSE, stringsAsFactors = FALSE)
    stopifnot(colnames(df)[1] == "Gene", !anyDuplicated(df[[1]]), !anyDuplicated(colnames(df)[-1]))
    df
}

wnt_genes <- c("CTNNB1", "AXIN1", "APC", "AXIN2", "RNF43", "ZNRF3")

mutation_qc <- list()

driver_qc <- list()

for (cohort in names(mutation_paths)) {
    cat("\n===== MUTATION AUDIT:", cohort, "=====")
    df <- read_mutation_matrix(mutation_paths[[cohort]])
    genes <- as.character(df[[1]])
    raw_values <- as.character(unlist(df[-1], use.names = FALSE))
    distinct_codes <- sort(unique(raw_values))
    binary <- all(!is.na(raw_values) & raw_values %in% c("0", "1"))
    cat("Genes:", nrow(df), "")
    cat("Patients:", ncol(df) - 1L, "")
    cat("Observed codes:", paste(head(distinct_codes, 15), collapse = ", "), "")
    cat("Strict binary 0/1:", binary, "")
    mutation_qc[[cohort]] <- data.frame(Cohort = cohort, N_Genes = nrow(df), N_Patients = ncol(df) - 
        1L, Binary_0_1 = binary, N_Missing_Cells = sum(is.na(raw_values)), N_Other_Codes = sum(!is.na(raw_values) & 
        !(raw_values %in% c("0", "1"))))
    if (!binary) {
        cat("WARNING: Matrix is not strictly binary.", "Do not interpret 1 as mutation before reviewing codes.")
        next
    }
    mut <- as.matrix(df[-1])
    storage.mode(mut) <- "integer"
    rownames(mut) <- genes
    colnames(mut) <- colnames(df)[-1]
    cat("Patients with at least one recorded variant:", sum(colSums(mut) > 0), "")
    cat("Median mutated genes per patient:", median(colSums(mut)), "")
    driver_table <- do.call(rbind, lapply(wnt_genes, function(g) {
        present <- g %in% rownames(mut)
        data.frame(Cohort = cohort, Gene = g, Recorded_In_Matrix = present, N_Mutated = if (present) {
            sum(mut[g, ] == 1)
        }
        else {
            NA_integer_
        }, N_Patients = ncol(mut), stringsAsFactors = FALSE)
    }))
    driver_qc[[cohort]] <- driver_table
    print(driver_table, row.names = FALSE)
    rm(df, mut, raw_values)
    invisible(gc(verbose = FALSE))
}

mutation_summary <- do.call(rbind, mutation_qc)

rownames(mutation_summary) <- NULL

driver_summary <- do.call(rbind, driver_qc)

if (!is.null(driver_summary)) {
    rownames(driver_summary) <- NULL
}

if (!is.null(driver_summary)) {
    cat("\n===== WNT DRIVER MUTATION SUMMARY =====")
    print(driver_summary, row.names = FALSE)
}

write.csv(mutation_summary, file.path(PROJECT_DIR, "08_reports", "Stage3_1B_mutation_coding_QC.csv"), 
    row.names = FALSE)

if (!is.null(driver_summary)) {
    write.csv(driver_summary, file.path(PROJECT_DIR, "08_reports", "Stage3_1B_WNT_driver_mutations.csv"), 
        row.names = FALSE)
}

PROJECT_DIR <- "HCC_OmicsAgent"

RAW_DIR <- "raw_data"

config <- readRDS(file.path(PROJECT_DIR, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))

mechanism_genes <- unique(c(unlist(config$gene_panels, use.names = FALSE), "CTNNB1"))


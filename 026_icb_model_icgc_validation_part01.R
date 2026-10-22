if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(data.table)

ICGC.dir <- "raw_data/ICGC_LICA_FR_FINAL"

ICGC.files <- c(expression = file.path(ICGC.dir, "ICGC_LICA_FR_expression_TPM_GeneSymbol.tsv.gz"), mutation = file.path(ICGC.dir, 
    "ICGC_LICA_FR_mutation_WES_GeneSymbol.tsv.gz"), clinical = file.path(ICGC.dir, "ICGC_LICA_FR_clinical_harmonized.tsv"))

stopifnot(all(file.exists(ICGC.files)))

ICGC.expr.raw <- fread(ICGC.files[["expression"]], data.table = FALSE, check.names = FALSE)

ICGC.mut.raw <- fread(ICGC.files[["mutation"]], data.table = FALSE, check.names = FALSE)

ICGC.clinical <- fread(ICGC.files[["clinical"]], data.table = FALSE, check.names = FALSE)

ICGC.make_matrix <- function(dat, label) {
    if (!"Gene" %in% names(dat)) {
        stop(label, " is missing the Gene column")
    }
    genes <- as.character(dat$Gene)
    sample_ids <- names(dat)[names(dat) != "Gene"]
    if (anyNA(genes) || any(!nzchar(trimws(genes)))) {
        stop(label, " contains missing or blank gene names")
    }
    if (anyDuplicated(genes)) {
        stop(label, " contains duplicate gene names; review before continuing")
    }
    if (anyDuplicated(sample_ids)) {
        stop(label, " contains duplicate patient column names")
    }
    mat <- as.matrix(dat[, sample_ids, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- genes
    if (any(!is.finite(mat))) {
        stop(label, " contains NA, Inf or nonnumeric values")
    }
    mat
}

ICGC.expr <- ICGC.make_matrix(ICGC.expr.raw, "ICGC expression")

ICGC.mut <- ICGC.make_matrix(ICGC.mut.raw, "ICGC mutation")

stopifnot(all(ICGC.expr >= 0), all(ICGC.mut %in% c(0, 1)))

stopifnot("Patient_ID" %in% names(ICGC.clinical))

ICGC.clinical$Patient_ID <- as.character(ICGC.clinical$Patient_ID)

stopifnot(!anyNA(ICGC.clinical$Patient_ID), all(nzchar(trimws(ICGC.clinical$Patient_ID))), !anyDuplicated(ICGC.clinical$Patient_ID))

ICGC.expr_ids <- colnames(ICGC.expr)

ICGC.mut_ids <- colnames(ICGC.mut)

ICGC.clin_ids <- ICGC.clinical$Patient_ID

ICGC.ID_QC <- data.frame(comparison = c("Expression vs mutation", "Expression vs clinical"), identical_patient_sets = c(setequal(ICGC.expr_ids, 
    ICGC.mut_ids), setequal(ICGC.expr_ids, ICGC.clin_ids)))

stopifnot(all(ICGC.ID_QC$identical_patient_sets))

ICGC.mut <- ICGC.mut[, ICGC.expr_ids, drop = FALSE]

ICGC.clinical <- ICGC.clinical[match(ICGC.expr_ids, ICGC.clinical$Patient_ID), , drop = FALSE]

rownames(ICGC.clinical) <- NULL

stopifnot(identical(colnames(ICGC.expr), colnames(ICGC.mut)), identical(colnames(ICGC.expr), ICGC.clinical$Patient_ID))

ICGC.required_clinical <- c("Patient_ID", "OS_time_months", "OS_status", "Gender", "Age", "Edmonson_grade", 
    "Vascular_invasion")

stopifnot(all(ICGC.required_clinical %in% names(ICGC.clinical)))

ICGC.to_numeric <- function(x, label) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (any(!is.na(x) & is.na(y))) {
        stop(label, " contains nonnumeric values")
    }
    if (any(!is.na(y) & !is.finite(y))) {
        stop(label, " contains nonfinite values")
    }
    y
}

for (nm in c("OS_time_months", "OS_status", "Age")) {
    ICGC.clinical[[nm]] <- ICGC.to_numeric(ICGC.clinical[[nm]], nm)
}

stopifnot(all(is.na(ICGC.clinical$OS_status) | ICGC.clinical$OS_status %in% c(0, 1)), all(is.na(ICGC.clinical$OS_time_months) | 
    ICGC.clinical$OS_time_months >= 0))

ICGC.OS_eligible <- with(ICGC.clinical, !is.na(OS_time_months) & OS_time_months > 0 & !is.na(OS_status))

ICGC.data_QC <- data.frame(item = c("Expression genes", "Mutation genes", "Aligned patients", "Patients with usable OS", 
    "Deaths among patients with usable OS"), value = c(nrow(ICGC.expr), nrow(ICGC.mut), ncol(ICGC.expr), 
    sum(ICGC.OS_eligible), sum(ICGC.clinical$OS_status[ICGC.OS_eligible] == 1)))

ICGC.driver_genes <- c("CTNNB1", "TP53", "AXIN1", "APC", "TERT")

ICGC.driver_mutation_summary <- do.call(rbind, lapply(ICGC.driver_genes, function(gene) {
    present <- gene %in% rownames(ICGC.mut)
    n_mut <- if (present) {
        sum(ICGC.mut[gene, ] > 0)
    }
    else {
        NA_integer_
    }
    data.frame(Gene = gene, present_in_matrix = present, n_patients = ncol(ICGC.mut), n_mut = n_mut, 
        frequency = n_mut/ncol(ICGC.mut))
}))

stopifnot("CTNNB1" %in% rownames(ICGC.mut))

ICGC.model_environment_QC <- data.frame(object = c("BA4_MODEL_BUNDLE", "run_ssgsea", "predict_BA4_ICBScore"), 
    available = c(exists("BA4_MODEL_BUNDLE"), exists("run_ssgsea", mode = "function"), exists("predict_BA4_ICBScore", 
        mode = "function")))

library(dplyr)

library(tibble)

library(GSVA)

library(BiocParallel)

library(glmnet)

ICGC.OS_excluded <- ICGC.clinical[!ICGC.OS_eligible, c("Patient_ID", "OS_time_months", "OS_status"), 
    drop = FALSE]

ICGC.bundle_required <- c("model", "lambda", "features", "background_genes", "gene_sets")

stopifnot(all(ICGC.bundle_required %in% names(BA4_MODEL_BUNDLE)), !any(vapply(BA4_MODEL_BUNDLE[ICGC.bundle_required], 
    is.null, logical(1))))

ICGC.background_genes <- BA4_MODEL_BUNDLE$background_genes

ICGC.BA4_gene_sets <- BA4_MODEL_BUNDLE$gene_sets

ICGC.BA4_feature_names <- BA4_MODEL_BUNDLE$features

stopifnot(length(ICGC.background_genes) > 0, !anyNA(ICGC.background_genes), !anyDuplicated(ICGC.background_genes), 
    length(ICGC.BA4_feature_names) == 4, !anyDuplicated(ICGC.BA4_feature_names), !anyDuplicated(names(ICGC.BA4_gene_sets)), 
    all(ICGC.BA4_feature_names %in% names(ICGC.BA4_gene_sets)))

ICGC.BA4_gene_sets <- ICGC.BA4_gene_sets[ICGC.BA4_feature_names]

ICGC.background_overlap <- intersect(ICGC.background_genes, rownames(ICGC.expr))

ICGC.background_missing <- setdiff(ICGC.background_genes, rownames(ICGC.expr))

ICGC.background_QC <- data.frame(training_background = length(ICGC.background_genes), ICGC_present = length(ICGC.background_overlap), 
    missing = length(ICGC.background_missing), coverage = length(ICGC.background_overlap)/length(ICGC.background_genes))

ICGC.log2_tpm <- log2(ICGC.expr + 1)

ICGC.background_present_ordered <- ICGC.background_genes[ICGC.background_genes %in% rownames(ICGC.log2_tpm)]

ICGC.log2_tpm_background <- ICGC.log2_tpm[ICGC.background_present_ordered, , drop = FALSE]

ICGC.BA4_expr_input <- ICGC.log2_tpm_background

stopifnot(nrow(ICGC.BA4_expr_input) > 0, all(is.finite(ICGC.BA4_expr_input)), identical(colnames(ICGC.BA4_expr_input), 
    ICGC.clinical$Patient_ID))

ICGC.BA4_gene_coverage <- bind_rows(lapply(names(ICGC.BA4_gene_sets), function(gs) {
    genes <- unique(ICGC.BA4_gene_sets[[gs]])
    present_expr <- intersect(genes, rownames(ICGC.expr))
    present_input <- intersect(genes, rownames(ICGC.BA4_expr_input))
    data.frame(gene_set = gs, total_genes = length(genes), present_in_expression = length(present_expr), 
        present_in_model_input = length(present_input), input_coverage = length(present_input)/length(genes), 
        missing_from_input = paste(setdiff(genes, present_input), collapse = ";"))
}))

ICGC.BA4_zero_gene_QC <- bind_rows(lapply(names(ICGC.BA4_gene_sets), function(gs) {
    genes <- intersect(unique(ICGC.BA4_gene_sets[[gs]]), rownames(ICGC.BA4_expr_input))
    x <- ICGC.BA4_expr_input[genes, , drop = FALSE]
    zero_all <- genes[rowSums(x > 0) == 0]
    constant <- if (length(genes) > 0) {
        genes[apply(x, 1, function(z) all(z == z[1]))]
    }
    else {
        character(0)
    }
    data.frame(gene_set = gs, n_genes = length(genes), n_all_zero = length(zero_all), all_zero_genes = paste(zero_all, 
        collapse = ";"), n_constant = length(constant), constant_genes = paste(constant, collapse = ";"))
}))

ICGC.BA4_gene_sets_available <- lapply(ICGC.BA4_gene_sets, function(gs) {
    intersect(unique(gs), rownames(ICGC.BA4_expr_input))
})

if (any(lengths(ICGC.BA4_gene_sets_available) < 2)) {
    stop("A BA4 gene set has fewer than two available genes; review coverage first.")
}

ICGC.ssgsea_BA4 <- run_ssgsea(expr_matrix = ICGC.BA4_expr_input, gene_sets = ICGC.BA4_gene_sets_available)

stopifnot(setequal(rownames(ICGC.ssgsea_BA4), ICGC.BA4_feature_names), setequal(colnames(ICGC.ssgsea_BA4), 
    colnames(ICGC.expr)))

ICGC.ssgsea_BA4 <- ICGC.ssgsea_BA4[ICGC.BA4_feature_names, colnames(ICGC.expr), drop = FALSE]

stopifnot(all(is.finite(ICGC.ssgsea_BA4)))

ICGC.BA4_features <- as.data.frame(t(ICGC.ssgsea_BA4), check.names = FALSE) %>% rownames_to_column("Patient_ID")

ICGC.BA4_features$BA4_ICBScore <- predict_BA4_ICBScore(newdata = ICGC.BA4_features, model = BA4_MODEL_BUNDLE$model, 
    features = BA4_MODEL_BUNDLE$features, lambda = BA4_MODEL_BUNDLE$lambda)

stopifnot(length(ICGC.BA4_features$BA4_ICBScore) == ncol(ICGC.expr), all(is.finite(ICGC.BA4_features$BA4_ICBScore)), 
    all(ICGC.BA4_features$BA4_ICBScore >= 0), all(ICGC.BA4_features$BA4_ICBScore <= 1))

ICGC.CTNNB1_status <- data.frame(Patient_ID = colnames(ICGC.mut), CTNNB1_mut = as.integer(ICGC.mut["CTNNB1", 
    ] > 0), stringsAsFactors = FALSE) %>% mutate(CTNNB1_status = factor(ifelse(CTNNB1_mut == 1, "CTNNB1_MUT", 
    "CTNNB1_WT"), levels = c("CTNNB1_WT", "CTNNB1_MUT")))

ICGC.analysis <- ICGC.BA4_features %>% left_join(ICGC.CTNNB1_status, by = "Patient_ID") %>% left_join(ICGC.clinical, 
    by = "Patient_ID")

stopifnot(nrow(ICGC.analysis) == ncol(ICGC.expr), !anyDuplicated(ICGC.analysis$Patient_ID), identical(ICGC.analysis$Patient_ID, 
    colnames(ICGC.expr)), !anyNA(ICGC.analysis$CTNNB1_mut))

ICGC.BA4_score_QC <- ICGC.analysis %>% group_by(CTNNB1_status) %>% summarise(n = n(), score_min = min(BA4_ICBScore), 
    score_Q1 = quantile(BA4_ICBScore, 0.25), score_median = median(BA4_ICBScore), score_mean = mean(BA4_ICBScore), 
    score_Q3 = quantile(BA4_ICBScore, 0.75), score_max = max(BA4_ICBScore), .groups = "drop")

library(dplyr)

ICGC.annotation_db <- org.Hs.eg.db::org.Hs.eg.db

ICGC.resolve_gene_ids <- function(gene_names, db) {
    gene_names <- unique(as.character(gene_names))
    result <- data.frame(input_name = gene_names, ENTREZID = "NA", mapping_type = "unmapped", n_entrez = 0L, 
        stringsAsFactors = FALSE)
    for (key_type in c("SYMBOL", "ALIAS")) {
        pending <- result$input_name[result$mapping_type == "unmapped"]
        valid_keys <- intersect(pending, AnnotationDbi::keys(db, keytype = key_type))
        if (length(valid_keys) == 0L) 
            next
        mapped <- AnnotationDbi::select(db, keys = valid_keys, columns = "ENTREZID", keytype = key_type)
        id_lists <- split(as.character(mapped$ENTREZID), mapped[[key_type]])
        for (gene in names(id_lists)) {
            ids <- unique(id_lists[[gene]])
            ids <- ids[!is.na(ids) & nzchar(ids)]
            if (length(ids) == 0L) 
                next
            i <- match(gene, result$input_name)
            result$n_entrez[i] <- length(ids)
            if (length(ids) == 1L) {
                result$ENTREZID[i] <- ids
                result$mapping_type[i] <- key_type
            }
            else {
                result$mapping_type[i] <- paste0(key_type, "_ambiguous")
            }
        }
    }
    result
}

ICGC.expr_annotation <- ICGC.resolve_gene_ids(rownames(ICGC.expr), ICGC.annotation_db)

ICGC.background_annotation <- ICGC.resolve_gene_ids(ICGC.background_genes, ICGC.annotation_db)

ICGC.annotation_versions <- data.frame(package = c("AnnotationDbi", "org.Hs.eg.db"), version = c(as.character(packageVersion("AnnotationDbi")), 
    as.character(packageVersion("org.Hs.eg.db"))))

ICGC.expr_id_counts <- ICGC.expr_annotation %>% filter(!is.na(ENTREZID)) %>% count(ENTREZID, name = "n_expression_names")

ICGC.background_id_counts <- ICGC.background_annotation %>% filter(!is.na(ENTREZID)) %>% count(ENTREZID, 
    name = "n_background_names")

ICGC.target_annotation <- ICGC.background_annotation %>% filter(input_name %in% ICGC.background_missing, 
    !is.na(ENTREZID)) %>% transmute(target_gene = input_name, ENTREZID, target_mapping_type = mapping_type)

ICGC.source_annotation <- ICGC.expr_annotation %>% filter(!is.na(ENTREZID)) %>% transmute(source_gene = input_name, 
    ENTREZID, source_mapping_type = mapping_type)

ICGC.symbol_mapping_audit <- as.data.frame(merge(ICGC.target_annotation, ICGC.source_annotation, by = "ENTREZID", 
    sort = FALSE)) %>% left_join(ICGC.expr_id_counts, by = "ENTREZID") %>% left_join(ICGC.background_id_counts, 
    by = "ENTREZID") %>% mutate(source_already_in_background = source_gene %in% ICGC.background_genes, 
    excluded_compound_name = source_gene == "RAB4B-EGLN2", eligible_one_to_one = n_expression_names == 
        1L & n_background_names == 1L & !source_already_in_background & !excluded_compound_name) %>% 
    arrange(target_gene, source_gene)

ICGC.symbol_mapping_candidates <- ICGC.symbol_mapping_audit %>% filter(eligible_one_to_one)

ICGC.mapping_QC <- data.frame(original_background_present = length(ICGC.background_overlap), candidate_recovered = n_distinct(ICGC.symbol_mapping_candidates$target_gene), 
    training_background = length(ICGC.background_genes))

ICGC.mapping_QC$projected_coverage <- with(ICGC.mapping_QC, (original_background_present + candidate_recovered)/training_background)

ICGC.BA4_missing_genes <- setdiff(unique(unlist(ICGC.BA4_gene_sets, use.names = FALSE)), rownames(ICGC.BA4_expr_input))

ICGC.BA4_mapping_audit <- data.frame(target_gene = ICGC.BA4_missing_genes) %>% left_join(ICGC.symbol_mapping_audit, 
    by = "target_gene")

library(dplyr)

library(tibble)

ICGC.symbol_mapping_used <- as.data.frame(ICGC.symbol_mapping_candidates)

stopifnot(nrow(ICGC.symbol_mapping_used) > 0, all(ICGC.symbol_mapping_used$eligible_one_to_one), !anyNA(ICGC.symbol_mapping_used$source_gene), 
    !anyNA(ICGC.symbol_mapping_used$target_gene), !anyDuplicated(ICGC.symbol_mapping_used$source_gene), 
    !anyDuplicated(ICGC.symbol_mapping_used$target_gene), all(ICGC.symbol_mapping_used$source_gene %in% 
        rownames(ICGC.expr)), !any(ICGC.symbol_mapping_used$target_gene %in% rownames(ICGC.expr)))

ICGC.expr_harmonized <- ICGC.expr

ICGC.rename_index <- match(ICGC.symbol_mapping_used$source_gene, rownames(ICGC.expr_harmonized))

ICGC.harmonized_names <- rownames(ICGC.expr_harmonized)

ICGC.harmonized_names[ICGC.rename_index] <- ICGC.symbol_mapping_used$target_gene

stopifnot(!anyDuplicated(ICGC.harmonized_names))

rownames(ICGC.expr_harmonized) <- ICGC.harmonized_names

stopifnot(identical(unname(ICGC.expr_harmonized), unname(ICGC.expr)), identical(colnames(ICGC.expr_harmonized), 
    ICGC.clinical$Patient_ID))

ICGC.background_overlap <- intersect(ICGC.background_genes, rownames(ICGC.expr_harmonized))

ICGC.background_missing <- setdiff(ICGC.background_genes, rownames(ICGC.expr_harmonized))

ICGC.background_QC <- data.frame(training_background = length(ICGC.background_genes), ICGC_present = length(ICGC.background_overlap), 
    missing = length(ICGC.background_missing), coverage = length(ICGC.background_overlap)/length(ICGC.background_genes))

ICGC.log2_tpm <- log2(ICGC.expr_harmonized + 1)

ICGC.background_present_ordered <- ICGC.background_genes[ICGC.background_genes %in% rownames(ICGC.log2_tpm)]

ICGC.log2_tpm_background <- ICGC.log2_tpm[ICGC.background_present_ordered, , drop = FALSE]

ICGC.BA4_expr_input <- ICGC.log2_tpm_background

stopifnot(all(is.finite(ICGC.BA4_expr_input)))

ICGC.BA4_gene_coverage <- bind_rows(lapply(names(ICGC.BA4_gene_sets), function(gs) {
    genes <- unique(ICGC.BA4_gene_sets[[gs]])
    present <- intersect(genes, rownames(ICGC.BA4_expr_input))
    data.frame(gene_set = gs, total_genes = length(genes), present_in_model_input = length(present), 
        input_coverage = length(present)/length(genes), missing_from_input = paste(setdiff(genes, present), 
            collapse = ";"))
}))

ICGC.BA4_zero_gene_QC <- bind_rows(lapply(names(ICGC.BA4_gene_sets), function(gs) {
    genes <- intersect(unique(ICGC.BA4_gene_sets[[gs]]), rownames(ICGC.BA4_expr_input))
    x <- ICGC.BA4_expr_input[genes, , drop = FALSE]
    constant <- if (length(genes) > 0L) {
        genes[apply(x, 1, function(z) all(z == z[1]))]
    }
    else {
        character(0)
    }
    zero_all <- genes[rowSums(x > 0) == 0]
    data.frame(gene_set = gs, n_genes = length(genes), n_all_zero = length(zero_all), all_zero_genes = paste(zero_all, 
        collapse = ";"), n_constant = length(constant), constant_genes = paste(constant, collapse = ";"))
}))

ICGC.BA4_gene_sets_available <- lapply(ICGC.BA4_gene_sets, function(gs) {
    intersect(unique(gs), rownames(ICGC.BA4_expr_input))
})

stopifnot(all(lengths(ICGC.BA4_gene_sets_available) >= 2L))

ICGC.ssgsea_BA4 <- run_ssgsea(expr_matrix = ICGC.BA4_expr_input, gene_sets = ICGC.BA4_gene_sets_available)

stopifnot(setequal(rownames(ICGC.ssgsea_BA4), ICGC.BA4_feature_names), setequal(colnames(ICGC.ssgsea_BA4), 
    colnames(ICGC.expr)))

ICGC.ssgsea_BA4 <- ICGC.ssgsea_BA4[ICGC.BA4_feature_names, colnames(ICGC.expr), drop = FALSE]

stopifnot(all(is.finite(ICGC.ssgsea_BA4)))

ICGC.BA4_features <- as.data.frame(t(ICGC.ssgsea_BA4), check.names = FALSE) %>% rownames_to_column("Patient_ID")

ICGC.BA4_features$BA4_ICBScore <- predict_BA4_ICBScore(newdata = ICGC.BA4_features, model = BA4_MODEL_BUNDLE$model, 
    features = BA4_MODEL_BUNDLE$features, lambda = BA4_MODEL_BUNDLE$lambda)

stopifnot(all(is.finite(ICGC.BA4_features$BA4_ICBScore)), all(ICGC.BA4_features$BA4_ICBScore >= 0), all(ICGC.BA4_features$BA4_ICBScore <= 
    1))

ICGC.analysis <- ICGC.BA4_features %>% left_join(ICGC.CTNNB1_status, by = "Patient_ID") %>% left_join(ICGC.clinical, 
    by = "Patient_ID")

stopifnot(nrow(ICGC.analysis) == ncol(ICGC.expr), !anyDuplicated(ICGC.analysis$Patient_ID), identical(ICGC.analysis$Patient_ID, 
    colnames(ICGC.expr)), !anyNA(ICGC.analysis$CTNNB1_mut))

ICGC.score_comparison <- ICGC.analysis_before_symbol_check %>% transmute(Patient_ID, score_before = BA4_ICBScore) %>% 
    inner_join(ICGC.analysis %>% transmute(Patient_ID, score_after = BA4_ICBScore), by = "Patient_ID") %>% 
    mutate(difference = score_after - score_before)

stopifnot(nrow(ICGC.score_comparison) == nrow(ICGC.analysis))

ICGC.score_comparison_QC <- with(ICGC.score_comparison, data.frame(n = length(score_before), Pearson_r = cor(score_before, 
    score_after), Spearman_rho = cor(score_before, score_after, method = "spearman"), median_change = median(difference), 
    median_absolute_change = median(abs(difference)), max_absolute_change = max(abs(difference))))

ICGC.BA4_score_QC <- ICGC.analysis %>% group_by(CTNNB1_status) %>% summarise(n = n(), score_min = min(BA4_ICBScore), 
    score_Q1 = quantile(BA4_ICBScore, 0.25), score_median = median(BA4_ICBScore), score_mean = mean(BA4_ICBScore), 
    score_Q3 = quantile(BA4_ICBScore, 0.75), score_max = max(BA4_ICBScore), .groups = "drop")

with(ICGC.score_comparison, {
    plot(score_before, score_after, xlab = "BA4 ICBScore before symbol harmonization", ylab = "BA4 ICBScore after symbol harmonization", 
        pch = 16, col = "#336699", asp = 1)
    abline(a = 0, b = 1, lty = 2, col = "grey50")
})

library(dplyr)

library(survival)

library(survminer)

library(ggplot2)

ICGC.surv <- ICGC.analysis %>% dplyr::select(Patient_ID, BA4_ICBScore, CTNNB1_status, OS_time_months, 
    OS_status) %>% mutate(OS_year = OS_time_months/12)

ICGC.surv_clean <- ICGC.surv %>% filter(is.finite(BA4_ICBScore), is.finite(OS_time_months), OS_time_months > 
    0, OS_status %in% c(0, 1))

ICGC.survival_QC <- data.frame(total_patients = nrow(ICGC.surv), included_patients = nrow(ICGC.surv_clean), 
    excluded_patients = nrow(ICGC.surv) - nrow(ICGC.surv_clean), deaths = sum(ICGC.surv_clean$OS_status == 
        1), censored = sum(ICGC.surv_clean$OS_status == 0))

ICGC.score_SD <- sd(ICGC.surv_clean$BA4_ICBScore)

stopifnot(is.finite(ICGC.score_SD), ICGC.score_SD > 0)

ICGC.surv_clean <- ICGC.surv_clean %>% mutate(ICBScore_per01 = BA4_ICBScore/0.1, ICBScore_z = as.numeric(scale(BA4_ICBScore)))

ICGC.cox_ICB_01 <- survival::coxph(Surv(OS_time_months, OS_status) ~ ICBScore_per01, data = ICGC.surv_clean, 
    ties = "efron", x = TRUE)

ICGC.cox_ICB_z <- survival::coxph(Surv(OS_time_months, OS_status) ~ ICBScore_z, data = ICGC.surv_clean, 
    ties = "efron", x = TRUE)

ICGC.ICBScore_median <- median(ICGC.surv_clean$BA4_ICBScore)

ICGC.surv_clean <- ICGC.surv_clean %>% mutate(ICBScore_group = factor(ifelse(BA4_ICBScore >= ICGC.ICBScore_median, 
    "High", "Low"), levels = c("Low", "High")))

ICGC.survival_group_QC <- ICGC.surv_clean %>% group_by(ICBScore_group) %>% summarise(n = n(), deaths = sum(OS_status == 
    1), censored = sum(OS_status == 0), score_min = min(BA4_ICBScore), score_max = max(BA4_ICBScore), 
    .groups = "drop")

stopifnot(all(table(ICGC.surv_clean$ICBScore_group) > 0))

ICGC.cox_ICB_group <- survival::coxph(Surv(OS_time_months, OS_status) ~ ICBScore_group, data = ICGC.surv_clean, 
    ties = "efron", x = TRUE)


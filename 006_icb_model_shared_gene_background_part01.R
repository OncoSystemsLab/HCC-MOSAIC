if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

datasets <- list(PRJEB34724 = PRJEB34724.tpm, GSE181946 = gse181946.tpm, GSE202069 = gse202069.tpm)

for (nm in names(datasets)) {
    x <- datasets[[nm]]
    cat("\n====================")
    cat(nm, "")
    cat("Genes:", nrow(x), "")
    cat("Samples:", ncol(x) - 1, "")
    cat("Duplicated Symbol:", sum(duplicated(x$Symbol)), "")
    cat("NA Symbol:", sum(is.na(x$Symbol)), "")
    cat("Empty Symbol:", sum(trimws(x$Symbol) == "", na.rm = TRUE), "")
    cat("NA expression:", sum(is.na(x[, -1])), "")
}

for (nm in names(datasets)) {
    x <- datasets[[nm]]
    expr <- as.matrix(x[, -1])
    storage.mode(expr) <- "numeric"
    sums <- colSums(expr, na.rm = TRUE)
    cat("", nm, "")
    print(summary(sums))
}

library(dplyr)

prepare_tpm <- function(df, dataset_name) {
    df <- as.data.frame(df, check.names = FALSE)
    if (!"Symbol" %in% colnames(df)) {
        stop(dataset_name, ": missing Symbol column")
    }
    df$Symbol <- trimws(as.character(df$Symbol))
    df <- df[!is.na(df$Symbol) & df$Symbol != "", , drop = FALSE]
    sample_cols <- setdiff(colnames(df), "Symbol")
    for (cc in sample_cols) {
        df[[cc]] <- as.numeric(df[[cc]])
    }
    if (anyDuplicated(df$Symbol) > 0) {
        message(dataset_name, ": found ", sum(duplicated(df$Symbol)), " duplicate symbols; summing expression by Symbol")
        df <- df %>% group_by(Symbol) %>% summarise(across(all_of(sample_cols), ~sum(.x, na.rm = TRUE)), 
            .groups = "drop")
    }
    mat <- as.matrix(df[, sample_cols, drop = FALSE])
    rownames(mat) <- df$Symbol
    storage.mode(mat) <- "numeric"
    if (anyNA(mat)) {
        stop(dataset_name, ": expression matrix contains NA values")
    }
    if (any(mat < 0)) {
        stop(dataset_name, ": TPM contains negative values")
    }
    if (anyDuplicated(rownames(mat)) > 0) {
        stop(dataset_name, ": duplicate symbols remain")
    }
    return(mat)
}

tpm_list <- list(PRJEB34724 = prepare_tpm(PRJEB34724.tpm, "PRJEB34724"), GSE181946 = prepare_tpm(gse181946.tpm, 
    "GSE181946"), GSE202069 = prepare_tpm(gse202069.tpm, "GSE202069"))

sapply(tpm_list, dim)

sapply(tpm_list, nrow)

common_genes <- Reduce(intersect, lapply(tpm_list, rownames))

length(common_genes)

common_genes <- sort(common_genes)

tpm_common <- lapply(tpm_list, function(x) {
    x[common_genes, , drop = FALSE]
})

sapply(tpm_common, nrow)

stopifnot(identical(rownames(tpm_common$PRJEB34724), rownames(tpm_common$GSE181946)), identical(rownames(tpm_common$PRJEB34724), 
    rownames(tpm_common$GSE202069)))

zero_summary <- data.frame(Dataset = names(tpm_common), All_zero_genes = sapply(tpm_common, function(x) {
    sum(rowSums(x > 0) == 0)
}))

zero_summary

expressed_each_cohort <- lapply(tpm_common, function(x) rowSums(x > 0) > 0)

all_cohort_expressed <- Reduce("&", expressed_each_cohort)

table(all_cohort_expressed)

background_genes <- common_genes[all_cohort_expressed]

tpm_background <- lapply(tpm_common, function(x) {
    x[background_genes, , drop = FALSE]
})

rm(tpm_list, tpm_common, zero_summary, x, expr, expressed_each_cohort)


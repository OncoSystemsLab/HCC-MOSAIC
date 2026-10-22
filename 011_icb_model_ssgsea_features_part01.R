if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

log2_tpm_background <- lapply(tpm_background, function(x) {
    log2(x + 1)
})

sapply(log2_tpm_background, dim)

stopifnot(all(sapply(log2_tpm_background, function(x) all(is.finite(x)))))

run_ssgsea <- function(expr_matrix, gene_sets) {
    stopifnot(is.matrix(expr_matrix))
    param <- ssgseaParam(exprData = expr_matrix, geneSets = gene_sets, minSize = 2, maxSize = 500, alpha = 0.25, 
        normalize = FALSE)
    es <- gsva(param, verbose = FALSE, BPPARAM = SerialParam())
    return(es)
}

ssgsea_list <- list()

for (nm in names(log2_tpm_background)) {
    cat("\nRunning ssGSEA:", nm, "")
    ssgsea_list[[nm]] <- run_ssgsea(expr_matrix = log2_tpm_background[[nm]], gene_sets = all_gene_sets)
    cat("Finished:", nm, " | Dimension:", paste(dim(ssgsea_list[[nm]]), collapse = " × "), "")
}

sapply(ssgsea_list, dim)

feature_order <- names(all_gene_sets)

ssgsea_list <- lapply(ssgsea_list, function(x) {
    missing_features <- setdiff(feature_order, rownames(x))
    if (length(missing_features) > 0) {
        stop("Missing ssGSEA features: ", paste(missing_features, collapse = ", "))
    }
    x[feature_order, , drop = FALSE]
})

stopifnot(identical(rownames(ssgsea_list$PRJEB34724), rownames(ssgsea_list$GSE181946)), identical(rownames(ssgsea_list$PRJEB34724), 
    rownames(ssgsea_list$GSE202069)))

score_list <- lapply(ssgsea_list, function(x) {
    score <- t(x)
    score <- as.data.frame(score, check.names = FALSE)
    return(score)
})

sapply(score_list, dim)

score_all <- bind_rows(lapply(names(score_list), function(nm) {
    x <- score_list[[nm]]
    x$sample_id <- rownames(x)
    x$cohort <- nm
    x
}))

score_all <- score_all %>% select(cohort, sample_id, everything())

dim(score_all)

table(score_all$cohort)

score_all <- score_all %>% mutate(patient_uid = paste(cohort, sample_id, sep = "__")) %>% select(patient_uid, 
    cohort, sample_id, everything())

anyDuplicated(score_all$patient_uid)

score_matrix <- score_all %>% select(all_of(feature_order))

sum(is.na(as.matrix(score_matrix)))

sum(!is.finite(as.matrix(score_matrix)))


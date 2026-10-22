if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

for (j in seq_len(ncol(X))) {
    if (registry$Scale_Valid[j]) {
        z_matrix[, j] <- (X[, j] - registry$Mean[j])/registry$SD[j]
    }
}

stopifnot(identical(dim(percentile_matrix), c(69L, 89L)), identical(dim(z_matrix), c(69L, 89L)), all(is.finite(percentile_matrix)), 
    all(percentile_matrix > 0), all(percentile_matrix < 100), all(is.finite(z_matrix[, registry$Scale_Valid, 
        drop = FALSE])))

if (any(registry$Scale_Valid)) {
    z_means <- colMeans(z_matrix[, registry$Scale_Valid, drop = FALSE])
    stopifnot(max(abs(z_means)) < 1e-10)
}

saveRDS(percentile_matrix, file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_89percentiles.rds"))

saveRDS(z_matrix, file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_89Zscores.rds"))

make_reference_df <- function(mat) {
    data.frame(patient_uid = rownames(mat), Cohort = reference_map$Cohort, mat, check.names = FALSE, 
        row.names = NULL)
}

write.csv(make_reference_df(percentile_matrix), file.path(PROJECT_DIR, "04_processed", "Stage2_69patients_89percentiles.csv"), 
    row.names = FALSE)

write.csv(make_reference_df(z_matrix), file.path(PROJECT_DIR, "04_processed", "Stage2_69patients_89Zscores.csv"), 
    row.names = FALSE, na = "")

example_patient <- "PRJEB34724__S02-T-RNA"

selected_features <- c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
    "AXIN2", "GLUL", "BAAT", "CYP7A1")

stopifnot(example_patient %in% rownames(X), all(selected_features %in% colnames(X)))

example_profile <- data.frame(Feature = selected_features, Type = registry$Type[match(selected_features, 
    registry$Feature)], Raw_Value = as.numeric(X[example_patient, selected_features]), Reference_Percentile = as.numeric(percentile_matrix[example_patient, 
    selected_features]), Reference_Z = as.numeric(z_matrix[example_patient, selected_features]))

passed <- (nrow(X) == 69L && ncol(X) == 89L && nrow(reference_registry) == 89L && all(is.finite(percentile_matrix)) && 
    all(percentile_matrix > 0 & percentile_matrix < 100))

if (passed) {
    cat("\nStage 2.3A: PASS")
} else {
    cat("\nStage 2.3A: CHECK REQUIRED")
}

PROJECT_DIR <- "HCC_OmicsAgent"

reference_bundle <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE, 
    stringsAsFactors = FALSE)

X <- reference_bundle$reference_values

registry <- reference_bundle$reference_registry

stopifnot(identical(dim(X), c(69L, 89L)), identical(colnames(X), registry$Feature), all(rownames(X) %in% 
    reference_map$patient_uid))

idx <- match(rownames(X), reference_map$patient_uid)

stopifnot(!anyNA(idx))

cohort_levels <- c("PRJEB34724", "GSE181946", "GSE202069")

cohort <- factor(reference_map$Cohort[idx], levels = cohort_levels)

stopifnot(!anyNA(cohort), all(table(cohort) > 1), all(registry$Scale_Valid))

Z <- sweep(X, 2, registry$Mean, FUN = "-")

Z <- sweep(Z, 2, registry$SD, FUN = "/")

stopifnot(all(is.finite(Z)))

feature_results <- lapply(seq_len(ncol(X)), function(j) {
    raw_values <- X[, j]
    test <- stats::kruskal.test(raw_values ~ cohort)
    H <- as.numeric(test$statistic)
    n <- length(raw_values)
    k <- nlevels(cohort)
    epsilon2 <- max(0, min(1, (H - k + 1)/(n - k)))
    raw_medians <- vapply(cohort_levels, function(cc) {
        median(raw_values[cohort == cc])
    }, numeric(1))
    z_medians <- vapply(cohort_levels, function(cc) {
        median(Z[cohort == cc, j])
    }, numeric(1))
    data.frame(Feature = colnames(X)[j], Type = registry$Type[j], Module = registry$Module[j], KW_H = H, 
        KW_P = as.numeric(test$p.value), Epsilon2 = epsilon2, Median_PRJEB34724 = raw_medians[1], Median_GSE181946 = raw_medians[2], 
        Median_GSE202069 = raw_medians[3], Median_Spread_Z = max(z_medians) - min(z_medians), stringsAsFactors = FALSE)
})

cohort_audit <- do.call(rbind, feature_results)

rownames(cohort_audit) <- NULL

cohort_audit$FDR_89 <- p.adjust(cohort_audit$KW_P, method = "BH")

cohort_audit <- cohort_audit[order(cohort_audit$FDR_89, -cohort_audit$Epsilon2), , drop = FALSE]

rownames(cohort_audit) <- NULL

focus_features <- c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
    "AXIN2", "GLUL", "BAAT", "CYP7A1")

focus_audit <- cohort_audit[match(focus_features, cohort_audit$Feature), , drop = FALSE]

stopifnot(!anyNA(focus_audit$Feature))

write.csv(cohort_audit, file.path(PROJECT_DIR, "08_reports", "Stage2_3B_89feature_cohort_audit.csv"), 
    row.names = FALSE)

write.csv(focus_audit, file.path(PROJECT_DIR, "08_reports", "Stage2_3B_CTNNB1_BA_Immune_focus.csv"), 
    row.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

signature_features <- registry$Feature[registry$Type == "ssGSEA"]

signature_X <- X[, signature_features, drop = FALSE]

stopifnot(nrow(signature_X) == 69, ncol(signature_X) == 26, all(apply(signature_X, 2, sd) > 0))

pca <- prcomp(signature_X, center = TRUE, scale. = TRUE)

variance <- 100 * pca$sdev^2/sum(pca$sdev^2)

pca_df <- data.frame(patient_uid = rownames(signature_X), Cohort = as.character(cohort), PC1 = pca$x[, 
    1], PC2 = pca$x[, 2], row.names = NULL)

pc1_test <- kruskal.test(PC1 ~ Cohort, data = pca_df)

pc2_test <- kruskal.test(PC2 ~ Cohort, data = pca_df)

figure_file <- file.path(PROJECT_DIR, "06_figures", "Stage2_3B_signature_PCA.png")

colors <- setNames(grDevices::hcl.colors(length(cohort_levels), "Dark 3"), cohort_levels)

png(filename = figure_file, width = 1200, height = 900, res = 150)

plot(pca_df$PC1, pca_df$PC2, pch = 19, col = colors[pca_df$Cohort], cex = 1.2, xlab = paste0("PC1 (", 
    round(variance[1], 1), "%)"), ylab = paste0("PC2 (", round(variance[2], 1), "%)"), main = "HCC Molecular Signatures - Cohort PCA")

abline(h = 0, v = 0, lty = 3, col = "gray75")

legend("topright", legend = cohort_levels, col = colors[cohort_levels], pch = 19, bty = "n")

dev.off()

write.csv(pca_df, file.path(PROJECT_DIR, "04_processed", "Stage2_3B_signature_PCA_coordinates.csv"), 
    row.names = FALSE)

saveRDS(pca, file.path(PROJECT_DIR, "02_reference", "Stage2_3B_signature_PCA.rds"))

stopifnot(nrow(pca_df) == 69, all(is.finite(pca_df$PC1)), all(is.finite(pca_df$PC2)), file.exists(figure_file))

predict_HCC_fingerprint_v1 <- function(tpm, sample_uid = NULL, project_dir = "HCC_OmicsAgent") {
    ba4 <- readRDS(file.path(project_dir, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))
    config <- readRDS(file.path(project_dir, "00_config", "Stage2_Molecular_Fingerprint_config_v1.rds"))
    reference <- readRDS(file.path(project_dir, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))
    cohort_qc <- read.csv(file.path(project_dir, "08_reports", "Stage2_3B_89feature_cohort_audit.csv"), 
        check.names = FALSE, stringsAsFactors = FALSE)
    tpm <- as.matrix(tpm)
    stopifnot(ncol(tpm) == 1L, !is.null(rownames(tpm)), !is.null(colnames(tpm)), !anyDuplicated(rownames(tpm)))
    storage.mode(tpm) <- "double"
    stopifnot(all(is.finite(tpm)), all(tpm >= 0))
    if (is.null(sample_uid)) {
        sample_uid <- colnames(tpm)[1]
    }
    stopifnot(length(sample_uid) == 1L, !is.na(sample_uid), nzchar(sample_uid))
    bg <- config$background_genes
    missing_bg <- setdiff(bg, rownames(tpm))
    if (length(missing_bg) > 0) {
        stop("Prediction requires the complete reference gene background.", "Currently missing ", length(missing_bg), " genes; examples: ", 
            paste(head(missing_bg, 10), collapse = ", "))
    }
    tpm_bg <- tpm[bg, , drop = FALSE]
    expr <- log2(tpm_bg + 1)
    signature_names <- names(config$signature_sets)
    pre <- config$preprocessing
    param <- GSVA::ssgseaParam(exprData = expr, geneSets = config$signature_sets, minSize = pre$ssgsea_minSize, 
        maxSize = pre$ssgsea_maxSize, alpha = pre$ssgsea_alpha, normalize = pre$ssgsea_normalize)
    es <- suppressWarnings(suppressMessages(GSVA::gsva(param, verbose = FALSE, BPPARAM = BiocParallel::SerialParam())))
    stopifnot(all(signature_names %in% rownames(es)))
    signature_values <- setNames(as.numeric(es[signature_names, 1]), signature_names)
    mechanism_genes <- unique(c(unlist(config$gene_panels, use.names = FALSE), "CTNNB1"))
    gene_values <- setNames(log2(as.numeric(tpm_bg[mechanism_genes, 1]) + 1), mechanism_genes)
    molecular_values <- c(signature_values, gene_values)
    ref_X <- reference$reference_values
    registry <- reference$reference_registry
    stopifnot(length(molecular_values) == 89L, identical(names(molecular_values), colnames(ref_X)), identical(colnames(ref_X), 
        registry$Feature), all(is.finite(molecular_values)))
    ba4_features <- ba4$features
    ba4_input <- matrix(as.numeric(signature_values[ba4_features]), nrow = 1, dimnames = list(sample_uid, 
        ba4_features))
    icb_score <- as.numeric(predict(ba4$model, newx = ba4_input, s = ba4$lambda, type = "response"))
    percentile <- vapply(seq_along(molecular_values), function(j) {
        r <- ref_X[, j]
        v <- molecular_values[j]
        100 * (sum(r < v) + 0.5 * sum(r == v))/length(r)
    }, numeric(1))
    z_score <- rep(NA_real_, length(molecular_values))
    valid <- registry$Scale_Valid
    z_score[valid] <- (molecular_values[valid] - registry$Mean[valid])/registry$SD[valid]
    outside_reference <- vapply(seq_along(molecular_values), function(j) {
        v <- molecular_values[j]
        r <- ref_X[, j]
        v < min(r) || v > max(r)
    }, logical(1))
    audit_idx <- match(names(molecular_values), cohort_qc$Feature)
    stopifnot(!anyNA(audit_idx))
    cohort_fdr <- cohort_qc$FDR_89[audit_idx]
    cohort_sensitive <- cohort_fdr < 0.05
    profile <- data.frame(Feature = names(molecular_values), Type = registry$Type, Module = registry$Module, 
        Raw_Value = unname(molecular_values), Reference_Percentile = unname(percentile), Reference_Z = unname(z_score), 
        Outside_Reference_Range = outside_reference, Cohort_FDR = cohort_fdr, Cohort_Associated = cohort_sensitive, 
        check.names = FALSE, row.names = NULL)
    summary <- data.frame(Sample_ID = sample_uid, ICBScore = icb_score, N_Signatures = length(signature_values), 
        N_Genes = length(gene_values), N_Features = nrow(profile), N_Cohort_Associated = sum(cohort_sensitive), 
        N_Outside_Reference = sum(outside_reference), Self_In_Reference = sample_uid %in% rownames(ref_X), 
        Reference_Type = reference$reference_type)
    list(summary = summary, profile = profile, signature_values = signature_values, gene_logTPM = gene_values, 
        input_qc = data.frame(Sample_ID = sample_uid, Input_Genes = nrow(tpm), Matched_Background = length(bg), 
            Missing_Background = 0L, Extra_Genes = nrow(tpm) - length(bg)))
}

PROJECT_DIR <- "HCC_OmicsAgent"

tpm_reference <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds"))

bg <- BA4_BUNDLE$background_genes

sample_id <- "S02-T-RNA"

sample_uid <- "PRJEB34724__S02-T-RNA"

one_tpm <- as.matrix(tpm_reference[["PRJEB34724"]][bg, sample_id, drop = FALSE])

fingerprint <- predict_HCC_fingerprint_v1(tpm = one_tpm, sample_uid = sample_uid, project_dir = PROJECT_DIR)

reference <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_Molecular_reference_69_v1.rds"))

ref_profile <- reference$reference_values[sample_uid, , drop = TRUE]

profile <- fingerprint$profile

feature_error <- max(abs(profile$Raw_Value - as.numeric(ref_profile[profile$Feature])))

reference_map <- read.csv(file.path(PROJECT_DIR, "08_reports", "Stage1_5A_cohort_inventory.csv"), check.names = FALSE)

expected_score <- reference_map$BA4_ReproducedScore[match(sample_uid, reference_map$patient_uid)]

score_error <- abs(fingerprint$summary$ICBScore - expected_score)

old_percentiles <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_89percentiles.rds"))

old_zscores <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage2_69patients_89Zscores.rds"))

percentile_error <- max(abs(profile$Reference_Percentile - as.numeric(old_percentiles[sample_uid, profile$Feature])))

z_error <- max(abs(profile$Reference_Z - as.numeric(old_zscores[sample_uid, profile$Feature])))

stopifnot(nrow(profile) == 89L, feature_error < 1e-08, score_error < 1e-08, percentile_error < 1e-08, 
    z_error < 1e-08)

focus_features <- c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_INTERFERON_GAMMA_RESPONSE", 
    "AXIN2", "GLUL", "BAAT", "SLC27A5", "CYP7A1", "CYP8B1", "NR1H4", "CD8A", "GZMB")

focus_profile <- profile[match(focus_features, profile$Feature), c("Feature", "Raw_Value", "Reference_Percentile", 
    "Reference_Z", "Cohort_Associated")]

safe_id <- gsub("[^A-Za-z0-9_-]", "_", sample_uid)

write.csv(profile, file.path(PROJECT_DIR, "05_results", paste0("Stage2_4_", safe_id, "_profile.csv")), 
    row.names = FALSE)

write.csv(fingerprint$summary, file.path(PROJECT_DIR, "05_results", paste0("Stage2_4_", safe_id, "_summary.csv")), 
    row.names = FALSE)

render_HCC_fingerprint_v1 <- function(fingerprint, project_dir = "HCC_OmicsAgent") {
    profile <- fingerprint$profile
    summary <- fingerprint$summary
    required_cols <- c("Feature", "Reference_Percentile", "Reference_Z", "Cohort_Associated", "Outside_Reference_Range")
    stopifnot(is.data.frame(profile), nrow(profile) == 89L, all(required_cols %in% colnames(profile)), 
        !anyDuplicated(profile$Feature), nrow(summary) == 1L, all(is.finite(profile$Reference_Percentile)), 
        all(profile$Reference_Percentile >= 0), all(profile$Reference_Percentile <= 100), all(is.finite(profile$Reference_Z)))
    uid <- as.character(summary$Sample_ID[1])
    safe_id <- gsub("[^A-Za-z0-9_-]", "_", uid)
    figure_dir <- file.path(project_dir, "06_figures")
    report_dir <- file.path(project_dir, "08_reports")
    dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
    display_modules <- list(`Driver and Wnt` = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "AXIN2", "GLUL"), 
        `Bile Acid Metabolism` = c("HALLMARK_BILE_ACID_METABOLISM", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", 
            "NR1H4"), `Immune Ecosystem` = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "CD8A", "GZMB", "LAG3"))
    display_labels <- c(HALLMARK_WNT_BETA_CATENIN_SIGNALING = "Wnt / beta-catenin", HALLMARK_BILE_ACID_METABOLISM = "Bile acid hallmark", 
        HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma response", AXIN2 = "AXIN2", GLUL = "GLUL", CYP7A1 = "CYP7A1", 
        CYP8B1 = "CYP8B1", BAAT = "BAAT", SLC27A5 = "SLC27A5", NR1H4 = "NR1H4", CD8A = "CD8A", GZMB = "GZMB", 
        LAG3 = "LAG3")
    display_features <- unlist(display_modules, use.names = FALSE)
    stopifnot(all(display_features %in% profile$Feature), all(display_features %in% names(display_labels)))
    figure_file <- file.path(figure_dir, paste0("Stage2_4_", safe_id, "_fingerprint.png"))
    module_colors <- c(`Driver and Wnt` = "#4E79A7", `Bile Acid Metabolism` = "#E39A43", `Immune Ecosystem` = "#4CA49A")
    grDevices::png(filename = figure_file, width = 1600, height = 1600, res = 160)
    old_par <- par(no.readonly = TRUE)
    par(mfrow = c(3, 1), mar = c(3.5, 13, 3.5, 4))
    for (module in names(display_modules)) {
        feature_names <- display_modules[[module]]
        idx <- match(feature_names, profile$Feature)
        df <- profile[idx, , drop = FALSE]
        p <- rev(df$Reference_Percentile)
        labels <- rev(paste0(unname(display_labels[feature_names]), ifelse(df$Cohort_Associated, " *", 
            "")))
        mids <- barplot(p, names.arg = labels, horiz = TRUE, las = 1, col = module_colors[module], border = NA, 
            xlim = c(0, 115), xlab = "Reference percentile", main = module, cex.names = 0.95, cex.main = 1.15, 
            axes = FALSE)
        axis(1, at = seq(0, 100, by = 20))
        abline(v = c(20, 80), lty = 3, col = "gray55")
        text(x = p + 2.5, y = mids, labels = paste0("P", sprintf("%.1f", p)), cex = 0.85, pos = 4, xpd = NA)
        box(bty = "l")
    }
    par(old_par)
    dev.off()
    focus_profile <- profile[match(display_features, profile$Feature), , drop = FALSE]
    focus_profile$Relative_Position <- ifelse(focus_profile$Reference_Percentile < 20, "Relatively low in the reference cohort", 
        ifelse(focus_profile$Reference_Percentile >= 80, "Relatively high in the reference cohort", "Within the middle range of the reference cohort"))
    get_p <- function(feature) {
        i <- match(feature, profile$Feature)
        if (is.na(i)) 
            stop("Missing feature: ", feature)
        profile$Reference_Percentile[i]
    }
    position <- function(p) {
        if (p < 20) 
            return("relatively low")
        if (p >= 80) 
            return("relatively high")
        return("within the middle range")
    }
    wnt_p <- get_p("HALLMARK_WNT_BETA_CATENIN_SIGNALING")
    ba_p <- get_p("HALLMARK_BILE_ACID_METABOLISM")
    ifng_p <- get_p("HALLMARK_INTERFERON_GAMMA_RESPONSE")
    baat_p <- get_p("BAAT")
    slc27a5_p <- get_p("SLC27A5")
    cyp7a1_p <- get_p("CYP7A1")
    cd8a_p <- get_p("CD8A")
    gzmb_p <- get_p("GZMB")
    interpretation <- c("HCC Molecular Fingerprint v1.0", "======================================", paste0("Sample: ", 
        uid), sprintf("BA4 ICBScore: %.6f", summary$ICBScore), "Reference: 69 HCC development-cohort samples", 
        "", "[1] Driver / Wnt-related transcription", sprintf("Wnt/beta-catenin: P%.1f (%s).", wnt_p, 
            position(wnt_p)), "This score describes Wnt-related transcriptional features, not CTNNB1 mutation status.", 
        "", "[2] Bile acid-related transcription", sprintf("Hallmark BA: P%.1f (%s).", ba_p, position(ba_p)), 
        sprintf("BAAT P%.1f; SLC27A5 P%.1f; CYP7A1 P%.1f.", baat_p, slc27a5_p, cyp7a1_p), "These results describe bile-acid-related gene expression; they do not measure bile acid concentration or synthesis flux.", 
        "", "[3] Immune-related transcription", sprintf("IFN-gamma P%.1f (%s).", ifng_p, position(ifng_p)), 
        sprintf("CD8A P%.1f; GZMB P%.1f.", cd8a_p, gzmb_p), "Bulk RNA expression reflects both cell composition and cell state and is not a direct measure of immune-cell cytotoxicity.", 
        "", "[4] Reference quality and limitations", paste0("Cohort-associated features in registry: ", 
            summary$N_Cohort_Associated), paste0("Features outside reference range: ", summary$N_Outside_Reference), 
        paste0("Sample included in reference: ", summary$Self_In_Reference), "An asterisk (*) marks a cohort-related distribution difference, not an individual patient abnormality.", 
        "The 20th/80th percentile thresholds are exploratory display categories, not clinical diagnostic cutoffs.", "ICBScore is the frozen BA4 model output, not a calibrated probability of clinical ICB response.", 
        "Mechanistic and clinical interpretations require independent WES, metabolomic or ICB outcome validation.")
    focus_file <- file.path(report_dir, paste0("Stage2_4_", safe_id, "_focus.csv"))
    text_file <- file.path(report_dir, paste0("Stage2_4_", safe_id, "_interpretation.txt"))
    write.csv(focus_profile, focus_file, row.names = FALSE, fileEncoding = "UTF-8")
    con <- file(text_file, open = "wt", encoding = "UTF-8")
    writeLines(enc2utf8(interpretation), con = con)
    close(con)
    list(sample_uid = uid, figure_file = figure_file, focus_file = focus_file, text_file = text_file, 
        focus_profile = focus_profile, interpretation = interpretation)
}

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("fingerprint"), exists("render_HCC_fingerprint_v1"))

report <- render_HCC_fingerprint_v1(fingerprint = fingerprint, project_dir = PROJECT_DIR)

stopifnot(file.exists(report$figure_file), file.exists(report$focus_file), file.exists(report$text_file), 
    nrow(report$focus_profile) == 13L)

library(ggplot2)


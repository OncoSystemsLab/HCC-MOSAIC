if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

library(ggplot2)

CPTAC.forest_bundle <- readRDS(paste0("ICB.model/", "CPTAC_validation_TPM_v19_159/", "CPTAC_159_validation_results.rds"))

CPTAC.forest_features <- CPTAC.forest_bundle$BA4_features

CPTAC.forest_mut <- CPTAC.forest_bundle$mutation_matrix

stopifnot(is.function(compare_CTNNB1_features), is.function(BA4.bootstrap_features), nrow(CPTAC.forest_features) == 
    159L, !anyDuplicated(CPTAC.forest_features$Patient_ID), all(CPTAC.forest_features$Patient_ID %in% 
    colnames(CPTAC.forest_mut)), "CTNNB1" %in% rownames(CPTAC.forest_mut), all(BA4_MODEL_BUNDLE$features %in% 
    names(CPTAC.forest_features)))

TCGA.feature_comparison <- compare_CTNNB1_features(TCGA.BA4_features, TCGA.mut, "TCGA")

ICGC.feature_comparison <- compare_CTNNB1_features(ICGC.BA4_features, ICGC.mut, "ICGC")

CPTAC.feature_comparison <- compare_CTNNB1_features(CPTAC.forest_features, CPTAC.forest_mut, "CPTAC")

BA4.feature_comparison <- bind_rows(TCGA.feature_comparison, ICGC.feature_comparison, CPTAC.feature_comparison)

TCGA.feature_bootstrap <- BA4.bootstrap_features(TCGA.BA4_features, TCGA.mut, "TCGA", B = 2000, seed = 2026)

ICGC.feature_bootstrap <- BA4.bootstrap_features(ICGC.BA4_features, ICGC.mut, "ICGC", B = 2000, seed = 2027)

CPTAC.feature_bootstrap <- BA4.bootstrap_features(CPTAC.forest_features, CPTAC.forest_mut, "CPTAC", B = 2000, 
    seed = 2028)

BA4.cohort_colors <- c(TCGA = "#D73027", ICGC = "#2878B5", CPTAC = "#009E73")

BA4.feature_labels <- c(HALLMARK_BILE_ACID_METABOLISM = "Bile acid metabolism", HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid metabolism", 
    HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "ROS pathway", ICI__HCC_INFLAMMATORY4 = "Inflammatory4")

BA4.cohort_offsets <- c(TCGA = 0.23, ICGC = 0, CPTAC = -0.23)

BA4.forest_data <- bind_rows(TCGA.feature_bootstrap, ICGC.feature_bootstrap, CPTAC.feature_bootstrap) %>% 
    left_join(BA4.feature_comparison %>% select(cohort, feature, P_value, FDR), by = c("cohort", "feature")) %>% 
    mutate(cohort = factor(cohort, levels = c("TCGA", "ICGC", "CPTAC")), feature_label = factor(unname(BA4.feature_labels[feature]), 
        levels = c("Inflammatory4", "ROS pathway", "Fatty acid metabolism", "Bile acid metabolism")), 
        plot_y = as.numeric(feature_label) + unname(BA4.cohort_offsets[as.character(cohort)]), result_label = sprintf("%.2f [%.2f, %.2f]; FDR = %.2g", 
            effect, CI_lower, CI_upper, FDR))

stopifnot(nrow(BA4.forest_data) == 12L, !anyDuplicated(BA4.forest_data[, c("cohort", "feature")]), !anyNA(BA4.forest_data$feature_label), 
    all(is.finite(BA4.forest_data$effect)), all(is.finite(BA4.forest_data$CI_lower)), all(is.finite(BA4.forest_data$CI_upper)), 
    all(is.finite(BA4.forest_data$FDR)))

BA4.p_feature_forest <- ggplot(BA4.forest_data, aes(y = plot_y, color = cohort)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey55") + geom_segment(aes(x = CI_lower, xend = CI_upper, yend = plot_y), 
    linewidth = 0.9) + geom_segment(aes(x = CI_lower, xend = CI_lower, y = plot_y - 0.035, yend = plot_y + 
    0.035), linewidth = 0.7) + geom_segment(aes(x = CI_upper, xend = CI_upper, y = plot_y - 0.035, yend = plot_y + 
    0.035), linewidth = 0.8) + geom_point(aes(x = effect, fill = cohort), shape = 21, color = "black", 
    stroke = 0.5, size = 4) + geom_text(aes(x = 1.05, label = result_label), hjust = 0, size = 3.6, show.legend = FALSE) + 
    scale_color_manual(values = BA4.cohort_colors, breaks = c("TCGA", "ICGC", "CPTAC"), name = "Cohort") + 
    scale_fill_manual(values = BA4.cohort_colors, guide = "none") + scale_y_continuous(breaks = 1:4, 
    labels = levels(BA4.forest_data$feature_label), limits = c(0.55, 4.45), expand = expansion(mult = 0)) + 
    scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 2.55)) + labs(title = "CTNNB1 mutation-associated ICI model features", 
    subtitle = paste0("Points: rank-biserial effect; ", "intervals: bootstrap 95% CI"), x = "Rank-biserial correlation", 
    y = NULL) + theme_classic(base_size = 13) + theme(axis.text.x = element_text(size = 13, color = "black"), 
    axis.text.y = element_text(size = 13, color = "black"), axis.title.x = element_text(size = 15), axis.ticks.y = element_blank(), 
    plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11), plot.caption = element_text(hjust = 0, 
        size = 9), legend.position = "top")

BA4.visual_dir <- paste0("ICB.model/", "BA4_visualization")

ggsave(file.path(BA4.visual_dir, "CTNNB1_four_feature_forest_three_cohorts.pdf"), plot = BA4.p_feature_forest, 
    width = 8, height = 6, bg = "white")

data.table::fwrite(BA4.forest_data, file.path(BA4.visual_dir, "CTNNB1_four_feature_effects_three_cohorts.tsv"), 
    sep = "\t")

library(dplyr)

library(ggplot2)

stopifnot(requireNamespace("AnnotationDbi", quietly = TRUE), requireNamespace("org.Hs.eg.db", quietly = TRUE), 
    requireNamespace("data.table", quietly = TRUE))

Combined.out_dir <- paste0("ICB.model/", "Combined_TCGA_ICGC_CPTAC")

dir.create(Combined.out_dir, recursive = TRUE, showWarnings = FALSE)

Combined.B <- 2000L

Combined.plot_nominal_only <- TRUE

Combined.CPTAC_bundle <- readRDS(paste0("ICB.model/", "CPTAC_validation_TPM_v19_159/", "CPTAC_159_validation_results.rds"))

Combined.features_list <- list(TCGA = TCGA.BA4_features, ICGC = ICGC.BA4_features, CPTAC = Combined.CPTAC_bundle$BA4_features)

Combined.mutation_list <- list(TCGA = TCGA.mut, ICGC = ICGC.mut, CPTAC = Combined.CPTAC_bundle$mutation_matrix)

Combined.patient_tables <- list()

Combined.aligned_mutations <- list()

Combined.patient_QC <- list()

for (cc in names(Combined.features_list)) {
    dat <- Combined.features_list[[cc]]
    mat <- as.matrix(Combined.mutation_list[[cc]])
    stopifnot(all(c("Patient_ID", "BA4_ICBScore") %in% names(dat)), !is.null(rownames(mat)), !is.null(colnames(mat)), 
        !anyDuplicated(rownames(mat)), !anyDuplicated(colnames(mat)), !anyDuplicated(dat$Patient_ID), 
        !anyNA(dat$Patient_ID), is.numeric(dat$BA4_ICBScore), all(is.finite(mat)), all(mat %in% c(0, 
            1)))
    dat <- dat %>% transmute(Patient_ID = as.character(Patient_ID), ICI_Score = BA4_ICBScore)
    eligible <- is.finite(dat$ICI_Score) & dat$Patient_ID %in% colnames(mat)
    Combined.patient_QC[[cc]] <- data.frame(cohort = cc, n_score_table = nrow(dat), n_paired_included = sum(eligible), 
        n_excluded = sum(!eligible))
    dat <- dat[eligible, , drop = FALSE]
    stopifnot(nrow(dat) > 0L)
    mat <- mat[, dat$Patient_ID, drop = FALSE]
    dat <- dat %>% mutate(cohort = cc, patient_uid = paste(cc, Patient_ID, sep = "::"))
    colnames(mat) <- dat$patient_uid
    Combined.patient_tables[[cc]] <- dat
    Combined.aligned_mutations[[cc]] <- mat
}

Combined.patients <- bind_rows(Combined.patient_tables)

Combined.patient_QC <- bind_rows(Combined.patient_QC)

stopifnot(!anyDuplicated(Combined.patients$patient_uid), all(Combined.patients$ICI_Score >= 0), all(Combined.patients$ICI_Score <= 
    1))

Combined.db <- org.Hs.eg.db::org.Hs.eg.db

Combined.all_original_genes <- unique(unlist(lapply(Combined.aligned_mutations, rownames), use.names = FALSE))

Combined.official_symbols <- AnnotationDbi::keys(Combined.db, keytype = "SYMBOL")

Combined.alias_keys <- intersect(setdiff(Combined.all_original_genes, Combined.official_symbols), AnnotationDbi::keys(Combined.db, 
    keytype = "ALIAS"))

Combined.gene_map <- data.frame(original_gene = Combined.all_original_genes, harmonized_gene = Combined.all_original_genes, 
    stringsAsFactors = FALSE)

if (length(Combined.alias_keys)) {
    Combined.alias_symbols <- AnnotationDbi::mapIds(Combined.db, keys = Combined.alias_keys, column = "SYMBOL", 
        keytype = "ALIAS", multiVals = "asNA")
    Combined.alias_symbols <- Combined.alias_symbols[!is.na(Combined.alias_symbols)]
    Combined.map_index <- match(names(Combined.alias_symbols), Combined.gene_map$original_gene)
    Combined.gene_map$harmonized_gene[Combined.map_index] <- unname(Combined.alias_symbols)
}

for (cc in names(Combined.aligned_mutations)) {
    mat <- Combined.aligned_mutations[[cc]]
    genes <- Combined.gene_map$harmonized_gene[match(rownames(mat), Combined.gene_map$original_gene)]
    mat <- rowsum(mat, group = genes, reorder = FALSE)
    mat <- (mat > 0) * 1L
    stopifnot(!anyDuplicated(rownames(mat)))
    Combined.aligned_mutations[[cc]] <- mat
}

Combined.all_genes <- sort(unique(unlist(lapply(Combined.aligned_mutations, rownames), use.names = FALSE)))

Combined.gene_availability <- data.frame(Gene = Combined.all_genes)

for (cc in names(Combined.aligned_mutations)) {
    Combined.gene_availability[[cc]] <- Combined.all_genes %in% rownames(Combined.aligned_mutations[[cc]])
}

Combined.gene_availability$available_in_all_cohorts <- rowSums(Combined.gene_availability[, names(Combined.aligned_mutations), 
    drop = FALSE]) == length(Combined.aligned_mutations)

Combined.common_genes <- Combined.gene_availability %>% filter(available_in_all_cohorts) %>% pull(Gene)

if (!length(Combined.common_genes)) {
    stop("No shared genes have known mutation status in all three cohorts.")
}

Combined.mut <- do.call(cbind, lapply(Combined.aligned_mutations, function(mat) {
    mat[Combined.common_genes, , drop = FALSE]
}))

stopifnot(identical(colnames(Combined.mut), Combined.patients$patient_uid), !anyNA(Combined.mut))

Combined.mutation_frequency <- data.frame(Gene = rownames(Combined.mut), n_MUT = as.integer(rowSums(Combined.mut))) %>% 
    mutate(n_total = ncol(Combined.mut), n_WT = n_total - n_MUT, mutation_frequency = n_MUT/n_total) %>% 
    arrange(desc(mutation_frequency), Gene)

Combined.selected_genes <- Combined.mutation_frequency %>% filter(mutation_frequency > 0.05) %>% pull(Gene)

if (!length(Combined.selected_genes)) {
    stop("No genes have mutation frequency strictly above 5%.")
}

Combined.rank_biserial <- function(x, y) {
    nx <- length(x)
    ny <- length(y)
    if (!nx || !ny) 
        return(NA_real_)
    r <- rank(c(x, y), ties.method = "average")
    U <- sum(r[seq_len(nx)]) - nx * (nx + 1)/2
    2 * U/(nx * ny) - 1
}

Combined.compare_gene <- function(gene, B = 2000L) {
    mut <- as.vector(Combined.mut[gene, ] == 1L)
    score <- Combined.patients$ICI_Score
    cohort <- Combined.patients$cohort
    x <- score[mut]
    y <- score[!mut]
    nx <- length(x)
    ny <- length(y)
    testable <- nx > 0L && ny > 0L
    effect <- p_value <- NA_real_
    ci <- c(NA_real_, NA_real_)
    if (testable) {
        effect <- Combined.rank_biserial(x, y)
        p_value <- if (length(unique(score)) == 1L) {
            1
        }
        else {
            wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
        }
        if (nx >= 2L && ny >= 2L) {
            mut_strata <- split(which(mut), cohort[mut])
            wt_strata <- split(which(!mut), cohort[!mut])
            resample_indices <- function(strata) {
                unlist(lapply(strata, function(idx) {
                  idx[sample.int(length(idx), size = length(idx), replace = TRUE)]
                }), use.names = FALSE)
            }
            boot_effect <- replicate(B, {
                im <- resample_indices(mut_strata)
                iw <- resample_indices(wt_strata)
                Combined.rank_biserial(score[im], score[iw])
            })
            ci <- unname(quantile(boot_effect, probs = c(0.025, 0.975)))
        }
    }
    data.frame(Gene = gene, n_MUT = nx, n_WT = ny, mutation_frequency = nx/(nx + ny), mean_MUT = if (nx) 
        mean(x)
    else NA_real_, mean_WT = if (ny) 
        mean(y)
    else NA_real_, median_MUT = if (nx) 
        median(x)
    else NA_real_, median_WT = if (ny) 
        median(y)
    else NA_real_, delta_mean = if (testable) 
        mean(x) - mean(y)
    else NA_real_, delta_median = if (testable) 
        median(x) - median(y)
    else NA_real_, rank_biserial = effect, CI_lower = ci[1], CI_upper = ci[2], P_value = p_value, testable = testable)
}

set.seed(2026)

Combined.gene_BA4_effects <- bind_rows(lapply(Combined.selected_genes, Combined.compare_gene, B = Combined.B))

Combined.gene_BA4_effects$FDR <- NA_real_

Combined.tested <- which(Combined.gene_BA4_effects$testable)

Combined.gene_BA4_effects$FDR[Combined.tested] <- p.adjust(Combined.gene_BA4_effects$P_value[Combined.tested], 
    method = "BH")

Combined.gene_BA4_effects <- Combined.gene_BA4_effects %>% mutate(Direction = case_when(!testable ~ "Not testable", 
    rank_biserial > 0 ~ "Higher in MUT", rank_biserial < 0 ~ "Lower in MUT", TRUE ~ "No rank difference")) %>% 
    arrange(FDR, desc(abs(rank_biserial)))

Combined.frequency_by_cohort <- bind_rows(lapply(names(Combined.aligned_mutations), function(cc) {
    mat <- Combined.aligned_mutations[[cc]][Combined.selected_genes, , drop = FALSE]
    data.frame(cohort = cc, Gene = rownames(mat), n_MUT = as.integer(rowSums(mat)), n_total = ncol(mat), 
        mutation_frequency = rowMeans(mat))
}))

Combined.effect_plot_data <- Combined.gene_BA4_effects %>% filter(testable)

if (Combined.plot_nominal_only) {
    Combined.effect_plot_data <- Combined.effect_plot_data %>% filter(P_value < 0.05)
}

if (!nrow(Combined.effect_plot_data)) {
    stop("No genes meet the current plot filters; consider setting plot_nominal_only to FALSE.")
}

Combined.effect_plot_data <- Combined.effect_plot_data %>% arrange(desc(rank_biserial)) %>% mutate(Gene_label = factor(Gene, 
    levels = Gene), Significance = factor(ifelse(FDR < 0.05, "FDR < 0.05", "FDR >= 0.05"), levels = c("FDR < 0.05", 
    "FDR >= 0.05")))

Combined.p_gene_BA4_effect <- ggplot(Combined.effect_plot_data, aes(x = rank_biserial, y = Gene_label)) + 
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") + geom_segment(aes(x = CI_lower, 
    xend = CI_upper, yend = Gene_label), color = "grey55", linewidth = 0.8, na.rm = TRUE) + geom_point(aes(fill = Significance, 
    size = -log10(pmax(P_value, .Machine$double.xmin))), shape = 21, color = "black", stroke = 0.6) + 
    scale_fill_manual(name = NULL, values = c(`FDR < 0.05` = "#D73027", `FDR >= 0.05` = "#2878B5")) + 
    scale_size_continuous(name = expression(-log[10](P)), range = c(2, 7)) + scale_x_continuous(limits = c(-1, 
    1), breaks = seq(-1, 1, 0.5)) + labs(x = "Rank-biserial correlation (MUT vs WT)", y = NULL, title = "Mutation-associated differences in ICI Score", 
    subtitle = sprintf("TCGA + ICGC + CPTAC | n = %d | Mutation frequency > 5%%", nrow(Combined.patients)), 
    caption = paste0("Pooled, unadjusted association; ", "95% CI: cohort-by-mutation-stratified bootstrap.", 
        "Positive values indicate higher ICI Score in MUT.")) + theme_classic(base_size = 13) + theme(axis.text.x = element_text(size = 14, 
    color = "black"), axis.text.y = element_text(size = 14, color = "black", face = "italic"), axis.title.x = element_text(size = 16), 
    legend.position = "top", plot.caption = element_text(hjust = 0, size = 9))

Combined.write <- function(dat, filename) {
    data.table::fwrite(dat, file.path(Combined.out_dir, filename), sep = "\t", na = "NA")
}

Combined.write(Combined.patient_QC, "01_patient_QC.tsv")

Combined.write(Combined.patients, "02_paired_patients_and_ICI_scores.tsv")

Combined.write(Combined.gene_map, "03_gene_symbol_mapping.tsv")

Combined.write(Combined.gene_availability, "04_gene_availability_by_cohort.tsv")

Combined.write(Combined.mutation_frequency, "05_pooled_mutation_frequency.tsv")

Combined.write(Combined.frequency_by_cohort, "06_selected_gene_frequency_by_cohort.tsv")

Combined.write(Combined.gene_BA4_effects, "07_pooled_gene_ICI_effects.tsv")

Combined.write(Combined.effect_plot_data, "08_pooled_effect_plot_data.tsv")

ggsave(file.path(Combined.out_dir, "09_pooled_gene_ICI_effect_forest.pdf"), plot = Combined.p_gene_BA4_effect, 
    width = 6, height = max(5, 0.28 * nrow(Combined.effect_plot_data) + 2.5), bg = "white", limitsize = FALSE)

saveRDS(list(patient_QC = Combined.patient_QC, patients = Combined.patients, gene_availability = Combined.gene_availability, 
    mutation_frequency = Combined.mutation_frequency, frequency_by_cohort = Combined.frequency_by_cohort, 
    effects = Combined.gene_BA4_effects, plot_data = Combined.effect_plot_data, bootstrap_B = Combined.B, 
    cohort_adjusted = FALSE), file.path(Combined.out_dir, "Combined_mutation_ICI_analysis.rds"))


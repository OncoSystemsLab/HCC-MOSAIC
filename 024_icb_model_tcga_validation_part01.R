if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(data.table)

library(dplyr)

library(tibble)

library(stringr)

TCGA_zip <- "raw_data/TCGA_LIHC_final.zip"

TCGA_dir <- "raw_data/TCGA_LIHC_final"

if (!dir.exists(TCGA_dir)) {
    dir.create(TCGA_dir, recursive = TRUE)
}

unzip(TCGA_zip, exdir = TCGA_dir)

TCGA.expr.raw <- fread(file.path(TCGA_dir, "TCGA_LIHC_expression_TPM.tsv.gz"), data.table = FALSE)

TCGA.expr <- as.data.frame(TCGA.expr.raw)

rownames(TCGA.expr) <- TCGA.expr$Gene

TCGA.expr$Gene <- NULL

TCGA.expr <- as.matrix(TCGA.expr)

storage.mode(TCGA.expr) <- "numeric"

TCGA.mut.raw <- fread(file.path(TCGA_dir, "TCGA_LIHC_mutation_matrix.tsv.gz"), data.table = FALSE)

TCGA.mut <- as.data.frame(TCGA.mut.raw)

rownames(TCGA.mut) <- TCGA.mut$Gene

TCGA.mut$Gene <- NULL

TCGA.mut <- as.matrix(TCGA.mut)

storage.mode(TCGA.mut) <- "numeric"

driver_genes <- c("CTNNB1", "TP53", "AXIN1", "APC", "TERT")

driver_genes %in% rownames(TCGA.mut)

TCGA.driver_mutation_summary <- data.frame(Gene = driver_genes, n_mut = rowSums(TCGA.mut[driver_genes, 
    , drop = FALSE] > 0))

TCGA.driver_mutation_summary$frequency <- TCGA.driver_mutation_summary$n_mut/ncol(TCGA.mut)

TCGA.driver_mutation_summary

TCGA.clinical <- fread(file.path(TCGA_dir, "TCGA_LIHC_clinical.tsv"), data.table = FALSE)

TCGA.clinical %>% select(Patient_ID, OS_time, OS_status, ajcc_pathologic_stage, tumor_grade, sex_at_birth, 
    age_at_diagnosis) %>% head()

TCGA.expr_ids <- colnames(TCGA.expr)

TCGA.mut_ids <- colnames(TCGA.mut)

TCGA.clin_ids <- TCGA.clinical$Patient_ID

stopifnot(setequal(TCGA.expr_ids, TCGA.mut_ids))

stopifnot(setequal(TCGA.expr_ids, TCGA.clin_ids))

TCGA.mut <- TCGA.mut[, TCGA.expr_ids, drop = FALSE]

TCGA.clinical <- TCGA.clinical[match(TCGA.expr_ids, TCGA.clinical$Patient_ID), , drop = FALSE]

stopifnot(identical(colnames(TCGA.expr), colnames(TCGA.mut)))

stopifnot(identical(colnames(TCGA.expr), TCGA.clinical$Patient_ID))

exists("BA4_MODEL_BUNDLE")

TCGA.background_genes <- BA4_MODEL_BUNDLE$background_genes

TCGA.BA4_gene_sets <- BA4_MODEL_BUNDLE$gene_sets

length(TCGA.background_genes)

names(TCGA.BA4_gene_sets)

TCGA.background_overlap <- intersect(TCGA.background_genes, rownames(TCGA.expr))

TCGA.background_missing <- setdiff(TCGA.background_genes, rownames(TCGA.expr))

TCGA.background_QC <- data.frame(training_background = length(TCGA.background_genes), TCGA_present = length(TCGA.background_overlap), 
    missing = length(TCGA.background_missing), coverage = length(TCGA.background_overlap)/length(TCGA.background_genes))

TCGA.background_QC

TCGA.BA4_gene_coverage <- bind_rows(lapply(names(TCGA.BA4_gene_sets), function(gs) {
    genes <- unique(TCGA.BA4_gene_sets[[gs]])
    present <- intersect(genes, rownames(TCGA.expr))
    missing <- setdiff(genes, rownames(TCGA.expr))
    data.frame(gene_set = gs, total_genes = length(genes), present_genes = length(present), missing_genes = length(missing), 
        coverage = length(present)/length(genes), missing_gene_names = paste(missing, collapse = ";"))
}))

TCGA.BA4_gene_coverage

TCGA.log2_tpm <- log2(TCGA.expr + 1)

TCGA.background_present_ordered <- TCGA.background_genes[TCGA.background_genes %in% rownames(TCGA.log2_tpm)]

TCGA.log2_tpm_background <- TCGA.log2_tpm[TCGA.background_present_ordered, , drop = FALSE]

TCGA.BA4_zero_gene_QC <- bind_rows(lapply(names(TCGA.BA4_gene_sets), function(gs) {
    genes <- intersect(TCGA.BA4_gene_sets[[gs]], rownames(TCGA.expr))
    x <- TCGA.expr[genes, , drop = FALSE]
    zero_all <- genes[rowSums(x > 0) == 0]
    data.frame(gene_set = gs, n_genes = length(genes), n_all_zero = length(zero_all), all_zero_genes = paste(zero_all, 
        collapse = ";"))
}))

TCGA.BA4_zero_gene_QC

TCGA.BA4_expr_input <- TCGA.log2_tpm_background

library(GSVA)

library(BiocParallel)

library(dplyr)

library(tibble)

library(glmnet)

TCGA.BA4_gene_sets_available <- lapply(TCGA.BA4_gene_sets, function(gs) {
    intersect(gs, rownames(TCGA.BA4_expr_input))
})

sapply(TCGA.BA4_gene_sets_available, length)

TCGA.ssgsea_BA4 <- run_ssgsea(expr_matrix = TCGA.BA4_expr_input, gene_sets = TCGA.BA4_gene_sets_available)

TCGA.BA4_features <- as.data.frame(t(TCGA.ssgsea_BA4))

TCGA.BA4_features <- TCGA.BA4_features %>% rownames_to_column("Patient_ID")

TCGA.BA4_features$BA4_ICBScore <- predict_BA4_ICBScore(newdata = TCGA.BA4_features)

TCGA.CTNNB1_status <- data.frame(Patient_ID = colnames(TCGA.mut), CTNNB1_mut = as.integer(TCGA.mut["CTNNB1", 
    , drop = TRUE] > 0))

TCGA.CTNNB1_status <- TCGA.CTNNB1_status %>% mutate(CTNNB1_status = ifelse(CTNNB1_mut == 1, "CTNNB1_MUT", 
    "CTNNB1_WT"))

TCGA.analysis <- TCGA.BA4_features %>% left_join(TCGA.CTNNB1_status, by = "Patient_ID") %>% left_join(TCGA.clinical, 
    by = "Patient_ID")

library(dplyr)

library(survival)

library(survminer)

library(ggplot2)

TCGA.surv <- TCGA.analysis %>% select(Patient_ID, BA4_ICBScore, OS_time, OS_status) %>% mutate(OS_time = as.numeric(OS_time), 
    OS_status = as.numeric(OS_status), OS_year = OS_time/365)

TCGA.surv_clean <- TCGA.surv %>% filter(!is.na(BA4_ICBScore), !is.na(OS_time), !is.na(OS_status), OS_time > 
    0, OS_status %in% c(0, 1))

TCGA.surv_clean <- TCGA.surv_clean %>% mutate(ICBScore_per01 = BA4_ICBScore/0.1)

cox_ICB_01 <- coxph(Surv(OS_time, OS_status) ~ ICBScore_per01, data = TCGA.surv_clean)

summary(cox_ICB_01)

TCGA.surv_clean <- TCGA.surv_clean %>% mutate(ICBScore_z = as.numeric(scale(BA4_ICBScore)))

cox_ICB_z <- coxph(Surv(OS_time, OS_status) ~ ICBScore_z, data = TCGA.surv_clean)

summary(cox_ICB_z)

extract_cox_result <- function(model, variable, description) {
    s <- summary(model)
    data.frame(variable = variable, description = description, HR = s$coefficients[variable, "exp(coef)"], 
        CI_lower = s$conf.int[variable, "lower .95"], CI_upper = s$conf.int[variable, "upper .95"], P_value = s$coefficients[variable, 
            "Pr(>|z|)"])
}

TCGA.ICBScore_cox_results <- bind_rows(extract_cox_result(cox_ICB_01, "ICBScore_per01", "Per 0.1 increase in BA4 ICBScore"), 
    extract_cox_result(cox_ICB_z, "ICBScore_z", "Per 1-SD increase in BA4 ICBScore"))

TCGA.ICBScore_cox_results

TCGA.ICBScore_median <- median(TCGA.surv_clean$BA4_ICBScore, na.rm = TRUE)

TCGA.ICBScore_median

TCGA.surv_clean <- TCGA.surv_clean %>% mutate(ICBScore_group = ifelse(BA4_ICBScore >= TCGA.ICBScore_median, 
    "High", "Low"), ICBScore_group = factor(ICBScore_group, levels = c("Low", "High")))

cox_ICB_group <- coxph(Surv(OS_time, OS_status) ~ ICBScore_group, data = TCGA.surv_clean)

summary(cox_ICB_group)

fit_ICB_KM <- survfit(Surv(OS_year, OS_status) ~ ICBScore_group, data = TCGA.surv_clean)

p_ICB_KM <- ggsurvplot(fit_ICB_KM, data = TCGA.surv_clean, risk.table = TRUE, pval = TRUE, conf.int = FALSE, 
    censor = TRUE, xlab = "Overall survival (years)", ylab = "Overall survival probability", legend.title = "BA4 ICBScore", 
    legend.labs = c("Low", "High"), risk.table.height = 0.25, break.time.by = 2, ggtheme = theme_classic(base_size = 13))

p_ICB_KM

library(dplyr)

library(tibble)

TCGA.ICB_median <- median(TCGA.BA4_features$BA4_ICBScore, na.rm = TRUE)

TCGA.ICB_median

TCGA.ICB_group <- data.frame(Patient_ID = TCGA.BA4_features$Patient_ID, BA4_ICBScore = TCGA.BA4_features$BA4_ICBScore) %>% 
    mutate(ICB_group = ifelse(BA4_ICBScore >= TCGA.ICB_median, "High", "Low"), ICB_group = factor(ICB_group, 
        levels = c("Low", "High")))

TCGA.mut_BA4 <- TCGA.mut[, TCGA.ICB_group$Patient_ID, drop = FALSE]

TCGA.mutation_frequency <- data.frame(Gene = rownames(TCGA.mut_BA4), n_mut = rowSums(TCGA.mut_BA4 > 0))

TCGA.mutation_frequency <- TCGA.mutation_frequency %>% mutate(mutation_frequency = n_mut/ncol(TCGA.mut_BA4))

TCGA.mutation_frequency %>% arrange(desc(mutation_frequency)) %>% head(30)

TCGA.common_mutation_genes <- TCGA.mutation_frequency %>% filter(mutation_frequency >= 0.05) %>% pull(Gene)

TCGA.mutation_frequency %>% filter(Gene %in% TCGA.common_mutation_genes) %>% arrange(desc(mutation_frequency))

TCGA.low_ids <- TCGA.ICB_group %>% filter(ICB_group == "Low") %>% pull(Patient_ID)

TCGA.high_ids <- TCGA.ICB_group %>% filter(ICB_group == "High") %>% pull(Patient_ID)

test_mutation_enrichment <- function(gene) {
    mut_low <- sum(TCGA.mut_BA4[gene, TCGA.low_ids] > 0)
    wt_low <- length(TCGA.low_ids) - mut_low
    mut_high <- sum(TCGA.mut_BA4[gene, TCGA.high_ids] > 0)
    wt_high <- length(TCGA.high_ids) - mut_high
    tab <- matrix(c(mut_low, wt_low, mut_high, wt_high), nrow = 2, byrow = TRUE)
    rownames(tab) <- c("Low", "High")
    colnames(tab) <- c("MUT", "WT")
    fisher_res <- fisher.test(tab)
    data.frame(Gene = gene, mut_low = mut_low, n_low = length(TCGA.low_ids), freq_low = mut_low/length(TCGA.low_ids), 
        mut_high = mut_high, n_high = length(TCGA.high_ids), freq_high = mut_high/length(TCGA.high_ids), 
        odds_ratio = unname(fisher_res$estimate), P_value = fisher_res$p.value)
}

TCGA.mutation_enrichment <- bind_rows(lapply(TCGA.common_mutation_genes, test_mutation_enrichment))

TCGA.mutation_enrichment <- TCGA.mutation_enrichment %>% mutate(FDR = p.adjust(P_value, method = "BH"), 
    delta_frequency = freq_low - freq_high, enrichment_direction = case_when(delta_frequency > 0 ~ "Low_enriched", 
        delta_frequency < 0 ~ "High_enriched", TRUE ~ "Equal"))

TCGA.mutation_enrichment_low <- TCGA.mutation_enrichment %>% arrange(FDR, desc(delta_frequency))

TCGA.low_enriched_mutations <- TCGA.mutation_enrichment %>% filter(enrichment_direction == "Low_enriched") %>% 
    arrange(FDR, desc(delta_frequency))

TCGA.low_enriched_mutations %>% head(30)

TCGA.mutation_enrichment %>% filter(Gene == "CTNNB1")

TCGA.mutation_enrichment_ranked <- TCGA.mutation_enrichment %>% arrange(P_value) %>% mutate(rank_P = row_number())

TCGA.mutation_enrichment_ranked %>% filter(Gene == "CTNNB1")

TCGA.mutation_enrichment_ranked_FDR <- TCGA.mutation_enrichment %>% arrange(FDR) %>% mutate(rank_FDR = row_number())

TCGA.mutation_enrichment_ranked_FDR %>% filter(Gene == "CTNNB1")

TCGA.mutation_enrichment %>% filter(Gene %in% c("CTNNB1", "TP53", "AXIN1", "APC", "TERT")) %>% arrange(P_value)

TCGA.score_vector <- TCGA.ICB_group$BA4_ICBScore

names(TCGA.score_vector) <- TCGA.ICB_group$Patient_ID

test_mutation_score <- function(gene) {
    mut_status <- TCGA.mut_BA4[gene, names(TCGA.score_vector)] > 0
    score_mut <- TCGA.score_vector[mut_status]
    score_wt <- TCGA.score_vector[!mut_status]
    wt <- wilcox.test(score_mut, score_wt, exact = FALSE)
    data.frame(Gene = gene, n_MUT = length(score_mut), n_WT = length(score_wt), mean_MUT = mean(score_mut), 
        mean_WT = mean(score_wt), median_MUT = median(score_mut), median_WT = median(score_wt), delta_median = median(score_mut) - 
            median(score_wt), P_value = wt$p.value)
}

TCGA.mutation_score_assoc <- bind_rows(lapply(TCGA.common_mutation_genes, test_mutation_score))

TCGA.mutation_score_assoc <- TCGA.mutation_score_assoc %>% mutate(FDR = p.adjust(P_value, method = "BH"), 
    direction = case_when(delta_median < 0 ~ "Lower_ICBScore_in_MUT", delta_median > 0 ~ "Higher_ICBScore_in_MUT", 
        TRUE ~ "Equal"))

TCGA.mutation_score_assoc_ranked <- TCGA.mutation_score_assoc %>% arrange(P_value)

TCGA.low_score_mutations <- TCGA.mutation_score_assoc %>% filter(direction == "Lower_ICBScore_in_MUT") %>% 
    arrange(P_value)

TCGA.mutation_score_assoc %>% filter(Gene == "CTNNB1")

TCGA.mutation_score_assoc_ranked %>% mutate(rank_P = row_number()) %>% filter(Gene == "CTNNB1")

TCGA.mutation_score_assoc %>% filter(Gene %in% c("CTNNB1", "TP53", "AXIN1", "APC")) %>% arrange(P_value)

stop("Install the required R packages before running this analysis.")

library(rstatix)

TCGA.CTNNB1_score_df <- data.frame(score = TCGA.score_vector, mutation = ifelse(TCGA.mut_BA4["CTNNB1", 
    names(TCGA.score_vector)] > 0, "MUT", "WT"))

TCGA.CTNNB1_effect <- TCGA.CTNNB1_score_df %>% wilcox_effsize(score ~ mutation)

TCGA.CTNNB1_effect

TCGA.mutation_integrated <- TCGA.mutation_score_assoc %>% select(Gene, n_MUT, n_WT, mean_MUT, mean_WT, 
    median_MUT, median_WT, delta_median, P_score = P_value, FDR_score = FDR, score_direction = direction) %>% 
    left_join(TCGA.mutation_enrichment %>% select(Gene, mut_low, freq_low, mut_high, freq_high, odds_ratio, 
        delta_frequency, P_group = P_value, FDR_group = FDR, group_direction = enrichment_direction), 
        by = "Gene")

TCGA.mutation_integrated <- TCGA.mutation_integrated %>% mutate(direction_consistent = score_direction == 
    "Lower_ICBScore_in_MUT" & group_direction == "Low_enriched", significant_continuous = FDR_score < 
    0.05, nominal_group = P_group < 0.05)

TCGA.low_ICB_mutation_candidates <- TCGA.mutation_integrated %>% filter(direction_consistent) %>% arrange(FDR_score, 
    P_group)

TCGA.low_ICB_mutation_candidates %>% head(20)

TCGA.high_confidence_low_ICB <- TCGA.mutation_integrated %>% filter(FDR_score < 0.05, P_group < 0.05, 
    delta_median < 0, delta_frequency > 0) %>% arrange(FDR_score)

TCGA.high_confidence_low_ICB

TCGA.mutation_integrated <- TCGA.mutation_integrated %>% arrange(FDR_score, P_score) %>% mutate(rank_continuous = row_number())

TCGA.mutation_integrated %>% filter(Gene == "CTNNB1")

TCGA.mutation_integrated %>% head(20)

TCGA.mutation_plot_data <- TCGA.mutation_score_assoc %>% mutate(minus_log10_FDR = -log10(pmax(FDR, 1e-300)), 
    label = ifelse(Gene %in% c("CTNNB1", "TP53", "AXIN1", "APC"), Gene, NA))

library(ggplot2)

library(ggrepel)

p_mut_score <- ggplot(TCGA.mutation_plot_data, aes(x = delta_median, y = minus_log10_FDR)) + geom_point(size = 2.5) + 
    geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) + 
    geom_text_repel(aes(label = label), na.rm = TRUE, size = 4, max.overlaps = Inf) + theme_classic(base_size = 13) + 
    labs(x = expression(Delta * " median BA4 ICBScore (MUT - WT)"), y = expression(-log[10] * "(FDR)"), 
        title = "Mutation-associated BA4 ICBScore changes in TCGA-LIHC")

p_mut_score

TCGA.CTNNB1_plot <- data.frame(Patient_ID = TCGA.BA4_features$Patient_ID, BA4_ICBScore = TCGA.BA4_features$BA4_ICBScore, 
    CTNNB1_status = ifelse(TCGA.mut["CTNNB1", TCGA.BA4_features$Patient_ID] > 0, "MUT", "WT"))

TCGA.CTNNB1_plot$CTNNB1_status <- factor(TCGA.CTNNB1_plot$CTNNB1_status, levels = c("WT", "MUT"))

p_CTNNB1_score <- ggplot(TCGA.CTNNB1_plot, aes(x = CTNNB1_status, y = BA4_ICBScore)) + geom_boxplot(width = 0.55, 
    outlier.shape = NA) + geom_jitter(width = 0.12, size = 1.5, alpha = 0.5) + theme_classic(base_size = 13) + 
    labs(x = "CTNNB1 status", y = "BA4 ICBScore", title = "BA4 ICBScore according to CTNNB1 mutation status")

p_CTNNB1_score

TCGA.CTNNB1_group_table <- TCGA.ICB_group %>% mutate(CTNNB1_status = ifelse(TCGA.mut["CTNNB1", Patient_ID] > 
    0, "MUT", "WT")) %>% count(ICB_group, CTNNB1_status) %>% group_by(ICB_group) %>% mutate(proportion = n/sum(n)) %>% 
    ungroup()

TCGA.CTNNB1_group_table

ggplot(TCGA.CTNNB1_group_table %>% filter(CTNNB1_status == "MUT"), aes(x = ICB_group, y = proportion)) + 
    geom_col(width = 0.6) + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + theme_classic(base_size = 13) + 
    labs(x = "BA4 ICBScore group", y = "CTNNB1 mutation frequency", title = "CTNNB1 mutation enrichment in the low-BA4 group")

library(dplyr)

library(tidyr)

library(ggplot2)

TCGA.CTNNB1_feature_data <- TCGA.BA4_features %>% mutate(CTNNB1_mut = as.integer(TCGA.mut["CTNNB1", Patient_ID] > 
    0), CTNNB1_status = ifelse(CTNNB1_mut == 1, "MUT", "WT"), CTNNB1_status = factor(CTNNB1_status, levels = c("WT", 
    "MUT")))

TCGA.CTNNB1_feature_summary <- TCGA.CTNNB1_feature_data %>% group_by(CTNNB1_status) %>% summarise(across(all_of(LOCKED_FEATURES), 
    list(mean = ~mean(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE), SD = ~sd(.x, na.rm = TRUE))), 
    .groups = "drop")

TCGA.CTNNB1_feature_summary

test_CTNNB1_feature <- function(feature) {
    x_mut <- TCGA.CTNNB1_feature_data[TCGA.CTNNB1_feature_data$CTNNB1_status == "MUT", feature, drop = TRUE]
    x_wt <- TCGA.CTNNB1_feature_data[TCGA.CTNNB1_feature_data$CTNNB1_status == "WT", feature, drop = TRUE]
    wt <- wilcox.test(x_mut, x_wt, exact = FALSE)
    data.frame(feature = feature, n_MUT = length(x_mut), n_WT = length(x_wt), mean_MUT = mean(x_mut), 
        mean_WT = mean(x_wt), median_MUT = median(x_mut), median_WT = median(x_wt), delta_mean = mean(x_mut) - 
            mean(x_wt), delta_median = median(x_mut) - median(x_wt), P_value = wt$p.value)
}

TCGA.CTNNB1_feature_test <- bind_rows(lapply(LOCKED_FEATURES, test_CTNNB1_feature))

TCGA.CTNNB1_feature_test <- TCGA.CTNNB1_feature_test %>% mutate(FDR = p.adjust(P_value, method = "BH"), 
    direction = case_when(delta_median > 0 ~ "Higher_in_CTNNB1_MUT", delta_median < 0 ~ "Lower_in_CTNNB1_MUT", 
        TRUE ~ "No_change")) %>% arrange(P_value)

TCGA.CTNNB1_feature_effect <- bind_rows(lapply(LOCKED_FEATURES, function(f) {
    x <- TCGA.CTNNB1_feature_data[[f]]
    mut <- TCGA.CTNNB1_feature_data$CTNNB1_status == "MUT"
    pooled_sd <- sd(x, na.rm = TRUE)
    data.frame(feature = f, standardized_mean_difference = (mean(x[mut], na.rm = TRUE) - mean(x[!mut], 
        na.rm = TRUE))/pooled_sd, standardized_median_difference = (median(x[mut], na.rm = TRUE) - median(x[!mut], 
        na.rm = TRUE))/pooled_sd)
}))

TCGA.CTNNB1_feature_test <- TCGA.CTNNB1_feature_test %>% left_join(TCGA.CTNNB1_feature_effect, by = "feature")

TCGA.CTNNB1_feature_test

TCGA.CTNNB1_feature_long <- TCGA.CTNNB1_feature_data %>% select(Patient_ID, CTNNB1_status, all_of(LOCKED_FEATURES)) %>% 
    pivot_longer(cols = all_of(LOCKED_FEATURES), names_to = "feature", values_to = "score")

ggplot(TCGA.CTNNB1_feature_long, aes(x = CTNNB1_status, y = score)) + geom_boxplot(width = 0.55, outlier.shape = NA) + 
    geom_jitter(width = 0.12, alpha = 0.25, size = 1) + facet_wrap(~feature, scales = "free_y", ncol = 2) + 
    theme_classic(base_size = 12) + labs(x = "CTNNB1 status", y = "ssGSEA score", title = "BA4 component features according to CTNNB1 mutation status")

ggplot(TCGA.CTNNB1_feature_data, aes(x = CTNNB1_status, y = HALLMARK_BILE_ACID_METABOLISM)) + geom_boxplot(width = 0.55, 
    outlier.shape = NA) + geom_jitter(width = 0.12, alpha = 0.35, size = 1.2) + theme_classic(base_size = 13) + 
    labs(x = "CTNNB1 status", y = "Bile acid metabolism ssGSEA score", title = "Bile acid metabolism according to CTNNB1 mutation status")

library(glmnet)

BA4.coef_matrix <- as.matrix(coef(FINAL_BA4_MODEL, s = lambda_final))

BA4.coef_matrix

BA4.component_coefficients <- data.frame(feature = LOCKED_FEATURES, coefficient = as.numeric(BA4.coef_matrix[LOCKED_FEATURES, 
    1]))

BA4.component_coefficients

TCGA.CTNNB1_component_decomposition <- TCGA.CTNNB1_feature_test %>% left_join(BA4.component_coefficients, 
    by = "feature") %>% mutate(contribution_eta_mean = coefficient * delta_mean, contribution_eta_median = coefficient * 
    delta_median)

TCGA.CTNNB1_component_decomposition %>% select(feature, mean_MUT, mean_WT, delta_mean, coefficient, contribution_eta_mean, 
    P_value, FDR)

TCGA.CTNNB1_component_decomposition <- TCGA.CTNNB1_component_decomposition %>% arrange(contribution_eta_mean)

TCGA.CTNNB1_component_decomposition %>% select(feature, delta_mean, coefficient, contribution_eta_mean, 
    standardized_mean_difference, FDR)

TCGA.CTNNB1_total_eta_difference <- sum(TCGA.CTNNB1_component_decomposition$contribution_eta_mean)

TCGA.CTNNB1_total_eta_difference

TCGA.CTNNB1_negative_components <- TCGA.CTNNB1_component_decomposition %>% filter(contribution_eta_mean < 
    0) %>% mutate(absolute_contribution = abs(contribution_eta_mean), relative_negative_contribution = absolute_contribution/sum(absolute_contribution))

TCGA.CTNNB1_negative_components %>% select(feature, contribution_eta_mean, relative_negative_contribution)

TCGA.CTNNB1_negative_components %>% mutate(contribution_percent = relative_negative_contribution * 100) %>% 
    select(feature, contribution_eta_mean, contribution_percent)

library(ggplot2)

TCGA.component_plot_data <- TCGA.CTNNB1_component_decomposition %>% mutate(feature_short = case_when(feature == 
    "ICI__HCC_INFLAMMATORY4" ~ "Inflammatory4", feature == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" ~ 
    "ROS", feature == "HALLMARK_BILE_ACID_METABOLISM" ~ "Bile acid", feature == "HALLMARK_FATTY_ACID_METABOLISM" ~ 
    "Fatty acid", TRUE ~ feature))

ggplot(TCGA.component_plot_data, aes(x = reorder(feature_short, contribution_eta_mean), y = contribution_eta_mean)) + 
    geom_col(width = 0.65) + geom_hline(yintercept = 0, linetype = 2) + coord_flip() + theme_classic(base_size = 13) + 
    labs(x = NULL, y = expression(Delta * " linear predictor contribution"), title = "Contribution of BA4 components to the CTNNB1-associated score shift")

library(dplyr)

library(tibble)

library(ggplot2)

library(ggrepel)

BA_genes <- TCGA.BA4_gene_sets_available[["HALLMARK_BILE_ACID_METABOLISM"]]

length(BA_genes)

TCGA.BA_expr <- TCGA.log2_tpm[BA_genes, TCGA.BA4_features$Patient_ID, drop = FALSE]

TCGA.CTNNB1_status_vector <- ifelse(TCGA.mut["CTNNB1", colnames(TCGA.BA_expr)] > 0, "MUT", "WT")

test_BA_gene <- function(gene) {
    x <- TCGA.BA_expr[gene, , drop = TRUE]
    x_MUT <- x[TCGA.CTNNB1_status_vector == "MUT"]
    x_WT <- x[TCGA.CTNNB1_status_vector == "WT"]
    wt <- wilcox.test(x_MUT, x_WT, exact = FALSE)
    data.frame(Gene = gene, n_MUT = length(x_MUT), n_WT = length(x_WT), mean_MUT = mean(x_MUT, na.rm = TRUE), 
        mean_WT = mean(x_WT, na.rm = TRUE), median_MUT = median(x_MUT, na.rm = TRUE), median_WT = median(x_WT, 
            na.rm = TRUE), delta_mean = mean(x_MUT, na.rm = TRUE) - mean(x_WT, na.rm = TRUE), delta_median = median(x_MUT, 
            na.rm = TRUE) - median(x_WT, na.rm = TRUE), P_value = wt$p.value)
}


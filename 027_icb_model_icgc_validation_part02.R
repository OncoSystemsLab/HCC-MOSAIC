if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

ICGC.extract_cox_result <- function(model, variable, description) {
    s <- summary(model)
    data.frame(variable = variable, description = description, n = model$n, events = model$nevent, HR = unname(s$coefficients[variable, 
        "exp(coef)"]), CI_lower = unname(s$conf.int[variable, "lower .95"]), CI_upper = unname(s$conf.int[variable, 
        "upper .95"]), P_value = unname(s$coefficients[variable, "Pr(>|z|)"]), row.names = NULL)
}

ICGC.ICBScore_cox_results <- bind_rows(ICGC.extract_cox_result(ICGC.cox_ICB_01, "ICBScore_per01", "Per 0.1 increase in BA4 ICBScore"), 
    ICGC.extract_cox_result(ICGC.cox_ICB_z, "ICBScore_z", "Per 1-SD increase in BA4 ICBScore"), ICGC.extract_cox_result(ICGC.cox_ICB_group, 
        "ICBScore_groupHigh", "High vs Low; ICGC median cutoff"))

ICGC.PH_continuous <- survival::cox.zph(ICGC.cox_ICB_01)

ICGC.PH_group <- survival::cox.zph(ICGC.cox_ICB_group)

ICGC.fit_ICB_KM <- survival::survfit(Surv(OS_year, OS_status) ~ ICBScore_group, data = ICGC.surv_clean)

ICGC.logrank <- survival::survdiff(Surv(OS_year, OS_status) ~ ICBScore_group, data = ICGC.surv_clean, 
    rho = 0)

ICGC.logrank_result <- data.frame(chisq = unname(ICGC.logrank$chisq), df = length(ICGC.logrank$n) - 1, 
    P_value = pchisq(ICGC.logrank$chisq, df = length(ICGC.logrank$n) - 1, lower.tail = FALSE))

ICGC.p_ICB_KM <- survminer::ggsurvplot(ICGC.fit_ICB_KM, data = ICGC.surv_clean, risk.table = TRUE, pval = TRUE, 
    conf.int = FALSE, censor = TRUE, xlab = "Overall survival (years)", ylab = "Overall survival probability", 
    legend.title = "BA4 ICBScore", legend.labs = c("Low", "High"), palette = c("#377EB8", "#E41A1C"), 
    risk.table.height = 0.25, break.time.by = 2, ggtheme = theme_classic(base_size = 13))

library(dplyr)

ICGC.ICB_median <- median(ICGC.BA4_features$BA4_ICBScore)

ICGC.ICB_group <- ICGC.BA4_features %>% dplyr::select(Patient_ID, BA4_ICBScore) %>% mutate(ICB_group = factor(ifelse(BA4_ICBScore >= 
    ICGC.ICB_median, "High", "Low"), levels = c("Low", "High")))

stopifnot(!anyDuplicated(ICGC.ICB_group$Patient_ID), !anyNA(ICGC.ICB_group$ICB_group), all(ICGC.ICB_group$Patient_ID %in% 
    colnames(ICGC.mut)))

ICGC.mut_BA4 <- ICGC.mut[, ICGC.ICB_group$Patient_ID, drop = FALSE]

stopifnot(identical(colnames(ICGC.mut_BA4), ICGC.ICB_group$Patient_ID), all(ICGC.mut_BA4 %in% c(0, 1)))

ICGC.low_ids <- ICGC.ICB_group$Patient_ID[ICGC.ICB_group$ICB_group == "Low"]

ICGC.high_ids <- ICGC.ICB_group$Patient_ID[ICGC.ICB_group$ICB_group == "High"]

stopifnot(length(ICGC.low_ids) > 0, length(ICGC.high_ids) > 0)

ICGC.mutation_frequency <- data.frame(Gene = rownames(ICGC.mut_BA4), n_mut = rowSums(ICGC.mut_BA4 > 0)) %>% 
    mutate(mutation_frequency = n_mut/ncol(ICGC.mut_BA4)) %>% arrange(desc(mutation_frequency), Gene)

ICGC.common_mutation_genes <- ICGC.mutation_frequency %>% filter(mutation_frequency >= 0.05) %>% pull(Gene)

stopifnot(length(ICGC.common_mutation_genes) > 0)

ICGC.test_mutation_enrichment <- function(gene) {
    stopifnot(gene %in% rownames(ICGC.mut_BA4))
    mut_low <- sum(ICGC.mut_BA4[gene, ICGC.low_ids] > 0)
    mut_high <- sum(ICGC.mut_BA4[gene, ICGC.high_ids] > 0)
    tab <- matrix(c(mut_low, length(ICGC.low_ids) - mut_low, mut_high, length(ICGC.high_ids) - mut_high), 
        nrow = 2, byrow = TRUE, dimnames = list(Group = c("Low", "High"), Status = c("MUT", "WT")))
    res <- fisher.test(tab, alternative = "two.sided", conf.level = 0.95)
    data.frame(Gene = gene, mut_low = mut_low, n_low = length(ICGC.low_ids), freq_low = mut_low/length(ICGC.low_ids), 
        mut_high = mut_high, n_high = length(ICGC.high_ids), freq_high = mut_high/length(ICGC.high_ids), 
        odds_ratio = unname(res$estimate), CI_lower = unname(res$conf.int[1]), CI_upper = unname(res$conf.int[2]), 
        P_value = res$p.value)
}

ICGC.mutation_enrichment <- bind_rows(lapply(ICGC.common_mutation_genes, ICGC.test_mutation_enrichment)) %>% 
    mutate(FDR = p.adjust(P_value, method = "BH"), delta_frequency = freq_low - freq_high, enrichment_direction = case_when(delta_frequency > 
        0 ~ "Low_enriched", delta_frequency < 0 ~ "High_enriched", TRUE ~ "Equal"), rank_P = min_rank(P_value), 
        rank_FDR = min_rank(FDR)) %>% arrange(FDR, P_value, desc(delta_frequency))

ICGC.low_enriched_mutations <- ICGC.mutation_enrichment %>% filter(enrichment_direction == "Low_enriched") %>% 
    arrange(FDR, P_value, desc(delta_frequency))

ICGC.significant_mutations <- ICGC.mutation_enrichment %>% filter(FDR < 0.05)

ICGC.CTNNB1_enrichment <- ICGC.mutation_enrichment %>% filter(Gene == "CTNNB1")

ICGC.driver_enrichment <- ICGC.mutation_enrichment %>% filter(Gene %in% ICGC.driver_genes) %>% arrange(P_value)

ICGC.driver_screening_QC <- data.frame(Gene = ICGC.driver_genes, present_in_matrix = ICGC.driver_genes %in% 
    rownames(ICGC.mut_BA4), included_in_5pct_screen = ICGC.driver_genes %in% ICGC.common_mutation_genes)

library(dplyr)

library(ggplot2)

ICGC.score_vector <- setNames(ICGC.ICB_group$BA4_ICBScore, ICGC.ICB_group$Patient_ID)

stopifnot(!anyDuplicated(names(ICGC.score_vector)), all(is.finite(ICGC.score_vector)), identical(names(ICGC.score_vector), 
    colnames(ICGC.mut_BA4)), all(ICGC.common_mutation_genes %in% rownames(ICGC.mut_BA4)))

ICGC.test_mutation_score <- function(gene) {
    mut_status <- ICGC.mut_BA4[gene, names(ICGC.score_vector)] > 0
    score_mut <- unname(ICGC.score_vector[mut_status])
    score_wt <- unname(ICGC.score_vector[!mut_status])
    n_mut <- length(score_mut)
    n_wt <- length(score_wt)
    testable <- n_mut > 0L && n_wt > 0L
    p_value <- NA_real_
    rank_biserial <- NA_real_
    if (testable) {
        pooled <- c(score_mut, score_wt)
        if (length(unique(pooled)) == 1L) {
            p_value <- 1
            rank_biserial <- 0
        }
        else {
            wt <- wilcox.test(score_mut, score_wt, alternative = "two.sided", exact = FALSE, correct = TRUE)
            p_value <- wt$p.value
            pooled_ranks <- rank(pooled, ties.method = "average")
            U_mut <- sum(pooled_ranks[seq_len(n_mut)]) - n_mut * (n_mut + 1)/2
            rank_biserial <- 2 * U_mut/(n_mut * n_wt) - 1
        }
    }
    data.frame(Gene = gene, n_MUT = n_mut, n_WT = n_wt, mean_MUT = if (n_mut > 0) 
        mean(score_mut)
    else NA_real_, mean_WT = if (n_wt > 0) 
        mean(score_wt)
    else NA_real_, median_MUT = if (n_mut > 0) 
        median(score_mut)
    else NA_real_, median_WT = if (n_wt > 0) 
        median(score_wt)
    else NA_real_, delta_median = if (testable) {
        median(score_mut) - median(score_wt)
    }
    else {
        NA_real_
    }, rank_biserial = rank_biserial, P_value = p_value, testable = testable)
}

ICGC.mutation_score_assoc <- bind_rows(lapply(ICGC.common_mutation_genes, ICGC.test_mutation_score)) %>% 
    mutate(FDR = p.adjust(P_value, method = "BH"), direction = case_when(is.na(delta_median) ~ "Not_tested", 
        delta_median < 0 ~ "Lower_ICBScore_in_MUT", delta_median > 0 ~ "Higher_ICBScore_in_MUT", TRUE ~ 
            "Equal_medians"), rank_P = min_rank(P_value), rank_FDR = min_rank(FDR))

ICGC.mutation_score_assoc_ranked <- ICGC.mutation_score_assoc %>% arrange(P_value, Gene)

ICGC.continuous_scan_QC <- data.frame(n_patients = length(ICGC.score_vector), n_screened_genes = nrow(ICGC.mutation_score_assoc), 
    n_tested_genes = sum(is.finite(ICGC.mutation_score_assoc$P_value)), n_nominal_P_lt_005 = sum(ICGC.mutation_score_assoc$P_value < 
        0.05, na.rm = TRUE), n_FDR_lt_005 = sum(ICGC.mutation_score_assoc$FDR < 0.05, na.rm = TRUE))

ICGC.low_score_mutations <- ICGC.mutation_score_assoc %>% filter(direction == "Lower_ICBScore_in_MUT") %>% 
    arrange(P_value, Gene)

ICGC.CTNNB1_score_assoc <- ICGC.mutation_score_assoc %>% filter(Gene == "CTNNB1")

ICGC.driver_score_assoc <- ICGC.mutation_score_assoc %>% filter(Gene %in% ICGC.driver_genes) %>% arrange(P_value, 
    Gene)

ICGC.CTNNB1_score_df <- data.frame(Patient_ID = names(ICGC.score_vector), score = unname(ICGC.score_vector), 
    mutation = factor(ifelse(ICGC.mut_BA4["CTNNB1", names(ICGC.score_vector)] > 0, "MUT", "WT"), levels = c("WT", 
        "MUT")))

ICGC.p_CTNNB1_score <- ggplot(ICGC.CTNNB1_score_df, aes(x = mutation, y = score, fill = mutation)) + 
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, 
    height = 0, seed = 2026), size = 1.8, alpha = 0.75) + scale_fill_manual(values = c(WT = "#377EB8", 
    MUT = "#E41A1C")) + labs(x = "CTNNB1 status", y = "BA4 ICBScore", title = "ICGC-LICA-FR", subtitle = sprintf("Wilcoxon P = %.3g; BH-FDR = %.3g; rank-biserial = %.3f", 
    ICGC.CTNNB1_score_assoc$P_value, ICGC.CTNNB1_score_assoc$FDR, ICGC.CTNNB1_score_assoc$rank_biserial)) + 
    theme_classic(base_size = 13) + theme(legend.position = "none")

library(dplyr)

library(ggplot2)

stopifnot(!anyDuplicated(ICGC.mutation_score_assoc$Gene), !anyDuplicated(ICGC.mutation_enrichment$Gene), 
    setequal(ICGC.mutation_score_assoc$Gene, ICGC.mutation_enrichment$Gene))

ICGC.mutation_integrated <- ICGC.mutation_score_assoc %>% dplyr::select(Gene, n_MUT, n_WT, mean_MUT, 
    mean_WT, median_MUT, median_WT, delta_median, rank_biserial, P_score = P_value, FDR_score = FDR, 
    score_direction = direction) %>% left_join(ICGC.mutation_enrichment %>% dplyr::select(Gene, mut_low, 
    n_low, freq_low, mut_high, n_high, freq_high, odds_ratio, CI_lower, CI_upper, delta_frequency, P_group = P_value, 
    FDR_group = FDR, group_direction = enrichment_direction), by = "Gene") %>% mutate(direction_consistent = delta_median < 
    0 & delta_frequency > 0, significant_continuous = FDR_score < 0.05, nominal_group = P_group < 0.05, 
    meets_original_criteria = direction_consistent & significant_continuous & nominal_group, both_FDR_significant = FDR_score < 
        0.05 & FDR_group < 0.05, rank_continuous = min_rank(P_score)) %>% arrange(FDR_score, P_score, 
    Gene)

ICGC.low_ICB_mutation_candidates <- ICGC.mutation_integrated %>% filter(direction_consistent) %>% arrange(FDR_score, 
    P_score, P_group)

ICGC.original_criteria_candidates <- ICGC.mutation_integrated %>% filter(meets_original_criteria)

ICGC.both_FDR_candidates <- ICGC.mutation_integrated %>% filter(direction_consistent, both_FDR_significant)

ICGC.integrated_QC <- ICGC.mutation_integrated %>% summarise(n_genes = n(), n_direction_consistent = sum(direction_consistent, 
    na.rm = TRUE), n_original_criteria = sum(meets_original_criteria, na.rm = TRUE), n_low_direction_both_FDR = sum(direction_consistent & 
    both_FDR_significant, na.rm = TRUE))

ICGC.CTNNB1_integrated <- ICGC.mutation_integrated %>% filter(Gene == "CTNNB1")

ICGC.mutation_plot_data <- ICGC.mutation_integrated %>% filter(is.finite(delta_median), is.finite(FDR_score)) %>% 
    mutate(minus_log10_FDR = -log10(pmax(FDR_score, 1e-300)), highlighted = Gene %in% ICGC.driver_genes)

ICGC.mutation_label_data <- ICGC.mutation_plot_data %>% filter(highlighted)

ICGC.p_mut_score <- ggplot(ICGC.mutation_plot_data, aes(x = delta_median, y = minus_log10_FDR)) + geom_point(aes(color = highlighted), 
    size = 2.2, alpha = 0.75) + scale_color_manual(values = c(`FALSE` = "grey65", `TRUE` = "#B2182B"), 
    guide = "none") + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = -log10(0.05), 
    linetype = 2, color = "grey40") + expand_limits(y = -log10(0.05) + 0.15) + labs(x = "Median BA4 ICBScore difference (MUT - WT)", 
    y = expression(-log[10](FDR)), title = "Mutation-associated BA4 ICBScore in ICGC-LICA-FR", subtitle = "Dashed horizontal line: BH-FDR = 0.05") + 
    theme_classic(base_size = 13)

if (requireNamespace("ggrepel", quietly = TRUE)) {
    ICGC.p_mut_score <- ICGC.p_mut_score + ggrepel::geom_text_repel(data = ICGC.mutation_label_data, 
        aes(label = Gene), seed = 2026, max.overlaps = Inf, size = 3.5)
} else {
    ICGC.p_mut_score <- ICGC.p_mut_score + geom_text(data = ICGC.mutation_label_data, aes(label = Gene), 
        vjust = -0.7, size = 3.5, check_overlap = TRUE)
}

ICGC.CTNNB1_group_df <- ICGC.ICB_group %>% mutate(CTNNB1_status = factor(ifelse(ICGC.mut_BA4["CTNNB1", 
    Patient_ID] > 0, "MUT", "WT"), levels = c("WT", "MUT")))

ICGC.CTNNB1_group_table <- as.data.frame(with(ICGC.CTNNB1_group_df, table(ICB_group, CTNNB1_status)), 
    stringsAsFactors = FALSE) %>% rename(n = Freq) %>% group_by(ICB_group) %>% mutate(group_n = sum(n), 
    proportion = n/group_n) %>% ungroup()

ICGC.CTNNB1_frequency_plot_data <- ICGC.CTNNB1_group_table %>% filter(CTNNB1_status == "MUT") %>% mutate(count_label = sprintf("%d/%d (%.1f%%)", 
    n, group_n, proportion * 100))

ICGC.p_CTNNB1_frequency <- ggplot(ICGC.CTNNB1_frequency_plot_data, aes(x = ICB_group, y = proportion, 
    fill = ICB_group)) + geom_col(width = 0.6) + geom_text(aes(label = count_label), vjust = -0.5, size = 4) + 
    scale_fill_manual(values = c(Low = "#377EB8", High = "#E41A1C"), guide = "none") + scale_y_continuous(limits = c(0, 
    1), labels = function(x) paste0(round(x * 100), "%"), expand = expansion(mult = c(0, 0.03))) + labs(x = "BA4 ICBScore group", 
    y = "CTNNB1 mutation frequency", title = "CTNNB1 mutation frequency in ICGC-LICA-FR", subtitle = sprintf("Fisher P = %.3g; BH-FDR = %.3g", 
        ICGC.CTNNB1_integrated$P_group, ICGC.CTNNB1_integrated$FDR_group)) + theme_classic(base_size = 13)

library(dplyr)

library(tidyr)

library(ggplot2)

ICGC.CTNNB1_feature_data <- ICGC.BA4_features %>% mutate(CTNNB1_mut = as.integer(ICGC.mut["CTNNB1", Patient_ID] > 
    0), CTNNB1_status = factor(ifelse(CTNNB1_mut == 1, "MUT", "WT"), levels = c("WT", "MUT")))

stopifnot(!anyDuplicated(ICGC.CTNNB1_feature_data$Patient_ID), all(ICGC.BA4_feature_names %in% names(ICGC.CTNNB1_feature_data)), 
    all(is.finite(as.matrix(ICGC.CTNNB1_feature_data[, ICGC.BA4_feature_names, drop = FALSE]))))

ICGC.CTNNB1_feature_summary <- ICGC.CTNNB1_feature_data %>% group_by(CTNNB1_status) %>% summarise(n = n(), 
    across(all_of(ICGC.BA4_feature_names), list(mean = ~mean(.x), median = ~median(.x), SD = ~sd(.x))), 
    .groups = "drop")

ICGC.test_CTNNB1_feature <- function(feature) {
    x <- ICGC.CTNNB1_feature_data[[feature]]
    mut <- ICGC.CTNNB1_feature_data$CTNNB1_status == "MUT"
    x_mut <- x[mut]
    x_wt <- x[!mut]
    n_mut <- length(x_mut)
    n_wt <- length(x_wt)
    stopifnot(n_mut > 0, n_wt > 0)
    cohort_sd <- sd(x)
    if (length(unique(x)) == 1L) {
        p_value <- 1
        rank_biserial <- 0
    }
    else {
        wt <- wilcox.test(x_mut, x_wt, alternative = "two.sided", exact = FALSE, correct = TRUE)
        p_value <- wt$p.value
        ranks <- rank(c(x_mut, x_wt), ties.method = "average")
        U_mut <- sum(ranks[seq_len(n_mut)]) - n_mut * (n_mut + 1)/2
        rank_biserial <- 2 * U_mut/(n_mut * n_wt) - 1
    }
    delta_mean <- mean(x_mut) - mean(x_wt)
    delta_median <- median(x_mut) - median(x_wt)
    data.frame(feature = feature, n_MUT = n_mut, n_WT = n_wt, mean_MUT = mean(x_mut), mean_WT = mean(x_wt), 
        median_MUT = median(x_mut), median_WT = median(x_wt), delta_mean = delta_mean, delta_median = delta_median, 
        cohort_SD = cohort_sd, mean_difference_per_cohort_SD = if (cohort_sd > 0) {
            delta_mean/cohort_sd
        }
        else {
            NA_real_
        }, median_difference_per_cohort_SD = if (cohort_sd > 0) {
            delta_median/cohort_sd
        }
        else {
            NA_real_
        }, rank_biserial = rank_biserial, P_value = p_value)
}

ICGC.CTNNB1_feature_test <- bind_rows(lapply(ICGC.BA4_feature_names, ICGC.test_CTNNB1_feature)) %>% mutate(FDR = p.adjust(P_value, 
    method = "BH"), direction = case_when(delta_median > 0 ~ "Higher_in_CTNNB1_MUT", delta_median < 0 ~ 
    "Lower_in_CTNNB1_MUT", TRUE ~ "Equal_medians")) %>% arrange(P_value)

ICGC.CTNNB1_feature_long <- ICGC.CTNNB1_feature_data %>% dplyr::select(Patient_ID, CTNNB1_status, all_of(ICGC.BA4_feature_names)) %>% 
    pivot_longer(cols = all_of(ICGC.BA4_feature_names), names_to = "feature", values_to = "score")

ICGC.feature_plot_labels <- setNames(paste0(gsub("_", " ", ICGC.CTNNB1_feature_test$feature, fixed = TRUE), 
    "", sprintf("P = %.3g; BH-FDR = %.3g", ICGC.CTNNB1_feature_test$P_value, ICGC.CTNNB1_feature_test$FDR)), 
    ICGC.CTNNB1_feature_test$feature)

ICGC.CTNNB1_feature_long$feature <- factor(ICGC.CTNNB1_feature_long$feature, levels = ICGC.BA4_feature_names)

ICGC.p_CTNNB1_features <- ggplot(ICGC.CTNNB1_feature_long, aes(x = CTNNB1_status, y = score, fill = CTNNB1_status)) + 
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, 
    height = 0, seed = 2026), alpha = 0.5, size = 1.2) + facet_wrap(~feature, scales = "free_y", ncol = 2, 
    labeller = labeller(feature = as_labeller(ICGC.feature_plot_labels))) + scale_fill_manual(values = c(WT = "#377EB8", 
    MUT = "#E41A1C")) + labs(x = "CTNNB1 status", y = "ssGSEA score", title = "BA4 component features in ICGC-LICA-FR") + 
    theme_classic(base_size = 12) + theme(legend.position = "none", strip.text = element_text(size = 9))

ICGC.BA_feature <- "HALLMARK_BILE_ACID_METABOLISM"

ICGC.BA_feature_result <- ICGC.CTNNB1_feature_test %>% filter(feature == ICGC.BA_feature)

ICGC.p_CTNNB1_BA <- ggplot(ICGC.CTNNB1_feature_data, aes(x = CTNNB1_status, y = .data[[ICGC.BA_feature]], 
    fill = CTNNB1_status)) + geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, 
    height = 0, seed = 2026), alpha = 0.5, size = 1.5) + scale_fill_manual(values = c(WT = "#377EB8", 
    MUT = "#E41A1C")) + labs(x = "CTNNB1 status", y = "Bile acid metabolism ssGSEA score", title = "ICGC-LICA-FR", 
    subtitle = sprintf("P = %.3g; BH-FDR = %.3g; rank-biserial = %.3f", ICGC.BA_feature_result$P_value, 
        ICGC.BA_feature_result$FDR, ICGC.BA_feature_result$rank_biserial)) + theme_classic(base_size = 13) + 
    theme(legend.position = "none")

library(glmnet)

library(dplyr)

library(ggplot2)

ICGC.BA4_coef_matrix <- as.matrix(coef(BA4_MODEL_BUNDLE$model, s = BA4_MODEL_BUNDLE$lambda))

stopifnot(ncol(ICGC.BA4_coef_matrix) == 1L, all(ICGC.BA4_feature_names %in% rownames(ICGC.BA4_coef_matrix)))

ICGC.BA4_component_coefficients <- data.frame(feature = ICGC.BA4_feature_names, coefficient = as.numeric(ICGC.BA4_coef_matrix[ICGC.BA4_feature_names, 
    1]))

ICGC.CTNNB1_component_decomposition <- ICGC.CTNNB1_feature_test %>% left_join(ICGC.BA4_component_coefficients, 
    by = "feature") %>% mutate(contribution_eta_mean = coefficient * delta_mean, contribution_eta_median = coefficient * 
    delta_median, contribution_direction = case_when(contribution_eta_mean < 0 ~ "Negative", contribution_eta_mean > 
    0 ~ "Positive", TRUE ~ "Zero")) %>% arrange(contribution_eta_mean)

stopifnot(nrow(ICGC.CTNNB1_component_decomposition) == 4L, all(is.finite(ICGC.CTNNB1_component_decomposition$contribution_eta_mean)))

ICGC.BA4_X <- as.matrix(ICGC.CTNNB1_feature_data[, ICGC.BA4_feature_names, drop = FALSE])

storage.mode(ICGC.BA4_X) <- "numeric"

ICGC.BA4_eta <- as.numeric(predict(BA4_MODEL_BUNDLE$model, newx = ICGC.BA4_X, s = BA4_MODEL_BUNDLE$lambda, 
    type = "link"))

ICGC.BA4_probability <- as.numeric(predict(BA4_MODEL_BUNDLE$model, newx = ICGC.BA4_X, s = BA4_MODEL_BUNDLE$lambda, 
    type = "response"))

stopifnot(all(is.finite(ICGC.BA4_eta)), isTRUE(all.equal(ICGC.BA4_probability, ICGC.CTNNB1_feature_data$BA4_ICBScore, 
    tolerance = 1e-08)))

ICGC.CTNNB1_eta_data <- ICGC.CTNNB1_feature_data %>% dplyr::select(Patient_ID, CTNNB1_status, BA4_ICBScore) %>% 
    mutate(eta = ICGC.BA4_eta)

ICGC.CTNNB1_eta_summary <- ICGC.CTNNB1_eta_data %>% group_by(CTNNB1_status) %>% summarise(n = n(), mean_eta = mean(eta), 
    median_eta = median(eta), mean_probability = mean(BA4_ICBScore), median_probability = median(BA4_ICBScore), 
    .groups = "drop")

ICGC.CTNNB1_total_eta_difference <- sum(ICGC.CTNNB1_component_decomposition$contribution_eta_mean)

ICGC.CTNNB1_direct_eta_difference <- with(ICGC.CTNNB1_eta_data, mean(eta[CTNNB1_status == "MUT"]) - mean(eta[CTNNB1_status == 
    "WT"]))

ICGC.decomposition_QC <- data.frame(sum_component_contributions = ICGC.CTNNB1_total_eta_difference, direct_mean_eta_difference = ICGC.CTNNB1_direct_eta_difference, 
    numerical_difference = ICGC.CTNNB1_total_eta_difference - ICGC.CTNNB1_direct_eta_difference, observed_mean_probability_difference = with(ICGC.CTNNB1_eta_data, 
        mean(BA4_ICBScore[CTNNB1_status == "MUT"]) - mean(BA4_ICBScore[CTNNB1_status == "WT"])))

stopifnot(isTRUE(all.equal(ICGC.CTNNB1_total_eta_difference, ICGC.CTNNB1_direct_eta_difference, tolerance = 1e-08)))

ICGC.CTNNB1_negative_components <- ICGC.CTNNB1_component_decomposition %>% filter(contribution_eta_mean < 
    0) %>% mutate(absolute_contribution = abs(contribution_eta_mean), relative_negative_contribution = absolute_contribution/sum(absolute_contribution), 
    contribution_percent = 100 * relative_negative_contribution)

ICGC.contribution_totals <- data.frame(negative_sum = sum(pmin(ICGC.CTNNB1_component_decomposition$contribution_eta_mean, 
    0)), positive_sum = sum(pmax(ICGC.CTNNB1_component_decomposition$contribution_eta_mean, 0)), net_eta_difference = ICGC.CTNNB1_total_eta_difference)

ICGC.component_plot_data <- ICGC.CTNNB1_component_decomposition %>% mutate(feature_short = case_when(feature == 
    "ICI__HCC_INFLAMMATORY4" ~ "Inflammatory4", feature == "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY" ~ 
    "ROS", feature == "HALLMARK_BILE_ACID_METABOLISM" ~ "Bile acid", feature == "HALLMARK_FATTY_ACID_METABOLISM" ~ 
    "Fatty acid", TRUE ~ feature))

ICGC.p_component_contribution <- ggplot(ICGC.component_plot_data, aes(x = reorder(feature_short, contribution_eta_mean), 
    y = contribution_eta_mean, fill = contribution_direction)) + geom_col(width = 0.65) + geom_hline(yintercept = 0, 
    linetype = 2) + coord_flip() + scale_fill_manual(values = c(Negative = "#377EB8", Positive = "#E41A1C", 
    Zero = "grey65"), guide = "none") + labs(x = NULL, y = "Contribution to mean linear predictor difference (MUT - WT)", 
    title = "BA4 component contributions in ICGC-LICA-FR", subtitle = sprintf("Net mean linear predictor difference = %.4f", 
        ICGC.CTNNB1_total_eta_difference)) + theme_classic(base_size = 13)

library(dplyr)

library(ggplot2)

ICGC.BA_genes <- unique(ICGC.BA4_gene_sets_available[["HALLMARK_BILE_ACID_METABOLISM"]])

ICGC.BA_expr <- ICGC.log2_tpm[ICGC.BA_genes, ICGC.BA4_features$Patient_ID, drop = FALSE]

ICGC.CTNNB1_status_vector <- setNames(ifelse(ICGC.mut["CTNNB1", colnames(ICGC.BA_expr)] > 0, "MUT", "WT"), 
    colnames(ICGC.BA_expr))

stopifnot(length(ICGC.BA_genes) > 0, all(is.finite(ICGC.BA_expr)), identical(colnames(ICGC.BA_expr), 
    names(ICGC.CTNNB1_status_vector)))


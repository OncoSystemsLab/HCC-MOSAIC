if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

CPTAC.compare <- function(values, is_mut, B = 0L) {
    stopifnot(length(values) == length(is_mut), all(is.finite(values)), !anyNA(is_mut))
    x <- values[is_mut]
    y <- values[!is_mut]
    nx <- length(x)
    ny <- length(y)
    testable <- nx > 0L && ny > 0L
    p <- effect <- NA_real_
    ci <- c(NA_real_, NA_real_)
    if (testable) {
        effect <- CPTAC.rank_biserial(x, y)
        p <- if (length(unique(c(x, y))) == 1L) {
            1
        }
        else {
            wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
        }
        if (B > 0L && nx >= 2L && ny >= 2L) {
            boot <- replicate(B, {
                CPTAC.rank_biserial(x[sample.int(nx, nx, replace = TRUE)], y[sample.int(ny, ny, replace = TRUE)])
            })
            ci <- unname(quantile(boot, c(0.025, 0.975)))
        }
    }
    cohort_sd <- sd(values)
    data.frame(n_MUT = nx, n_WT = ny, mean_MUT = if (nx) 
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
    else NA_real_, mean_difference_per_cohort_SD = if (testable && cohort_sd > 0) {
        (mean(x) - mean(y))/cohort_sd
    }
    else {
        NA_real_
    }, rank_biserial = effect, CI_lower = ci[1], CI_upper = ci[2], P_value = p, testable = testable)
}

CPTAC.adjust_BH <- function(dat) {
    dat$FDR <- NA_real_
    i <- which(is.finite(dat$P_value))
    dat$FDR[i] <- p.adjust(dat$P_value[i], method = "BH")
    dat
}

CPTAC.surv <- CPTAC.analysis %>% filter(is.finite(OS_time_native), OS_time_native > 0, OS_status %in% 
    c(0, 1), is.finite(BA4_ICBScore))

CPTAC.survival_QC <- data.frame(total = nrow(CPTAC.analysis), included = nrow(CPTAC.surv), excluded = nrow(CPTAC.analysis) - 
    nrow(CPTAC.surv), deaths = sum(CPTAC.surv$OS_status == 1), time_unit_assumed = CPTAC.OS_unit)

CPTAC.write(CPTAC.survival_QC, "08_survival_QC.tsv")

CPTAC.write(CPTAC.analysis %>% filter(!Patient_ID %in% CPTAC.surv$Patient_ID), "08_survival_excluded.tsv")

stopifnot(nrow(CPTAC.surv) > 0, sum(CPTAC.surv$OS_status == 1) > 0, sd(CPTAC.surv$BA4_ICBScore) > 0, 
    CPTAC.OS_unit %in% c("months", "days", "years"))

CPTAC.surv <- CPTAC.surv %>% mutate(OS_year = OS_time_native/switch(CPTAC.OS_unit, months = 12, days = 365.25, 
    years = 1), ICBScore_per01 = BA4_ICBScore/0.1, ICBScore_z = as.numeric(scale(BA4_ICBScore)), ICBScore_group = factor(ifelse(BA4_ICBScore >= 
    median(BA4_ICBScore), "High", "Low"), levels = c("Low", "High")))

stopifnot(nlevels(droplevels(CPTAC.surv$ICBScore_group)) == 2L)

CPTAC.cox_models <- list(per_0.1_score = coxph(Surv(OS_year, OS_status) ~ ICBScore_per01, data = CPTAC.surv, 
    ties = "efron", x = TRUE), per_1_SD = coxph(Surv(OS_year, OS_status) ~ ICBScore_z, data = CPTAC.surv, 
    ties = "efron", x = TRUE), High_vs_Low = coxph(Surv(OS_year, OS_status) ~ ICBScore_group, data = CPTAC.surv, 
    ties = "efron", x = TRUE))

CPTAC.cox_results <- bind_rows(lapply(names(CPTAC.cox_models), function(nm) {
    fit <- CPTAC.cox_models[[nm]]
    z <- summary(fit)
    data.frame(Model = nm, term = rownames(z$coefficients), HR = z$conf.int[, "exp(coef)"], CI_lower = z$conf.int[, 
        "lower .95"], CI_upper = z$conf.int[, "upper .95"], P_value = z$coefficients[, "Pr(>|z|)"], n = fit$n, 
        events = fit$nevent)
}))

CPTAC.write(CPTAC.cox_results, "08_Cox_results.tsv")

CPTAC.PH_tests <- lapply(CPTAC.cox_models, survival::cox.zph)

CPTAC.PH_results <- bind_rows(lapply(names(CPTAC.PH_tests), function(nm) {
    as.data.frame(CPTAC.PH_tests[[nm]]$table) %>% rownames_to_column("term") %>% mutate(Model = nm)
}))

CPTAC.write(CPTAC.PH_results, "08_Cox_PH_tests.tsv")

CPTAC.KM_fit <- survfit(Surv(OS_year, OS_status) ~ ICBScore_group, data = CPTAC.surv)

CPTAC.KM_plot <- survminer::ggsurvplot(CPTAC.KM_fit, data = CPTAC.surv, risk.table = TRUE, pval = TRUE, 
    conf.int = FALSE, palette = c("#2878B5", "#D73027"), legend.labs = c("Low", "High"), legend.title = "BA4 Score", 
    xlab = "Overall survival (years)", ylab = "Overall survival probability", title = "CPTAC / CHCC-HBV", 
    ggtheme = theme_classic(base_size = 13))

pdf(file.path(CPTAC.out_dir, "08_BA4_OS_KM.pdf"), width = 8, height = 7)

dev.off()

CPTAC.ICB_group <- CPTAC.BA4_features %>% mutate(ICB_group = factor(ifelse(BA4_ICBScore >= median(BA4_ICBScore), 
    "High", "Low"), levels = c("Low", "High")))

CPTAC.mutation_frequency <- data.frame(Gene = rownames(CPTAC.mut), n_MUT = as.integer(rowSums(CPTAC.mut))) %>% 
    mutate(n_WT = ncol(CPTAC.mut) - n_MUT, mutation_frequency = n_MUT/ncol(CPTAC.mut)) %>% arrange(desc(mutation_frequency), 
    Gene)

CPTAC.common_mutation_genes <- CPTAC.mutation_frequency %>% filter(mutation_frequency > CPTAC.mutation_threshold) %>% 
    pull(Gene)

if (!length(CPTAC.common_mutation_genes)) {
    stop("No genes meet the mutation-frequency threshold.")
}

CPTAC.write(CPTAC.mutation_frequency, "09_all_mutation_frequencies.tsv")

CPTAC.is_low <- CPTAC.ICB_group$ICB_group == "Low"

CPTAC.mutation_enrichment <- bind_rows(lapply(CPTAC.common_mutation_genes, function(gene) {
    m <- CPTAC.mut[gene, ] == 1
    a <- sum(m & CPTAC.is_low)
    b <- sum(!m & CPTAC.is_low)
    c <- sum(m & !CPTAC.is_low)
    d <- sum(!m & !CPTAC.is_low)
    tab <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, dimnames = list(c("Low", "High"), c("MUT", "WT")))
    if (sum(m) == length(m)) {
        return(data.frame(Gene = gene, mut_low = a, n_low = a + b, mut_high = c, n_high = c + d, freq_low = a/(a + 
            b), freq_high = c/(c + d), odds_ratio = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, 
            P_value = NA_real_))
    }
    ft <- fisher.test(tab, alternative = "two.sided")
    data.frame(Gene = gene, mut_low = a, n_low = a + b, mut_high = c, n_high = c + d, freq_low = a/(a + 
        b), freq_high = c/(c + d), odds_ratio = unname(ft$estimate), CI_lower = ft$conf.int[1], CI_upper = ft$conf.int[2], 
        P_value = ft$p.value)
})) %>% CPTAC.adjust_BH() %>% mutate(delta_frequency = freq_low - freq_high, enrichment_direction = case_when(delta_frequency > 
    0 ~ "Low_enriched", delta_frequency < 0 ~ "High_enriched", TRUE ~ "Equal")) %>% arrange(FDR, P_value)

CPTAC.write(CPTAC.mutation_enrichment, "09_mutation_HighLow_Fisher.tsv")

set.seed(2026)

CPTAC.gene_BA4_effects <- bind_rows(lapply(CPTAC.common_mutation_genes, function(gene) {
    cbind(data.frame(Gene = gene, mutation_frequency = mean(CPTAC.mut[gene, ] == 1)), CPTAC.compare(CPTAC.BA4_features$BA4_ICBScore, 
        CPTAC.mut[gene, ] == 1, B = CPTAC.bootstrap_B))
})) %>% CPTAC.adjust_BH() %>% mutate(Direction = case_when(!testable ~ "Not testable", rank_biserial > 
    0 ~ "Higher in MUT", rank_biserial < 0 ~ "Lower in MUT", TRUE ~ "No rank difference")) %>% arrange(FDR, 
    desc(abs(rank_biserial)))

CPTAC.mutation_score_assoc <- CPTAC.gene_BA4_effects

CPTAC.write(CPTAC.gene_BA4_effects, "10_gene_BA4_effects.tsv")

CPTAC.effect_plot_data <- CPTAC.gene_BA4_effects %>% filter(testable) %>% arrange(rank_biserial) %>% 
    mutate(Gene_label = factor(Gene, levels = Gene), Significance = if_else(FDR < 0.05, "FDR < 0.05", 
        "FDR >= 0.05"))

CPTAC.p_gene_BA4_effect <- ggplot(CPTAC.effect_plot_data, aes(x = rank_biserial, y = Gene_label)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey60") + geom_segment(aes(x = CI_lower, xend = CI_upper, yend = Gene_label), 
    color = "grey55", linewidth = 0.6, na.rm = TRUE) + geom_point(aes(fill = Significance, size = -log10(pmax(P_value, 
    .Machine$double.xmin))), shape = 21, color = "black", stroke = 0.6) + scale_fill_manual(name = NULL, 
    values = c(`FDR < 0.05` = "#D73027", `FDR >= 0.05` = "#2878B5")) + scale_size_continuous(name = expression(-log[10](P)), 
    range = c(2, 7)) + scale_x_continuous(limits = c(-1, 1)) + labs(x = "Rank-biserial correlation (MUT vs WT)", 
    y = NULL, title = "Mutation-associated differences in BA4 Score", subtitle = "CPTAC / CHCC-HBV | Mutation frequency > 5%", 
    caption = "Lines: bootstrap 95% CI; positive values: higher score in MUT.") + theme_classic(base_size = 13) + 
    theme(axis.text.x = element_text(size = 14, color = "black"), axis.text.y = element_text(size = 14, 
        color = "black", face = "italic"), axis.title.x = element_text(size = 16), legend.position = "top")

CPTAC.save_plot(CPTAC.p_gene_BA4_effect, "10_gene_BA4_effect_forest.pdf", width = 10, height = max(5, 
    0.28 * nrow(CPTAC.effect_plot_data) + 2))

CPTAC.write(CPTAC.effect_plot_data, "10_gene_BA4_effect_plot_data.tsv")

CPTAC.mutation_integrated <- CPTAC.gene_BA4_effects %>% select(Gene, mutation_frequency, n_MUT, n_WT, 
    delta_mean, delta_median, rank_biserial, P_score = P_value, FDR_score = FDR) %>% left_join(CPTAC.mutation_enrichment %>% 
    select(Gene, freq_low, freq_high, delta_frequency, odds_ratio, P_group = P_value, FDR_group = FDR), 
    by = "Gene") %>% mutate(direction_consistent = delta_median < 0 & delta_frequency > 0, meets_original_criteria = direction_consistent & 
    FDR_score < 0.05 & P_group < 0.05, both_FDR_significant = FDR_score < 0.05 & FDR_group < 0.05) %>% 
    arrange(FDR_score, P_score)

CPTAC.write(CPTAC.mutation_integrated, "11_mutation_integrated.tsv")

CPTAC.write(CPTAC.mutation_integrated %>% filter(meets_original_criteria), "11_original_criteria_candidates.tsv")

CPTAC.is_CTNNB1_mut <- CPTAC.analysis$CTNNB1_mut == 1L

stopifnot(any(CPTAC.is_CTNNB1_mut), any(!CPTAC.is_CTNNB1_mut))

CPTAC.CTNNB1_feature_test <- bind_rows(lapply(CPTAC.BA4_feature_names, function(f) {
    cbind(data.frame(feature = f), CPTAC.compare(CPTAC.analysis[[f]], CPTAC.is_CTNNB1_mut))
})) %>% CPTAC.adjust_BH()

CPTAC.CTNNB1_score_test <- CPTAC.compare(CPTAC.analysis$BA4_ICBScore, CPTAC.is_CTNNB1_mut, B = CPTAC.bootstrap_B)

CPTAC.write(CPTAC.CTNNB1_feature_test, "12_CTNNB1_four_features.tsv")

CPTAC.write(CPTAC.CTNNB1_score_test, "12_CTNNB1_BA4_score.tsv")

CPTAC.feature_long <- CPTAC.analysis %>% select(Patient_ID, CTNNB1_status, all_of(CPTAC.BA4_feature_names)) %>% 
    pivot_longer(all_of(CPTAC.BA4_feature_names), names_to = "feature", values_to = "ssGSEA_score")

CPTAC.p_features <- ggplot(CPTAC.feature_long, aes(CTNNB1_status, ssGSEA_score, fill = CTNNB1_status)) + 
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, 
    height = 0, seed = 2026), size = 1.2, alpha = 0.5) + facet_wrap(~feature, scales = "free_y", ncol = 2) + 
    scale_fill_manual(values = c(WT = "#2878B5", MUT = "#D73027")) + labs(x = "CTNNB1 status", y = "ssGSEA score", 
    title = "CTNNB1 mutation and BA4 component features") + theme_classic(base_size = 12) + theme(legend.position = "none")

CPTAC.save_plot(CPTAC.p_features, "12_CTNNB1_four_features.pdf", width = 11, height = 7)

CPTAC.p_CTNNB1_score <- ggplot(CPTAC.analysis, aes(CTNNB1_status, BA4_ICBScore, fill = CTNNB1_status)) + 
    geom_boxplot(outlier.shape = NA, width = 0.55, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, 
    height = 0, seed = 2026), alpha = 0.5, size = 1.5) + scale_fill_manual(values = c(WT = "#2878B5", 
    MUT = "#D73027")) + labs(x = "CTNNB1 status", y = "BA4 Score", subtitle = sprintf("Wilcoxon P = %.3g; rank-biserial = %.3f", 
    CPTAC.CTNNB1_score_test$P_value, CPTAC.CTNNB1_score_test$rank_biserial)) + theme_classic(base_size = 13) + 
    theme(legend.position = "none")

CPTAC.save_plot(CPTAC.p_CTNNB1_score, "12_CTNNB1_BA4_score.pdf")

CPTAC.coef_matrix <- as.matrix(coef(BA4_MODEL_BUNDLE$model, s = BA4_MODEL_BUNDLE$lambda))

CPTAC.component_coefficients <- data.frame(feature = CPTAC.BA4_feature_names, coefficient = as.numeric(CPTAC.coef_matrix[CPTAC.BA4_feature_names, 
    1]))

CPTAC.CTNNB1_component_decomposition <- CPTAC.CTNNB1_feature_test %>% left_join(CPTAC.component_coefficients, 
    by = "feature") %>% mutate(contribution_eta_mean = coefficient * delta_mean, contribution_direction = if_else(contribution_eta_mean >= 
    0, "Positive", "Negative")) %>% arrange(contribution_eta_mean)

CPTAC.direct_eta_difference <- mean(CPTAC.analysis$BA4_eta[CPTAC.is_CTNNB1_mut]) - mean(CPTAC.analysis$BA4_eta[!CPTAC.is_CTNNB1_mut])

CPTAC.decomposition_QC <- data.frame(direct_delta_mean_eta = CPTAC.direct_eta_difference, sum_component_contributions = sum(CPTAC.CTNNB1_component_decomposition$contribution_eta_mean))

stopifnot(isTRUE(all.equal(CPTAC.decomposition_QC$direct_delta_mean_eta, CPTAC.decomposition_QC$sum_component_contributions, 
    tolerance = 1e-08)))

CPTAC.write(CPTAC.CTNNB1_component_decomposition, "13_CTNNB1_component_decomposition.tsv")

CPTAC.write(CPTAC.decomposition_QC, "13_decomposition_QC.tsv")

CPTAC.p_decomposition <- ggplot(CPTAC.CTNNB1_component_decomposition, aes(x = contribution_eta_mean, 
    y = reorder(feature, contribution_eta_mean), fill = contribution_direction)) + geom_col(width = 0.65) + 
    geom_vline(xintercept = 0, linetype = "dashed") + scale_fill_manual(values = c(Positive = "#D73027", 
    Negative = "#2878B5")) + labs(x = "Contribution to mean eta difference (MUT - WT)", y = NULL, fill = NULL, 
    title = "CTNNB1-associated BA4 model decomposition") + theme_classic(base_size = 12)

CPTAC.save_plot(CPTAC.p_decomposition, "13_CTNNB1_component_decomposition.pdf", width = 10)

CPTAC.BA_feature <- "HALLMARK_BILE_ACID_METABOLISM"

stopifnot(CPTAC.BA_feature %in% names(CPTAC.BA4_gene_sets_available))

CPTAC.BA_genes <- CPTAC.BA4_gene_sets_available[[CPTAC.BA_feature]]

CPTAC.CTNNB1_BA_gene_test <- bind_rows(lapply(CPTAC.BA_genes, function(gene) {
    cbind(data.frame(Gene = gene), CPTAC.compare(as.numeric(CPTAC.expr[gene, CPTAC.patient_ids]), CPTAC.is_CTNNB1_mut))
})) %>% CPTAC.adjust_BH() %>% arrange(FDR, P_value)

CPTAC.write(CPTAC.CTNNB1_BA_gene_test, "14_CTNNB1_BA_pathway_genes.tsv")

CPTAC.write(CPTAC.CTNNB1_BA_gene_test %>% filter(FDR < 0.05), "14_CTNNB1_BA_significant_genes.tsv")

CPTAC.BA_volcano_data <- CPTAC.CTNNB1_BA_gene_test %>% mutate(significance = case_when(FDR < 0.05 & delta_mean > 
    0 ~ "Higher in MUT", FDR < 0.05 & delta_mean < 0 ~ "Lower in MUT", TRUE ~ "Not significant"), minus_log10_FDR = -log10(pmax(FDR, 
    .Machine$double.xmin)))

CPTAC.p_BA_volcano <- ggplot(CPTAC.BA_volcano_data, aes(delta_mean, minus_log10_FDR, color = significance)) + 
    geom_point(size = 2, alpha = 0.8) + geom_vline(xintercept = 0, linetype = "dashed") + geom_hline(yintercept = -log10(0.05), 
    linetype = "dashed") + scale_color_manual(values = c(`Higher in MUT` = "#D73027", `Lower in MUT` = "#2878B5", 
    `Not significant` = "grey70")) + labs(x = paste0("Mean ", CPTAC.expression_label, " difference (MUT - WT)"), 
    y = expression(-log[10](FDR)), color = NULL, title = "CTNNB1-associated bile acid gene expression") + 
    theme_classic(base_size = 12)

CPTAC.save_plot(CPTAC.p_BA_volcano, "14_CTNNB1_BA_gene_volcano.pdf", width = 10)

stopifnot("BAAT" %in% rownames(CPTAC.expr))

CPTAC.BAAT_data <- CPTAC.analysis %>% mutate(BAAT = as.numeric(CPTAC.expr["BAAT", Patient_ID]))

CPTAC.BAAT_test <- cbind(data.frame(Gene = "BAAT", included_in_frozen_BA_set = "BAAT" %in% CPTAC.BA_genes), 
    CPTAC.compare(CPTAC.BAAT_data$BAAT, CPTAC.is_CTNNB1_mut, B = CPTAC.bootstrap_B))

CPTAC.BA_plus_BAAT_test <- bind_rows(CPTAC.CTNNB1_BA_gene_test %>% filter(Gene != "BAAT") %>% select(-FDR), 
    CPTAC.BAAT_test %>% select(-included_in_frozen_BA_set)) %>% CPTAC.adjust_BH() %>% arrange(FDR, P_value)

CPTAC.write(CPTAC.BAAT_test, "15_BAAT_separate_test.tsv")

CPTAC.write(CPTAC.BA_plus_BAAT_test, "15_BA_pathway_plus_BAAT_tests.tsv")

CPTAC.BAAT_correlation_results <- bind_rows(lapply(c(CPTAC.BA_feature, "BA4_ICBScore"), function(target) {
    x <- CPTAC.BAAT_data$BAAT
    y <- CPTAC.BAAT_data[[target]]
    if (sd(x) == 0 || sd(y) == 0) {
        return(data.frame(target = target, n = length(x), rho = NA_real_, P_value = NA_real_))
    }
    z <- cor.test(x, y, method = "spearman", exact = FALSE)
    data.frame(target = target, n = length(x), rho = unname(z$estimate), P_value = z$p.value)
})) %>% CPTAC.adjust_BH()

CPTAC.write(CPTAC.BAAT_correlation_results, "15_BAAT_Spearman_correlations.tsv")

CPTAC.p_BAAT <- ggplot(CPTAC.BAAT_data, aes(CTNNB1_status, BAAT, fill = CTNNB1_status)) + geom_boxplot(width = 0.55, 
    outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, height = 0, 
    seed = 2026), alpha = 0.5, size = 1.5) + scale_fill_manual(values = c(WT = "#2878B5", MUT = "#D73027")) + 
    labs(x = "CTNNB1 status", y = paste("BAAT", CPTAC.expression_label), title = "BAAT expression in CPTAC / CHCC-HBV", 
        subtitle = sprintf("Wilcoxon P = %.3g; rank-biserial = %.3f", CPTAC.BAAT_test$P_value, CPTAC.BAAT_test$rank_biserial)) + 
    theme_classic(base_size = 13) + theme(legend.position = "none")

CPTAC.save_plot(CPTAC.p_BAAT, "15_BAAT_CTNNB1_boxplot.pdf")

CPTAC.settings <- list(cohort = "CPTAC / CHCC-HBV", n_patients = length(CPTAC.patient_ids), expression_scale = CPTAC.expression_label, 
    expression_file = CPTAC.count_file, mutation_source = c(published_table = CPTAC.xlsx_file, T777 = CPTAC.wes_variants_file), 
    mutation_source_md5 = CPTAC.wes_md5, mutation_sample_qc = CPTAC.wes_qc, excluded_patient = character(0), 
    mutation_frequency_rule = "> 0.05", variant_functions = CPTAC.variant_functions, OS_time_unit_assumed = CPTAC.OS_unit, 
    bootstrap_B = CPTAC.bootstrap_B, ssgsea = list(minSize = 2, maxSize = 500, alpha = 0.25, normalize = FALSE), 
    model_refitted = FALSE)

saveRDS(list(settings = CPTAC.settings, patient_ids = CPTAC.patient_ids, pairing_QC = CPTAC.selected_mapping, 
    background_QC = CPTAC.background_QC, gene_coverage = CPTAC.BA4_gene_coverage, symbol_mapping = CPTAC.symbol_mapping_used, 
    BA4_features = CPTAC.BA4_features, analysis = CPTAC.analysis, mutation_matrix = CPTAC.mut, mutation_frequency = CPTAC.mutation_frequency, 
    mutation_enrichment = CPTAC.mutation_enrichment, mutation_score_effects = CPTAC.gene_BA4_effects, 
    mutation_integrated = CPTAC.mutation_integrated, cox_models = CPTAC.cox_models, PH_tests = CPTAC.PH_tests, 
    CTNNB1_feature_test = CPTAC.CTNNB1_feature_test, component_decomposition = CPTAC.CTNNB1_component_decomposition, 
    BA_gene_test = CPTAC.CTNNB1_BA_gene_test, BAAT_test = CPTAC.BAAT_test, BAAT_correlations = CPTAC.BAAT_correlation_results), 
    file.path(CPTAC.out_dir, "CPTAC_159_validation_results.rds"))

capture.output(sessionInfo(), file = file.path(CPTAC.out_dir, "sessionInfo.txt"))


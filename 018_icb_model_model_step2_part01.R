if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

library(ggplot2)

library(pROC)

stopifnot(nrow(BA4_OOCO_predictions) == 69)

stopifnot(anyDuplicated(BA4_OOCO_predictions$patient_uid) == 0)

stopifnot(all(BA4_OOCO_predictions$response %in% c(0L, 1L)))

stopifnot(all(is.finite(BA4_OOCO_predictions$OOCO_ICBScore)))

BA4_OOCO_predictions <- BA4_OOCO_predictions %>% mutate(response_label = factor(ifelse(response == 1, 
    "Responder", "Non-responder"), levels = c("Non-responder", "Responder")))

cohorts_final <- c("GSE181946", "GSE202069", "PRJEB34724")

BA4_auc_ci_list <- list()

for (cc in cohorts_final) {
    df <- BA4_OOCO_predictions %>% filter(cohort == cc)
    roc_obj <- pROC::roc(response = df$response, predictor = df$OOCO_ICBScore, levels = c(0, 1), direction = "<", 
        quiet = TRUE)
    ci_delong <- pROC::ci.auc(roc_obj, method = "delong")
    set.seed(20260922)
    ci_boot <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
    BA4_auc_ci_list[[cc]] <- data.frame(cohort = cc, n = nrow(df), responders = sum(df$response == 1), 
        nonresponders = sum(df$response == 0), AUC = as.numeric(pROC::auc(roc_obj)), DeLong_low = as.numeric(ci_delong[1]), 
        DeLong_high = as.numeric(ci_delong[3]), Bootstrap_low = as.numeric(ci_boot[1]), Bootstrap_high = as.numeric(ci_boot[3]))
}

BA4_auc_ci <- bind_rows(BA4_auc_ci_list)

BA4_auc_ci

BA4_auc_ci_display <- BA4_auc_ci %>% mutate(AUC_CI_bootstrap = paste0(sprintf("%.3f", AUC), " (", sprintf("%.3f", 
    Bootstrap_low), "–", sprintf("%.3f", Bootstrap_high), ")"), AUC_CI_DeLong = paste0(sprintf("%.3f", 
    AUC), " (", sprintf("%.3f", DeLong_low), "–", sprintf("%.3f", DeLong_high), ")")) %>% select(cohort, 
    n, responders, nonresponders, AUC_CI_bootstrap, AUC_CI_DeLong)

BA4_auc_ci_display

BA4_score_summary <- BA4_OOCO_predictions %>% group_by(cohort, response_label) %>% summarise(n = n(), 
    mean_score = mean(OOCO_ICBScore), sd_score = sd(OOCO_ICBScore), median_score = median(OOCO_ICBScore), 
    Q1 = quantile(OOCO_ICBScore, 0.25), Q3 = quantile(OOCO_ICBScore, 0.75), .groups = "drop")

BA4_score_summary

BA4_score_wilcox <- bind_rows(lapply(cohorts_final, function(cc) {
    df <- BA4_OOCO_predictions %>% filter(cohort == cc)
    wt <- wilcox.test(OOCO_ICBScore ~ response, data = df, exact = FALSE)
    data.frame(cohort = cc, median_NR = median(df$OOCO_ICBScore[df$response == 0]), median_R = median(df$OOCO_ICBScore[df$response == 
        1]), delta_median = median(df$OOCO_ICBScore[df$response == 1]) - median(df$OOCO_ICBScore[df$response == 
        0]), Wilcoxon_P = wt$p.value)
}))

BA4_score_wilcox

p_score <- ggplot(BA4_OOCO_predictions, aes(x = response_label, y = OOCO_ICBScore)) + geom_boxplot(width = 0.55, 
    outlier.shape = NA) + geom_jitter(width = 0.12, height = 0, size = 2, alpha = 0.75) + facet_wrap(~cohort, 
    scales = "free_y", nrow = 1) + labs(x = NULL, y = "OOCO BA4 ICB Score", title = "Out-of-cohort prediction of anti-PD-1 response") + 
    theme_classic(base_size = 13) + theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 12))

p_score

ggsave(filename = "STEP2_BA4_OOCO_score_distribution.pdf", plot = p_score, width = 10, height = 4.5)

BA4_roc_list <- list()

BA4_roc_curve_list <- list()

for (cc in cohorts_final) {
    df <- BA4_OOCO_predictions %>% filter(cohort == cc)
    roc_obj <- pROC::roc(response = df$response, predictor = df$OOCO_ICBScore, levels = c(0, 1), direction = "<", 
        quiet = TRUE)
    BA4_roc_list[[cc]] <- roc_obj
    roc_coords <- pROC::coords(roc_obj, x = "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
    roc_df <- data.frame(cohort = cc, FPR = 1 - roc_coords$specificity, TPR = roc_coords$sensitivity)
    BA4_roc_curve_list[[cc]] <- roc_df
}

BA4_roc_curve_df <- bind_rows(BA4_roc_curve_list)

roc_legend_labels <- setNames(paste0(BA4_auc_ci$cohort, "  AUC=", sprintf("%.3f", BA4_auc_ci$AUC)), BA4_auc_ci$cohort)

p_roc <- ggplot(BA4_roc_curve_df, aes(x = FPR, y = TPR, color = cohort)) + geom_step(linewidth = 1) + 
    geom_abline(intercept = 0, slope = 1, linetype = 2) + coord_equal() + scale_x_continuous(limits = c(0, 
    1)) + scale_y_continuous(limits = c(0, 1)) + scale_color_discrete(labels = roc_legend_labels) + labs(x = "False positive rate", 
    y = "True positive rate", color = "Cohort", title = "Locked BA4 model: leave-one-cohort-out ROC") + 
    theme_classic(base_size = 13)

p_roc

ggsave(filename = "STEP2_BA4_OOCO_ROC.pdf", plot = p_roc, width = 7, height = 6)

make_pr_curve <- function(y, score, cohort_name) {
    ord <- order(score, decreasing = TRUE)
    y_sorted <- y[ord]
    TP <- cumsum(y_sorted == 1)
    FP <- cumsum(y_sorted == 0)
    total_positive <- sum(y == 1)
    recall <- TP/total_positive
    precision <- TP/(TP + FP)
    data.frame(cohort = cohort_name, recall = c(0, recall), precision = c(1, precision))
}

BA4_pr_curve_df <- bind_rows(lapply(cohorts_final, function(cc) {
    df <- BA4_OOCO_predictions %>% filter(cohort == cc)
    make_pr_curve(y = df$response, score = df$OOCO_ICBScore, cohort_name = cc)
}))

prevalence_df <- BA4_OOCO_predictions %>% group_by(cohort) %>% summarise(prevalence = mean(response), 
    .groups = "drop")

p_pr <- ggplot(BA4_pr_curve_df, aes(x = recall, y = precision)) + geom_step(linewidth = 1) + geom_hline(data = prevalence_df, 
    aes(yintercept = prevalence), linetype = 2) + facet_wrap(~cohort, nrow = 1) + scale_x_continuous(limits = c(0, 
    1)) + scale_y_continuous(limits = c(0, 1)) + labs(x = "Recall", y = "Precision", title = "Locked BA4 model: out-of-cohort precision-recall curves") + 
    theme_classic(base_size = 13)

p_pr

ggsave(filename = "STEP2_BA4_OOCO_PR_curves.pdf", plot = p_pr, width = 10, height = 4.5)

BA4_auc_ci_plot <- BA4_auc_ci %>% mutate(cohort = factor(cohort, levels = rev(cohorts_final)))

p_auc_ci <- ggplot(BA4_auc_ci_plot, aes(x = AUC, y = cohort)) + geom_vline(xintercept = 0.5, linetype = 2) + 
    geom_errorbarh(aes(xmin = Bootstrap_low, xmax = Bootstrap_high), height = 0.15) + geom_point(size = 3) + 
    scale_x_continuous(limits = c(0, 1)) + labs(x = "ROC-AUC (bootstrap 95% CI)", y = NULL, title = "Out-of-cohort discrimination of the locked BA4 model") + 
    theme_classic(base_size = 13)

p_auc_ci

ggsave(filename = "STEP2_BA4_cohort_AUC_95CI.pdf", plot = p_auc_ci, width = 7, height = 4.5)

BA4_cohort_score_summary <- BA4_OOCO_predictions %>% group_by(cohort) %>% summarise(mean = mean(OOCO_ICBScore), 
    sd = sd(OOCO_ICBScore), median = median(OOCO_ICBScore), Q1 = quantile(OOCO_ICBScore, 0.25), Q3 = quantile(OOCO_ICBScore, 
        0.75), min = min(OOCO_ICBScore), max = max(OOCO_ICBScore), .groups = "drop")

BA4_cohort_score_summary

write.csv(BA4_auc_ci, "STEP2_BA4_cohort_AUC_95CI.csv", row.names = FALSE)

write.csv(BA4_score_summary, "STEP2_BA4_score_summary_by_response.csv", row.names = FALSE)

write.csv(BA4_score_wilcox, "STEP2_BA4_score_Wilcoxon.csv", row.names = FALSE)

write.csv(BA4_cohort_score_summary, "STEP2_BA4_cohort_score_scale_summary.csv", row.names = FALSE)

saveRDS(BA4_roc_list, "STEP2_BA4_ROC_objects.rds")

BA4_auc_ci

BA4_score_summary

BA4_score_wilcox

BA4_cohort_score_summary


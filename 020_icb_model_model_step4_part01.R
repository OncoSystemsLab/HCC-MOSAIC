if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

library(ggplot2)

boot_data <- BA4_OOCO_predictions %>% select(patient_uid, cohort, sample_id, response, OOCO_ICBScore)

stopifnot(nrow(boot_data) == 69)

stopifnot(anyDuplicated(boot_data$patient_uid) == 0)

table(boot_data$cohort, boot_data$response)

stratified_bootstrap_sample <- function(data, seed) {
    set.seed(seed)
    boot_df <- data %>% group_by(cohort, response) %>% slice_sample(prop = 1, replace = TRUE) %>% ungroup()
    return(boot_df)
}

test_boot <- stratified_bootstrap_sample(data = boot_data, seed = 20260924)

nrow(test_boot)

table(test_boot$cohort, test_boot$response)

run_one_performance_bootstrap <- function(iteration, data, seed_base = 20260924) {
    boot_df <- stratified_bootstrap_sample(data = data, seed = seed_base + iteration)
    cohorts <- c("GSE181946", "GSE202069", "PRJEB34724")
    cohort_metrics <- list()
    for (cc in cohorts) {
        df_cc <- boot_df %>% filter(cohort == cc)
        m <- calculate_metrics(y = df_cc$response, score = df_cc$OOCO_ICBScore)
        m$cohort <- cc
        cohort_metrics[[cc]] <- m
    }
    cohort_metrics <- bind_rows(cohort_metrics)
    pooled_m <- calculate_metrics(y = boot_df$response, score = boot_df$OOCO_ICBScore)
    data.frame(iteration = iteration, Macro_AUC = mean(cohort_metrics$AUC), Macro_AP = mean(cohort_metrics$AP), 
        Macro_APnorm = mean(cohort_metrics$AP_norm), Min_AUC = min(cohort_metrics$AUC), SD_AUC = sd(cohort_metrics$AUC), 
        Macro_Brier = mean(cohort_metrics$Brier), Macro_Brier_skill = mean(cohort_metrics$Brier_skill), 
        Pooled_AUC = pooled_m$AUC, Pooled_AP = pooled_m$AP, Pooled_APnorm = pooled_m$AP_norm, Pooled_Brier = pooled_m$Brier, 
        AUC_GSE181946 = cohort_metrics$AUC[cohort_metrics$cohort == "GSE181946"], AUC_GSE202069 = cohort_metrics$AUC[cohort_metrics$cohort == 
            "GSE202069"], AUC_PRJEB34724 = cohort_metrics$AUC[cohort_metrics$cohort == "PRJEB34724"])
}

bootstrap_test_1 <- run_one_performance_bootstrap(iteration = 1, data = boot_data)

bootstrap_test_1

B_boot <- 2000

BA4_bootstrap_list <- vector("list", B_boot)

for (b in seq_len(B_boot)) {
    BA4_bootstrap_list[[b]] <- run_one_performance_bootstrap(iteration = b, data = boot_data, seed_base = 20260924)
    if (b%%200 == 0) {
        cat("Bootstrap:", b, "/", B_boot, "")
    }
}

BA4_bootstrap_results <- bind_rows(BA4_bootstrap_list)

dim(BA4_bootstrap_results)

summary(BA4_bootstrap_results$Macro_AUC)

observed_metrics_boot <- c(Macro_AUC = observed_Macro_AUC, Macro_AP = BA4_macro_metrics$Macro_AP[1], 
    Macro_APnorm = observed_Macro_APnorm, Min_AUC = observed_Min_AUC, SD_AUC = BA4_macro_metrics$SD_AUC[1], 
    Macro_Brier = observed_Macro_Brier, Macro_Brier_skill = BA4_macro_metrics$Macro_Brier_skill[1], Pooled_AUC = observed_Pooled_AUC, 
    Pooled_AP = BA4_pooled_metrics$AP[1], Pooled_APnorm = BA4_pooled_metrics$AP_norm[1], Pooled_Brier = BA4_pooled_metrics$Brier[1])

observed_metrics_boot

bootstrap_ci <- function(x, conf = 0.95) {
    alpha <- 1 - conf
    q <- quantile(x, probs = c(alpha/2, 0.5, 1 - alpha/2), na.rm = TRUE, names = FALSE)
    names(q) <- c("low", "median", "high")
    return(q)
}

primary_metrics <- c("Macro_AUC", "Macro_APnorm", "Min_AUC", "Macro_Brier", "Macro_Brier_skill")

BA4_bootstrap_primary_summary <- bind_rows(lapply(primary_metrics, function(metric_name) {
    x <- BA4_bootstrap_results[[metric_name]]
    ci <- bootstrap_ci(x)
    data.frame(Metric = metric_name, Observed = observed_metrics_boot[metric_name], Bootstrap_mean = mean(x, 
        na.rm = TRUE), Bootstrap_median = median(x, na.rm = TRUE), CI_low = unname(ci["low"]), CI_high = unname(ci["high"]), 
        Bootstrap_SD = sd(x, na.rm = TRUE))
}))

BA4_bootstrap_primary_summary

secondary_metrics <- c("Pooled_AUC", "Pooled_AP", "Pooled_APnorm", "Pooled_Brier")

BA4_bootstrap_secondary_summary <- bind_rows(lapply(secondary_metrics, function(metric_name) {
    x <- BA4_bootstrap_results[[metric_name]]
    ci <- bootstrap_ci(x)
    data.frame(Metric = metric_name, Observed = observed_metrics_boot[metric_name], Bootstrap_mean = mean(x, 
        na.rm = TRUE), Bootstrap_median = median(x, na.rm = TRUE), CI_low = unname(ci["low"]), CI_high = unname(ci["high"]), 
        Bootstrap_SD = sd(x, na.rm = TRUE))
}))

BA4_bootstrap_secondary_summary

bootstrap_ci(BA4_bootstrap_results$Macro_AUC)

quantile(BA4_bootstrap_results$Macro_AUC, probs = c(0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975))

mean(BA4_bootstrap_results$Macro_AUC > 0.5)

mean(BA4_bootstrap_results$Min_AUC > 0.5)

BA4_bootstrap_stability_simple <- data.frame(Metric = c("P(Macro_AUC > 0.5)", "P(Min_AUC > 0.5)", "P(Macro_APnorm > 0)"), 
    Value = c(mean(BA4_bootstrap_results$Macro_AUC > 0.5), mean(BA4_bootstrap_results$Min_AUC > 0.5), 
        mean(BA4_bootstrap_results$Macro_APnorm > 0)))

BA4_bootstrap_stability_simple

cohort_auc_variables <- c(GSE181946 = "AUC_GSE181946", GSE202069 = "AUC_GSE202069", PRJEB34724 = "AUC_PRJEB34724")

BA4_bootstrap_cohort_auc <- bind_rows(lapply(names(cohort_auc_variables), function(cc) {
    varname <- cohort_auc_variables[[cc]]
    x <- BA4_bootstrap_results[[varname]]
    ci <- bootstrap_ci(x)
    observed_auc <- BA4_outer_metrics$AUC[BA4_outer_metrics$cohort == cc]
    data.frame(cohort = cc, Observed_AUC = observed_auc, Bootstrap_mean = mean(x, na.rm = TRUE), CI_low = unname(ci["low"]), 
        CI_high = unname(ci["high"]))
}))

BA4_bootstrap_cohort_auc

p_boot_auc <- ggplot(BA4_bootstrap_results, aes(x = Macro_AUC)) + geom_histogram(bins = 40) + geom_vline(xintercept = observed_Macro_AUC, 
    linetype = 2, linewidth = 1) + labs(x = "Bootstrap macro ROC-AUC", y = "Bootstrap samples", title = "Sampling uncertainty of the locked BA4 model") + 
    theme_classic(base_size = 13)

p_boot_auc

p_boot_minauc <- ggplot(BA4_bootstrap_results, aes(x = Min_AUC)) + geom_histogram(bins = 40) + geom_vline(xintercept = observed_Min_AUC, 
    linetype = 2, linewidth = 1) + geom_vline(xintercept = 0.5, linetype = 3) + labs(x = "Bootstrap minimum cohort ROC-AUC", 
    y = "Bootstrap samples", title = "Worst-cohort robustness of the locked BA4 model") + theme_classic(base_size = 13)

p_boot_minauc

write.csv(BA4_bootstrap_results, "STEP4_BA4_bootstrap_2000_raw.csv", row.names = FALSE)

write.csv(BA4_bootstrap_primary_summary, "STEP4_BA4_bootstrap_primary_summary.csv", row.names = FALSE)

write.csv(BA4_bootstrap_secondary_summary, "STEP4_BA4_bootstrap_secondary_summary.csv", row.names = FALSE)

write.csv(BA4_bootstrap_cohort_auc, "STEP4_BA4_bootstrap_cohort_AUC.csv", row.names = FALSE)

write.csv(BA4_bootstrap_stability_simple, "STEP4_BA4_bootstrap_stability_summary.csv", row.names = FALSE)

saveRDS(BA4_bootstrap_results, "STEP4_BA4_bootstrap_2000_results.rds")

ggsave("STEP4_BA4_bootstrap_MacroAUC.pdf", p_boot_auc, width = 7, height = 5)

ggsave("STEP4_BA4_bootstrap_MinAUC.pdf", p_boot_minauc, width = 7, height = 5)

BA4_bootstrap_primary_summary

BA4_bootstrap_cohort_auc

BA4_bootstrap_stability_simple

summary(BA4_bootstrap_results$Macro_AUC)


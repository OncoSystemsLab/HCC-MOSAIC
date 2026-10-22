if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(glmnet)

library(dplyr)

cohorts_final <- c("GSE181946", "GSE202069", "PRJEB34724")

BA4_reference_coefficients <- bind_rows(lapply(cohorts_final, function(cc) {
    lambda_use <- BA4_outer_metrics$selected_lambda[BA4_outer_metrics$cohort == cc]
    fit <- BA4_IECV$models[[cc]]
    cf <- as.matrix(coef(fit, s = lambda_use))
    data.frame(heldout_cohort = cc, feature = rownames(cf), coefficient = as.numeric(cf[, 1]))
}))

BA4_reference_coefficients

BA4_reference_coefficients %>% filter(feature != "(Intercept)")

model_data_BA4 <- model_data_BA4 %>% mutate(cohort = trimws(as.character(cohort)), sample_id = trimws(as.character(sample_id)), 
    patient_uid = paste(cohort, sample_id, sep = "__"))

BA4_OOCO_predictions <- BA4_OOCO_predictions %>% mutate(cohort = trimws(as.character(cohort)), sample_id = trimws(as.character(sample_id)), 
    patient_uid = paste(cohort, sample_id, sep = "__"))

stopifnot(anyDuplicated(model_data_BA4$patient_uid) == 0)

stopifnot(anyDuplicated(BA4_OOCO_predictions$patient_uid) == 0)

stopifnot(setequal(model_data_BA4$patient_uid, BA4_OOCO_predictions$patient_uid))

for (cc in cohorts_final) {
    ids_model <- model_data_BA4 %>% filter(cohort == cc) %>% pull(patient_uid)
    ids_pred <- BA4_OOCO_predictions %>% filter(cohort == cc) %>% pull(patient_uid)
    cat("", cc, "")
    cat("model N =", length(ids_model), "")
    cat("prediction N =", length(ids_pred), "")
    cat("same patient set =", setequal(ids_model, ids_pred), "")
}

bootstrap_training_data <- function(data, seed) {
    set.seed(seed)
    strata <- interaction(data$cohort, data$response, drop = TRUE)
    idx_list <- split(seq_len(nrow(data)), strata)
    boot_idx <- unlist(lapply(idx_list, function(idx) {
        sample(idx, size = length(idx), replace = TRUE)
    }), use.names = FALSE)
    boot_data <- data[boot_idx, , drop = FALSE]
    return(boot_data)
}

test_train <- model_data_BA4 %>% filter(cohort != "GSE181946")

table(test_train$cohort, test_train$response)

test_train_boot <- bootstrap_training_data(data = test_train, seed = 20260925)

table(test_train_boot$cohort, test_train_boot$response)

table(test_train$cohort, test_train$response)

anyDuplicated(test_train_boot$patient_uid)

run_one_training_bootstrap <- function(iteration, heldout_name, data, features, lambda_use, reference_predictions, 
    seed_base = 20260925) {
    train_df <- data %>% filter(.data$cohort != .env$heldout_name)
    test_df <- data %>% filter(.data$cohort == .env$heldout_name) %>% mutate(cohort = trimws(as.character(cohort)), 
        sample_id = trimws(as.character(sample_id)), patient_uid = paste(cohort, sample_id, sep = "__")) %>% 
        arrange(patient_uid)
    boot_train <- bootstrap_training_data(data = train_df, seed = seed_base + iteration)
    X_train <- boot_train %>% select(all_of(features)) %>% as.matrix()
    storage.mode(X_train) <- "numeric"
    y_train <- as.integer(boot_train$response)
    cohort_train <- as.character(boot_train$cohort)
    X_test <- test_df %>% select(all_of(features)) %>% as.matrix()
    storage.mode(X_test) <- "numeric"
    y_test <- as.integer(test_df$response)
    weights_train <- make_balanced_weights(y = y_train, cohort = cohort_train)
    fit <- glmnet(x = X_train, y = y_train, family = "binomial", alpha = 0, lambda = lambda_use, weights = weights_train, 
        standardize = TRUE, intercept = TRUE, maxit = 1e+05)
    pred <- as.numeric(predict(fit, newx = X_test, type = "response", s = lambda_use))
    ref_df <- reference_predictions %>% filter(.data$cohort == .env$heldout_name) %>% transmute(cohort = trimws(as.character(cohort)), 
        sample_id = trimws(as.character(sample_id)), patient_uid = paste(cohort, sample_id, sep = "__"), 
        reference_score = OOCO_ICBScore)
    stopifnot(anyDuplicated(test_df$patient_uid) == 0)
    stopifnot(anyDuplicated(ref_df$patient_uid) == 0)
    stopifnot(nrow(test_df) == nrow(ref_df))
    stopifnot(setequal(test_df$patient_uid, ref_df$patient_uid))
    match_idx <- match(test_df$patient_uid, ref_df$patient_uid)
    stopifnot(!anyNA(match_idx))
    ref_score <- ref_df$reference_score[match_idx]
    m <- calculate_metrics(y = y_test, score = pred)
    pearson_r <- suppressWarnings(cor(pred, ref_score, method = "pearson"))
    spearman_r <- suppressWarnings(cor(pred, ref_score, method = "spearman"))
    score_rmse <- sqrt(mean((pred - ref_score)^2))
    cf <- as.matrix(coef(fit, s = lambda_use))
    coef_df <- data.frame(iteration = iteration, heldout_cohort = heldout_name, feature = rownames(cf), 
        coefficient = as.numeric(cf[, 1]))
    prediction_df <- data.frame(iteration = iteration, heldout_cohort = heldout_name, patient_uid = test_df$patient_uid, 
        response = test_df$response, reference_score = ref_score, bootstrap_score = pred)
    summary_df <- data.frame(iteration = iteration, heldout_cohort = heldout_name, AUC = m$AUC, AP = m$AP, 
        AP_norm = m$AP_norm, Brier = m$Brier, Pearson_r = pearson_r, Spearman_r = spearman_r, Score_RMSE = score_rmse)
    return(list(summary = summary_df, coefficients = coef_df, predictions = prediction_df))
}

test_heldout <- "GSE181946"

test_lambda <- BA4_outer_metrics$selected_lambda[BA4_outer_metrics$cohort == test_heldout]

training_boot_test <- run_one_training_bootstrap(iteration = 1, heldout_name = test_heldout, data = model_data_BA4, 
    features = LOCKED_FEATURES, lambda_use = test_lambda, reference_predictions = BA4_OOCO_predictions, 
    seed_base = 20260925)

training_boot_test$summary

training_boot_test$coefficients

B_stability <- 1000

training_stability_summary <- list()

training_stability_coef <- list()

training_stability_pred <- list()

counter <- 1

for (cc_i in seq_along(cohorts_final)) {
    cc <- cohorts_final[cc_i]
    lambda_use <- BA4_outer_metrics$selected_lambda[BA4_outer_metrics$cohort == cc]
    cat("\n====================================", "Training stability held-out:", cc, "", "Lambda:", 
        lambda_use, "", "====================================")
    for (b in seq_len(B_stability)) {
        tmp <- run_one_training_bootstrap(iteration = b, heldout_name = cc, data = model_data_BA4, features = LOCKED_FEATURES, 
            lambda_use = lambda_use, reference_predictions = BA4_OOCO_predictions, seed_base = 20260925 + 
                cc_i * 1e+05)
        training_stability_summary[[counter]] <- tmp$summary
        training_stability_coef[[counter]] <- tmp$coefficients
        training_stability_pred[[counter]] <- tmp$predictions
        counter <- counter + 1
        if (b%%100 == 0) {
            cat(cc, ":", b, "/", B_stability, "")
        }
    }
}

BA4_training_stability_summary <- bind_rows(training_stability_summary)

BA4_training_stability_coef <- bind_rows(training_stability_coef)

BA4_training_stability_predictions <- bind_rows(training_stability_pred)

dim(BA4_training_stability_summary)

dim(BA4_training_stability_coef)

dim(BA4_training_stability_predictions)

BA4_prediction_stability <- BA4_training_stability_summary %>% group_by(heldout_cohort) %>% summarise(n_boot = n(), 
    AUC_mean = mean(AUC, na.rm = TRUE), AUC_median = median(AUC, na.rm = TRUE), AUC_low = quantile(AUC, 
        0.025, na.rm = TRUE), AUC_high = quantile(AUC, 0.975, na.rm = TRUE), Spearman_median = median(Spearman_r, 
        na.rm = TRUE), Spearman_Q1 = quantile(Spearman_r, 0.25, na.rm = TRUE), Spearman_Q3 = quantile(Spearman_r, 
        0.75, na.rm = TRUE), Pearson_median = median(Pearson_r, na.rm = TRUE), RMSE_median = median(Score_RMSE, 
        na.rm = TRUE), pct_AUC_gt_05 = mean(AUC > 0.5, na.rm = TRUE), pct_Spearman_gt_08 = mean(Spearman_r > 
        0.8, na.rm = TRUE), .groups = "drop")

BA4_prediction_stability

BA4_coef_only <- BA4_training_stability_coef %>% filter(feature != "(Intercept)")

BA4_coefficient_stability <- BA4_coef_only %>% group_by(heldout_cohort, feature) %>% summarise(n_boot = n(), 
    median_coef = median(coefficient, na.rm = TRUE), mean_coef = mean(coefficient, na.rm = TRUE), CI_low = quantile(coefficient, 
        0.025, na.rm = TRUE), CI_high = quantile(coefficient, 0.975, na.rm = TRUE), pct_positive = mean(coefficient > 
        0, na.rm = TRUE), pct_negative = mean(coefficient < 0, na.rm = TRUE), sign_stability = pmax(pct_positive, 
        pct_negative), .groups = "drop")

BA4_coefficient_stability

BA4_patient_prediction_stability <- BA4_training_stability_predictions %>% group_by(heldout_cohort, patient_uid, 
    response) %>% summarise(reference_score = first(reference_score), bootstrap_mean = mean(bootstrap_score, 
    na.rm = TRUE), bootstrap_median = median(bootstrap_score, na.rm = TRUE), bootstrap_SD = sd(bootstrap_score, 
    na.rm = TRUE), CI_low = quantile(bootstrap_score, 0.025, na.rm = TRUE), CI_high = quantile(bootstrap_score, 
    0.975, na.rm = TRUE), .groups = "drop")

summary(BA4_training_stability_summary$Spearman_r)

BA4_training_stability_summary %>% group_by(heldout_cohort) %>% summarise(min = min(Spearman_r, na.rm = TRUE), 
    Q1 = quantile(Spearman_r, 0.25, na.rm = TRUE), median = median(Spearman_r, na.rm = TRUE), Q3 = quantile(Spearman_r, 
        0.75, na.rm = TRUE), max = max(Spearman_r, na.rm = TRUE), .groups = "drop")

p_coef_stability <- ggplot(BA4_coef_only, aes(x = coefficient, y = feature)) + geom_vline(xintercept = 0, 
    linetype = 2) + geom_boxplot(outlier.shape = NA) + facet_wrap(~heldout_cohort, scales = "free_x") + 
    labs(x = "Ridge coefficient across bootstrap training sets", y = NULL, title = "Coefficient stability of the locked BA4 model") + 
    theme_classic(base_size = 12)

p_coef_stability

p_prediction_stability <- ggplot(BA4_training_stability_summary, aes(x = Spearman_r)) + geom_histogram(bins = 40) + 
    facet_wrap(~heldout_cohort, nrow = 1) + geom_vline(xintercept = 0.8, linetype = 2) + labs(x = "Spearman correlation with original OOCO prediction", 
    y = "Bootstrap training models", title = "Prediction stability under training-set resampling") + 
    theme_classic(base_size = 12)

p_prediction_stability

write.csv(BA4_prediction_stability, "STEP5_BA4_prediction_stability.csv", row.names = FALSE)

write.csv(BA4_coefficient_stability, "STEP5_BA4_coefficient_stability.csv", row.names = FALSE)

write.csv(BA4_patient_prediction_stability, "STEP5_BA4_patient_prediction_stability.csv", row.names = FALSE)

saveRDS(BA4_training_stability_summary, "STEP5_BA4_training_stability_summary.rds")

saveRDS(BA4_training_stability_coef, "STEP5_BA4_training_stability_coefficients.rds")

saveRDS(BA4_training_stability_predictions, "STEP5_BA4_training_stability_predictions.rds")

ggsave("STEP5_BA4_coefficient_stability.pdf", p_coef_stability, width = 11, height = 5)

ggsave("STEP5_BA4_prediction_stability.pdf", p_prediction_stability, width = 10, height = 4.5)

BA4_prediction_stability

BA4_coefficient_stability

summary(BA4_training_stability_summary$Spearman_r)

BA4_patient_prediction_stability %>% summarise(median_patient_SD = median(bootstrap_SD), max_patient_SD = max(bootstrap_SD), 
    median_CI_width = median(CI_high - CI_low))


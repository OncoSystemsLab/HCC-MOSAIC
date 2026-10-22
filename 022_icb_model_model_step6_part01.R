if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(glmnet)

library(dplyr)

outer_lambdas <- BA4_outer_metrics$selected_lambda

outer_lambdas

lambda_final <- exp(median(log(outer_lambdas)))

lambda_final

BA4_final_lambda_info <- data.frame(cohort = BA4_outer_metrics$cohort, outer_lambda = BA4_outer_metrics$selected_lambda)

BA4_final_lambda_info

X_final <- model_data_BA4 %>% select(all_of(LOCKED_FEATURES)) %>% as.matrix()

storage.mode(X_final) <- "numeric"

y_final <- as.integer(model_data_BA4$response)

cohort_final <- as.character(model_data_BA4$cohort)

dim(X_final)

table(y_final)

table(cohort_final, y_final)

stopifnot(nrow(X_final) == 69)

stopifnot(ncol(X_final) == 4)

stopifnot(all(is.finite(X_final)))

stopifnot(identical(colnames(X_final), LOCKED_FEATURES))

weights_final <- make_balanced_weights(y = y_final, cohort = cohort_final)

summary(weights_final)

tapply(weights_final, cohort_final, sum)

data.frame(cohort = cohort_final, response = y_final, weight = weights_final) %>% group_by(cohort, response) %>% 
    summarise(n = n(), total_weight = sum(weight), .groups = "drop")

FINAL_BA4_MODEL <- glmnet(x = X_final, y = y_final, family = "binomial", alpha = 0, lambda = lambda_final, 
    weights = weights_final, standardize = TRUE, intercept = TRUE, maxit = 1e+05)

final_coef_matrix <- as.matrix(coef(FINAL_BA4_MODEL, s = lambda_final))

BA4_final_coefficients <- data.frame(feature = rownames(final_coef_matrix), coefficient = as.numeric(final_coef_matrix[, 
    1]))

BA4_final_coefficients

BA4_final_coefficients %>% mutate(direction = case_when(coefficient > 0 ~ "Responder-high", coefficient < 
    0 ~ "Non-responder-high", TRUE ~ "Zero"))

feature_sd <- apply(X_final, 2, sd)

BA4_standardized_effects <- BA4_final_coefficients %>% filter(feature != "(Intercept)") %>% mutate(feature_SD = feature_sd[feature], 
    standardized_coefficient = coefficient * feature_SD, direction = ifelse(standardized_coefficient > 
        0, "Responder-high", "Non-responder-high"))

BA4_standardized_effects

BA4_final_fullfit_score <- as.numeric(predict(FINAL_BA4_MODEL, newx = X_final, type = "response", s = lambda_final))

BA4_final_training_scores <- model_data_BA4 %>% select(patient_uid, cohort, sample_id, response) %>% 
    mutate(BA4_FinalScore = BA4_final_fullfit_score)

score_compare <- BA4_final_training_scores %>% inner_join(BA4_OOCO_predictions %>% select(patient_uid, 
    OOCO_ICBScore), by = "patient_uid")

nrow(score_compare)

cor(score_compare$BA4_FinalScore, score_compare$OOCO_ICBScore, method = "spearman")

cor(score_compare$BA4_FinalScore, score_compare$OOCO_ICBScore, method = "pearson")

BA4_final_apparent_metrics <- calculate_metrics(y = y_final, score = BA4_final_fullfit_score)

BA4_final_apparent_metrics

beta0 <- BA4_final_coefficients$coefficient[BA4_final_coefficients$feature == "(Intercept)"]

beta_vec <- BA4_final_coefficients %>% filter(feature != "(Intercept)")

beta_numeric <- beta_vec$coefficient

names(beta_numeric) <- beta_vec$feature

manual_eta <- beta0 + X_final[, names(beta_numeric), drop = FALSE] %*% beta_numeric

manual_score <- as.numeric(plogis(manual_eta))

max(abs(manual_score - BA4_final_fullfit_score))

predict_BA4_ICBScore <- function(newdata, model = FINAL_BA4_MODEL, features = LOCKED_FEATURES, lambda = lambda_final) {
    missing_features <- setdiff(features, colnames(newdata))
    if (length(missing_features) > 0) {
        stop("Missing BA4 features: ", paste(missing_features, collapse = ", "))
    }
    X_new <- newdata[, features, drop = FALSE] %>% as.matrix()
    storage.mode(X_new) <- "numeric"
    if (any(!is.finite(X_new))) {
        stop("newdata contains NA/Inf.")
    }
    score <- as.numeric(predict(model, newx = X_new, type = "response", s = lambda))
    return(score)
}

test_score <- predict_BA4_ICBScore(model_data_BA4)

max(abs(test_score - BA4_final_fullfit_score))

BA4_feature_reference <- bind_rows(lapply(LOCKED_FEATURES, function(f) {
    x <- model_data_BA4[[f]]
    data.frame(feature = f, mean = mean(x), SD = sd(x), median = median(x), Q1 = quantile(x, 0.25), Q3 = quantile(x, 
        0.75), min = min(x), max = max(x))
}))

BA4_feature_reference

BA4_gene_sets <- all_gene_sets[LOCKED_FEATURES]

names(BA4_gene_sets)

background_genes

BA4_MODEL_BUNDLE <- list(model_name = "HCC_BA4_ICBScore_v1", model_type = "L2-penalized logistic regression", 
    model = FINAL_BA4_MODEL, lambda = lambda_final, alpha = 0, features = LOCKED_FEATURES, coefficients = BA4_final_coefficients, 
    standardized_effects = BA4_standardized_effects, gene_sets = BA4_gene_sets, background_genes = background_genes, 
    feature_reference = BA4_feature_reference, preprocessing = list(expression_input = "TPM", transformation = "log2(TPM + 1)", 
        enrichment = "ssGSEA", ssgsea_minSize = 2, ssgsea_maxSize = 500, ssgsea_alpha = 0.25, ssgsea_normalize = FALSE), 
    development_cohorts = c("PRJEB34724", "GSE181946", "GSE202069"), n_training = 69, score_definition = paste("Continuous responder-like ICB score.", 
        "Higher score indicates more responder-like", "transcriptional state.", "Not an absolute response probability."))

saveRDS(BA4_MODEL_BUNDLE, file = "FINAL_HCC_BA4_ICBScore_v1_bundle.rds")

write.csv(BA4_final_coefficients, "FINAL_BA4_coefficients.csv", row.names = FALSE)

write.csv(BA4_standardized_effects, "FINAL_BA4_standardized_effects.csv", row.names = FALSE)

write.csv(BA4_feature_reference, "FINAL_BA4_training_feature_reference.csv", row.names = FALSE)

write.csv(BA4_final_training_scores, "FINAL_BA4_training_scores_69patients.csv", row.names = FALSE)

write.csv(BA4_final_lambda_info, "FINAL_BA4_lambda_source.csv", row.names = FALSE)

capture.output(sessionInfo(), file = "FINAL_BA4_R_sessionInfo.txt")

BA4_validation_overview <- data.frame(item = c("Development cohorts", "Patients", "Responders", "Non-responders", 
    "Locked features", "IECV Macro AUC", "IECV Min AUC", "IECV Macro APnorm", "Permutation P Macro AUC", 
    "Bootstrap P(Macro AUC > 0.5)", "Bootstrap P(Min AUC > 0.5)", "Median training-bootstrap Spearman", 
    "BA coefficient direction", "FA coefficient direction", "ROS coefficient direction", "Inflammatory4 coefficient direction"), 
    value = c("3", "69", "21", "48", "4", sprintf("%.3f", observed_Macro_AUC), sprintf("%.3f", observed_Min_AUC), 
        sprintf("%.3f", observed_Macro_APnorm), sprintf("%.4f", perm_p_Macro_AUC), sprintf("%.3f", BA4_bootstrap_stability_simple$Value[BA4_bootstrap_stability_simple$Metric == 
            "P(Macro_AUC > 0.5)"]), sprintf("%.3f", BA4_bootstrap_stability_simple$Value[BA4_bootstrap_stability_simple$Metric == 
            "P(Min_AUC > 0.5)"]), sprintf("%.3f", median(BA4_training_stability_summary$Spearman_r, na.rm = TRUE)), 
        "Negative", "Negative", "Positive", "Positive"))

BA4_validation_overview


if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(glmnet)

library(dplyr)

library(tibble)

LOCKED_FEATURES <- c("ICI__HCC_INFLAMMATORY4", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_BILE_ACID_METABOLISM", 
    "HALLMARK_FATTY_ACID_METABOLISM")

FINAL_MODEL_NAME <- "BA4_Ridge"

stopifnot(all(LOCKED_FEATURES %in% colnames(model_data)))

LOCKED_FEATURES

model_data_BA4 <- model_data %>% select(patient_uid, cohort, sample_id, response, all_of(LOCKED_FEATURES))

dim(model_data_BA4)

stopifnot(nrow(model_data_BA4) == 69)

stopifnot(anyDuplicated(model_data_BA4$patient_uid) == 0)

X_check <- model_data_BA4 %>% select(all_of(LOCKED_FEATURES)) %>% as.matrix()

storage.mode(X_check) <- "numeric"

stopifnot(all(is.finite(X_check)))

stopifnot(all(model_data_BA4$response %in% c(0L, 1L)))

make_balanced_weights <- function(y, cohort) {
    y <- as.integer(y)
    cohort <- as.character(cohort)
    cohorts <- unique(cohort)
    w <- numeric(length(y))
    for (cc in cohorts) {
        for (yy in c(0L, 1L)) {
            idx <- which(cohort == cc & y == yy)
            if (length(idx) > 0) {
                w[idx] <- 1/(length(cohorts) * 2 * length(idx))
            }
        }
    }
    w <- w/mean(w)
    return(w)
}

test_weights <- make_balanced_weights(y = model_data_BA4$response, cohort = model_data_BA4$cohort)

summary(test_weights)

tapply(test_weights, model_data_BA4$cohort, sum)

make_stratified_foldid <- function(y, cohort, K = 3, seed = 123) {
    set.seed(seed)
    y <- as.integer(y)
    cohort <- as.character(cohort)
    foldid <- integer(length(y))
    strata <- interaction(cohort, y, drop = TRUE)
    for (s in levels(strata)) {
        idx <- which(strata == s)
        idx <- sample(idx)
        foldid[idx] <- rep(seq_len(K), length.out = length(idx))
    }
    return(foldid)
}

test_foldid <- make_stratified_foldid(y = model_data_BA4$response, cohort = model_data_BA4$cohort, K = 3, 
    seed = 20260922)

table(test_foldid)

table(model_data_BA4$cohort, model_data_BA4$response, test_foldid)

lambda_grid_final <- 10^seq(2, -4, length.out = 80)

range(lambda_grid_final)

length(lambda_grid_final)

choose_lambda_final <- function(X, y, cohort, lambda_grid, repeats = 50, K = 3, seed = 20260922) {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    y <- as.integer(y)
    cohort <- as.character(cohort)
    lambda_values <- numeric(repeats)
    weights_fit <- make_balanced_weights(y = y, cohort = cohort)
    for (r in seq_len(repeats)) {
        foldid <- make_stratified_foldid(y = y, cohort = cohort, K = K, seed = seed + r)
        cvfit <- cv.glmnet(x = X, y = y, family = "binomial", alpha = 0, lambda = lambda_grid, weights = weights_fit, 
            foldid = foldid, type.measure = "deviance", standardize = TRUE, intercept = TRUE, maxit = 1e+05)
        lambda_values[r] <- cvfit$lambda.1se
    }
    lambda_final <- exp(median(log(lambda_values)))
    return(list(lambda = lambda_final, lambda_distribution = lambda_values, lambda_min = min(lambda_values), 
        lambda_max = max(lambda_values), lambda_median = median(lambda_values)))
}

calc_auc <- function(y, score) {
    y <- as.integer(y)
    n_pos <- sum(y == 1)
    n_neg <- sum(y == 0)
    if (n_pos == 0 || n_neg == 0) {
        return(NA_real_)
    }
    r <- rank(score, ties.method = "average")
    auc <- (sum(r[y == 1]) - n_pos * (n_pos + 1)/2)/(n_pos * n_neg)
    return(as.numeric(auc))
}

calc_ap <- function(y, score) {
    y <- as.integer(y)
    if (sum(y == 1) == 0) {
        return(NA_real_)
    }
    ord <- order(score, decreasing = TRUE)
    y_sorted <- y[ord]
    cumulative_tp <- cumsum(y_sorted == 1)
    precision_at_k <- cumulative_tp/seq_along(y_sorted)
    AP <- mean(precision_at_k[y_sorted == 1])
    return(as.numeric(AP))
}

calculate_metrics <- function(y, score) {
    y <- as.integer(y)
    score <- as.numeric(score)
    prevalence <- mean(y)
    AUC <- calc_auc(y, score)
    AP <- calc_ap(y, score)
    AP_norm <- if (prevalence < 1) {
        (AP - prevalence)/(1 - prevalence)
    }
    else {
        NA_real_
    }
    Brier <- mean((score - y)^2)
    Brier_null <- prevalence * (1 - prevalence)
    Brier_skill <- if (Brier_null > 0) {
        1 - Brier/Brier_null
    }
    else {
        NA_real_
    }
    data.frame(n = length(y), responders = sum(y == 1), nonresponders = sum(y == 0), prevalence = prevalence, 
        AUC = AUC, AP = AP, AP_norm = AP_norm, Brier = Brier, Brier_skill = Brier_skill)
}

run_locked_BA4_IECV <- function(data, features, lambda_grid, inner_repeats = 50, seed = 20260922) {
    X_all <- data %>% select(all_of(features)) %>% as.matrix()
    storage.mode(X_all) <- "numeric"
    y_all <- as.integer(data$response)
    cohort_all <- as.character(data$cohort)
    cohorts <- sort(unique(cohort_all))
    prediction_list <- list()
    metric_list <- list()
    lambda_list <- list()
    model_list <- list()
    for (outer_i in seq_along(cohorts)) {
        heldout <- cohorts[outer_i]
        cat("\n========================================", "FINAL BA4 OUTER HELD-OUT:", heldout, "", "========================================")
        test_idx <- cohort_all == heldout
        train_idx <- !test_idx
        X_train <- X_all[train_idx, , drop = FALSE]
        X_test <- X_all[test_idx, , drop = FALSE]
        y_train <- y_all[train_idx]
        y_test <- y_all[test_idx]
        cohort_train <- cohort_all[train_idx]
        cat("Training N:", length(y_train), "| R:", sum(y_train == 1), "| NR:", sum(y_train == 0), "")
        cat("Test N:", length(y_test), "| R:", sum(y_test == 1), "| NR:", sum(y_test == 0), "")
        lambda_result <- choose_lambda_final(X = X_train, y = y_train, cohort = cohort_train, lambda_grid = lambda_grid, 
            repeats = inner_repeats, K = 3, seed = seed + outer_i * 1000)
        best_lambda <- lambda_result$lambda
        cat("Selected lambda:", signif(best_lambda, 6), "")
        cat("Lambda range across repeats:", signif(lambda_result$lambda_min, 5), "-", signif(lambda_result$lambda_max, 
            5), "")
        weights_train <- make_balanced_weights(y = y_train, cohort = cohort_train)
        fit <- glmnet(x = X_train, y = y_train, family = "binomial", alpha = 0, lambda = best_lambda, 
            weights = weights_train, standardize = TRUE, intercept = TRUE, maxit = 1e+05)
        pred <- as.numeric(predict(fit, newx = X_test, type = "response", s = best_lambda))
        metrics <- calculate_metrics(y = y_test, score = pred)
        metrics$cohort <- heldout
        metrics$selected_lambda <- best_lambda
        metric_list[[heldout]] <- metrics
        pred_df <- data[test_idx, c("patient_uid", "cohort", "sample_id", "response")]
        pred_df$OOCO_ICBScore <- pred
        pred_df$heldout_cohort <- heldout
        pred_df$selected_lambda <- best_lambda
        prediction_list[[heldout]] <- pred_df
        lambda_list[[heldout]] <- lambda_result
        model_list[[heldout]] <- fit
    }
    return(list(predictions = bind_rows(prediction_list), metrics = bind_rows(metric_list), lambda_results = lambda_list, 
        models = model_list))
}

set.seed(20260922)

BA4_IECV <- run_locked_BA4_IECV(data = model_data_BA4, features = LOCKED_FEATURES, lambda_grid = lambda_grid_final, 
    inner_repeats = 50, seed = 20260922)

BA4_outer_metrics <- BA4_IECV$metrics %>% select(cohort, n, responders, nonresponders, prevalence, AUC, 
    AP, AP_norm, Brier, Brier_skill, selected_lambda)

BA4_outer_metrics

BA4_macro_metrics <- BA4_outer_metrics %>% summarise(Macro_AUC = mean(AUC), SD_AUC = sd(AUC), Min_AUC = min(AUC), 
    Macro_AP = mean(AP), Macro_AP_norm = mean(AP_norm), Macro_Brier = mean(Brier), Macro_Brier_skill = mean(Brier_skill))

BA4_macro_metrics

BA4_OOCO_predictions <- BA4_IECV$predictions %>% arrange(cohort, sample_id)

nrow(BA4_OOCO_predictions)

anyDuplicated(BA4_OOCO_predictions$patient_uid)

setdiff(model_data_BA4$patient_uid, BA4_OOCO_predictions$patient_uid)

summary(BA4_OOCO_predictions$OOCO_ICBScore)

stopifnot(all(BA4_OOCO_predictions$OOCO_ICBScore >= 0 & BA4_OOCO_predictions$OOCO_ICBScore <= 1))

BA4_pooled_metrics <- calculate_metrics(y = BA4_OOCO_predictions$response, score = BA4_OOCO_predictions$OOCO_ICBScore)

BA4_pooled_metrics

BA4_lambda_summary <- bind_rows(lapply(names(BA4_IECV$lambda_results), function(cc) {
    x <- BA4_IECV$lambda_results[[cc]]$lambda_distribution
    data.frame(cohort_heldout = cc, lambda_median = median(x), lambda_mean = mean(x), lambda_sd = sd(x), 
        lambda_min = min(x), lambda_max = max(x), lambda_Q1 = quantile(x, 0.25), lambda_Q3 = quantile(x, 
            0.75))
}))

BA4_lambda_summary

BA4_lambda_boundary <- bind_rows(lapply(names(BA4_IECV$lambda_results), function(cc) {
    x <- BA4_IECV$lambda_results[[cc]]$lambda_distribution
    data.frame(cohort_heldout = cc, pct_at_min_grid = mean(x == min(lambda_grid_final)), pct_at_max_grid = mean(x == 
        max(lambda_grid_final)))
}))

BA4_lambda_boundary

saveRDS(BA4_IECV, file = "STEP1_BA4_locked_IECV_full_results.rds")

write.csv(BA4_OOCO_predictions, file = "STEP1_BA4_OOCO_predictions_69patients.csv", row.names = FALSE)

write.csv(BA4_outer_metrics, file = "STEP1_BA4_outer_metrics.csv", row.names = FALSE)

write.csv(BA4_macro_metrics, file = "STEP1_BA4_macro_metrics.csv", row.names = FALSE)

write.csv(BA4_lambda_summary, file = "STEP1_BA4_lambda_stability.csv", row.names = FALSE)

saveRDS(LOCKED_FEATURES, file = "STEP1_BA4_locked_features.rds")

BA4_outer_metrics

BA4_macro_metrics

BA4_pooled_metrics

BA4_lambda_summary

BA4_lambda_boundary

summary(BA4_OOCO_predictions$OOCO_ICBScore)


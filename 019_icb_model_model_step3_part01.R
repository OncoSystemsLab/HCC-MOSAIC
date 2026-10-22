if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

observed_Macro_AUC <- BA4_macro_metrics$Macro_AUC[1]

observed_Macro_APnorm <- BA4_macro_metrics$Macro_AP_norm[1]

observed_Min_AUC <- BA4_macro_metrics$Min_AUC[1]

observed_Macro_Brier <- BA4_macro_metrics$Macro_Brier[1]

observed_Pooled_AUC <- BA4_pooled_metrics$AUC[1]

data.frame(Metric = c("Macro_AUC", "Macro_APnorm", "Min_AUC", "Macro_Brier", "Pooled_AUC"), Observed = c(observed_Macro_AUC, 
    observed_Macro_APnorm, observed_Min_AUC, observed_Macro_Brier, observed_Pooled_AUC))

permute_response_within_cohort <- function(data, seed) {
    set.seed(seed)
    permuted_data <- data %>% group_by(cohort) %>% mutate(response = sample(response, size = n(), replace = FALSE)) %>% 
        ungroup()
    return(permuted_data)
}

original_count <- model_data_BA4 %>% count(cohort, response) %>% arrange(cohort, response) %>% as.data.frame()

test_permuted_data <- permute_response_within_cohort(data = model_data_BA4, seed = 20260923)

permuted_count <- test_permuted_data %>% count(cohort, response) %>% arrange(cohort, response) %>% as.data.frame()

stopifnot(all(original_count$cohort == permuted_count$cohort), all(original_count$response == permuted_count$response), 
    all(original_count$n == permuted_count$n))

n_labels_changed <- sum(model_data_BA4$response != test_permuted_data$response)

stopifnot(n_labels_changed > 0)

run_one_BA4_permutation <- function(iteration, data, features, lambda_grid, inner_repeats = 50, seed_base = 20260923) {
    perm_data <- permute_response_within_cohort(data = data, seed = seed_base + iteration)
    perm_fit <- NULL
    invisible(capture.output(perm_fit <- suppressWarnings(run_locked_BA4_IECV(data = perm_data, features = features, 
        lambda_grid = lambda_grid, inner_repeats = inner_repeats, seed = seed_base + iteration * 1000))))
    m <- perm_fit$metrics
    pred <- perm_fit$predictions
    pooled_m <- calculate_metrics(y = pred$response, score = pred$OOCO_ICBScore)
    data.frame(iteration = iteration, Macro_AUC = mean(m$AUC, na.rm = TRUE), Macro_AP = mean(m$AP, na.rm = TRUE), 
        Macro_APnorm = mean(m$AP_norm, na.rm = TRUE), Min_AUC = min(m$AUC, na.rm = TRUE), Macro_Brier = mean(m$Brier, 
            na.rm = TRUE), Pooled_AUC = pooled_m$AUC, AUC_GSE181946 = m$AUC[m$cohort == "GSE181946"], 
        AUC_GSE202069 = m$AUC[m$cohort == "GSE202069"], AUC_PRJEB34724 = m$AUC[m$cohort == "PRJEB34724"], 
        status = "OK", error_message = "NA")
}

perm_test_1 <- run_one_BA4_permutation(iteration = 1, data = model_data_BA4, features = LOCKED_FEATURES, 
    lambda_grid = lambda_grid_final, inner_repeats = 50, seed_base = 20260923)

perm_test_1

perm_smoke <- vector("list", 10)

for (b in 1:10) {
    cat("Smoke permutation:", b, "/10")
    perm_smoke[[b]] <- run_one_BA4_permutation(iteration = b, data = model_data_BA4, features = LOCKED_FEATURES, 
        lambda_grid = lambda_grid_final, inner_repeats = 50, seed_base = 20260923)
}

perm_smoke_df <- bind_rows(perm_smoke)

perm_smoke_df

summary(perm_smoke_df$Macro_AUC)

make_perm_error_row <- function(iteration, message) {
    data.frame(iteration = iteration, Macro_AUC = NA_real_, Macro_AP = NA_real_, Macro_APnorm = NA_real_, 
        Min_AUC = NA_real_, Macro_Brier = NA_real_, Pooled_AUC = NA_real_, AUC_GSE181946 = NA_real_, 
        AUC_GSE202069 = NA_real_, AUC_PRJEB34724 = NA_real_, status = "ERROR", error_message = as.character(message))
}

B_perm <- 1000

perm_checkpoint_file <- "STEP3_BA4_permutation_checkpoint.rds"

if (file.exists(perm_checkpoint_file)) {
    cat("Existing checkpoint found.")
    permutation_results <- readRDS(perm_checkpoint_file)
    if (length(permutation_results) < B_perm) {
        length(permutation_results) <- B_perm
    }
} else {
    permutation_results <- vector("list", B_perm)
}

completed_perm <- which(!vapply(permutation_results, is.null, logical(1)))

length(completed_perm)

for (b in setdiff(seq_len(B_perm), completed_perm)) {
    cat("Permutation:", b, "/", B_perm, "")
    permutation_results[[b]] <- tryCatch({
        run_one_BA4_permutation(iteration = b, data = model_data_BA4, features = LOCKED_FEATURES, lambda_grid = lambda_grid_final, 
            inner_repeats = 50, seed_base = 20260923)
    }, error = function(e) {
        make_perm_error_row(iteration = b, message = conditionMessage(e))
    })
    if (b%%10 == 0 || b == B_perm) {
        saveRDS(permutation_results, perm_checkpoint_file)
        cat("Checkpoint saved:", b, "/", B_perm, "")
    }
}

BA4_permutation_results <- bind_rows(permutation_results)

table(BA4_permutation_results$status)

BA4_perm_ok <- BA4_permutation_results %>% filter(status == "OK")

nrow(BA4_perm_ok)

summary(BA4_perm_ok$Macro_AUC)

quantile(BA4_perm_ok$Macro_AUC, probs = c(0.025, 0.5, 0.95, 0.975, 0.99))

summary(BA4_perm_ok$Macro_APnorm)

summary(BA4_perm_ok$Macro_Brier)

perm_p_Macro_AUC <- (1 + sum(BA4_perm_ok$Macro_AUC >= observed_Macro_AUC))/(nrow(BA4_perm_ok) + 1)

perm_p_Macro_APnorm <- (1 + sum(BA4_perm_ok$Macro_APnorm >= observed_Macro_APnorm))/(nrow(BA4_perm_ok) + 
    1)

perm_p_Min_AUC <- (1 + sum(BA4_perm_ok$Min_AUC >= observed_Min_AUC))/(nrow(BA4_perm_ok) + 1)

perm_p_Macro_Brier <- (1 + sum(BA4_perm_ok$Macro_Brier <= observed_Macro_Brier))/(nrow(BA4_perm_ok) + 
    1)

perm_p_Pooled_AUC <- (1 + sum(BA4_perm_ok$Pooled_AUC >= observed_Pooled_AUC))/(nrow(BA4_perm_ok) + 1)

BA4_permutation_summary <- data.frame(Metric = c("Macro_AUC", "Macro_APnorm", "Min_AUC", "Macro_Brier", 
    "Pooled_AUC"), Observed = c(observed_Macro_AUC, observed_Macro_APnorm, observed_Min_AUC, observed_Macro_Brier, 
    observed_Pooled_AUC), Null_mean = c(mean(BA4_perm_ok$Macro_AUC), mean(BA4_perm_ok$Macro_APnorm), 
    mean(BA4_perm_ok$Min_AUC), mean(BA4_perm_ok$Macro_Brier), mean(BA4_perm_ok$Pooled_AUC)), Null_SD = c(sd(BA4_perm_ok$Macro_AUC), 
    sd(BA4_perm_ok$Macro_APnorm), sd(BA4_perm_ok$Min_AUC), sd(BA4_perm_ok$Macro_Brier), sd(BA4_perm_ok$Pooled_AUC)), 
    Permutation_P = c(perm_p_Macro_AUC, perm_p_Macro_APnorm, perm_p_Min_AUC, perm_p_Macro_Brier, perm_p_Pooled_AUC))

BA4_permutation_summary

BA4_perm_auc_quantiles <- quantile(BA4_perm_ok$Macro_AUC, probs = c(0.025, 0.5, 0.95, 0.975, 0.99))

BA4_perm_auc_quantiles

p_perm_auc <- ggplot(BA4_perm_ok, aes(x = Macro_AUC)) + geom_histogram(bins = 40) + geom_vline(xintercept = observed_Macro_AUC, 
    linetype = 2, linewidth = 1) + labs(x = "Permuted macro ROC-AUC", y = "Number of permutations", title = "Permutation test of the locked BA4 model") + 
    theme_classic(base_size = 13)

p_perm_auc

ggsave("STEP3_BA4_permutation_MacroAUC.pdf", p_perm_auc, width = 7, height = 5)

p_perm_ap <- ggplot(BA4_perm_ok, aes(x = Macro_APnorm)) + geom_histogram(bins = 40) + geom_vline(xintercept = observed_Macro_APnorm, 
    linetype = 2, linewidth = 1) + labs(x = "Permuted macro normalized AP", y = "Number of permutations", 
    title = "Permutation test of BA4 precision-recall performance") + theme_classic(base_size = 13)

p_perm_ap

write.csv(BA4_permutation_results, "STEP3_BA4_permutation_1000_raw.csv", row.names = FALSE)

write.csv(BA4_permutation_summary, "STEP3_BA4_permutation_summary.csv", row.names = FALSE)

saveRDS(BA4_permutation_results, "STEP3_BA4_permutation_1000_results.rds")

table(BA4_permutation_results$status)

BA4_permutation_summary

BA4_perm_auc_quantiles

summary(BA4_perm_ok$Macro_AUC)


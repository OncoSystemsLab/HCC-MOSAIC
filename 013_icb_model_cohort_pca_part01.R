if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

X_qc <- score_all %>% select(all_of(feature_order)) %>% as.matrix()

X_qc_scaled <- scale(X_qc)

pca_result <- prcomp(X_qc_scaled, center = FALSE, scale. = FALSE)

pca_df <- data.frame(patient_uid = score_all$patient_uid, cohort = score_all$cohort, sample_id = score_all$sample_id, 
    PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2])

pve <- (pca_result$sdev^2/sum(pca_result$sdev^2)) * 100

ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort)) + geom_point(size = 3, alpha = 0.8) + labs(x = paste0("PC1 (", 
    round(pve[1], 1), "%)"), y = paste0("PC2 (", round(pve[2], 1), "%)"), title = "PCA of Hallmark + ICI ssGSEA scores") + 
    theme_classic(base_size = 13)


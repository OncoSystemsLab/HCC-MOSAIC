if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

library(tidyr)

library(ggplot2)

library(patchwork)

library(pROC)

library(scales)

library(forcats)

library(stringr)

library(cowplot)

library(ggrepel)

library(ComplexHeatmap)

library(circlize)

library(grid)

library(gridExtra)

library(RColorBrewer)

library(DiagrammeR)

theme_nature <- function(base_size = 12) {
    theme_classic(base_size = base_size) + theme(axis.title = element_text(color = "black", face = "bold"), 
        axis.text = element_text(color = "black"), plot.title = element_text(face = "bold", hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), legend.title = element_text(face = "bold"), legend.text = element_text(color = "black"), 
        strip.background = element_rect(fill = "white", color = "black"), strip.text = element_text(face = "bold"), 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5))
}

cohort_cols <- c(GSE181946 = "#D95F02", GSE202069 = "#1B9E77", PRJEB34724 = "#7570B3")

response_cols <- c(`Non-responder` = "#4D4D4D", Responder = "#D7301F")

feature_cols <- c(ICI__HCC_INFLAMMATORY4 = "#B2182B", HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "#EF8A62", 
    HALLMARK_BILE_ACID_METABOLISM = "#2166AC", HALLMARK_FATTY_ACID_METABOLISM = "#4393C3")

pretty_feature <- c(ICI__HCC_INFLAMMATORY4 = "HCC inflammatory 4", HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "Reactive oxygen species", 
    HALLMARK_BILE_ACID_METABOLISM = "Bile acid metabolism", HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid metabolism")

get_auc_ci <- function(df, score_col = "OOCO_ICBScore", y_col = "response", n_boot = 2000) {
    roc_obj <- pROC::roc(df[[y_col]], df[[score_col]], quiet = TRUE, direction = "<")
    ci_obj <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot)
    data.frame(AUC = as.numeric(pROC::auc(roc_obj)), CI_low = as.numeric(ci_obj[1]), CI_high = as.numeric(ci_obj[3]))
}

format_p <- function(p) {
    ifelse(p < 0.001, "P < 0.001", paste0("P = ", sprintf("%.3f", p)))
}

graph_text <- "\ndigraph BA4_workflow {\n\ngraph [layout = dot, rankdir = TB, fontsize = 20]\nnode [shape = box, style = 'rounded,filled', color = '#2C3E50',\n      fillcolor = '#F8F9FA', fontname = Helvetica, fontsize = 18]\nedge [color = '#2C3E50', penwidth = 1.5]\n\nA [label = 'Three pretreatment HCC anti-PD-1 cohorts\\nPRJEB34724 (n=35), GSE181946 (n=17), GSE202069 (n=17)\\nTotal n=69']\nB [label = 'RNA-seq preprocessing\\nTPM conversion -> log2(TPM+1)\\nCommon background genes: 21,308']\nC [label = 'ssGSEA feature scoring\\n50 Hallmark + 7 ICI signatures\\nTotal 57 candidate features']\nD [label = 'Feature discovery\\n127 ICI combinations evaluated\\nHallmark incremental addition across anchor models']\nE [label = 'Locked BA4 feature set\\n1) HCC inflammatory4\\n2) ROS pathway\\n3) Bile acid metabolism\\n4) Fatty acid metabolism']\nF [label = 'Model fitting\\nRidge logistic regression\\nalpha = 0, family = binomial\\nBalanced cohort-by-response weights']\nG [label = 'Nested internal-external cross-validation\\nLeave-one-cohort-out IECV\\nInner repeated CV for lambda selection']\nH [label = 'Robustness evaluation\\nWithin-cohort permutation (n=1000)\\nStratified bootstrap (n=2000)\\nTraining-set bootstrap stability (n=1000/fold)']\nI [label = 'Final locked BA4 model\\nlambda = 1.060026\\nBA4 score = responder-like ICB score']\n\nA -> B -> C -> D -> E -> F -> G -> H -> I\n}"

grViz(graph_text)

rsvg::rsvg_pdf(charToRaw(DiagrammeRsvg::export_svg(grViz(graph_text))), "Figure1_BA4_workflow.pdf")

rsvg::rsvg_png(charToRaw(DiagrammeRsvg::export_svg(grViz(graph_text))), "Figure1_BA4_workflow.png", width = 3200, 
    height = 4200)

roc_plot_df <- bind_rows(lapply(split(BA4_OOCO_predictions, BA4_OOCO_predictions$cohort), function(df) {
    roc_obj <- roc(df$response, df$OOCO_ICBScore, quiet = TRUE, direction = "<")
    coords_df <- coords(roc_obj, x = "all", ret = c("specificity", "sensitivity"))
    data.frame(cohort = unique(df$cohort), FPR = 1 - coords_df$specificity, TPR = coords_df$sensitivity, 
        AUC = as.numeric(auc(roc_obj)))
}))

auc_label_df <- BA4_outer_metrics %>% transmute(cohort, label = paste0(cohort, "  AUC=", sprintf("%.3f", 
    AUC)))

p_roc <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = cohort)) + geom_path(linewidth = 1.1) + geom_abline(slope = 1, 
    intercept = 0, linetype = 2, color = "black") + scale_color_manual(values = cohort_cols, labels = auc_label_df$label) + 
    coord_equal() + labs(x = "False positive rate", y = "True positive rate", color = "Cohort", title = "Locked BA4 model: leave-one-cohort-out ROC") + 
    theme_nature(base_size = 13)

p_roc

pr_curve_df <- bind_rows(lapply(split(BA4_OOCO_predictions, BA4_OOCO_predictions$cohort), function(df) {
    y <- df$response
    score <- df$OOCO_ICBScore
    ord <- order(score, decreasing = TRUE)
    y_ord <- y[ord]
    tp <- cumsum(y_ord == 1)
    fp <- cumsum(y_ord == 0)
    fn <- sum(y_ord == 1) - tp
    precision <- tp/(tp + fp)
    recall <- tp/(tp + fn)
    data.frame(cohort = unique(df$cohort), Recall = c(0, recall), Precision = c(1, precision), prevalence = mean(y))
}))

prev_df <- BA4_OOCO_predictions %>% group_by(cohort) %>% summarise(prevalence = mean(response), .groups = "drop")

p_pr <- ggplot(pr_curve_df, aes(x = Recall, y = Precision, color = cohort)) + geom_path(linewidth = 1.1) + 
    geom_hline(data = prev_df, aes(yintercept = prevalence, color = cohort), linetype = 2, alpha = 0.7, 
        show.legend = FALSE) + scale_color_manual(values = cohort_cols) + facet_wrap(~cohort, nrow = 1) + 
    labs(x = "Recall", y = "Precision", title = "Locked BA4 model: out-of-cohort precision-recall curves") + 
    theme_nature(base_size = 13)

p_pr

auc_ci_df <- bind_rows(lapply(split(BA4_OOCO_predictions, BA4_OOCO_predictions$cohort), function(df) {
    out <- get_auc_ci(df, n_boot = 2000)
    out$cohort <- unique(df$cohort)
    out
})) %>% left_join(BA4_outer_metrics %>% select(cohort, n, responders, nonresponders), by = "cohort") %>% 
    mutate(cohort = factor(cohort, levels = rev(c("GSE181946", "GSE202069", "PRJEB34724"))))

p_auc_forest <- ggplot(auc_ci_df, aes(x = AUC, y = cohort, color = cohort)) + geom_vline(xintercept = 0.5, 
    linetype = 2) + geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.15, linewidth = 1) + 
    geom_point(size = 3) + scale_color_manual(values = cohort_cols, guide = "none") + scale_x_continuous(limits = c(0, 
    1.02)) + labs(x = "ROC-AUC (bootstrap 95% CI)", y = NULL, title = "Out-of-cohort discrimination of the locked BA4 model") + 
    theme_nature(base_size = 13)

p_auc_forest

plot_df <- BA4_OOCO_predictions %>% mutate(response_label = ifelse(response == 1, "Responder", "Non-responder"))

wilcox_df <- plot_df %>% group_by(cohort) %>% summarise(p = wilcox.test(OOCO_ICBScore ~ response_label)$p.value, 
    y = max(OOCO_ICBScore) * 1.05, label = format_p(p), .groups = "drop")

p_score_box <- ggplot(plot_df, aes(x = response_label, y = OOCO_ICBScore, fill = response_label)) + geom_boxplot(width = 0.65, 
    outlier.shape = NA, alpha = 0.35, color = "black") + geom_jitter(width = 0.12, size = 2.2, alpha = 0.85) + 
    geom_text(data = wilcox_df, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 4) + 
    scale_fill_manual(values = response_cols) + facet_wrap(~cohort, nrow = 1, scales = "free_y") + labs(x = NULL, 
    y = "OOCO BA4 ICB Score", title = "Out-of-cohort prediction of anti-PD-1 response") + theme_nature(base_size = 13) + 
    theme(legend.position = "none")

p_score_box

fig2 <- (p_roc | p_auc_forest)/(p_pr | p_score_box) + plot_annotation(tag_levels = "A")

fig2

ggsave("Figure2_BA4_overall_performance.tiff", fig2, width = 14, height = 10, dpi = 300, compression = "lzw")

waterfall_df <- BA4_OOCO_predictions %>% mutate(response_label = ifelse(response == 1, "Responder", "Non-responder")) %>% 
    arrange(OOCO_ICBScore) %>% mutate(patient_order = row_number())

p_waterfall <- ggplot(waterfall_df, aes(x = patient_order, y = OOCO_ICBScore, fill = response_label)) + 
    geom_col(width = 0.9, color = "white", linewidth = 0.15) + scale_fill_manual(values = response_cols) + 
    labs(x = "Patients (ordered by BA4 ICB Score)", y = "OOCO BA4 ICB Score", fill = NULL, title = "Patient-level distribution of the BA4 ICB Score") + 
    theme_nature(base_size = 13) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    legend.position = "top")

p_waterfall

density_df <- BA4_OOCO_predictions %>% mutate(response_label = ifelse(response == 1, "Responder", "Non-responder"))

p_density <- ggplot(density_df, aes(x = OOCO_ICBScore, fill = response_label, color = response_label)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 18, alpha = 0.25, position = "identity") + geom_density(linewidth = 1.1, 
    alpha = 0.15) + scale_fill_manual(values = response_cols) + scale_color_manual(values = response_cols) + 
    labs(x = "OOCO BA4 ICB Score", y = "Density", title = "Distribution of BA4 ICB Score by response status") + 
    theme_nature(base_size = 13)

p_density

heat_df <- model_data_BA4 %>% select(patient_uid, cohort, sample_id, response, all_of(LOCKED_FEATURES)) %>% 
    left_join(BA4_OOCO_predictions %>% select(patient_uid, OOCO_ICBScore), by = "patient_uid") %>% arrange(OOCO_ICBScore)

mat_heat <- heat_df %>% select(all_of(LOCKED_FEATURES), OOCO_ICBScore) %>% as.matrix()

rownames(mat_heat) <- heat_df$patient_uid

mat_heat_z <- scale(mat_heat)

mat_heat_z <- t(mat_heat_z)

ha_top <- HeatmapAnnotation(Cohort = heat_df$cohort, Response = ifelse(heat_df$response == 1, "Responder", 
    "Non-responder"), col = list(Cohort = cohort_cols, Response = response_cols), annotation_name_side = "left")

ht <- Heatmap(mat_heat_z, name = "Z-score", top_annotation = ha_top, col = colorRamp2(c(-2, 0, 2), c("#2166AC", 
    "white", "#B2182B")), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = FALSE, row_labels = c(pretty_feature[LOCKED_FEATURES], 
    "BA4 ICB Score"), column_title = "69 HCC patients ordered by BA4 score", row_names_gp = gpar(fontsize = 10, 
    fontface = "bold"))

draw(ht)

pdf("Figure3C_patient_feature_heatmap.pdf", width = 14, height = 5)

draw(ht)

dev.off()

fig3_ab <- p_waterfall/p_density + plot_annotation(tag_levels = "A")

fig3_ab

ggsave("Figure3_AB_BA4_score_landscape.tiff", fig3_ab, width = 12, height = 9, dpi = 300, compression = "lzw")

coef_df <- BA4_standardized_effects %>% mutate(feature_label = pretty_feature[feature], feature_label = factor(feature_label, 
    levels = feature_label[order(standardized_coefficient)]))

p_coef <- ggplot(coef_df, aes(x = standardized_coefficient, y = feature_label, fill = direction)) + geom_col(width = 0.7, 
    color = "black") + geom_vline(xintercept = 0, linetype = 2) + scale_fill_manual(values = c(`Responder-high` = "#D7301F", 
    `Non-responder-high` = "#2166AC")) + labs(x = "Standardized coefficient", y = NULL, fill = NULL, 
    title = "Feature contributions in the final locked BA4 model") + theme_nature(base_size = 13)

p_coef

feature_long <- model_data_BA4 %>% select(patient_uid, response, all_of(LOCKED_FEATURES)) %>% mutate(response_label = ifelse(response == 
    1, "Responder", "Non-responder")) %>% pivot_longer(cols = all_of(LOCKED_FEATURES), names_to = "feature", 
    values_to = "value") %>% mutate(feature_label = pretty_feature[feature])

p_feature_box <- ggplot(feature_long, aes(x = response_label, y = value, fill = response_label)) + geom_boxplot(outlier.shape = NA, 
    alpha = 0.35, width = 0.65) + geom_jitter(width = 0.12, size = 1.6, alpha = 0.75) + scale_fill_manual(values = response_cols) + 
    facet_wrap(~feature_label, scales = "free_y", ncol = 2) + labs(x = NULL, y = "ssGSEA score", title = "Distribution of BA4 model features by response status") + 
    theme_nature(base_size = 12) + theme(legend.position = "none")

p_feature_box

feature_score_df <- model_data_BA4 %>% select(patient_uid, all_of(LOCKED_FEATURES)) %>% left_join(BA4_OOCO_predictions %>% 
    select(patient_uid, OOCO_ICBScore), by = "patient_uid") %>% pivot_longer(cols = all_of(LOCKED_FEATURES), 
    names_to = "feature", values_to = "value") %>% mutate(feature_label = pretty_feature[feature])

cor_stat_df <- feature_score_df %>% group_by(feature_label) %>% summarise(rho = cor(value, OOCO_ICBScore, 
    method = "spearman"), p = cor.test(value, OOCO_ICBScore, method = "spearman")$p.value, x = min(value), 
    y = max(OOCO_ICBScore), label = paste0("rho = ", sprintf("%.2f", rho), "", format_p(p)), .groups = "drop")

p_feature_corr <- ggplot(feature_score_df, aes(x = value, y = OOCO_ICBScore)) + geom_point(alpha = 0.75, 
    size = 2, color = "#2C3E50") + geom_smooth(method = "lm", se = TRUE, color = "#B2182B", linewidth = 0.8) + 
    geom_text(data = cor_stat_df, aes(x = x, y = y, label = label), hjust = 0, vjust = 1, size = 3.6, 
        inherit.aes = FALSE) + facet_wrap(~feature_label, scales = "free_x", ncol = 2) + labs(x = "Feature ssGSEA score", 
    y = "OOCO BA4 ICB Score", title = "Association between BA4 model features and BA4 score") + theme_nature(base_size = 12)

p_feature_corr

cor_mat <- cor(model_data_BA4 %>% select(all_of(LOCKED_FEATURES)) %>% as.matrix(), method = "spearman")

cor_df <- as.data.frame(as.table(cor_mat)) %>% rename(Feature1 = Var1, Feature2 = Var2, Correlation = Freq) %>% 
    mutate(Feature1 = pretty_feature[Feature1], Feature2 = pretty_feature[Feature2])

p_cor_heat <- ggplot(cor_df, aes(x = Feature1, y = Feature2, fill = Correlation)) + geom_tile(color = "white") + 
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 4) + scale_fill_gradient2(low = "#B2182B", 
    mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1)) + coord_fixed() + labs(x = NULL, 
    y = NULL, fill = "Spearman\nrho", title = "Correlation structure of locked BA4 model features") + 
    theme_minimal(base_size = 12) + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(face = "bold", 
    hjust = 0.5))

p_cor_heat

fig4 <- (p_coef | p_cor_heat)/(p_feature_box | p_feature_corr) + plot_annotation(tag_levels = "A")

fig4

ggsave("Figure4_BA4_feature_contribution.tiff", fig4, width = 15, height = 12, dpi = 300, compression = "lzw")

perm_df <- BA4_perm_ok

p_perm <- ggplot(perm_df, aes(x = Macro_AUC)) + geom_histogram(bins = 35, fill = "grey60", color = "white") + 
    geom_vline(xintercept = observed_Macro_AUC, color = "black", linewidth = 1.2) + geom_vline(xintercept = BA4_perm_auc_quantiles["97.5%"], 
    linetype = 3, linewidth = 0.9) + labs(x = "Permutation macro ROC-AUC", y = "Permutation runs", title = "Permutation test of the locked BA4 model", 
    subtitle = paste0("Observed Macro AUC = ", sprintf("%.3f", observed_Macro_AUC), "; empirical P = ", 
        sprintf("%.4f", BA4_permutation_summary$Permutation_P[1]))) + theme_nature(base_size = 13)

p_perm

p_boot_macro <- ggplot(BA4_bootstrap_results, aes(x = Macro_AUC)) + geom_histogram(bins = 35, fill = "grey50", 
    color = "white") + geom_vline(xintercept = observed_Macro_AUC, linetype = 2, linewidth = 1.2) + labs(x = "Bootstrap macro ROC-AUC", 
    y = "Bootstrap samples", title = "Sampling uncertainty of the locked BA4 model") + theme_nature(base_size = 13)

p_boot_macro

p_boot_min <- ggplot(BA4_bootstrap_results, aes(x = Min_AUC)) + geom_histogram(bins = 35, fill = "grey50", 
    color = "white") + geom_vline(xintercept = observed_Min_AUC, linetype = 2, linewidth = 1.2) + geom_vline(xintercept = 0.5, 
    linetype = 3, linewidth = 0.9) + labs(x = "Bootstrap minimum cohort ROC-AUC", y = "Bootstrap samples", 
    title = "Worst-cohort robustness of the locked BA4 model") + theme_nature(base_size = 13)

p_boot_min

p_stability <- ggplot(BA4_training_stability_summary, aes(x = Spearman_r)) + geom_histogram(bins = 40, 
    fill = "grey45", color = "white") + geom_vline(xintercept = 0.8, linetype = 2) + facet_wrap(~heldout_cohort, 
    nrow = 1) + labs(x = "Spearman correlation with original OOCO prediction", y = "Bootstrap training models", 
    title = "Prediction stability under training-set resampling") + theme_nature(base_size = 13)

p_stability

fig5 <- (p_perm | p_boot_macro)/(p_boot_min | p_stability) + plot_annotation(tag_levels = "A")

fig5

ggsave("Figure5_BA4_robustness.tiff", fig5, width = 14, height = 10, dpi = 300, compression = "lzw")

pca_df$cohort <- factor(pca_df$cohort, levels = c("GSE181946", "GSE202069", "PRJEB34724"))

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = cohort)) + geom_point(size = 3, alpha = 0.85) + 
    stat_ellipse(level = 0.68, linewidth = 0.8, alpha = 0.5) + scale_color_manual(values = cohort_cols) + 
    labs(x = paste0("PC1 (", round(pve[1], 1), "%)"), y = paste0("PC2 (", round(pve[2], 1), "%)"), title = "PCA of Hallmark and ICI ssGSEA features") + 
    theme_nature(base_size = 13)

p_pca

p_compare <- ggplot(score_compare, aes(x = OOCO_ICBScore, y = BA4_FinalScore, color = cohort)) + geom_point(size = 2.6, 
    alpha = 0.9) + geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.8) + scale_color_manual(values = cohort_cols) + 
    labs(x = "OOCO BA4 ICB Score", y = "Final full-fit BA4 score", title = "Consistency between OOCO and final full-fit BA4 scores", 
        subtitle = paste0("Spearman rho = ", sprintf("%.3f", cor(score_compare$BA4_FinalScore, score_compare$OOCO_ICBScore, 
            method = "spearman")))) + theme_nature(base_size = 13)

p_compare

coef_stab_df <- BA4_coefficient_stability %>% mutate(feature_label = pretty_feature[feature])

BA4_boot_coef_long <- BA4_boot_coef_long %>% mutate(coefficient_x1000 = coefficient * 1000, feature_label = factor(feature_label, 
    levels = c("HCC inflammatory 4", "Reactive oxygen species", "Bile acid metabolism", "Fatty acid metabolism")))

p_coef_stab <- ggplot(BA4_boot_coef_long, aes(x = coefficient_x1000, y = feature_label)) + geom_vline(xintercept = 0, 
    linetype = 2, linewidth = 0.6, color = "grey40") + geom_boxplot(aes(fill = feature), width = 0.65, 
    outlier.shape = NA, linewidth = 0.5, alpha = 0.8) + scale_fill_manual(values = feature_cols, guide = "none") + 
    facet_wrap(~heldout_cohort, scales = "free_x", nrow = 1) + labs(x = expression(paste("Ridge coefficient (", 
    beta, " × 10"^3, ")")), y = NULL, title = "Coefficient stability under training-set resampling") + 
    theme_nature(base_size = 12) + theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", 
    hjust = 0.5))

p_coef_stab

ggsave("Supplementary_BA4_coefficient_stability.tiff", p_coef_stab, width = 12, height = 5, dpi = 600, 
    compression = "lzw")

ggsave("Supplementary_BA4_coefficient_stability.pdf", p_coef_stab, width = 12, height = 5)

patient_uncertainty_df <- BA4_patient_prediction_stability %>% left_join(BA4_OOCO_predictions %>% select(patient_uid, 
    cohort, OOCO_ICBScore), by = "patient_uid") %>% mutate(response_label = ifelse(response == 1, "Responder", 
    "Non-responder"))

p_patient_sd <- ggplot(patient_uncertainty_df, aes(x = OOCO_ICBScore, y = bootstrap_SD, color = cohort, 
    shape = response_label)) + geom_point(size = 2.7, alpha = 0.85) + scale_color_manual(values = cohort_cols) + 
    scale_shape_manual(values = c(`Non-responder` = 16, Responder = 17)) + labs(x = "OOCO BA4 ICB Score", 
    y = "Bootstrap SD of predicted score", color = "Cohort", shape = "Response", title = "Patient-level prediction uncertainty") + 
    theme_nature(base_size = 13)

p_patient_sd

patient_uncertainty_plot_df <- patient_uncertainty_df %>% arrange(OOCO_ICBScore) %>% mutate(patient_order = row_number())

p_patient_ci <- ggplot(patient_uncertainty_plot_df, aes(x = patient_order, y = bootstrap_median)) + geom_linerange(aes(ymin = CI_low, 
    ymax = CI_high, color = response_label), linewidth = 0.7, alpha = 0.75) + geom_point(aes(color = response_label), 
    size = 2) + geom_point(aes(y = OOCO_ICBScore), shape = 1, size = 2.4, stroke = 0.7, color = "black") + 
    scale_color_manual(values = response_cols) + labs(x = "Patients ordered by OOCO BA4 ICB Score", y = "BA4 ICB Score", 
    color = NULL, title = "Patient-level prediction stability", subtitle = "Points and bars show bootstrap median and 95% interval; open circles indicate reference OOCO score") + 
    theme_nature(base_size = 12) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
    legend.position = "top")

p_patient_ci

ggsave("Supplementary_patient_prediction_uncertainty.tiff", p_patient_ci, width = 13, height = 6, dpi = 600, 
    compression = "lzw")

p_ici_scan <- ggplot(ICI_subset_summary, aes(x = n_features, y = Macro_AP_norm)) + geom_jitter(width = 0.15, 
    alpha = 0.5, color = "#2C3E50") + stat_summary(fun = max, geom = "line", aes(group = 1), color = "#D7301F", 
    linewidth = 1) + stat_summary(fun = max, geom = "point", size = 3, color = "#D7301F") + labs(x = "Number of ICI signatures", 
    y = "Macro normalized AP", title = "Performance scan across all 127 ICI signature combinations") + 
    theme_nature(base_size = 13)

p_ici_scan

hallmark_plot_df <- hallmark_increment_summary %>% mutate(highlight = hallmark == "HALLMARK_BILE_ACID_METABOLISM")

p_hallmark_inc <- ggplot(hallmark_plot_df, aes(x = mean_Delta_AUC, y = mean_Delta_APnorm)) + geom_point(aes(color = highlight), 
    size = 3, alpha = 0.85) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + 
    geom_text_repel(data = hallmark_plot_df %>% filter(highlight | mean_Delta_APnorm > 0.05), aes(label = hallmark), 
        size = 3.2, max.overlaps = 20) + scale_color_manual(values = c(`FALSE` = "grey60", `TRUE` = "#D7301F"), 
    guide = "none") + labs(x = "Mean Delta AUC", y = "Mean Delta normalized AP", title = "Incremental contribution of Hallmark features across BA4 anchor models") + 
    theme_nature(base_size = 13)

p_hallmark_inc


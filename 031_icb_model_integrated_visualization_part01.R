if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

library(ggplot2)

ICGC.enriched_top5_data <- ICGC.mutation_enrichment %>% filter(is.finite(P_value), delta_frequency != 
    0) %>% mutate(Group = factor(ifelse(delta_frequency > 0, "ICI Low", "ICI High"), levels = c("ICI High", 
    "ICI Low")), mutation_percent = ifelse(Group == "ICI High", freq_high, freq_low) * 100, significance = -log10(pmax(P_value, 
    1e-300))) %>% group_by(Group) %>% arrange(P_value, desc(mutation_percent), Gene, .by_group = TRUE) %>% 
    slice_head(n = 5) %>% ungroup() %>% arrange(Group, P_value, desc(mutation_percent), Gene)

ICGC.enriched_top5_data$Gene <- factor(ICGC.enriched_top5_data$Gene, levels = rev(ICGC.enriched_top5_data$Gene))

ICGC.group_colors <- c(`ICI High` = "#D95F02", `ICI Low` = "#1B9E77")

ICGC.p_enriched_top5 <- ggplot(ICGC.enriched_top5_data, aes(x = mutation_percent, y = Gene, color = Group)) + 
    geom_segment(aes(x = 0, xend = mutation_percent, yend = Gene), linewidth = 1, alpha = 0.8) + geom_point(aes(size = significance, 
    fill = Group), shape = 21, color = "black", stroke = 0.8) + geom_text(aes(label = sprintf("%.1f%%", 
    mutation_percent)), hjust = -0.65, size = 3.8, show.legend = FALSE) + scale_color_manual(values = ICGC.group_colors, 
    name = "Enriched group") + scale_fill_manual(values = ICGC.group_colors, guide = "none") + scale_size_continuous(name = expression(-log[10](P)), 
    range = c(3, 9)) + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = function(x) paste0(x, 
    "%"), expand = expansion(mult = c(0, 0.03))) + scale_y_discrete(expand = expansion(add = 0.7)) + 
    labs(title = "Top mutation enrichment candidates", subtitle = "ICGC-LICA-FR | Top 5 per group ranked by Fisher P", 
        x = "Mutation frequency within enriched group", y = NULL, caption = paste0("Overall mutation frequency ≥5%; point size: −log10(Fisher P).", 
            "Selection does not require statistical significance.")) + theme_classic(base_size = 13) + 
    theme(axis.text.y = element_text(color = "black", face = "italic", size = 11), plot.title = element_text(face = "bold", 
        size = 15), plot.subtitle = element_text(size = 10), plot.caption = element_text(hjust = 0, size = 9), 
        legend.position = "right")

if (all(c("ICI High", "ICI Low") %in% as.character(ICGC.enriched_top5_data$Group))) {
    ICGC.p_enriched_top5 <- ICGC.p_enriched_top5 + geom_hline(yintercept = sum(ICGC.enriched_top5_data$Group == 
        "ICI Low") + 0.5, linetype = "dashed", color = "grey70")
}

ICGC.visual_dir <- "ICB.model/ICGC_visualization"

dir.create(ICGC.visual_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(ICGC.visual_dir, "ICI_enriched_top5_lollipop.pdf"), plot = ICGC.p_enriched_top5, width = 9, 
    height = 6, bg = "white")

library(dplyr)

library(ggplot2)

TCGA.enriched_top5_data <- TCGA.mutation_enrichment_ranked_FDR %>% filter(is.finite(P_value), delta_frequency != 
    0) %>% mutate(Group = factor(ifelse(delta_frequency > 0, "ICI Low", "ICI High"), levels = c("ICI High", 
    "ICI Low")), mutation_percent = ifelse(Group == "ICI High", freq_high, freq_low) * 100, significance = -log10(pmax(P_value, 
    1e-300))) %>% group_by(Group) %>% arrange(P_value, desc(mutation_percent), Gene, .by_group = TRUE) %>% 
    slice_head(n = 5) %>% ungroup() %>% arrange(Group, P_value, desc(mutation_percent), Gene)

stopifnot(nrow(TCGA.enriched_top5_data) > 0, !anyDuplicated(TCGA.enriched_top5_data$Gene))

TCGA.enriched_top5_data$Gene <- factor(TCGA.enriched_top5_data$Gene, levels = rev(TCGA.enriched_top5_data$Gene))

TCGA.group_colors <- c(`ICI High` = "#D95F02", `ICI Low` = "#1B9E77")

TCGA.enriched_top5_data <- TCGA.mutation_enrichment %>% filter(is.finite(P_value), delta_frequency != 
    0) %>% mutate(Group = factor(ifelse(delta_frequency > 0, "ICI Low", "ICI High"), levels = c("ICI High", 
    "ICI Low")), mutation_percent = ifelse(Group == "ICI High", freq_high, freq_low) * 100, significance = -log10(pmax(P_value, 
    1e-300))) %>% group_by(Group) %>% arrange(P_value, desc(mutation_percent), Gene, .by_group = TRUE) %>% 
    slice_head(n = 5) %>% ungroup() %>% arrange(Group, desc(mutation_percent), P_value, Gene)

TCGA.enriched_top5_data$Gene <- factor(TCGA.enriched_top5_data$Gene, levels = rev(TCGA.enriched_top5_data$Gene))

TCGA.p_enriched_top5 <- ggplot(TCGA.enriched_top5_data, aes(x = mutation_percent, y = Gene, color = Group)) + 
    geom_segment(aes(x = 0, xend = mutation_percent, yend = Gene), linewidth = 1, alpha = 0.8) + geom_point(aes(size = significance, 
    fill = Group), shape = 21, color = "black", stroke = 0.8) + geom_text(aes(label = sprintf("%.1f%%", 
    mutation_percent)), hjust = -0.65, size = 3.8, show.legend = FALSE) + scale_color_manual(values = TCGA.group_colors, 
    name = "Enriched group") + scale_fill_manual(values = TCGA.group_colors, guide = "none") + scale_size_continuous(name = expression(-log[10](P)), 
    range = c(3, 9)) + scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = function(x) paste0(x, 
    "%"), expand = expansion(mult = c(0, 0.03))) + scale_y_discrete(expand = expansion(add = 0.7)) + 
    labs(title = "Top mutation enrichment candidates", subtitle = "TCGA-LIHC | Top 5 per group ranked by Fisher P", 
        x = "Mutation frequency", y = NULL) + theme_classic(base_size = 13) + theme(axis.text.y = element_text(color = "black", 
    face = "italic", size = 11), plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 10), 
    plot.caption = element_text(hjust = 0, size = 9), legend.position = "right")

if (all(c("ICI High", "ICI Low") %in% as.character(TCGA.enriched_top5_data$Group))) {
    TCGA.p_enriched_top5 <- TCGA.p_enriched_top5 + geom_hline(yintercept = sum(TCGA.enriched_top5_data$Group == 
        "ICI Low") + 0.5, linetype = "dashed", color = "grey70")
}

TCGA.visual_dir <- "ICB.model/TCGA_visualization"

dir.create(TCGA.visual_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(TCGA.visual_dir, "ICI_enriched_top5_lollipop.pdf"), plot = TCGA.p_enriched_top5, width = 5, 
    height = 4.5, bg = "white")

library(dplyr)

library(ggplot2)

ICGC.BA39_genes <- ICGC.CTNNB1_BA_gene_test %>% filter(FDR < 0.05) %>% pull(Gene)

ICGC.BA40_genes <- union(ICGC.BA39_genes, "BAAT")

stopifnot(length(ICGC.BA39_genes) == 39, length(ICGC.BA40_genes) == 40, "BAAT" %in% rownames(ICGC.log2_tpm))

ICGC.BAAT_ids <- ICGC.BA4_features$Patient_ID

ICGC.BAAT_expression <- as.numeric(ICGC.log2_tpm["BAAT", ICGC.BAAT_ids])

ICGC.BAAT_is_mut <- ICGC.mut["CTNNB1", ICGC.BAAT_ids] > 0

stopifnot(all(is.finite(ICGC.BAAT_expression)), any(ICGC.BAAT_is_mut), any(!ICGC.BAAT_is_mut))

ICGC.BAAT_volcano_test <- data.frame(Gene = "BAAT", delta_mean = mean(ICGC.BAAT_expression[ICGC.BAAT_is_mut]) - 
    mean(ICGC.BAAT_expression[!ICGC.BAAT_is_mut]), P_value = wilcox.test(ICGC.BAAT_expression[ICGC.BAAT_is_mut], 
    ICGC.BAAT_expression[!ICGC.BAAT_is_mut], alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value)

ICGC.BA113_results <- bind_rows(ICGC.CTNNB1_BA_gene_test %>% dplyr::select(Gene, delta_mean, P_value) %>% 
    filter(Gene != "BAAT"), ICGC.BAAT_volcano_test) %>% mutate(FDR_113 = p.adjust(P_value, method = "BH"))

stopifnot(nrow(ICGC.BA113_results) == 113, !anyDuplicated(ICGC.BA113_results$Gene))

ICGC.BA40_volcano_data <- ICGC.BA113_results %>% filter(Gene %in% ICGC.BA40_genes) %>% mutate(minus_log10_FDR = -log10(pmax(FDR_113, 
    1e-300)), Direction = factor(case_when(FDR_113 < 0.05 & delta_mean > 0 ~ "Higher in MUT", FDR_113 < 
    0.05 & delta_mean < 0 ~ "Lower in MUT", TRUE ~ "Not significant"), levels = c("Higher in MUT", "Lower in MUT", 
    "Not significant")))

ICGC.p_BA40_volcano <- ggplot(ICGC.BA40_volcano_data, aes(x = delta_mean, y = minus_log10_FDR)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey65") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
    color = "grey65") + geom_point(aes(fill = Direction), shape = 21, color = "black", stroke = 0.35, 
    size = 3.5, alpha = 0.9) + scale_fill_manual(values = c(`Higher in MUT` = "#D73027", `Lower in MUT` = "#2878B5", 
    `Not significant` = "#BDBDBD"), drop = FALSE, name = NULL) + scale_y_continuous(expand = expansion(mult = c(0.03, 
    0.15))) + labs(title = "CTNNB1-associated bile acid metabolism genes", subtitle = "ICGC-LICA-FR | 39 previously significant genes + BAAT", 
    x = "Mean log2(TPM + 1) difference", y = expression(-log[10](FDR))) + theme_classic(base_size = 13) + 
    theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 11), 
        plot.caption = element_text(hjust = 0, size = 9), legend.position = "top")

ICGC.BA40_label_data <- ICGC.BA40_volcano_data %>% filter(Gene %in% c("CYP8B1", "BAAT", "CYP27A1", "SLC27A5"))

if (requireNamespace("ggrepel", quietly = TRUE)) {
    ICGC.p_BA40_volcano <- ICGC.p_BA40_volcano + ggrepel::geom_text_repel(data = ICGC.BA40_label_data, 
        aes(label = Gene), color = "#D73027", fontface = "bold", size = 4, seed = 2026, max.overlaps = Inf, 
        box.padding = 0.5, point.padding = 0.4, min.segment.length = 0, show.legend = FALSE)
} else {
    ICGC.p_BA40_volcano <- ICGC.p_BA40_volcano + geom_text(data = ICGC.BA40_label_data, aes(label = Gene), 
        color = "#D73027", fontface = "bold", size = 8, vjust = -0.8, show.legend = FALSE)
}

ICGC.visual_dir <- "ICB.model/ICGC_visualization"

dir.create(ICGC.visual_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(ICGC.visual_dir, "CTNNB1_BA40_volcano.pdf"), plot = ICGC.p_BA40_volcano, width = 6, 
    height = 6, bg = "white")

data.table::fwrite(ICGC.BA40_volcano_data, file.path(ICGC.visual_dir, "CTNNB1_BA40_volcano_data.tsv"), 
    sep = "\t")

library(dplyr)

library(ggplot2)

TCGA.BA_volcano_genes <- union(TCGA.CTNNB1_BA_gene_test$Gene[TCGA.CTNNB1_BA_gene_test$FDR < 0.05], "BAAT")

TCGA.BAAT_ids <- TCGA.BA4_features$Patient_ID

TCGA.BAAT_expression <- as.numeric(TCGA.log2_tpm["BAAT", TCGA.BAAT_ids])

TCGA.BAAT_is_mut <- TCGA.mut["CTNNB1", TCGA.BAAT_ids] > 0

stopifnot(all(is.finite(TCGA.BAAT_expression)), !anyNA(TCGA.BAAT_is_mut), any(TCGA.BAAT_is_mut), any(!TCGA.BAAT_is_mut))

TCGA.BAAT_volcano_test <- data.frame(Gene = "BAAT", delta_mean = mean(TCGA.BAAT_expression[TCGA.BAAT_is_mut]) - 
    mean(TCGA.BAAT_expression[!TCGA.BAAT_is_mut]), P_value = wilcox.test(TCGA.BAAT_expression[TCGA.BAAT_is_mut], 
    TCGA.BAAT_expression[!TCGA.BAAT_is_mut], alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value)

TCGA.BA_volcano_all <- bind_rows(TCGA.CTNNB1_BA_gene_test %>% dplyr::select(Gene, delta_mean, P_value) %>% 
    filter(Gene != "BAAT"), TCGA.BAAT_volcano_test) %>% mutate(FDR_combined = p.adjust(P_value, method = "BH"))

stopifnot(!anyDuplicated(TCGA.BA_volcano_all$Gene))

TCGA.BA_volcano_data <- TCGA.BA_volcano_all %>% filter(Gene %in% TCGA.BA_volcano_genes) %>% mutate(minus_log10_FDR = -log10(pmax(FDR_combined, 
    1e-300)), Direction = factor(case_when(FDR_combined < 0.05 & delta_mean > 0 ~ "Higher in MUT", FDR_combined < 
    0.05 & delta_mean < 0 ~ "Lower in MUT", TRUE ~ "Not significant"), levels = c("Higher in MUT", "Lower in MUT", 
    "Not significant")))

TCGA.BA_label_data <- TCGA.BA_volcano_data %>% filter(Gene %in% c("CYP8B1", "BAAT", "CYP27A1", "SLC27A5"))

TCGA.p_BA_volcano <- ggplot(TCGA.BA_volcano_data, aes(x = delta_mean, y = minus_log10_FDR)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey65") + geom_hline(yintercept = -log10(0.05), linetype = "dashed", 
    color = "grey65") + geom_point(aes(fill = Direction), shape = 21, color = "black", stroke = 0.35, 
    size = 3.5, alpha = 0.9) + scale_fill_manual(values = c(`Higher in MUT` = "#D73027", `Lower in MUT` = "#2878B5", 
    `Not significant` = "#BDBDBD"), drop = FALSE, name = NULL) + xlim(-1, 2) + scale_y_continuous(expand = expansion(mult = c(0.03, 
    0.15))) + labs(title = "CTNNB1-associated bile acid metabolism genes", subtitle = paste0("TCGA-LIHC | Previously significant genes + BAAT (n = ", 
    nrow(TCGA.BA_volcano_data), ")"), x = "Mean log2(TPM + 1) difference", y = expression(-log[10](FDR))) + 
    theme_classic(base_size = 13) + theme(plot.title = element_text(face = "bold", size = 15), plot.subtitle = element_text(size = 11), 
    plot.caption = element_text(hjust = 0, size = 9), legend.position = "top")

TCGA.BA_label_all <- TCGA.BA_volcano_data %>% mutate(label = ifelse(Gene %in% c("CYP8B1", "BAAT", "CYP27A1", 
    "SLC27A5"), Gene, ""))

TCGA.p_BA_volcano <- TCGA.p_BA_volcano + ggrepel::geom_text_repel(data = TCGA.BA_label_all, aes(label = label), 
    color = "#D73027", fontface = "bold", size = 4, point.size = 3.5, point.padding = grid::unit(0.6, 
        "lines"), box.padding = grid::unit(0.7, "lines"), force = 5, force_pull = 0.2, max.overlaps = Inf, 
    max.time = 5, max.iter = 50000, min.segment.length = 0, segment.color = "grey45", segment.size = 0.4, 
    seed = 2026, show.legend = FALSE)

TCGA.visual_dir <- "ICB.model/TCGA_visualization"

dir.create(TCGA.visual_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(TCGA.visual_dir, "CTNNB1_BA_volcano.pdf"), plot = TCGA.p_BA_volcano, width = 6, height = 6, 
    bg = "white")

data.table::fwrite(TCGA.BA_volcano_data, file.path(TCGA.visual_dir, "CTNNB1_BA_volcano_data.tsv"), sep = "\t")

library(dplyr)

compare_CTNNB1_features <- function(features, mutation, cohort) {
    feature_names <- BA4_MODEL_BUNDLE$features
    ids <- features$Patient_ID
    stopifnot(!anyDuplicated(ids), all(ids %in% colnames(mutation)), "CTNNB1" %in% rownames(mutation), 
        all(feature_names %in% names(features)))
    mut <- as.numeric(mutation["CTNNB1", ids]) > 0
    stopifnot(!anyNA(mut), any(mut), any(!mut))
    result <- bind_rows(lapply(feature_names, function(f) {
        x <- features[[f]]
        stopifnot(is.numeric(x), all(is.finite(x)))
        x_mut <- x[mut]
        x_wt <- x[!mut]
        n_mut <- length(x_mut)
        n_wt <- length(x_wt)
        ranks <- rank(c(x_mut, x_wt), ties.method = "average")
        U <- sum(ranks[seq_len(n_mut)]) - n_mut * (n_mut + 1)/2
        data.frame(cohort = cohort, feature = f, n_MUT = n_mut, n_WT = n_wt, median_MUT = median(x_mut), 
            median_WT = median(x_wt), delta_median = median(x_mut) - median(x_wt), rank_biserial = 2 * 
                U/(n_mut * n_wt) - 1, P_value = if (length(unique(x)) == 1L) {
                1
            }
            else {
                wilcox.test(x_mut, x_wt, exact = FALSE, correct = TRUE)$p.value
            })
    }))
    result %>% mutate(FDR = p.adjust(P_value, method = "BH"), rank_positive_effect = min_rank(desc(rank_biserial)), 
        rank_P = min_rank(P_value)) %>% arrange(desc(rank_biserial))
}

TCGA.feature_comparison <- compare_CTNNB1_features(TCGA.BA4_features, TCGA.mut, "TCGA")

ICGC.feature_comparison <- compare_CTNNB1_features(ICGC.BA4_features, ICGC.mut, "ICGC")

BA4.feature_comparison <- bind_rows(TCGA.feature_comparison, ICGC.feature_comparison)

BA4.feature_consistency <- TCGA.feature_comparison %>% dplyr::select(feature, effect_TCGA = rank_biserial, 
    FDR_TCGA = FDR) %>% inner_join(ICGC.feature_comparison %>% dplyr::select(feature, effect_ICGC = rank_biserial, 
    FDR_ICGC = FDR), by = "feature") %>% mutate(same_direction = sign(effect_TCGA) == sign(effect_ICGC), 
    positive_in_both = effect_TCGA > 0 & effect_ICGC > 0, FDR_significant_in_both = FDR_TCGA < 0.05 & 
        FDR_ICGC < 0.05)

BA4.bile_acid_comparison <- BA4.feature_comparison %>% filter(feature == "HALLMARK_BILE_ACID_METABOLISM") %>% 
    dplyr::select(cohort, n_MUT, n_WT, rank_biserial, P_value, FDR, rank_positive_effect, rank_P)

library(dplyr)

library(ggplot2)

BA4.bootstrap_features <- function(features, mutation, cohort, B = 2000, seed = 2026) {
    feature_names <- BA4_MODEL_BUNDLE$features
    ids <- features$Patient_ID
    X <- as.matrix(features[, feature_names, drop = FALSE])
    mut <- as.numeric(mutation["CTNNB1", ids]) > 0
    stopifnot(!anyDuplicated(ids), all(is.finite(X)), !anyNA(mut), sum(mut) >= 2, sum(!mut) >= 2)
    X_mut <- X[mut, , drop = FALSE]
    X_wt <- X[!mut, , drop = FALSE]
    effect_function <- function(a, b) {
        n_a <- length(a)
        n_b <- length(b)
        U <- sum(rank(c(a, b), ties.method = "average")[seq_len(n_a)]) - n_a * (n_a + 1)/2
        2 * U/(n_a * n_b) - 1
    }
    observed <- vapply(seq_along(feature_names), function(j) effect_function(X_mut[, j], X_wt[, j]), 
        numeric(1))
    set.seed(seed)
    boot_effects <- replicate(B, {
        i_mut <- sample.int(nrow(X_mut), nrow(X_mut), replace = TRUE)
        i_wt <- sample.int(nrow(X_wt), nrow(X_wt), replace = TRUE)
        vapply(seq_along(feature_names), function(j) {
            effect_function(X_mut[i_mut, j], X_wt[i_wt, j])
        }, numeric(1))
    })
    CI <- t(apply(boot_effects, 1, quantile, probs = c(0.025, 0.975), names = FALSE))
    data.frame(cohort = cohort, feature = feature_names, effect = observed, CI_lower = CI[, 1], CI_upper = CI[, 
        2], n_MUT = nrow(X_mut), n_WT = nrow(X_wt))
}

TCGA.feature_bootstrap <- BA4.bootstrap_features(TCGA.BA4_features, TCGA.mut, "TCGA", seed = 2026)

ICGC.feature_bootstrap <- BA4.bootstrap_features(ICGC.BA4_features, ICGC.mut, "ICGC", seed = 2027)

BA4.forest_data <- bind_rows(TCGA.feature_bootstrap, ICGC.feature_bootstrap) %>% left_join(BA4.feature_comparison %>% 
    dplyr::select(cohort, feature, FDR), by = c("cohort", "feature")) %>% mutate(cohort = factor(cohort, 
    levels = c("TCGA", "ICGC")), feature_label = recode(feature, HALLMARK_BILE_ACID_METABOLISM = "Bile acid metabolism", 
    HALLMARK_FATTY_ACID_METABOLISM = "Fatty acid metabolism", HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY = "ROS pathway", 
    ICI__HCC_INFLAMMATORY4 = "Inflammatory4"), feature_label = factor(feature_label, levels = c("Inflammatory4", 
    "ROS pathway", "Fatty acid metabolism", "Bile acid metabolism")), plot_y = as.numeric(feature_label) + 
    ifelse(cohort == "TCGA", 0.13, -0.13), result_label = sprintf("%.2f [%.2f, %.2f]; FDR = %.2g", effect, 
    CI_lower, CI_upper, FDR))

BA4.p_feature_forest <- ggplot(BA4.forest_data, aes(y = plot_y, color = cohort)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey55") + geom_segment(aes(x = CI_lower, xend = CI_upper, yend = plot_y), 
    linewidth = 0.9) + geom_segment(aes(x = CI_lower, xend = CI_lower, y = plot_y - 0.04, yend = plot_y + 
    0.04), linewidth = 0.7) + geom_segment(aes(x = CI_upper, xend = CI_upper, y = plot_y - 0.04, yend = plot_y + 
    0.04), linewidth = 0.7) + geom_point(aes(x = effect, fill = cohort), shape = 21, color = "black", 
    stroke = 0.35, size = 3.5) + geom_text(aes(x = 1.05, label = result_label), hjust = 0, size = 3.5, 
    show.legend = FALSE) + scale_color_manual(values = c(TCGA = "#D73027", ICGC = "#2878B5"), name = "Cohort") + 
    scale_fill_manual(values = c(TCGA = "#D73027", ICGC = "#2878B5"), guide = "none") + scale_y_continuous(breaks = 1:4, 
    labels = levels(BA4.forest_data$feature_label), limits = c(0.6, 4.4), expand = expansion(mult = 0)) + 
    scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 2.1)) + labs(title = "CTNNB1 mutation-associated ICI model features", 
    subtitle = "Points: rank-biserial effect; intervals: bootstrap 95% CI", x = "Rank-biserial correlation", 
    y = NULL, caption = paste0("Right-hand labels: effect [95% CI]; within-cohort four-feature BH-FDR.", 
        "Confidence intervals: 2,000 stratified patient bootstrap resamples.")) + theme_classic(base_size = 13) + 
    theme(axis.text.y = element_text(color = "black", size = 12), axis.ticks.y = element_blank(), plot.title = element_text(face = "bold", 
        size = 15), plot.subtitle = element_text(size = 11), plot.caption = element_text(hjust = 0, size = 9), 
        legend.position = "top")

BA4.visual_dir <- "ICB.model/BA4_comparison"

dir.create(BA4.visual_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(BA4.visual_dir, "CTNNB1_four_feature_forest.pdf"), plot = BA4.p_feature_forest, width = 8, 
    height = 5, bg = "white")

data.table::fwrite(BA4.forest_data, file.path(BA4.visual_dir, "CTNNB1_four_feature_effects.tsv"), sep = "\t")

library(dplyr)

library(ggplot2)

library(pROC)

BA4.roc_data <- BA4_OOCO_predictions %>% dplyr::select(cohort, response, OOCO_ICBScore) %>% mutate(cohort = as.character(cohort), 
    response = as.numeric(as.character(response)))

BA4.cohort_order <- c("GSE181946", "GSE202069", "PRJEB34724")

stopifnot(setequal(unique(BA4.roc_data$cohort), BA4.cohort_order), !anyNA(BA4.roc_data$response), all(BA4.roc_data$response %in% 
    c(0, 1)), all(is.finite(BA4.roc_data$OOCO_ICBScore)))

BA4.cohort_colors <- c(GSE181946 = "#D95F02", GSE202069 = "#1B9E77", PRJEB34724 = "#7570B3")

set.seed(2026)

BA4.roc_objects <- setNames(lapply(BA4.cohort_order, function(cc) {
    dat <- BA4.roc_data %>% filter(cohort == cc)
    stopifnot(length(unique(dat$response)) == 2)
    pROC::roc(response = dat$response, predictor = dat$OOCO_ICBScore, levels = c(0, 1), direction = "<", 
        quiet = TRUE)
}), BA4.cohort_order)

BA4.roc_summary <- bind_rows(lapply(BA4.cohort_order, function(cc) {
    roc_obj <- BA4.roc_objects[[cc]]
    ci <- pROC::ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
    data.frame(cohort = cc, n = sum(BA4.roc_data$cohort == cc), AUC = as.numeric(pROC::auc(roc_obj)), 
        CI_lower = as.numeric(ci[1]), CI_upper = as.numeric(ci[3]))
}))

BA4.roc_curve_data <- bind_rows(lapply(BA4.cohort_order, function(cc) {
    roc_obj <- BA4.roc_objects[[cc]]
    data.frame(cohort = cc, FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities) %>% arrange(FPR, 
        TPR)
}))

BA4.roc_legend <- setNames(sprintf("%s\nAUC %.3f (95%% CI %.3f–%.3f)", BA4.roc_summary$cohort, BA4.roc_summary$AUC, 
    BA4.roc_summary$CI_lower, BA4.roc_summary$CI_upper), BA4.roc_summary$cohort)

BA4.roc_legend <- setNames(sprintf("%s  AUC = %.3f", BA4.roc_summary$cohort, BA4.roc_summary$AUC), BA4.roc_summary$cohort)

BA4.roc_legend_order <- BA4.roc_summary %>% arrange(desc(AUC)) %>% pull(cohort)

BA4.p_OOCO_ROC <- BA4.p_OOCO_ROC + scale_color_manual(values = BA4.cohort_colors, breaks = BA4.cohort_order, 
    labels = BA4.roc_legend, name = NULL) + theme(panel.border = element_rect(color = "black", fill = NA, 
    linewidth = 0.7), axis.line = element_blank(), legend.position = c(0.97, 0.04), legend.justification = c(1, 
    0), legend.background = element_rect(fill = "white", color = NA), legend.text = element_text(size = 10), 
    legend.key.width = grid::unit(0.8, "cm"), legend.key.height = grid::unit(0.5, "cm"), legend.margin = margin(5, 
        5, 5, 5))

BA4.p_OOCO_ROC <- BA4.p_OOCO_ROC + scale_color_manual(values = BA4.cohort_colors, breaks = BA4.roc_legend_order, 
    labels = BA4.roc_legend[BA4.roc_legend_order], name = NULL) + guides(color = guide_legend(ncol = 1, 
    reverse = FALSE))

BA4.p_OOCO_ROC <- BA4.p_OOCO_ROC + theme(text = element_text(size = 16), plot.title = element_text(size = 21, 
    face = "bold"), plot.subtitle = element_text(size = 16), axis.title = element_text(size = 18), axis.text = element_text(size = 15, 
    color = "black"), legend.text = element_text(size = 14), legend.key.width = grid::unit(1, "cm"), 
    legend.key.height = grid::unit(0.65, "cm"))

BA4.roc_dir <- "ICB.model/BA4_visualization"

dir.create(BA4.roc_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(BA4.roc_dir, "BA4_leave_one_cohort_out_ROC.pdf"), plot = BA4.p_OOCO_ROC, width = 9, 
    height = 6, bg = "white")

library(dplyr)

library(ggplot2)

BA4.make_PR <- function(cc) {
    dat <- BA4_OOCO_predictions %>% filter(cohort == cc) %>% transmute(response = as.numeric(as.character(response)), 
        score = OOCO_ICBScore)
    stopifnot(!anyNA(dat$response), all(dat$response %in% c(0, 1)), all(is.finite(dat$score)), length(unique(dat$response)) == 
        2)
    n_positive <- sum(dat$response == 1)
    prevalence <- mean(dat$response)
    curve <- dat %>% group_by(score) %>% summarise(positives = sum(response == 1), negatives = sum(response == 
        0), .groups = "drop") %>% arrange(desc(score)) %>% mutate(TP = cumsum(positives), FP = cumsum(negatives), 
        Recall = TP/n_positive, Precision = TP/(TP + FP), delta_recall = Recall - lag(Recall, default = 0))
    AP <- sum(curve$delta_recall * curve$Precision)
    list(curve = bind_rows(data.frame(Recall = 0, Precision = 1), curve %>% dplyr::select(Recall, Precision)) %>% 
        mutate(cohort = cc), metrics = data.frame(cohort = cc, n = nrow(dat), responders = n_positive, 
        prevalence = prevalence, AP = AP, AP_lift = AP/prevalence))
}

BA4.PR_results <- lapply(BA4.cohort_order, BA4.make_PR)

BA4.PR_curve <- bind_rows(lapply(BA4.PR_results, function(x) x$curve))

BA4.PR_metrics <- bind_rows(lapply(BA4.PR_results, function(x) x$metrics)) %>% arrange(desc(AP))

BA4.PR_display <- BA4.PR_curve %>% group_by(cohort, Recall) %>% summarise(Precision = max(Precision), 
    .groups = "drop") %>% arrange(cohort, Recall)

BA4.PR_check <- BA4.PR_display %>% group_by(cohort) %>% summarise(AP_display = sum((Recall - lag(Recall, 
    default = 0)) * Precision), .groups = "drop") %>% left_join(BA4.PR_metrics %>% dplyr::select(cohort, 
    AP), by = "cohort")


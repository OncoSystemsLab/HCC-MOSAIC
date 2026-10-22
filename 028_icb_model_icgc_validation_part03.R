if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

ICGC.test_BA_gene <- function(gene) {
    x <- unname(ICGC.BA_expr[gene, ])
    mut <- ICGC.CTNNB1_status_vector == "MUT"
    x_mut <- x[mut]
    x_wt <- x[!mut]
    n_mut <- length(x_mut)
    n_wt <- length(x_wt)
    stopifnot(n_mut > 0, n_wt > 0)
    cohort_sd <- sd(x)
    if (length(unique(x)) == 1L) {
        p_value <- 1
        rank_biserial <- 0
    }
    else {
        p_value <- wilcox.test(x_mut, x_wt, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
        ranks <- rank(c(x_mut, x_wt), ties.method = "average")
        U_mut <- sum(ranks[seq_len(n_mut)]) - n_mut * (n_mut + 1)/2
        rank_biserial <- 2 * U_mut/(n_mut * n_wt) - 1
    }
    delta_mean <- mean(x_mut) - mean(x_wt)
    data.frame(Gene = gene, n_MUT = n_mut, n_WT = n_wt, mean_MUT = mean(x_mut), mean_WT = mean(x_wt), 
        median_MUT = median(x_mut), median_WT = median(x_wt), delta_mean = delta_mean, delta_median = median(x_mut) - 
            median(x_wt), mean_difference_per_cohort_SD = if (cohort_sd > 0) {
            delta_mean/cohort_sd
        }
        else {
            NA_real_
        }, rank_biserial = rank_biserial, P_value = p_value)
}

ICGC.CTNNB1_BA_gene_test <- bind_rows(lapply(ICGC.BA_genes, ICGC.test_BA_gene)) %>% mutate(FDR = p.adjust(P_value, 
    method = "BH"), direction = case_when(delta_mean > 0 ~ "Higher_in_CTNNB1_MUT", delta_mean < 0 ~ "Lower_in_CTNNB1_MUT", 
    TRUE ~ "Equal_means")) %>% arrange(FDR, P_value, desc(abs(mean_difference_per_cohort_SD)))

ICGC.CTNNB1_BA_significant <- ICGC.CTNNB1_BA_gene_test %>% filter(FDR < 0.05)

ICGC.CTNNB1_BA_up <- ICGC.CTNNB1_BA_significant %>% filter(delta_mean > 0)

ICGC.CTNNB1_BA_down <- ICGC.CTNNB1_BA_significant %>% filter(delta_mean < 0)

ICGC.BA_gene_direction_summary <- data.frame(total_BA_genes = length(ICGC.BA_genes), significant = nrow(ICGC.CTNNB1_BA_significant), 
    higher_in_MUT = nrow(ICGC.CTNNB1_BA_up), lower_in_MUT = nrow(ICGC.CTNNB1_BA_down))

ICGC.BA_key_genes <- c("CYP7A1", "CYP8B1", "CYP27A1", "HSD3B7", "AKR1D1", "SLC27A5", "BAAT")

ICGC.BA_key_gene_results <- data.frame(Gene = ICGC.BA_key_genes) %>% mutate(present_in_expression = Gene %in% 
    rownames(ICGC.log2_tpm), included_in_BA_test = Gene %in% ICGC.BA_genes) %>% left_join(ICGC.CTNNB1_BA_gene_test, 
    by = "Gene")

ICGC.BA_volcano_data <- ICGC.CTNNB1_BA_gene_test %>% mutate(minus_log10_FDR = -log10(pmax(FDR, 1e-300)), 
    plot_group = case_when(FDR < 0.05 & delta_mean > 0 ~ "Higher in MUT", FDR < 0.05 & delta_mean < 0 ~ 
        "Lower in MUT", TRUE ~ "Not significant"), label = ifelse(Gene %in% ICGC.BA_key_genes | (FDR < 
        0.01 & abs(mean_difference_per_cohort_SD) > 0.4), Gene, "NA"))

ICGC.p_BA_volcano <- ggplot(ICGC.BA_volcano_data, aes(x = delta_mean, y = minus_log10_FDR)) + geom_point(aes(color = plot_group), 
    size = 2.3, alpha = 0.8) + scale_color_manual(values = c(`Higher in MUT` = "#E41A1C", `Lower in MUT` = "#377EB8", 
    `Not significant` = "grey65")) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = -log10(0.05), 
    linetype = 2) + labs(x = "Mean log2(TPM + 1) difference (MUT - WT)", y = expression(-log[10](FDR)), 
    color = NULL, title = "CTNNB1-associated bile acid pathway genes", subtitle = "ICGC-LICA-FR") + theme_classic(base_size = 13)

ICGC.BA_label_data <- ICGC.BA_volcano_data %>% filter(!is.na(label))

if (requireNamespace("ggrepel", quietly = TRUE)) {
    ICGC.p_BA_volcano <- ICGC.p_BA_volcano + ggrepel::geom_text_repel(data = ICGC.BA_label_data, aes(label = label), 
        seed = 2026, size = 3.3, max.overlaps = Inf)
} else {
    ICGC.p_BA_volcano <- ICGC.p_BA_volcano + geom_text(data = ICGC.BA_label_data, aes(label = label), 
        vjust = -0.7, size = 3.3, check_overlap = TRUE)
}

stopifnot("BAAT" %in% rownames(ICGC.log2_tpm))

ICGC.BAAT_data <- ICGC.BA4_features %>% dplyr::select(Patient_ID, BA4_ICBScore, HALLMARK_BILE_ACID_METABOLISM) %>% 
    mutate(BAAT = as.numeric(ICGC.log2_tpm["BAAT", Patient_ID]), CTNNB1_status = factor(ifelse(ICGC.mut["CTNNB1", 
        Patient_ID] > 0, "MUT", "WT"), levels = c("WT", "MUT")))

ICGC.BAAT_summary <- ICGC.BAAT_data %>% group_by(CTNNB1_status) %>% summarise(n = n(), mean = mean(BAAT), 
    median = median(BAAT), SD = sd(BAAT), .groups = "drop")

ICGC.BAAT_test <- ICGC.CTNNB1_BA_gene_test %>% filter(Gene == "BAAT")

ICGC.p_BAAT <- ggplot(ICGC.BAAT_data, aes(x = CTNNB1_status, y = BAAT, fill = CTNNB1_status)) + geom_boxplot(width = 0.55, 
    outlier.shape = NA, alpha = 0.65) + geom_point(position = position_jitter(width = 0.12, height = 0, 
    seed = 2026), alpha = 0.5, size = 1.5) + scale_fill_manual(values = c(WT = "#377EB8", MUT = "#E41A1C")) + 
    labs(x = "CTNNB1 status", y = "BAAT log2(TPM + 1)", title = "BAAT expression in ICGC-LICA-FR") + 
    theme_classic(base_size = 13) + theme(legend.position = "none")

ICGC.BAAT_correlation_results <- bind_rows(lapply(c("HALLMARK_BILE_ACID_METABOLISM", "BA4_ICBScore"), 
    function(target) {
        x <- ICGC.BAAT_data$BAAT
        y <- ICGC.BAAT_data[[target]]
        valid <- is.finite(x) & is.finite(y)
        x <- x[valid]
        y <- y[valid]
        testable <- length(x) >= 3 && length(unique(x)) > 1 && length(unique(y)) > 1
        res <- if (testable) {
            cor.test(x, y, method = "spearman", exact = FALSE)
        }
        else {
            NULL
        }
        data.frame(target = target, n = length(x), Spearman_rho = if (testable) {
            unname(res$estimate)
        }
        else {
            NA_real_
        }, P_value = if (testable) {
            res$p.value
        }
        else {
            NA_real_
        })
    })) %>% mutate(FDR_two_correlations = p.adjust(P_value, method = "BH"))

ICGC.BAAT_mut <- with(ICGC.BAAT_data, BAAT[CTNNB1_status == "MUT"])

ICGC.BAAT_wt <- with(ICGC.BAAT_data, BAAT[CTNNB1_status == "WT"])

ICGC.BAAT_wilcox <- wilcox.test(ICGC.BAAT_mut, ICGC.BAAT_wt, alternative = "two.sided", exact = FALSE, 
    correct = TRUE)

ICGC.BAAT_n_mut <- length(ICGC.BAAT_mut)

ICGC.BAAT_n_wt <- length(ICGC.BAAT_wt)

ICGC.BAAT_ranks <- rank(c(ICGC.BAAT_mut, ICGC.BAAT_wt), ties.method = "average")

ICGC.BAAT_U <- sum(ICGC.BAAT_ranks[seq_len(ICGC.BAAT_n_mut)]) - ICGC.BAAT_n_mut * (ICGC.BAAT_n_mut + 
    1)/2

ICGC.BAAT_test <- data.frame(Gene = "BAAT", analysis = "Separate single-gene analysis", included_in_frozen_BA_set = "BAAT" %in% 
    ICGC.BA_genes, n_MUT = ICGC.BAAT_n_mut, n_WT = ICGC.BAAT_n_wt, mean_MUT = mean(ICGC.BAAT_mut), mean_WT = mean(ICGC.BAAT_wt), 
    median_MUT = median(ICGC.BAAT_mut), median_WT = median(ICGC.BAAT_wt), delta_mean = mean(ICGC.BAAT_mut) - 
        mean(ICGC.BAAT_wt), delta_median = median(ICGC.BAAT_mut) - median(ICGC.BAAT_wt), rank_biserial = 2 * 
        ICGC.BAAT_U/(ICGC.BAAT_n_mut * ICGC.BAAT_n_wt) - 1, P_value = ICGC.BAAT_wilcox$p.value)

ICGC.p_BAAT <- ICGC.p_BAAT + labs(subtitle = sprintf("Separate BAAT analysis: Wilcoxon P = %.3g", ICGC.BAAT_test$P_value))


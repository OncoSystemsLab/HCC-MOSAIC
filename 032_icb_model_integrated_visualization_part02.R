if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

stopifnot(all(abs(BA4.PR_check$AP_display - BA4.PR_check$AP) < 1e-10))

BA4.PR_legend <- setNames(sprintf("%s  AP = %.3f | Responders = %.1f%%", BA4.PR_metrics$cohort, BA4.PR_metrics$AP, 
    100 * BA4.PR_metrics$prevalence), BA4.PR_metrics$cohort)

BA4.p_OOCO_PR <- ggplot(BA4.PR_curve, aes(x = Recall, y = Precision, color = cohort)) + geom_hline(data = BA4.PR_metrics, 
    aes(yintercept = prevalence, color = cohort), linetype = "dashed", linewidth = 0.6, alpha = 0.5, 
    show.legend = FALSE) + geom_step(data = BA4.PR_display, direction = "vh", linewidth = 1.1) + scale_color_manual(values = BA4.cohort_colors, 
    breaks = BA4.PR_metrics$cohort, labels = BA4.PR_legend[BA4.PR_metrics$cohort], name = NULL) + scale_x_continuous(breaks = seq(0, 
    1, 0.2), expand = expansion(mult = 0)) + scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = 0)) + 
    coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1)) + labs(title = "BA4 model precision–recall curves", 
    subtitle = "Leave-one-cohort-out predictions", x = "Recall (sensitivity)", y = "Precision (positive predictive value)") + 
    theme_classic(base_size = 16) + theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7), 
    axis.line = element_blank(), plot.title = element_text(size = 21, face = "bold"), plot.subtitle = element_text(size = 16), 
    axis.title = element_text(size = 18), axis.text = element_text(size = 15, color = "black"), legend.position = "bottom", 
    legend.text = element_text(size = 12), legend.key.width = grid::unit(1, "cm"), plot.caption = element_text(size = 10, 
        hjust = 0)) + guides(color = guide_legend(ncol = 1))

BA4.p_OOCO_PR <- BA4.p_OOCO_PR + theme(legend.position = c(0.5, 0.02), legend.justification = c(0.5, 
    0), legend.background = element_rect(fill = "white", color = NA), legend.text = element_text(size = 11), 
    legend.key.width = grid::unit(0.8, "cm"), legend.key.height = grid::unit(0.45, "cm"), legend.margin = margin(4, 
        6, 4, 6)) + guides(color = guide_legend(ncol = 1))

dir.create(BA4.roc_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(BA4.roc_dir, "BA4_leave_one_cohort_out_PR.pdf"), plot = BA4.p_OOCO_PR, width = 9, height = 6, 
    bg = "white")

data.table::fwrite(BA4.PR_metrics, file.path(BA4.roc_dir, "BA4_leave_one_cohort_out_PR_metrics.tsv"), 
    sep = "\t")

library(dplyr)

library(ggplot2)

BA4.density_data <- BA4_OOCO_predictions %>% transmute(score = OOCO_ICBScore, response = as.numeric(as.character(response)))

stopifnot(all(is.finite(BA4.density_data$score)), !anyNA(BA4.density_data$response), all(BA4.density_data$response %in% 
    c(0, 1)), length(unique(BA4.density_data$response)) == 2)

BA4.density_data <- BA4.density_data %>% mutate(Group = factor(ifelse(response == 1, "Responder", "Non-responder"), 
    levels = c("Non-responder", "Responder")))

BA4.score_R <- BA4.density_data$score[BA4.density_data$response == 1]

BA4.score_NR <- BA4.density_data$score[BA4.density_data$response == 0]

BA4.density_stats <- data.frame(n_nonresponder = length(BA4.score_NR), n_responder = length(BA4.score_R), 
    median_nonresponder = median(BA4.score_NR), median_responder = median(BA4.score_R), MD = median(BA4.score_R) - 
        median(BA4.score_NR), P_value = wilcox.test(BA4.score_R, BA4.score_NR, alternative = "two.sided", 
        exact = FALSE, correct = TRUE)$p.value)

BA4.density_binwidth <- 0.05

BA4.density_breaks <- seq(floor(min(BA4.density_data$score)/BA4.density_binwidth) * BA4.density_binwidth - 
    BA4.density_binwidth, ceiling(max(BA4.density_data$score)/BA4.density_binwidth) * BA4.density_binwidth + 
    BA4.density_binwidth, by = BA4.density_binwidth)

BA4.density_outline <- bind_rows(lapply(levels(BA4.density_data$Group), function(g) {
    h <- hist(BA4.density_data$score[BA4.density_data$Group == g], breaks = BA4.density_breaks, plot = FALSE)
    data.frame(score = c(h$breaks[1], as.vector(rbind(head(h$breaks, -1), tail(h$breaks, -1))), tail(h$breaks, 
        1)), density = c(0, rep(h$density, each = 2), 0), Group = g)
}))

BA4.density_labels <- c(`Non-responder` = sprintf("Non-responder (n = %d)", length(BA4.score_NR)), Responder = sprintf("Responder (n = %d)", 
    length(BA4.score_R)))

BA4.p_score_density <- ggplot(BA4.density_outline, aes(x = score, y = density, color = Group)) + geom_path(linewidth = 1.2) + 
    scale_color_manual(values = c(`Non-responder` = "#2166AC", Responder = "#D73027"), breaks = c("Non-responder", 
        "Responder"), labels = BA4.density_labels, name = NULL) + scale_x_continuous(breaks = seq(0, 
    1, 0.1), expand = expansion(mult = c(0.03, 0.03))) + scale_y_continuous(expand = expansion(mult = c(0, 
    0.28))) + labs(title = sprintf("MD = %.3f     P = %.3g", BA4.density_stats$MD, BA4.density_stats$P_value), 
    subtitle = "Pooled leave-one-cohort-out predictions", x = "ICI Score", y = "Probability density", 
    caption = paste0("MD: median score difference (Responder − Non-responder).", "P: two-sided Wilcoxon test; pooled, unadjusted for cohort.")) + 
    theme_classic(base_size = 16) + theme(plot.title = element_text(size = 19), plot.subtitle = element_text(size = 12), 
    axis.title = element_text(size = 18), axis.text = element_text(size = 14, color = "black"), legend.position = c(0.03, 
        0.98), legend.justification = c(0, 1), legend.background = element_rect(fill = "white", color = NA), 
    legend.text = element_text(size = 13), legend.key.height = grid::unit(0.5, "cm"), plot.caption = element_text(size = 10, 
        hjust = 0))

BA4.density_dir <- "ICB.model/BA4_visualization"

dir.create(BA4.density_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(BA4.density_dir, "BA4_score_density_outline.pdf"), plot = BA4.p_score_density, width = 6, 
    height = 6, bg = "white")

library(data.table)

library(dplyr)

library(readxl)

PRJEB.wes_dir <- paste0("ICB.model/PRJEB_capture_CDS80/RNA_WES_overlap_29/")

PRJEB.mut_raw <- fread(file.path(PRJEB.wes_dir, "Gene_Mutation.tsv"), data.table = FALSE, check.names = FALSE)

PRJEB.clin_original <- read_excel("raw_data/PRJEB34724/PRJEB34724_clin.xlsx", sheet = "Sheet1")

stopifnot("gene" %in% names(PRJEB.mut_raw), !anyNA(PRJEB.mut_raw$gene), !anyDuplicated(PRJEB.mut_raw$gene), 
    all(c("sample", "label") %in% names(PRJEB.clin_original)))

PRJEB.mut <- as.matrix(PRJEB.mut_raw[, names(PRJEB.mut_raw) != "gene", drop = FALSE])

stopifnot(is.numeric(PRJEB.mut))

rownames(PRJEB.mut) <- PRJEB.mut_raw$gene

stopifnot(!anyDuplicated(colnames(PRJEB.mut)), all(is.na(PRJEB.mut) | PRJEB.mut %in% c(0, 1)), "CTNNB1" %in% 
    rownames(PRJEB.mut))

PRJEB.clin_match <- PRJEB.clin_original %>% transmute(RNA_sample = trimws(as.character(sample)), Patient_ID = sub("-T-RNA$", 
    "", RNA_sample), label_original = trimws(as.character(label)), response_group = case_when(tolower(label_original) == 
    "non-responder" ~ "Non-responder", tolower(label_original) == "responder" ~ "Responder", TRUE ~ "NA"))

stopifnot(!anyNA(PRJEB.clin_match$Patient_ID), !anyDuplicated(PRJEB.clin_match$Patient_ID))

PRJEB.match_QC <- data.frame(clinical_patients = nrow(PRJEB.clin_match), mutation_patients = ncol(PRJEB.mut), 
    matched_patients = sum(colnames(PRJEB.mut) %in% PRJEB.clin_match$Patient_ID))

stopifnot(all(colnames(PRJEB.mut) %in% PRJEB.clin_match$Patient_ID))

library(dplyr)

library(ggplot2)

PRJEB.plot_genes <- c("CTNNB1", "TP53")

stopifnot(all(PRJEB.plot_genes %in% rownames(PRJEB.mut)))

PRJEB.gene_mutation_results <- setNames(vector("list", 2), PRJEB.plot_genes)

for (gene in PRJEB.plot_genes) {
    analysis_data <- data.frame(Patient_ID = colnames(PRJEB.mut), mutation = as.numeric(PRJEB.mut[gene, 
        ])) %>% left_join(PRJEB.clin_match, by = "Patient_ID") %>% mutate(response_group = factor(response_group, 
        levels = c("Non-responder", "Responder"))) %>% filter(!is.na(response_group), !is.na(mutation))
    tab <- with(analysis_data, table(response_group, factor(mutation, levels = c(1, 0), labels = c("MUT", 
        "Negative"))))
    stopifnot(all(rowSums(tab) > 0))
    test <- fisher.test(tab, alternative = "two.sided")
    frequency <- analysis_data %>% group_by(response_group) %>% summarise(n = n(), n_MUT = sum(mutation == 
        1), n_negative = sum(mutation == 0), .groups = "drop") %>% rowwise() %>% mutate(Gene = gene, 
        mutation_percent = 100 * n_MUT/n, CI_lower = 100 * binom.test(n_MUT, n)$conf.int[1], CI_upper = 100 * 
            binom.test(n_MUT, n)$conf.int[2], label = sprintf("%d/%d (%.1f%%)", n_MUT, n, mutation_percent)) %>% 
        ungroup()
    PRJEB.gene_mutation_results[[gene]] <- list(analysis_data = analysis_data, frequency = frequency, 
        contingency_table = tab, fisher_result = data.frame(Gene = gene, comparison = "Non-responder vs Responder", 
            odds_ratio = unname(test$estimate), CI_lower = test$conf.int[1], CI_upper = test$conf.int[2], 
            P_value = test$p.value))
}

library(dplyr)

library(ggplot2)

PRJEB.two_gene_plot_data <- PRJEB.two_gene_plot_data %>% mutate(Gene = factor(Gene, levels = c("CTNNB1", 
    "TP53")), response_group = factor(response_group, levels = c("Non-responder", "Responder")))

PRJEB.response_dodge <- position_dodge(width = 0.55)

PRJEB.p_two_gene_mutation_rate <- ggplot(PRJEB.two_gene_plot_data, aes(x = Gene, y = mutation_percent, 
    group = response_group, fill = response_group)) + geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper, 
    color = response_group), position = PRJEB.response_dodge, width = 0.12, linewidth = 0.9, show.legend = FALSE) + 
    geom_point(position = PRJEB.response_dodge, shape = 21, size = 5, color = "black", stroke = 0.6) + 
    geom_text(aes(y = CI_upper + ifelse(response_group == "Non-responder", 8, 5), label = label, color = response_group), 
        position = PRJEB.two_gene_dodge, size = 3.5, fontface = "bold", show.legend = FALSE) + scale_fill_manual(name = "Response", 
    values = c(`Non-responder` = "#377EB8", Responder = "#D73027")) + scale_color_manual(values = c(`Non-responder` = "#377EB8", 
    Responder = "#D73027"), guide = "none") + scale_y_continuous(breaks = seq(0, 100, 20), labels = function(x) paste0(x, 
    "%"), expand = expansion(mult = 0)) + coord_cartesian(ylim = c(-3, 113)) + scale_x_discrete(expand = expansion(add = 0.7)) + 
    labs(title = "CTNNB1 and TP53 mutation frequency", subtitle = paste0("PRJEB34724 | Two-sided Fisher test", 
        PRJEB.two_gene_P_label), x = NULL, y = "Mutation frequency") + theme_classic(base_size = 15) + 
    theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11), 
        axis.text.x = element_text(face = "italic", color = "black", size = 15), axis.text.y = element_text(color = "black", 
            size = 13), axis.title.y = element_text(size = 15), legend.title = element_text(size = 12), 
        legend.text = element_text(size = 12), plot.caption = element_text(hjust = 0, size = 9))

ggsave(filename = paste0("ICB.model/BA4_visualization/", "CTNNB1_TP53_mutation_rate_by_gene.pdf"), plot = PRJEB.p_two_gene_mutation_rate, 
    width = 7, height = 4.5, bg = "white")

library(dplyr)

library(ggplot2)

ICGC.effect_patients <- ICGC.BA4_features %>% transmute(Patient_ID = as.character(Patient_ID), BA4_ICBScore = BA4_ICBScore)

stopifnot(nrow(ICGC.effect_patients) > 0L, !anyNA(ICGC.effect_patients$Patient_ID), !anyDuplicated(ICGC.effect_patients$Patient_ID), 
    is.numeric(ICGC.effect_patients$BA4_ICBScore), all(is.finite(ICGC.effect_patients$BA4_ICBScore)), 
    all(ICGC.effect_patients$Patient_ID %in% colnames(ICGC.mut)), !is.null(rownames(ICGC.mut)), !anyDuplicated(rownames(ICGC.mut)), 
    !anyDuplicated(colnames(ICGC.mut)))

ICGC.effect_mut <- ICGC.mut[, ICGC.effect_patients$Patient_ID, drop = FALSE]

stopifnot(identical(colnames(ICGC.effect_mut), ICGC.effect_patients$Patient_ID), !anyNA(ICGC.effect_mut), 
    all(ICGC.effect_mut %in% c(0, 1)))

ICGC.effect_score <- ICGC.effect_patients$BA4_ICBScore

ICGC.effect_frequency <- data.frame(Gene = rownames(ICGC.effect_mut), n_MUT = as.integer(rowSums(ICGC.effect_mut == 
    1)), stringsAsFactors = FALSE) %>% mutate(n_WT = ncol(ICGC.effect_mut) - n_MUT, mutation_frequency = n_MUT/ncol(ICGC.effect_mut)) %>% 
    filter(mutation_frequency > 0.05) %>% arrange(Gene)

if (nrow(ICGC.effect_frequency) == 0L) {
    stop("No genes have mutation frequency above 5% in this analysis population.")
}

ICGC.rank_biserial <- function(x, y) {
    nx <- length(x)
    ny <- length(y)
    if (nx == 0L || ny == 0L) 
        return(NA_real_)
    ranks <- rank(c(x, y), ties.method = "average")
    U_mut <- sum(ranks[seq_len(nx)]) - nx * (nx + 1)/2
    2 * U_mut/(nx * ny) - 1
}

ICGC.compare_gene_BA4 <- function(gene, B = 2000L) {
    is_mut <- as.vector(ICGC.effect_mut[gene, ] == 1)
    x <- ICGC.effect_score[is_mut]
    y <- ICGC.effect_score[!is_mut]
    nx <- length(x)
    ny <- length(y)
    testable <- nx > 0L && ny > 0L
    effect <- p_value <- NA_real_
    ci <- c(NA_real_, NA_real_)
    if (testable) {
        effect <- ICGC.rank_biserial(x, y)
        p_value <- if (length(unique(c(x, y))) == 1L) {
            1
        }
        else {
            wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
        }
        if (nx >= 2L && ny >= 2L) {
            boot_effect <- replicate(B, {
                xb <- x[sample.int(nx, nx, replace = TRUE)]
                yb <- y[sample.int(ny, ny, replace = TRUE)]
                ICGC.rank_biserial(xb, yb)
            })
            ci <- unname(quantile(boot_effect, probs = c(0.025, 0.975), names = FALSE))
        }
    }
    data.frame(Gene = gene, n_MUT = nx, n_WT = ny, mutation_frequency = nx/(nx + ny), mean_MUT = if (nx > 
        0L) 
        mean(x)
    else NA_real_, mean_WT = if (ny > 0L) 
        mean(y)
    else NA_real_, median_MUT = if (nx > 0L) 
        median(x)
    else NA_real_, median_WT = if (ny > 0L) 
        median(y)
    else NA_real_, delta_mean = if (testable) 
        mean(x) - mean(y)
    else NA_real_, delta_median = if (testable) 
        median(x) - median(y)
    else NA_real_, rank_biserial = effect, CI_lower = ci[1], CI_upper = ci[2], P_value = p_value, testable = testable)
}

set.seed(2026)

ICGC.gene_BA4_effects <- bind_rows(lapply(ICGC.effect_frequency$Gene, ICGC.compare_gene_BA4, B = 2000L))

ICGC.gene_BA4_effects$FDR <- NA_real_

ICGC.effect_tested <- which(ICGC.gene_BA4_effects$testable)

ICGC.gene_BA4_effects$FDR[ICGC.effect_tested] <- p.adjust(ICGC.gene_BA4_effects$P_value[ICGC.effect_tested], 
    method = "BH")

ICGC.gene_BA4_effects <- ICGC.gene_BA4_effects %>% mutate(Direction = case_when(!testable ~ "Not testable", 
    rank_biserial > 0 ~ "Higher in MUT", rank_biserial < 0 ~ "Lower in MUT", TRUE ~ "No rank difference")) %>% 
    arrange(FDR, desc(abs(rank_biserial)))

ICGC.effect_plot_data <- ICGC.gene_BA4_effects %>% filter(testable) %>% arrange(rank_biserial) %>% mutate(Gene_label = sprintf("%s (MUT=%d, WT=%d)", 
    Gene, n_MUT, n_WT), Gene_label = factor(Gene_label, levels = Gene_label), Significance = if_else(FDR < 
    0.05, "FDR < 0.05", "FDR >= 0.05"))

ICGC.effect_plot_data = ICGC.effect_plot_data[ICGC.effect_plot_data$P_value < 0.05, ]

ICGC.p_gene_BA4_effect <- ggplot(ICGC.effect_plot_data, aes(x = rank_biserial, y = Gene_label)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey60") + geom_segment(aes(x = CI_lower, xend = CI_upper, yend = Gene_label), 
    color = "grey55", linewidth = 0.6, na.rm = TRUE) + geom_point(aes(fill = Significance, size = -log10(pmax(P_value, 
    .Machine$double.xmin))), shape = 21, color = "black", stroke = 0.6) + scale_size_continuous(name = expression(-log[10](P_value)), 
    range = c(2, 7)) + scale_fill_manual(name = NULL, values = c(`FDR < 0.05` = "#D73027", `FDR >= 0.05` = "#2878B5")) + 
    scale_x_continuous(limits = c(-1, 1)) + labs(x = "Rank-biserial correlation (MUT vs WT)", y = NULL, 
    color = NULL, title = "Mutation-associated differences in ICI Score", subtitle = "ICGC | Mutation frequency > 5% | Bootstrap 95% CI", 
    caption = "Positive values indicate higher ICI Score in MUT.") + theme_classic(base_size = 12) + 
    theme(legend.position = "top")

ICGC.p_gene_BA4_effect <- ICGC.p_gene_BA4_effect + scale_y_discrete(labels = setNames(ICGC.effect_plot_data$Gene, 
    ICGC.effect_plot_data$Gene_label))

ICGC.p_gene_BA4_effect <- ICGC.p_gene_BA4_effect + theme(axis.text.x = element_text(size = 12, color = "black"), 
    axis.text.y = element_text(size = 12, color = "black", face = "italic"), axis.title.x = element_text(size = 14), 
    axis.title.y = element_text(size = 14))

ICGC.effect_dir <- file.path("ICB.model", "ICGC_visualization")

dir.create(ICGC.effect_dir, recursive = TRUE, showWarnings = FALSE)

write.table(ICGC.gene_BA4_effects, file.path(ICGC.effect_dir, "ICGC_gene_BA4_effects.tsv"), sep = "\t", 
    quote = FALSE, row.names = FALSE, na = "NA")

ggsave(filename = file.path(ICGC.effect_dir, "ICGC_gene_BA4_effect_forest.pdf"), plot = ICGC.p_gene_BA4_effect, 
    width = 9, height = max(4, 0.28 * nrow(ICGC.effect_plot_data) + 2), limitsize = FALSE)

library(dplyr)

library(ggplot2)

TCGA.effect_patients <- TCGA.BA4_features %>% transmute(Patient_ID = as.character(Patient_ID), BA4_ICBScore = BA4_ICBScore)

stopifnot(nrow(TCGA.effect_patients) > 0L, !anyNA(TCGA.effect_patients$Patient_ID), !anyDuplicated(TCGA.effect_patients$Patient_ID), 
    is.numeric(TCGA.effect_patients$BA4_ICBScore), all(is.finite(TCGA.effect_patients$BA4_ICBScore)), 
    all(TCGA.effect_patients$Patient_ID %in% colnames(TCGA.mut)), !is.null(rownames(TCGA.mut)), !anyDuplicated(rownames(TCGA.mut)), 
    !anyDuplicated(colnames(TCGA.mut)))

TCGA.effect_mut <- TCGA.mut[, TCGA.effect_patients$Patient_ID, drop = FALSE]

stopifnot(identical(colnames(TCGA.effect_mut), TCGA.effect_patients$Patient_ID), !anyNA(TCGA.effect_mut), 
    all(TCGA.effect_mut %in% c(0, 1)))

TCGA.effect_score <- TCGA.effect_patients$BA4_ICBScore

TCGA.effect_frequency <- data.frame(Gene = rownames(TCGA.effect_mut), n_MUT = as.integer(rowSums(TCGA.effect_mut == 
    1)), stringsAsFactors = FALSE) %>% mutate(n_WT = ncol(TCGA.effect_mut) - n_MUT, mutation_frequency = n_MUT/ncol(TCGA.effect_mut)) %>% 
    filter(mutation_frequency > 0.05) %>% arrange(Gene)

if (nrow(TCGA.effect_frequency) == 0L) {
    stop("No genes have mutation frequency above 5% in this analysis population.")
}

TCGA.rank_biserial <- function(x, y) {
    nx <- length(x)
    ny <- length(y)
    if (nx == 0L || ny == 0L) 
        return(NA_real_)
    pooled_ranks <- rank(c(x, y), ties.method = "average")
    U_mut <- sum(pooled_ranks[seq_len(nx)]) - nx * (nx + 1)/2
    2 * U_mut/(nx * ny) - 1
}

TCGA.compare_gene_BA4 <- function(gene, B = 2000L) {
    is_mut <- as.vector(TCGA.effect_mut[gene, ] == 1)
    x <- TCGA.effect_score[is_mut]
    y <- TCGA.effect_score[!is_mut]
    nx <- length(x)
    ny <- length(y)
    testable <- nx > 0L && ny > 0L
    effect <- NA_real_
    p_value <- NA_real_
    ci <- c(NA_real_, NA_real_)
    if (testable) {
        effect <- TCGA.rank_biserial(x, y)
        p_value <- if (length(unique(c(x, y))) == 1L) {
            1
        }
        else {
            wilcox.test(x, y, alternative = "two.sided", exact = FALSE, correct = TRUE)$p.value
        }
        if (nx >= 2L && ny >= 2L) {
            boot_effect <- replicate(B, {
                xb <- x[sample.int(nx, nx, replace = TRUE)]
                yb <- y[sample.int(ny, ny, replace = TRUE)]
                TCGA.rank_biserial(xb, yb)
            })
            ci <- unname(quantile(boot_effect, probs = c(0.025, 0.975), names = FALSE))
        }
    }
    data.frame(Gene = gene, n_MUT = nx, n_WT = ny, mutation_frequency = nx/(nx + ny), mean_MUT = if (nx > 
        0L) 
        mean(x)
    else NA_real_, mean_WT = if (ny > 0L) 
        mean(y)
    else NA_real_, median_MUT = if (nx > 0L) 
        median(x)
    else NA_real_, median_WT = if (ny > 0L) 
        median(y)
    else NA_real_, delta_mean = if (testable) 
        mean(x) - mean(y)
    else NA_real_, delta_median = if (testable) 
        median(x) - median(y)
    else NA_real_, rank_biserial = effect, CI_lower = ci[1], CI_upper = ci[2], P_value = p_value, testable = testable)
}

set.seed(2026)

TCGA.gene_BA4_effects <- bind_rows(lapply(TCGA.effect_frequency$Gene, TCGA.compare_gene_BA4, B = 2000L))

TCGA.gene_BA4_effects$FDR <- NA_real_

TCGA.effect_tested <- which(TCGA.gene_BA4_effects$testable)

TCGA.gene_BA4_effects$FDR[TCGA.effect_tested] <- p.adjust(TCGA.gene_BA4_effects$P_value[TCGA.effect_tested], 
    method = "BH")

TCGA.gene_BA4_effects <- TCGA.gene_BA4_effects %>% mutate(Direction = case_when(!testable ~ "Not testable", 
    rank_biserial > 0 ~ "Higher in MUT", rank_biserial < 0 ~ "Lower in MUT", TRUE ~ "No rank difference")) %>% 
    arrange(FDR, desc(abs(rank_biserial)))

TCGA.effect_plot_data <- TCGA.gene_BA4_effects %>% filter(testable) %>% arrange(rank_biserial) %>% mutate(Gene_label = sprintf("%s (MUT=%d, WT=%d)", 
    Gene, n_MUT, n_WT), Gene_label = factor(Gene_label, levels = Gene_label), Significance = if_else(FDR < 
    0.05, "FDR < 0.05", "FDR >= 0.05"))

TCGA.effect_plot_data = TCGA.effect_plot_data[TCGA.effect_plot_data$P_value < 0.05, ]

TCGA.p_gene_BA4_effect <- ggplot(TCGA.effect_plot_data, aes(x = rank_biserial, y = Gene_label)) + geom_vline(xintercept = 0, 
    linetype = "dashed", color = "grey60") + geom_segment(aes(x = CI_lower, xend = CI_upper, yend = Gene_label), 
    color = "grey55", linewidth = 0.6, na.rm = TRUE) + geom_point(aes(fill = Significance, size = -log10(pmax(P_value, 
    .Machine$double.xmin))), shape = 21, color = "black", stroke = 0.6) + scale_size_continuous(name = expression(-log[10](P_value)), 
    range = c(2, 7)) + scale_fill_manual(name = NULL, values = c(`FDR < 0.05` = "#D73027", `FDR >= 0.05` = "#2878B5")) + 
    scale_x_continuous(limits = c(-1, 1)) + labs(x = "Rank-biserial correlation (MUT vs WT)", y = NULL, 
    color = NULL, title = "Mutation-associated differences in ICI Score", subtitle = "TCGA-LIHC | Mutation frequency > 5% | Bootstrap 95% CI", 
    caption = "Positive values indicate higher ICI Score in MUT.") + theme_classic(base_size = 8) + theme(legend.position = "top")

TCGA.p_gene_BA4_effect <- TCGA.p_gene_BA4_effect + scale_y_discrete(labels = setNames(TCGA.effect_plot_data$Gene, 
    TCGA.effect_plot_data$Gene_label))


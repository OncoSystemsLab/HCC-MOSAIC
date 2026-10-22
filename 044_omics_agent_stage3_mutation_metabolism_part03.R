if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

write.csv(adjusted_cross, file.path(outdir, "Stage3_2C_cross_cohort_adjusted.csv"), row.names = FALSE)

write.csv(focus_adjusted, file.path(outdir, "Stage3_2C_adjusted_focus_genes.csv"), row.names = FALSE)

stopifnot(nrow(adjusted_cross) == 63L, all(is.finite(adjusted_cross$Beta_Adjusted_TCGA)), all(is.finite(adjusted_cross$Beta_Adjusted_ICGC)))

library(ggplot2)

PROJECT_DIR <- "HCC_OmicsAgent"

input_file <- file.path(PROJECT_DIR, "05_results", "Stage3_2C_cross_cohort_adjusted.csv")

stopifnot(file.exists(input_file))

adjusted_cross <- read.csv(input_file, stringsAsFactors = FALSE, check.names = FALSE)

gene_modules <- list(`Wnt Signaling` = c("NOTUM", "AXIN2", "GLUL", "TBX3"), `Bile Acid / Hepatic Metabolism` = c("CYP7A1", 
    "CYP8B1", "BAAT", "SLC27A5", "SLC10A1", "ABCB11", "ABCC2", "ABCC4", "SLCO1B3", "CYP27A1", "NR1H4", 
    "NR0B2"), `Immune Ecosystem` = c("CD8A", "GZMB", "IFNG", "NKG7", "HAVCR2"))

plot_genes <- unlist(gene_modules, use.names = FALSE)

stopifnot(!anyDuplicated(plot_genes), all(plot_genes %in% adjusted_cross$Gene))

gene_module_map <- data.frame(Gene = plot_genes, Module = rep(names(gene_modules), lengths(gene_modules)), 
    stringsAsFactors = FALSE)

plot_list <- lapply(c("TCGA", "ICGC"), function(cohort) {
    data.frame(Gene = adjusted_cross$Gene, Cohort = cohort, Beta = adjusted_cross[[paste0("Beta_Adjusted_", 
        cohort)]], CI_Lower = adjusted_cross[[paste0("CI_Lower_", cohort)]], CI_Upper = adjusted_cross[[paste0("CI_Upper_", 
        cohort)]], FDR = adjusted_cross[[paste0("FDR_63_", cohort)]], stringsAsFactors = FALSE)
})

plot_df <- do.call(rbind, plot_list)

plot_df <- merge(plot_df, gene_module_map, by = "Gene")

plot_df$Module <- factor(plot_df$Module, levels = names(gene_modules))

plot_df$Cohort <- factor(plot_df$Cohort, levels = c("TCGA", "ICGC"))

plot_df$Gene <- factor(plot_df$Gene, levels = rev(plot_genes))

plot_df$Significance <- ifelse(plot_df$FDR < 0.05, "FDR < 0.05", "Not significant")

plot_df$Significance <- factor(plot_df$Significance, levels = c("FDR < 0.05", "Not significant"))

stopifnot(nrow(plot_df) == length(plot_genes) * 2L, all(is.finite(plot_df$Beta)), all(is.finite(plot_df$CI_Lower)), 
    all(is.finite(plot_df$CI_Upper)))

p <- ggplot(plot_df, aes(y = Gene)) + geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, 
    color = "#8693A1") + geom_segment(aes(x = CI_Lower, xend = CI_Upper, yend = Gene, color = Significance), 
    linewidth = 0.9) + geom_point(aes(x = Beta, color = Significance), size = 2.5) + facet_grid(Module ~ 
    Cohort, scales = "free_y", space = "free_y") + scale_color_manual(values = c(`FDR < 0.05` = "#176B91", 
    `Not significant` = "#A7B0BA"), drop = FALSE) + scale_y_discrete(drop = TRUE) + scale_x_continuous(expand = expansion(mult = c(0.04, 
    0.07))) + labs(title = paste("CTNNB1 mutation-associated", "transcriptional changes"), subtitle = paste("Age- and sex-adjusted linear models", 
    "| TCGA-LIHC (n=361)", "| ICGC LICA-FR (n=96)"), x = "Adjusted coefficient (log2[TPM + 1])", y = NULL, 
    color = NULL, caption = paste("Points: adjusted coefficients;", "horizontal lines: 95% HC3 robust confidence intervals.", 
        "\nFDR: Benjamini-Hochberg correction", "across 63 genes within each cohort.", "\nAssociations do not establish causality", 
        "or bile acid concentrations.")) + theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(), strip.background = element_rect(fill = "#EDF2F7", color = NA), 
    strip.text = element_text(face = "bold", size = 10, color = "#24364B"), axis.text.y = element_text(face = "bold", 
        color = "#253245", size = 9), legend.position = "bottom", legend.text = element_text(size = 10), 
    plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 10, margin = margin(b = 13)), 
    plot.caption = element_text(hjust = 0, size = 8.5, lineheight = 1.25, margin = margin(t = 12)), panel.spacing.y = grid::unit(1.1, 
        "lines"), plot.margin = margin(15, 25, 15, 15))

figure_dir <- file.path(PROJECT_DIR, "06_figures")

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(figure_dir, "Stage3_2C_CTNNB1_adjusted_forest.png")

pdf_file <- file.path(figure_dir, "Stage3_2C_CTNNB1_adjusted_forest.pdf")

ggsave(filename = png_file, plot = p, width = 12.5, height = 11, dpi = 300, bg = "white")

ggsave(filename = pdf_file, plot = p, width = 12.5, height = 11, device = "pdf")

write.csv(plot_df, file.path(PROJECT_DIR, "05_results", "Stage3_2C_forest_plot_data.csv"), row.names = FALSE)

stopifnot(file.exists(png_file), file.exists(pdf_file), nrow(plot_df) == 42L)

library(ggplot2)

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("p"), inherits(p, "ggplot"))

p_fixed <- p + facet_grid(Module ~ Cohort, scales = "free_y", space = "free_y", labeller = labeller(Module = c(`Wnt Signaling` = "Wnt", 
    `Bile Acid / Hepatic Metabolism` = "Bile acid / liver", `Immune Ecosystem` = "Immune"))) + theme(strip.text.y.right = element_text(angle = 270, 
    face = "bold", size = 10))

ggsave(file.path(PROJECT_DIR, "06_figures", "Stage3_2C_CTNNB1_adjusted_forest_v1_1.png"), p_fixed, width = 12.5, 
    height = 11, dpi = 300, bg = "white")

PROJECT_DIR <- "HCC_OmicsAgent"

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

clinical_fields <- list(TCGA_LIHC = list(Age = "age_at_index", Categorical = c("sex_at_birth", "tumor_grade", 
    "ajcc_pathologic_stage", "ishak_fibrosis_score")), ICGC_LICA_FR = list(Age = "Age", Categorical = c("Gender", 
    "Edmonson_grade", "Hepatitis_B", "Hepatitis_C", "Alcohol_intake", "Fibrosis_stage", "Vascular_invasion")))

clean_clinical <- function(x) {
    x <- trimws(as.character(x))
    missing_codes <- c("", "NA", "N/A", "UNKNOWN", "NOT REPORTED", "NOT AVAILABLE", "[NOT REPORTED]", 
        "[NOT AVAILABLE]", "--")
    x[is.na(x) | toupper(x) %in% missing_codes] <- "NA"
    x
}

age_results <- list()

category_results <- list()

for (cohort in names(joint_data)) {
    joint <- joint_data[[cohort]]
    clinical <- joint$clinical
    ids <- joint$patient_ids
    fields <- clinical_fields[[cohort]]
    stopifnot(identical(clinical$Patient_ID, ids), fields$Age %in% names(clinical), all(fields$Categorical %in% 
        names(clinical)))
    mutation <- as.integer(joint$wnt_mutation["CTNNB1", ids])
    group <- ifelse(mutation == 1L, "CTNNB1_Mut", "No_Recorded_Mutation")
    group <- factor(group, levels = c("No_Recorded_Mutation", "CTNNB1_Mut"))
    cat("\n====================================")
    cat("COHORT:", cohort, "")
    cat("====================================")
    print(table(group))
    age <- suppressWarnings(as.numeric(clean_clinical(clinical[[fields$Age]])))
    for (g in levels(group)) {
        v <- age[group == g]
        v <- v[is.finite(v)]
        one <- data.frame(Cohort = cohort, Group = g, N_Available = length(v), Age_Median = median(v), 
            Age_Q1 = as.numeric(quantile(v, 0.25)), Age_Q3 = as.numeric(quantile(v, 0.75)))
        age_results[[length(age_results) + 1L]] <- one
    }
    cat("\nAGE BY MUTATION GROUP:")
    print(do.call(rbind, tail(age_results, 2)), row.names = FALSE)
    for (variable in fields$Categorical) {
        values <- clean_clinical(clinical[[variable]])
        values[is.na(values)] <- "<MISSING>"
        cat("\nVARIABLE:", variable, "")
        tab <- table(Value = values, Group = group)
        print(tab)
        for (g in levels(group)) {
            group_n <- sum(group == g)
            one <- data.frame(Cohort = cohort, Variable = variable, Category = rownames(tab), Group = g, 
                N = as.integer(tab[, g]), Group_Total = group_n, Percentage = round(100 * as.integer(tab[, 
                  g])/group_n, 2), stringsAsFactors = FALSE)
            category_results[[length(category_results) + 1L]] <- one
        }
    }
}

age_balance <- do.call(rbind, age_results)

category_balance <- do.call(rbind, category_results)

rownames(age_balance) <- NULL

rownames(category_balance) <- NULL

report_dir <- file.path(PROJECT_DIR, "08_reports")

dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(age_balance, file.path(report_dir, "Stage3_2D_age_balance.csv"), row.names = FALSE)

write.csv(category_balance, file.path(report_dir, "Stage3_2D_clinical_category_balance.csv"), row.names = FALSE)

stopifnot(nrow(age_balance) == 4L, nrow(category_balance) > 0L)

PROJECT_DIR <- "HCC_OmicsAgent"

joint_data <- list(TCGA_LIHC = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_TCGA_LIHC_joint_v1.rds")), 
    ICGC_LICA_FR = readRDS(file.path(PROJECT_DIR, "04_processed", "Stage3_2A_ICGC_LICA_FR_joint_v1.rds")))

clean_value <- function(x) {
    x <- trimws(as.character(x))
    x[toupper(x) %in% c("", "NA", "N/A", "UNKNOWN", "NOT REPORTED", "NOT AVAILABLE", "[NOT REPORTED]", 
        "[NOT AVAILABLE]")] <- "NA"
    x
}

yes_no_factor <- function(x) {
    x <- tolower(clean_value(x))
    z <- ifelse(x == "yes", "Yes", ifelse(x == "no", "No", "NA"))
    factor(z, levels = c("No", "Yes"))
}

make_clinical_dataset <- function(joint) {
    cl <- joint$clinical
    ids <- joint$patient_ids
    stopifnot(identical(cl$Patient_ID, ids))
    mutation <- as.integer(joint$wnt_mutation["CTNNB1", ids])
    if (joint$cohort == "TCGA_LIHC") {
        age <- as.numeric(cl$age_at_index)
        sex_raw <- clean_value(cl$sex_at_birth)
    }
    else {
        age <- as.numeric(cl$Age)
        sex_raw <- clean_value(cl$Gender)
    }
    sex_raw <- toupper(sex_raw)
    sex <- factor(ifelse(sex_raw %in% c("M", "MALE"), "Male", ifelse(sex_raw %in% c("F", "FEMALE"), "Female", 
        "NA")), levels = c("Female", "Male"))
    dat <- data.frame(Patient_ID = ids, Mutation = mutation, Age = age, Sex = sex)
    if (joint$cohort == "TCGA_LIHC") {
        grade <- clean_value(cl$tumor_grade)
        dat$Grade <- factor(ifelse(grade %in% c("G1", "G2"), "Low", ifelse(grade %in% c("G3", "G4"), 
            "High", "NA")), levels = c("Low", "High"))
        stage <- clean_value(cl$ajcc_pathologic_stage)
        early <- c("Stage I", "Stage II")
        advanced <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV", "Stage IVA", 
            "Stage IVB")
        dat$Stage <- factor(ifelse(stage %in% early, "Early", ifelse(stage %in% advanced, "Advanced", 
            "NA")), levels = c("Early", "Advanced"))
    }
    else {
        grade <- clean_value(cl$Edmonson_grade)
        dat$Grade <- factor(ifelse(grade == "I-II", "Low", ifelse(grade == "III-IV", "High", "NA")), 
            levels = c("Low", "High"))
        dat$Alcohol <- yes_no_factor(cl$Alcohol_intake)
        dat$HBV <- yes_no_factor(cl$Hepatitis_B)
        dat$HCV <- yes_no_factor(cl$Hepatitis_C)
    }
    stopifnot(all(is.finite(dat$Age)), all(dat$Age > 0 & dat$Age < 120), !anyNA(dat$Sex), all(dat$Mutation %in% 
        c(0, 1)))
    dat
}

model_specs <- list(TCGA_LIHC = list(M0_AgeSex = Y ~ Mutation + Age + Sex, M1_Grade = Y ~ Mutation + 
    Age + Sex + Grade, M2_Stage = Y ~ Mutation + Age + Sex + Stage, M3_GradeStage = Y ~ Mutation + Age + 
    Sex + Grade + Stage), ICGC_LICA_FR = list(M0_AgeSex = Y ~ Mutation + Age + Sex, M1_Etiology = Y ~ 
    Mutation + Age + Sex + Alcohol + HBV + HCV, M2_Grade = Y ~ Mutation + Age + Sex + Grade, M3_EtiologyGrade = Y ~ 
    Mutation + Age + Sex + Alcohol + HBV + HCV + Grade))

base_formula <- Y ~ Mutation + Age + Sex

fit_hc3 <- function(gene, expr, covariates, formula) {
    df <- covariates
    df$Y <- as.numeric(expr[gene, df$Patient_ID])
    fit <- lm(formula, data = df, na.action = na.fail)
    X <- model.matrix(fit)
    if (fit$rank != ncol(X)) {
        stop("Regression rank deficiency: ", gene)
    }
    inv_xtx <- solve(crossprod(X))
    h <- hatvalues(fit)
    e_hc3 <- residuals(fit)/pmax(1 - h, 1e-10)
    V <- inv_xtx %*% crossprod(X * e_hc3) %*% inv_xtx
    beta <- unname(coef(fit)["Mutation"])
    se <- sqrt(V["Mutation", "Mutation"])
    if (!is.finite(se) || se <= 0) {
        stop("Invalid robust SE: ", gene)
    }
    df_res <- df.residual(fit)
    pval <- 2 * pt(-abs(beta/se), df = df_res)
    tcrit <- qt(0.975, df_res)
    data.frame(Gene = gene, Beta = beta, SE_HC3 = se, CI_Lower = beta - tcrit * se, CI_Upper = beta + 
        tcrit * se, P_Value = pval, stringsAsFactors = FALSE)
}

all_results <- list()

model_qc <- list()

for (cohort in names(joint_data)) {
    joint <- joint_data[[cohort]]
    cov <- make_clinical_dataset(joint)
    expr <- joint$mechanism_logTPM
    stopifnot(nrow(expr) == 63L, identical(colnames(expr), cov$Patient_ID))
    for (model_name in names(model_specs[[cohort]])) {
        formula <- model_specs[[cohort]][[model_name]]
        required <- setdiff(all.vars(formula), "Y")
        valid <- complete.cases(cov[, required, drop = FALSE])
        dat <- droplevels(cov[valid, , drop = FALSE])
        n_mut <- sum(dat$Mutation == 1)
        n_ref <- sum(dat$Mutation == 0)
        stopifnot(n_mut >= 10L, n_ref >= 10L)
        res <- do.call(rbind, lapply(rownames(expr), fit_hc3, expr = expr, covariates = dat, formula = formula))
        same_n_base <- do.call(rbind, lapply(rownames(expr), fit_hc3, expr = expr, covariates = dat, 
            formula = base_formula))
        stopifnot(identical(res$Gene, same_n_base$Gene))
        res$Beta_M0_SameN <- same_n_base$Beta
        res$Magnitude_Change_Pct <- ifelse(abs(res$Beta_M0_SameN) >= 0.1, 100 * (abs(res$Beta)/abs(res$Beta_M0_SameN) - 
            1), NA_real_)
        res$FDR_63 <- p.adjust(res$P_Value, method = "BH")
        res$Cohort <- cohort
        res$Model <- model_name
        res$N <- nrow(dat)
        res$N_Mutation <- n_mut
        res$N_Reference <- n_ref
        key <- paste(cohort, model_name, sep = "__")
        all_results[[key]] <- res
        model_qc[[key]] <- data.frame(Cohort = cohort, Model = model_name, N = nrow(dat), N_Mutation = n_mut, 
            N_Reference = n_ref, N_Excluded = nrow(cov) - nrow(dat), FDR_Significant = sum(res$FDR_63 < 
                0.05))
        cat(cohort, "|", model_name, "| N =", nrow(dat), "| Mut =", n_mut, "| FDR < 0.05 =", sum(res$FDR_63 < 
            0.05), "")
    }
}

sensitivity_results <- do.call(rbind, all_results)

sensitivity_qc <- do.call(rbind, model_qc)

rownames(sensitivity_results) <- NULL

rownames(sensitivity_qc) <- NULL

previous <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage3_2C_cross_cohort_adjusted.csv"), stringsAsFactors = FALSE)

for (cohort in names(joint_data)) {
    x <- subset(sensitivity_results, Cohort == cohort & Model == "M0_AgeSex")
    col <- if (cohort == "TCGA_LIHC") {
        "Beta_Adjusted_TCGA"
    }
    else {
        "Beta_Adjusted_ICGC"
    }
    old_beta <- previous[match(x$Gene, previous$Gene), col]
    error <- max(abs(x$Beta - old_beta))
    cat("\nM0 reproduction error:", cohort, "=", error, "")
    stopifnot(is.finite(error), error < 1e-08)
}

tcga_m3 <- subset(sensitivity_results, Cohort == "TCGA_LIHC" & Model == "M3_GradeStage")

icgc_m3 <- subset(sensitivity_results, Cohort == "ICGC_LICA_FR" & Model == "M3_EtiologyGrade")

keep <- c("Gene", "Beta", "P_Value", "FDR_63", "Magnitude_Change_Pct")

cross_m3 <- merge(tcga_m3[, keep], icgc_m3[, keep], by = "Gene", suffixes = c("_TCGA", "_ICGC"))

cross_m3$Same_Direction <- with(cross_m3, sign(Beta_TCGA) == sign(Beta_ICGC) & Beta_TCGA != 0 & Beta_ICGC != 
    0)

cross_m3$Strict_Replication <- with(cross_m3, Same_Direction & FDR_63_TCGA < 0.05 & FDR_63_ICGC < 0.05)

focus_genes <- c("AXIN2", "GLUL", "CYP8B1", "BAAT", "ABCB11", "ABCC2", "ABCC4", "SLCO1B3", "NR0B2", "HAVCR2", 
    "CD8A")

focus_m3 <- cross_m3[match(focus_genes, cross_m3$Gene), , drop = FALSE]

stopifnot(!anyNA(focus_m3$Gene))

outdir <- file.path(PROJECT_DIR, "05_results")

report_dir <- file.path(PROJECT_DIR, "08_reports")

write.csv(sensitivity_results, file.path(outdir, "Stage3_2D_sensitivity_all_models.csv"), row.names = FALSE)

write.csv(cross_m3, file.path(outdir, "Stage3_2D_cross_cohort_M3.csv"), row.names = FALSE)

write.csv(focus_m3, file.path(outdir, "Stage3_2D_focus_M3.csv"), row.names = FALSE)

write.csv(sensitivity_qc, file.path(report_dir, "Stage3_2D_sensitivity_model_QC.csv"), row.names = FALSE)

stopifnot(nrow(sensitivity_results) == 63L * 8L, nrow(sensitivity_qc) == 8L, nrow(cross_m3) == 63L, all(is.finite(sensitivity_results$Beta)), 
    all(is.finite(sensitivity_results$P_Value)))

library(ggplot2)

PROJECT_DIR <- "HCC_OmicsAgent"

sensitivity <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage3_2D_sensitivity_all_models.csv"), 
    check.names = FALSE)

m0_cross <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage3_2C_cross_cohort_adjusted.csv"))

m3_cross <- read.csv(file.path(PROJECT_DIR, "05_results", "Stage3_2D_cross_cohort_M3.csv"))

m0_genes <- m0_cross$Gene[m0_cross$Strict_Replication]

m3_genes <- m3_cross$Gene[m3_cross$Strict_Replication]

stable_genes <- intersect(m0_genes, m3_genes)

gene_modules <- list(Wnt = c("AXIN2", "GLUL"), `Bile acid` = c("CYP8B1", "BAAT", "ABCB11", "ABCC2", "ABCC4", 
    "SLCO1B3", "NR0B2"), Immune = c("HAVCR2", "CD8A"))

focus_genes <- unlist(gene_modules, use.names = FALSE)

module_map <- data.frame(Gene = focus_genes, Module = rep(names(gene_modules), lengths(gene_modules)))

plot_df <- subset(sensitivity, Gene %in% focus_genes & Model %in% c("M3_GradeStage", "M3_EtiologyGrade"))

stopifnot(nrow(plot_df) == 22L, all(focus_genes %in% plot_df$Gene))

plot_df <- merge(plot_df, module_map, by = "Gene", sort = FALSE)

plot_df$Cohort <- factor(ifelse(plot_df$Cohort == "TCGA_LIHC", "TCGA", "ICGC"), levels = c("TCGA", "ICGC"))

plot_df$Module <- factor(plot_df$Module, levels = names(gene_modules))

plot_df$Gene <- factor(plot_df$Gene, levels = rev(focus_genes))

plot_df$Delta_Beta <- with(plot_df, Beta - Beta_M0_SameN)

plot_df$Direction_Preserved <- with(plot_df, sign(Beta) == sign(Beta_M0_SameN))

p_stability <- ggplot(plot_df, aes(y = Gene)) + geom_vline(xintercept = 0, color = "#9EA8B3", linetype = "dashed") + 
    geom_segment(aes(x = Beta_M0_SameN, xend = Beta, yend = Gene), color = "#A0B0BE", linewidth = 1.2) + 
    geom_point(aes(x = Beta_M0_SameN), shape = 21, fill = "white", color = "#64748B", size = 3, stroke = 1) + 
    geom_point(aes(x = Beta), shape = 16, color = "#176B91", size = 3) + facet_grid(Module ~ Cohort, 
    scales = "free_y", space = "free_y") + labs(title = "Stability of CTNNB1-associated effects", subtitle = paste("M0 vs cohort-specific M3", 
    "| Same-patient comparison"), x = "Adjusted coefficient (log2[TPM + 1])", y = NULL, caption = paste("Open circles: age/sex-adjusted M0;", 
    "filled circles: M3.", "\nM0 is refitted on the exact", "complete-case subset used by M3.", "\nTCGA M3 adds grade and stage;", 
    "ICGC M3 adds etiology and grade.")) + theme_minimal(base_size = 11) + theme(panel.grid.major.y = element_blank(), 
    panel.grid.minor = element_blank(), strip.background = element_rect(fill = "#EDF2F7", color = NA), 
    strip.text = element_text(face = "bold", size = 10), axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold", 
        size = 15), plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2), plot.margin = margin(14, 
        25, 14, 14))

write.csv(plot_df, file.path(PROJECT_DIR, "05_results", "Stage3_2D_M0_M3_same_patient_comparison.csv"), 
    row.names = FALSE)

write.csv(data.frame(Gene = stable_genes), file.path(PROJECT_DIR, "05_results", "Stage3_2D_stable_replication_genes.csv"), 
    row.names = FALSE)

ggsave(file.path(PROJECT_DIR, "06_figures", "Stage3_2D_CTNNB1_model_stability.png"), plot = p_stability, 
    width = 12, height = 8.5, dpi = 300, bg = "white")


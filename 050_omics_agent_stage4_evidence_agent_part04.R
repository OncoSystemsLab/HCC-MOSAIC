if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

export_summary <- do.call(rbind, lapply(names(exported), function(nm) {
    x <- exported[[nm]]
    data.frame(Test = nm, Export_Status = x$status, N_Files = length(x$files), Report_Path = x$report_file, 
        stringsAsFactors = FALSE)
}))

write.csv(export_summary, file.path(PROJECT_DIR, "08_reports", "Stage4_3A_report_export_summary.csv"), 
    row.names = FALSE)

stopifnot(nrow(export_summary) == 7L, all(export_summary$Export_Status == "success"))

library(ggplot2)

PROJECT_DIR <- "HCC_OmicsAgent"

previous <- readRDS(file.path(PROJECT_DIR, "08_reports", "Stage4_2B_auto_execution_results.rds"))

stopifnot(all(c("Full_RNA", "Cohort") %in% names(previous)), identical(previous$Full_RNA$status, "success"), 
    identical(previous$Cohort$status, "success"))

preview_dir <- file.path(PROJECT_DIR, "08_reports", "Stage4_3B_preview")

dir.create(preview_dir, recursive = TRUE, showWarnings = FALSE)

stage2 <- previous$Full_RNA$results$stage2

stopifnot(is.list(stage2), all(c("prediction", "report", "output_files") %in% names(stage2)))

inspect_structure <- function(x, path = "stage2", depth = 0L, max_depth = 3L) {
    shape <- if (is.null(dim(x))) {
        ""
    }
    else {
        paste(dim(x), collapse = "x")
    }
    current <- data.frame(Path = path, Class = paste(class(x), collapse = "/"), Length = length(x), Shape = shape, 
        stringsAsFactors = FALSE)
    if (depth >= max_depth || !is.list(x) || length(x) == 0L) {
        return(current)
    }
    nm <- names(x)
    if (is.null(nm)) {
        nm <- paste0("item_", seq_along(x))
    }
    nm[is.na(nm) | !nzchar(nm)] <- paste0("item_", which(is.na(nm) | !nzchar(nm)))
    children <- lapply(seq_along(x), function(i) {
        inspect_structure(x[[i]], path = paste0(path, "$", nm[i]), depth = depth + 1L, max_depth = max_depth)
    })
    do.call(rbind, c(list(current), children))
}

schema <- inspect_structure(stage2)

rownames(schema) <- NULL

write.csv(schema, file.path(preview_dir, "Stage4_3B_Stage2_schema.csv"), row.names = FALSE)

candidate_rows <- grepl("icb|score|prob|response", schema$Path, ignore.case = TRUE)

functional <- previous$Full_RNA$results$functional

stopifnot(is.data.frame(functional$scores), is.data.frame(functional$percentiles), nrow(functional$scores) == 
    1L, nrow(functional$percentiles) == 1L)

modules <- setdiff(names(functional$scores), "Patient_ID")

patient_df <- data.frame(Module = modules, Score = as.numeric(functional$scores[1, modules, drop = TRUE]), 
    Percentile = as.numeric(functional$percentiles[1, modules, drop = TRUE]), stringsAsFactors = FALSE)

stopifnot(nrow(patient_df) == 4L, all(is.finite(patient_df$Score)), all(is.finite(patient_df$Percentile)))

patient_df$Module <- factor(patient_df$Module, levels = rev(modules))

p_patient <- ggplot(patient_df, aes(x = Percentile, y = Module)) + geom_col(width = 0.58, fill = "#2877A8") + 
    geom_vline(xintercept = 50, linetype = "dashed", color = "#8C98A4") + geom_text(aes(label = sprintf("P%.1f", 
    Percentile)), hjust = -0.12, size = 3.6) + scale_x_continuous(limits = c(0, 108), breaks = seq(0, 
    100, 25), expand = c(0, 0)) + labs(title = "HCC functional expression profile", subtitle = paste0("Sample: ", 
    functional$scores$Patient_ID[1]), x = "Reference percentile", y = NULL, caption = paste("Relative to 69 development-cohort patients.", 
    "Not a clinical reference range or direct", "measurement of pathway activity.")) + theme_minimal(base_size = 12) + 
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", size = 15), plot.caption = element_text(size = 9, hjust = 0), 
        plot.margin = margin(15, 28, 15, 15))

patient_plot_file <- file.path(preview_dir, "Stage4_3B_patient_functional_profile.png")

ggsave(patient_plot_file, plot = p_patient, width = 10, height = 5, dpi = 300, bg = "white")

association <- previous$Cohort$results$association

required_fields <- c("Feature", "Beta_Adjusted", "CI_Lower", "CI_Upper", "FDR_BH")

stopifnot(is.data.frame(association), all(required_fields %in% names(association)), nrow(association) == 
    4L)

association$Feature <- factor(association$Feature, levels = rev(association$Feature))

association$FDR_Significant <- association$FDR_BH < 0.05

p_cohort <- ggplot(association, aes(x = Beta_Adjusted, y = Feature)) + geom_vline(xintercept = 0, color = "#98A5B3", 
    linetype = "dashed") + geom_segment(aes(x = CI_Lower, xend = CI_Upper, yend = Feature), linewidth = 0.9, 
    color = "#64748B") + geom_point(aes(color = FDR_Significant), size = 3.4) + scale_color_manual(values = c(`TRUE` = "#147C79", 
    `FALSE` = "#94A3B8"), labels = c(`FALSE` = "FDR >= 0.05", `TRUE` = "FDR < 0.05"), name = "Significance") + 
    labs(title = "CTNNB1-associated functional expression", subtitle = paste("TCGA-LIHC | Age/sex-adjusted", 
        "HC3 robust regression"), x = "Adjusted coefficient with 95% CI", y = NULL, caption = paste("Fixed 69-patient reference scores.", 
        "Observational association, not evidence", "of causality or ICB treatment response.")) + theme_minimal(base_size = 12) + 
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), legend.position = "bottom", 
        plot.title = element_text(face = "bold", size = 15), plot.caption = element_text(hjust = 0, size = 9), 
        plot.margin = margin(15, 20, 15, 15))

cohort_plot_file <- file.path(preview_dir, "Stage4_3B_TCGA_association_forest.png")

ggsave(cohort_plot_file, plot = p_cohort, width = 10, height = 5, dpi = 300, bg = "white")

stopifnot(file.exists(patient_plot_file), file.exists(cohort_plot_file), file.exists(file.path(preview_dir, 
    "Stage4_3B_Stage2_schema.csv")))

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "Stage4_3B2_report_acceptance.R"))

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "Stage4_4_launch_acceptance.R"))

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_RunPipeline_v1.R"))

result <- run_HCC_pipeline_v1(input = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"), project_dir = PROJECT_DIR, 
    sample_uid = "PRJEB34724__S02-T-RNA", run_label = "HCC_Patient_Analysis")

result$analysis_status

result$export_status

result$report_file


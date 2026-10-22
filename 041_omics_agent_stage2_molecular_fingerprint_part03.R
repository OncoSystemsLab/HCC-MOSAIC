if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

render_HCC_fingerprint_v1_1 <- function(fingerprint, project_dir = "HCC_OmicsAgent") {
    profile <- fingerprint$profile
    sm <- fingerprint$summary
    stopifnot(nrow(profile) == 89, nrow(sm) == 1, !anyDuplicated(profile$Feature))
    display_modules <- list(`Driver / Wnt` = c("HALLMARK_WNT_BETA_CATENIN_SIGNALING", "AXIN2", "GLUL"), 
        `Bile Acid Metabolism` = c("HALLMARK_BILE_ACID_METABOLISM", "CYP7A1", "CYP8B1", "BAAT", "SLC27A5", 
            "NR1H4"), `Immune Ecosystem` = c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "CD8A", "GZMB", "LAG3"))
    display_labels <- c(HALLMARK_WNT_BETA_CATENIN_SIGNALING = "Wnt / beta-catenin", HALLMARK_BILE_ACID_METABOLISM = "Bile acid hallmark", 
        HALLMARK_INTERFERON_GAMMA_RESPONSE = "IFN-gamma response", AXIN2 = "AXIN2", GLUL = "GLUL", CYP7A1 = "CYP7A1", 
        CYP8B1 = "CYP8B1", BAAT = "BAAT", SLC27A5 = "SLC27A5", NR1H4 = "NR1H4", CD8A = "CD8A", GZMB = "GZMB", 
        LAG3 = "LAG3")
    plot_df <- data.frame(Module = rep(names(display_modules), lengths(display_modules)), Feature = unlist(display_modules, 
        use.names = FALSE), stringsAsFactors = FALSE)
    idx <- match(plot_df$Feature, profile$Feature)
    stopifnot(!anyNA(idx))
    plot_df$Percentile <- profile$Reference_Percentile[idx]
    plot_df$Cohort_Flag <- profile$Cohort_Associated[idx]
    plot_df$Display_Label <- paste0(unname(display_labels[plot_df$Feature]), ifelse(plot_df$Cohort_Flag, 
        " *", ""))
    plot_df$Percentile_Label <- sprintf("P%.1f", plot_df$Percentile)
    plot_df$Module <- factor(plot_df$Module, levels = names(display_modules))
    plot_df$Display_Label <- factor(plot_df$Display_Label, levels = rev(unique(plot_df$Display_Label)))
    module_colors <- c(`Driver / Wnt` = "#507BA8", `Bile Acid Metabolism` = "#DE9846", `Immune Ecosystem` = "#4DA49B")
    uid <- as.character(sm$Sample_ID[1])
    subtitle_text <- sprintf("Sample: %s     |     BA4 ICBScore: %.4f     |     Reference: 69 HCC development samples", 
        uid, sm$ICBScore[1])
    p <- ggplot(plot_df, aes(x = Percentile, y = Display_Label, fill = Module)) + geom_vline(xintercept = c(20, 
        80), linetype = "dashed", color = "grey70", linewidth = 0.4) + geom_col(width = 0.72, orientation = "y") + 
        geom_text(aes(label = Percentile_Label), nudge_x = 2, hjust = 0, size = 3.5, color = "#253245") + 
        facet_grid(rows = vars(Module), scales = "free_y", space = "free_y", switch = "y") + scale_fill_manual(values = module_colors, 
        guide = "none") + scale_x_continuous(limits = c(0, 112), breaks = seq(0, 100, 20), expand = expansion(mult = c(0, 
        0))) + labs(title = "HCC Molecular Fingerprint", subtitle = subtitle_text, x = "Reference percentile", 
        y = NULL, caption = paste("* Cohort-associated feature (BH-FDR < 0.05).", "P20/P80 are exploratory display cutoffs, not clinical thresholds.", 
            "\nICBScore is not a clinically calibrated response probability.", "RNA signatures do not establish CTNNB1 mutation or bile acid concentration.")) + 
        theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 20, color = "#253245"), 
        plot.subtitle = element_text(size = 10, color = "#596579", margin = margin(b = 18)), plot.caption = element_text(size = 9, 
            color = "#687382", hjust = 0, margin = margin(t = 15)), strip.placement = "outside", strip.background = element_rect(fill = "#EDF2F7", 
            color = NA), strip.text.y.left = element_text(angle = 0, face = "bold", size = 10, color = "#253245", 
            margin = margin(8, 8, 8, 8)), panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major.x = element_line(color = "#EBEFF3", linewidth = 0.3), panel.spacing.y = grid::unit(0.5, 
            "cm"), axis.text.y = element_text(size = 10, color = "#344054"), axis.text.x = element_text(size = 10, 
            color = "#596579"), plot.margin = margin(20, 25, 20, 20))
    safe_id <- gsub("[^A-Za-z0-9_-]", "_", uid)
    figure_file <- file.path(project_dir, "06_figures", paste0("Stage2_4_", safe_id, "_fingerprint_v1_1.png"))
    dir.create(dirname(figure_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(filename = figure_file, plot = p, width = 13, height = 8.5, units = "in", dpi = 200, bg = "white")
    list(plot = p, figure_file = figure_file, plotted_features = nrow(plot_df))
}

PROJECT_DIR <- "HCC_OmicsAgent"

report_v1_1 <- render_HCC_fingerprint_v1_1(fingerprint = fingerprint, project_dir = PROJECT_DIR)

stopifnot(file.exists(report_v1_1$figure_file), report_v1_1$plotted_features == 13L)

library(ggplot2)

render_HCC_fingerprint_v1_2 <- function(fingerprint, project_dir = "HCC_OmicsAgent") {
    old_report <- render_HCC_fingerprint_v1_1(fingerprint = fingerprint, project_dir = project_dir)
    uid <- fingerprint$summary$Sample_ID[1]
    score <- fingerprint$summary$ICBScore[1]
    subtitle <- sprintf(paste0("Sample: %s    |    BA4 ICBScore: %.4f", "Exploratory reference: 69 HCC development samples"), 
        uid, score)
    notes <- paste("* = cohort-associated feature (BH-FDR < 0.05).", "P20/P80 are visual reference guides, not clinical thresholds.", 
        "ICBScore is not a clinically calibrated response probability.", "RNA does not establish CTNNB1 mutation or bile acid concentration.", 
        sep = "")
    p <- old_report$plot + facet_wrap(vars(Module), ncol = 1, scales = "free_y", strip.position = "top") + 
        labs(subtitle = subtitle, caption = notes) + theme(strip.background = element_rect(fill = "#EDF2F7", 
        color = NA), strip.text = element_text(hjust = 0, face = "bold", size = 11, color = "#253245", 
        margin = margin(7, 10, 7, 10)), plot.subtitle = element_text(size = 10, lineheight = 1.4, margin = margin(b = 14)), 
        plot.caption = element_text(size = 9, hjust = 0, lineheight = 1.3, margin = margin(t = 16)), 
        plot.margin = margin(16, 24, 16, 16))
    safe_id <- gsub("[^A-Za-z0-9_-]", "_", uid)
    output_file <- file.path(project_dir, "06_figures", paste0("Stage2_4_", safe_id, "_fingerprint_v1_2.png"))
    ggsave(filename = output_file, plot = p, width = 12.5, height = 9.5, units = "in", dpi = 200, bg = "white")
    stopifnot(file.exists(output_file))
    list(plot = p, figure_file = output_file, plotted_features = old_report$plotted_features)
}

report_v1_2 <- render_HCC_fingerprint_v1_2(fingerprint)

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(exists("predict_HCC_fingerprint_v1"), exists("render_HCC_fingerprint_v1_1"), exists("render_HCC_fingerprint_v1_2"))

read_HCC_TPM_single <- function(input_file) {
    if (!file.exists(input_file)) {
        stop("File does not exist: ", input_file)
    }
    ext <- tolower(tools::file_ext(input_file))
    sep <- switch(ext, csv = ",", tsv = "\t", txt = "\t", stop("Only CSV, TSV or TXT files are supported"))
    df <- read.table(input_file, header = TRUE, sep = sep, check.names = FALSE, stringsAsFactors = FALSE, 
        fileEncoding = "UTF-8-BOM", comment.char = "")
    if (ncol(df) != 2L) {
        stop("Two columns are required: GeneSymbol and patient TPM")
    }
    genes <- trimws(as.character(df[[1]]))
    values <- suppressWarnings(as.numeric(as.character(df[[2]])))
    if (anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
        stop("Gene names contain missing, blank or duplicate values")
    }
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
        stop("TPM contains invalid values")
    }
    matrix(values, ncol = 1L, dimnames = list(genes, colnames(df)[2]))
}

run_HCC_OmicsAgent <- function(input_file, sample_uid = NULL, project_dir = PROJECT_DIR) {
    suppressPackageStartupMessages({
        library(ggplot2)
        library(glmnet)
    })
    tpm <- read_HCC_TPM_single(input_file)
    if (is.null(sample_uid)) {
        sample_uid <- colnames(tpm)[1]
    }
    cat("\n===== HCC-OmicsAgent START =====")
    cat("Sample:", sample_uid, "")
    cat("Input genes:", nrow(tpm), "")
    result <- predict_HCC_fingerprint_v1(tpm = tpm, sample_uid = sample_uid, project_dir = project_dir)
    report <- render_HCC_fingerprint_v1_2(fingerprint = result, project_dir = project_dir)
    output_dir <- file.path(project_dir, "05_results")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    safe_id <- gsub("[^A-Za-z0-9_-]", "_", sample_uid)
    summary_file <- file.path(output_dir, paste0(safe_id, "_summary.csv"))
    profile_file <- file.path(output_dir, paste0(safe_id, "_profile.csv"))
    qc_file <- file.path(output_dir, paste0(safe_id, "_input_QC.csv"))
    write.csv(result$summary, summary_file, row.names = FALSE)
    write.csv(result$profile, profile_file, row.names = FALSE)
    write.csv(result$input_qc, qc_file, row.names = FALSE)
    stopifnot(nrow(result$profile) == 89L, file.exists(report$figure_file), file.exists(summary_file), 
        file.exists(profile_file), file.exists(qc_file))
    cat("\n===== AGENT OUTPUT =====")
    print(result$summary, row.names = FALSE)
    cat("Figure:", report$figure_file, "")
    cat("Profile:", profile_file, "")
    cat("\nHCC-OmicsAgent: PASS")
    list(prediction = result, report = report, output_files = c(summary = summary_file, profile = profile_file, 
        qc = qc_file, figure = report$figure_file))
}

PROJECT_DIR <- "HCC_OmicsAgent"

tpm_reference <- readRDS(file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds"))

ba4 <- readRDS(file.path(PROJECT_DIR, "01_models", "FINAL_HCC_BA4_ICBScore_v1_bundle.rds"))

bg <- ba4$background_genes

test_expression <- tpm_reference[["PRJEB34724"]][bg, "S02-T-RNA"]

test_df <- data.frame(GeneSymbol = bg, check.names = FALSE)

test_df[["S02-T-RNA"]] <- sprintf("%.17g", as.numeric(test_expression))

test_file <- file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv")

write.csv(test_df, test_file, row.names = FALSE)

stopifnot(file.exists(test_file), nrow(test_df) == 21308L)

result <- run_HCC_OmicsAgent(input_file = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"), sample_uid = "PRJEB34724__S02-T-RNA", 
    project_dir = PROJECT_DIR)

expected_score <- 0.57168779

stopifnot(abs(result$prediction$summary$ICBScore - expected_score) < 1e-06, nrow(result$prediction$profile) == 
    89L)

agent_dir <- file.path(PROJECT_DIR, "07_agent")

dir.create(agent_dir, recursive = TRUE, showWarnings = FALSE)

function_names <- c("predict_HCC_fingerprint_v1", "render_HCC_fingerprint_v1_1", "render_HCC_fingerprint_v1_2", 
    "read_HCC_TPM_single", "run_HCC_OmicsAgent")

stopifnot(all(vapply(function_names, function(nm) {
    exists(nm, envir = .GlobalEnv, inherits = FALSE) && is.function(get(nm, envir = .GlobalEnv))
}, logical(1))))

agent_script <- file.path(agent_dir, "HCC_OmicsAgent_v1.R")

dump(list = function_names, file = agent_script, envir = .GlobalEnv)

invisible(parse(file = agent_script))

test_env <- new.env(parent = globalenv())

sys.source(agent_script, envir = test_env)

stopifnot(all(vapply(function_names, function(nm) {
    exists(nm, envir = test_env, inherits = FALSE) && is.function(get(nm, envir = test_env))
}, logical(1))))

PROJECT_DIR <- "HCC_OmicsAgent"

agent_script <- file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v1.R")

function_names <- c("predict_HCC_fingerprint_v1", "render_HCC_fingerprint_v1_1", "render_HCC_fingerprint_v1_2", 
    "read_HCC_TPM_single", "run_HCC_OmicsAgent")

stopifnot(file.exists(agent_script), all(vapply(function_names, function(nm) is.function(get(nm, envir = .GlobalEnv)), 
    logical(1))))

agent_fn <- run_HCC_OmicsAgent

agent_formals <- formals(agent_fn)

agent_formals$project_dir <- PROJECT_DIR

formals(agent_fn) <- agent_formals

run_HCC_OmicsAgent <- agent_fn

stopifnot(identical(eval(formals(run_HCC_OmicsAgent)$project_dir, envir = baseenv()), PROJECT_DIR))

dump(list = function_names, file = agent_script, envir = .GlobalEnv)

invisible(parse(file = agent_script))

PROJECT_DIR <- "HCC_OmicsAgent"

agent_script <- file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v1.R")

test_file <- file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv")

cold_script <- file.path(PROJECT_DIR, "07_agent", "Stage2_5D_cold_start_test.R")

stopifnot(file.exists(agent_script), file.exists(test_file))

quote_r <- function(x) {
    paste(deparse(x), collapse = "")
}

script_lines <- c("", paste0("source(", quote_r(agent_script), ")"), paste0("res <- run_HCC_OmicsAgent(", 
    "input_file = ", quote_r(test_file), ",", "sample_uid = 'PRJEB34724__S02-T-RNA')"), "stopifnot(", 
    "  nrow(res$prediction$profile) == 89L,", "  abs(res$prediction$summary$ICBScore -", "      0.57168779) < 1e-6,", 
    "  file.exists(res$report$figure_file),", "  all(file.exists(res$output_files))", ")", "cat('Stage 2.5D: PASS\\n')")

writeLines(script_lines, con = cold_script, useBytes = TRUE)

rscript_candidates <- c(file.path(R.home("bin"), "Rscript.exe"), file.path(R.home(), "bin", "Rscript.exe"), 
    unname(Sys.which("Rscript")))

rscript_candidates <- rscript_candidates[nzchar(rscript_candidates) & file.exists(rscript_candidates)]

if (length(rscript_candidates) == 0L) {
    stop("The Rscript executable was not found")
}

rscript_exe <- rscript_candidates[1]

log_output <- suppressWarnings(system2(command = rscript_exe, args = shQuote(cold_script), stdout = TRUE, 
    stderr = TRUE))

exit_status <- attr(log_output, "status")

if (is.null(exit_status)) {
    exit_status <- 0L
}

writeLines(log_output, file.path(PROJECT_DIR, "08_reports", "Stage2_5D_cold_start_log.txt"))

stopifnot(as.integer(exit_status) == 0L, any(grepl("Stage 2.5D: PASS", log_output, fixed = TRUE)))

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "07_agent", "HCC_OmicsAgent_v1.R"))

input_file <- file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv")

stopifnot(file.exists(input_file))

original <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(ncol(original) == 2L, nrow(original) == 21308L, "CD274" %in% original[[1]])

test_dir <- tempfile("HCC_Agent_QC_")

dir.create(test_dir, recursive = TRUE)

bad_cases <- list()

bad_cases$Missing_Core_Gene <- original[original[[1]] != "CD274", , drop = FALSE]

duplicate_df <- original

duplicate_df[2, 1] <- duplicate_df[1, 1]

bad_cases$Duplicate_Gene <- duplicate_df

negative_df <- original

negative_df[1, 2] <- -1

bad_cases$Negative_TPM <- negative_df

invalid_df <- original

invalid_df[2, 2] <- "NOT_A_NUMBER"

bad_cases$NonNumeric_TPM <- invalid_df

multi_df <- original

multi_df[["Patient_002"]] <- multi_df[[2]]

bad_cases$Multiple_Patients <- multi_df

case_paths <- vapply(names(bad_cases), function(nm) {
    path <- file.path(test_dir, paste0(nm, ".csv"))
    write.csv(bad_cases[[nm]], path, row.names = FALSE)
    path
}, character(1))

case_paths <- c(case_paths, Missing_File = file.path(test_dir, "nonexistent.csv"))

expected_errors <- c(Missing_Core_Gene = "Prediction requires the complete reference gene background", Duplicate_Gene = "Gene names contain missing, blank or duplicate values", 
    Negative_TPM = "TPM contains invalid values", NonNumeric_TPM = "TPM contains invalid values", Multiple_Patients = "Two columns are required", 
    Missing_File = "File does not exist")

stopifnot(identical(names(case_paths), names(expected_errors)))

qc_results <- lapply(names(case_paths), function(case_name) {
    path <- case_paths[[case_name]]
    error_message <- tryCatch({
        invisible(capture.output(suppressWarnings(run_HCC_OmicsAgent(input_file = path, sample_uid = "QC_Test", 
            project_dir = PROJECT_DIR))))
        "NA"
    }, error = function(e) {
        conditionMessage(e)
    })
    was_rejected <- !is.na(error_message)
    expected_message <- grepl(expected_errors[[case_name]], error_message, fixed = TRUE)
    data.frame(Test = case_name, Rejected = was_rejected, Expected_Error = expected_message, PASS = was_rejected && 
        expected_message, Error_Message = error_message, stringsAsFactors = FALSE)
})

qc_results <- do.call(rbind, qc_results)

rownames(qc_results) <- NULL

qc_file <- file.path(PROJECT_DIR, "08_reports", "Stage2_5E_negative_input_QC.csv")

write.csv(qc_results, qc_file, row.names = FALSE)

stopifnot(nrow(qc_results) == 6L, all(qc_results$PASS))


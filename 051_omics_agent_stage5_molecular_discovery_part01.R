if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1A_HCC_reference_data_audit.R")

stopifnot(file.exists(script_file))

source(script_file)

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_1A_reference_audit", "Stage5_1A_reference_audit.rds")

stopifnot(file.exists(output_file))

audit <- readRDS(output_file)

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1A2_ICGC_gene_id_diagnosis.R")

stopifnot(file.exists(script_file))

source(script_file)

diagnostic_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_1A_reference_audit", "Stage5_1A2_ICGC_identifier_diagnostics.csv")

stopifnot(file.exists(diagnostic_file))

diagnostic <- read.csv(diagnostic_file, check.names = FALSE)

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1A3_ICGC_gene_symbol_repair_validation.R")

stopifnot(file.exists(script_file))

output_file <- file.path(PROJECT_DIR, "04_processed", "Stage5_1A3_ICGC_LICA_FR_full_TPM_GeneSymbol_clean_v1.rds")

if (file.exists(output_file)) {
    stop("The output file already exists; review its version before overwriting.")
}

source(script_file)

PROJECT_DIR <- "HCC_OmicsAgent"

if (!requireNamespace("matrixStats", quietly = TRUE)) {
    stop("Install the required R packages before running this analysis.")
}

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B1_build_TCGA_global_reference.R")

stopifnot(file.exists(script_file))

output_file <- file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds")

if (file.exists(output_file)) {
    stop("The target RDS already exists; review it before overwriting.")
}

source(script_file)

stopifnot(file.exists(output_file))

ref <- readRDS(output_file)

ref$version

ref$status

dim(ref$reference_logTPM)

sum(ref$gene_statistics$Screenable_Exploratory)

stopifnot(version_pass, count_pass)

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B2_build_HCC_GeneSetRegistry_v1.R")

output_file <- file.path(PROJECT_DIR, "00_config", "HCC_GeneSetRegistry_v1.rds")

stopifnot(file.exists(script_file))

stopifnot(!file.exists(output_file))

source(script_file)

stopifnot(file.exists(output_file))

registry <- readRDS(output_file)

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(requireNamespace("matrixStats", quietly = TRUE))

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B3A_build_fixed_universe_preflight.R")

output_file <- file.path(PROJECT_DIR, "02_reference", "HCC_FixedUniverse_TCGA_ICGC_v1.rds")

input_files <- c(file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds"), file.path(PROJECT_DIR, 
    "00_config", "HCC_GeneSetRegistry_v1.rds"), file.path(PROJECT_DIR, "04_processed", "Stage5_1A3_ICGC_LICA_FR_full_TPM_GeneSymbol_clean_v1.rds"), 
    file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"))

stopifnot(file.exists(script_file))

stopifnot(all(file.exists(input_files)))

if (file.exists(output_file)) {
    message("The fixed-background file already exists; skipping reconstruction and checking it directly.")
} else {
    source(script_file)
}

stopifnot(file.exists(output_file))

PROJECT_DIR <- "HCC_OmicsAgent"

if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Install the required R packages before running this analysis.")
}

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B3B_score_calibrate_HCC_programs_v1.R")

stopifnot(file.exists(script_file))

input_files <- c(TCGA = file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds"), Registry = file.path(PROJECT_DIR, 
    "00_config", "HCC_GeneSetRegistry_v1.rds"), Fixed = file.path(PROJECT_DIR, "02_reference", "HCC_FixedUniverse_TCGA_ICGC_v1.rds"), 
    ICGC = file.path(PROJECT_DIR, "04_processed", "Stage5_1A3_ICGC_LICA_FR_full_TPM_GeneSymbol_clean_v1.rds"))

stopifnot(all(file.exists(input_files)))

output_file <- file.path(PROJECT_DIR, "02_reference", "HCC_PathwayScore_Reference_TCGA_ICGC_v1.rds")

if (file.exists(output_file)) {
    stop("The output file already exists; repeated computation and overwriting are not allowed.")
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(file.exists(output_file))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B3C_mask_shift_sensitivity_v1.R")

input_files <- c(TCGA = file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds"), ICGC = file.path(PROJECT_DIR, 
    "04_processed", "Stage5_1A3_ICGC_LICA_FR_full_TPM_GeneSymbol_clean_v1.rds"), Fixed = file.path(PROJECT_DIR, 
    "02_reference", "HCC_FixedUniverse_TCGA_ICGC_v1.rds"), ReferenceScores = file.path(PROJECT_DIR, "02_reference", 
    "HCC_PathwayScore_Reference_TCGA_ICGC_v1.rds"), PatientRNA = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"))

stopifnot(requireNamespace("Matrix", quietly = TRUE))

stopifnot(file.exists(script_file))

stopifnot(all(file.exists(input_files)))

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_1B3C_mask_shift_sensitivity", "Stage5_1B3C_reference_mask_sensitivity_v1.rds")

if (file.exists(output_file)) {
    stop("Stage 5.1B.3C output already exists; review it before overwriting.")
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(file.exists(output_file))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1B3D_HCC_program_transport_triage_v1.R")

input_files <- c(ReferenceScores = file.path(PROJECT_DIR, "02_reference", "HCC_PathwayScore_Reference_TCGA_ICGC_v1.rds"), 
    MaskSensitivity = file.path(PROJECT_DIR, "08_reports", "Stage5_1B3C_mask_shift_sensitivity", "Stage5_1B3C_reference_mask_sensitivity_v1.rds"))

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_1B3D_program_transportability", "Stage5_1B3D_program_transportability_v1.rds")

stopifnot(file.exists(script_file))

stopifnot(all(file.exists(input_files)))

if (file.exists(output_file)) {
    stop("Stage 5.1B.3D file already exists; review its version first.")
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(file.exists(output_file))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

PROJECT_DIR <- "HCC_OmicsAgent"

if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Install the required R packages before running this analysis.")
}

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_2A0_ICB69_source_and_label_audit_v1.R")

stopifnot(file.exists(script_file))

RAW_DIR <- file.path(dirname(PROJECT_DIR), "raw_data")

input_files <- c(file.path(RAW_DIR, "PRJEB34724", "PRJEB34724_clin.xlsx"), file.path(RAW_DIR, "PRJEB34724", 
    "gene_TPM.tsv"), file.path(RAW_DIR, "PRJEB34724", "gene_annotation.tsv"), file.path(RAW_DIR, "GSE181946", 
    "GSE181946.xlsx"), file.path(RAW_DIR, "GSE181946", "Human.GRCh38.p13.annot.tsv"), file.path(RAW_DIR, 
    "GSE202069", "GSE202069_clin.xlsx"), file.path(RAW_DIR, "GSE202069", "GSE202069_raw_counts_GRCh38.p13_NCBI.tsv"), 
    file.path(RAW_DIR, "GSE202069", "Human.GRCh38.p13.annot.tsv"))

if (!all(file.exists(input_files))) {
    print(input_files[!file.exists(input_files)])
    stop("Source cohort files are incomplete; review the paths first.")
}

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_2A0_ICB69_audit", "Stage5_2A0_audit_result.rds")

if (file.exists(output_file)) {
    stop("The audit RDS already exists; review the previous results first.")
}

source(script_file)

stopifnot(file.exists(output_file))

audit_icb <- readRDS(output_file)

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_2A1_reuse_and_score_ICB69_v1.R")

input_files <- c(ICI_TPM = file.path(PROJECT_DIR, "02_reference", "Stage1_training_3cohort_TPM_background.rds"), 
    Legacy57 = file.path(PROJECT_DIR, "02_reference", "Stage1_training_57features.rds"), GeneSetRegistry = file.path(PROJECT_DIR, 
        "00_config", "HCC_GeneSetRegistry_v1.rds"), LabelAudit = file.path(PROJECT_DIR, "08_reports", 
        "Stage5_2A0_ICB69_audit", "Stage5_2A0_audit_result.rds"))

stopifnot(requireNamespace("Matrix", quietly = TRUE), file.exists(script_file), all(file.exists(input_files)))

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_2A1_reused_ICB_global_scoring", "Stage5_2A1_ICB69_global_program_scores_v1.rds")

if (file.exists(output_file)) {
    stop("Stage 5.2A.1 output already exists; review it before overwriting.")
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(file.exists(output_file))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2A2_ICB69_response_program_association_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2A3_ICB69_robustness_evidence_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2A4_freeze_ICB_ResponseProgram_Atlas_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2B1_build_HCC_mutation_program_bridge_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2B2_HCC_genomewide_mutation_program_association_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2B3_HCC_mutation_program_HC3_evidence_triage_v1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2B4_freeze_HCC_DriverProgram_Atlas_v1_1.R"))

source(file.path(PROJECT_DIR, "07_agent", "Stage5_2C1_build_HCC_MolecularEvidence_Registry_v1.R"))

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(requireNamespace("Matrix", quietly = TRUE))

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1C1_TCGA_masked_LOO_method_validation_v1.R")

input_files <- c(TCGA = file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds"), Fixed = file.path(PROJECT_DIR, 
    "02_reference", "HCC_FixedUniverse_TCGA_ICGC_v1.rds"), Scores = file.path(PROJECT_DIR, "02_reference", 
    "HCC_PathwayScore_Reference_TCGA_ICGC_v1.rds"), Mask = file.path(PROJECT_DIR, "08_reports", "Stage5_1B3C_mask_shift_sensitivity", 
    "Stage5_1B3C_reference_mask_sensitivity_v1.rds"), Triage = file.path(PROJECT_DIR, "08_reports", "Stage5_1B3D_program_transportability", 
    "Stage5_1B3D_program_transportability_v1.rds"), Patient = file.path(PROJECT_DIR, "03_input", "Patient_TPM.csv"))

stopifnot(file.exists(script_file), all(file.exists(input_files)))

output_file <- file.path(PROJECT_DIR, "08_reports", "Stage5_1C1_masked_LOO_method_validation", "Stage5_1C1_TCGA_masked_LOO_method_validation_v1.rds")

if (file.exists(output_file)) {
    stop("Stage 5.1C.1 already exists; review the previous results first.")
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(file.exists(output_file))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

PROJECT_DIR <- "HCC_OmicsAgent"

script_file <- file.path(PROJECT_DIR, "07_agent", "Stage5_1C2_HCC_heldout_scanner_and_input_audit_v1.R")

input_files <- c(TCGA = file.path(PROJECT_DIR, "02_reference", "HCC_GlobalReference_TCGA_v1.rds"), Fixed = file.path(PROJECT_DIR, 
    "02_reference", "HCC_FixedUniverse_TCGA_ICGC_v1.rds"), Scores = file.path(PROJECT_DIR, "02_reference", 
    "HCC_PathwayScore_Reference_TCGA_ICGC_v1.rds"), Triage = file.path(PROJECT_DIR, "08_reports", "Stage5_1B3D_program_transportability", 
    "Stage5_1B3D_program_transportability_v1.rds"), MethodQC = file.path(PROJECT_DIR, "08_reports", "Stage5_1C1_masked_LOO_method_validation", 
    "Stage5_1C1_TCGA_masked_LOO_method_validation_v1.rds"), Patient = file.path(PROJECT_DIR, "03_input", 
    "Patient_TPM.csv"))

stopifnot(file.exists(script_file), all(file.exists(input_files)))

output_dir <- file.path(PROJECT_DIR, "08_reports", "Stage5_1C2_heldout_scanner")

output_file <- file.path(output_dir, "Stage5_1C2_heldout_core342_demo_v1.rds")

expected_files <- file.path(output_dir, c("Stage5_1C2_heldout_core342_demo_v1.rds", "Stage5_1C2_TCGA_heldout_core342_scores.csv", 
    "Stage5_1C2_TCGA_heldout_top20.csv", "Stage5_1C2_TCGA_heldout_arithmetic_QC.csv", "Stage5_1C2_S02_input_provenance_QC.csv", 
    "Stage5_1C2_source_manifest.csv", "Stage5_1C2_internal_demonstration_report.md", "Stage5_1C2_TCGA_heldout_landscape.png"))

if (any(file.exists(expected_files))) {
    stop("Stage 5.1C.2 output files already exist; review them before overwriting: ", paste(basename(expected_files[file.exists(expected_files)]), 
        collapse = ", "))
}

md5_before <- tools::md5sum(input_files)

source(script_file)

stopifnot(all(file.exists(expected_files)))

md5_after <- tools::md5sum(input_files)

stopifnot(identical(md5_before, md5_after))

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(requireNamespace("DBI", quietly = TRUE), requireNamespace("RSQLite", quietly = TRUE))

atlas_dir <- file.path(PROJECT_DIR, "10_reference_atlas", "output")

db_v1 <- file.path(atlas_dir, "hcc_reference_atlas_v1.sqlite")

db_v11 <- file.path(atlas_dir, "hcc_reference_atlas_v1_1_staging.sqlite")

api_script <- file.path(PROJECT_DIR, "07_agent", "HCC_MolecularEvidence_Query_v1.R")

stopifnot(file.exists(db_v1), file.exists(api_script), !file.exists(db_v11))

source(api_script)

registry <- load_HCC_MolecularEvidence_Registry_v1(PROJECT_DIR, verify_sources = TRUE)

index <- registry$program_index

stopifnot(identical(registry$version, "HCC_MolecularEvidence_Registry_v1"), nrow(index) == 9420L, !anyDuplicated(index$Program_ID), 
    sum(index$ICB_Association_Available) == 5273L, sum(index$Driver_Association_Available) == 3637L, 
    sum(index$ICB_Association_Available & index$Driver_Association_Available) == 3578L)

build_atlas_v11_index <- function() {
    md5_v1 <- unname(tools::md5sum(db_v1))
    md5_registry <- unname(tools::md5sum(file.path(PROJECT_DIR, "02_reference", "HCC_MolecularEvidence_Registry_v1.rds")))
    con1 <- DBI::dbConnect(RSQLite::SQLite(), db_v1)
    on.exit(if (DBI::dbIsValid(con1)) DBI::dbDisconnect(con1), add = TRUE)
    stopifnot(identical(DBI::dbGetQuery(con1, "PRAGMA quick_check")[[1]][1], "ok"))
    core_tables <- c("samples", "cohorts", "gene_reference_stats", "feature_registry", "feature_values", 
        "mutation_events", "mutation_gene_summary", "ici_response_labels")
    stopifnot(all(core_tables %in% DBI::dbListTables(con1)))
    count_rows <- function(con, table_name) {
        sql <- paste0("SELECT COUNT(*) AS n FROM ", as.character(DBI::dbQuoteIdentifier(con, table_name)))
        as.numeric(DBI::dbGetQuery(con, sql)$n[1])
    }
    counts_before <- vapply(core_tables, function(x) count_rows(con1, x), numeric(1))
    target_sql <- as.character(DBI::dbQuoteString(con1, normalizePath(db_v11, winslash = "/", mustWork = FALSE)))
    DBI::dbExecute(con1, paste("VACUUM INTO", target_sql))
    DBI::dbDisconnect(con1)
    con2 <- DBI::dbConnect(RSQLite::SQLite(), db_v11)
    on.exit(if (DBI::dbIsValid(con2)) DBI::dbDisconnect(con2), add = TRUE)
    stopifnot(!DBI::dbExistsTable(con2, "molecular_program_index"))
    DBI::dbWriteTable(con2, "molecular_program_index", index, row.names = FALSE)
    DBI::dbExecute(con2, paste("CREATE UNIQUE INDEX idx_molecular_program_id", "ON molecular_program_index (Program_ID)"))
    manifest <- data.frame(Atlas_Extension_Version = "HCC_Reference_Atlas_v1.1_STAGING", Molecular_Registry_Version = registry$version, 
        Original_Atlas_MD5 = md5_v1, Molecular_Registry_MD5 = md5_registry, N_Programs = nrow(index), 
        Patient_Abnormality_Approved = FALSE, Clinical_ICI_Validation_Approved = FALSE, Created_At = as.character(Sys.time()), 
        stringsAsFactors = FALSE)
    DBI::dbWriteTable(con2, "atlas_v1_1_manifest", manifest, row.names = FALSE)
    counts_after <- vapply(core_tables, function(x) count_rows(con2, x), numeric(1))
    stopifnot(identical(counts_before, counts_after), count_rows(con2, "molecular_program_index") == 
        9420, identical(DBI::dbGetQuery(con2, "PRAGMA quick_check")[[1]][1], "ok"))
    stopifnot(identical(unname(tools::md5sum(db_v1)), md5_v1))
    cat("\n===== ATLAS V1.1 STAGING SUMMARY =====")
    cat("New database:", db_v11, "")
    cat("Programs indexed:", nrow(index), "")
    cat("Original tables unchanged: TRUE")
    cat("Original Atlas v1 MD5 unchanged: TRUE")
    cat("Patient clinical scoring approved: FALSE")
    cat("\nAtlas v1.1 Step 01: PASS")
    invisible(data.frame(Table = core_tables, N_Rows = counts_after))
}

atlas_v11_result <- build_atlas_v11_index()

PROJECT_DIR <- "HCC_OmicsAgent"

stopifnot(requireNamespace("DBI", quietly = TRUE), requireNamespace("RSQLite", quietly = TRUE))

source(file.path(PROJECT_DIR, "10_reference_atlas", "08_extend_v11_program_reference.R"))

atlas_v11_step02 <- build_atlas_v11_program_reference(PROJECT_DIR)

source(file.path(PROJECT_DIR, "10_reference_atlas", "09_query_v11_program_reference.R"))

ex <- hcc_v11_program_profile(program_id = "Hallmark__HALLMARK_ANGIOGENESIS", cohort = "TCGA_LIHC", PROJECT_DIR = PROJECT_DIR)

stopifnot(nrow(ex$distribution) == 1L, nrow(ex$patient_reference_scores) == 361L)

ex_icgc <- hcc_v11_program_profile(program_id = "Hallmark__HALLMARK_ANGIOGENESIS", cohort = "ICGC_LICA_FR", 
    PROJECT_DIR = PROJECT_DIR)

stopifnot(nrow(ex_icgc$patient_reference_scores) == 96L)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "10_reference_atlas", "11_accept_v11_molecular_evidence.R"))

atlas_v11_step03 <- accept_atlas_v11_step03(PROJECT_DIR)

PROJECT_DIR <- "HCC_OmicsAgent"

source(file.path(PROJECT_DIR, "10_reference_atlas", "12_accept_v11_release_candidate.R"))

atlas_v11_step04 <- accept_atlas_v11_step04(PROJECT_DIR)


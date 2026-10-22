if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(dplyr)

prepare_clin <- function(clin_df, cohort_name) {
    clin_df <- as.data.frame(clin_df, check.names = FALSE)
    stopifnot(all(c("sample", "label") %in% colnames(clin_df)))
    clin_df <- clin_df %>% transmute(cohort = cohort_name, sample_id = trimws(as.character(sample)), 
        label_original = trimws(as.character(label)), label_std = tolower(trimws(as.character(label))))
    clin_df <- clin_df %>% mutate(response = case_when(label_std == "responder" ~ 1L, label_std == "non-responder" ~ 
        0L, TRUE ~ NA_integer_), patient_uid = paste(cohort, sample_id, sep = "__"))
    return(clin_df)
}

PRJEB34724.meta <- prepare_clin(PRJEB34724.clin, "PRJEB34724")

GSE181946.meta <- prepare_clin(GSE181946.clin, "GSE181946")

GSE202069.meta <- prepare_clin(GSE202069.clin, "GSE202069")

metadata_all <- bind_rows(PRJEB34724.meta, GSE181946.meta, GSE202069.meta)

dim(metadata_all)

metadata_all %>% filter(is.na(response))

metadata_all %>% count(cohort, label_original)

metadata_summary <- metadata_all %>% count(cohort, response) %>% mutate(response_label = ifelse(response == 
    1, "Responder", "Non-responder"))

metadata_summary

response_summary <- metadata_all %>% group_by(cohort) %>% summarise(n = n(), responders = sum(response == 
    1), nonresponders = sum(response == 0), response_rate = mean(response), .groups = "drop")

response_summary


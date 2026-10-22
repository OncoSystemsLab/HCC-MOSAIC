if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

write.csv(score_all, file = "HCC_3cohort_ssGSEA_57features.csv", row.names = FALSE)

saveRDS(score_list, file = "HCC_3cohort_ssGSEA_score_list.rds")

saveRDS(score_all, file = "HCC_3cohort_ssGSEA_57features.rds")


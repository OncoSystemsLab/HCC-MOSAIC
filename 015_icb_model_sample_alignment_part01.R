if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

setdiff(score_all$patient_uid, metadata_all$patient_uid)

setdiff(metadata_all$patient_uid, score_all$patient_uid)

anyDuplicated(metadata_all$patient_uid)

anyDuplicated(score_all$patient_uid)

model_data <- score_all %>% inner_join(metadata_all %>% select(patient_uid, response, label_original), 
    by = "patient_uid")

model_data <- model_data %>% arrange(cohort, sample_id)

table(model_data$cohort, model_data$response)

model_data %>% group_by(cohort) %>% summarise(n = n(), R = sum(response == 1), NR = sum(response == 0), 
    response_rate = mean(response), .groups = "drop")


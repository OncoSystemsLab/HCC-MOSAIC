if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

saveRDS(all_gene_sets, file = "HCC_ICB_57_gene_sets.rds")

write.csv(hallmark_coverage, file = "Hallmark_gene_set_coverage.csv", row.names = FALSE)

write.csv(ici_coverage, file = "ICI_gene_set_coverage.csv", row.names = FALSE)

writeLines(background_genes, con = "common_background_20268_genes.txt")


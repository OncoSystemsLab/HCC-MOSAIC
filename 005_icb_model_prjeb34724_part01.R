if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

PRJEB34724.ann <- read.table("raw_data/PRJEB34724/gene_annotation.tsv", header = TRUE, sep = "\t", quote = "", 
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

PRJEB34724.bulk <- read.table("raw_data/PRJEB34724/gene_TPM.tsv", header = TRUE, sep = "\t", quote = "", 
    comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

PRJEB34724.clin <- read_xlsx("raw_data/PRJEB34724/PRJEB34724_clin.xlsx", sheet = "Sheet1")

{
    library(dplyr)
    PRJEB34724.annotation <- PRJEB34724.ann %>% select(gene_id, gene_symbol, gene_type) %>% filter(!is.na(gene_id), 
        gene_id != "", !is.na(gene_symbol), gene_symbol != "")
    PRJEB34724.tpm.annotated <- PRJEB34724.bulk %>% left_join(PRJEB34724.annotation, by = "gene_id")
    head(PRJEB34724.tpm.annotated)
    sum(is.na(PRJEB34724.tpm.annotated$gene_symbol))
    PRJEB34724.tpm <- PRJEB34724.tpm.annotated %>% filter(!is.na(gene_symbol), gene_symbol != "") %>% 
        select(gene_symbol, where(is.numeric)) %>% group_by(gene_symbol) %>% summarise(across(everything(), 
        sum, na.rm = TRUE), .groups = "drop")
    head(PRJEB34724.tpm)
    sum(duplicated(PRJEB34724.tpm$gene_symbol))
    colSums(PRJEB34724.bulk[, -1])
    colSums(PRJEB34724.tpm[, -1])
    colnames(PRJEB34724.tpm)[1] = "Symbol"
}


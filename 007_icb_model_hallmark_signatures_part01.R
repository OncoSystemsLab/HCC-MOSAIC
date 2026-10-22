if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

library(msigdbr)

library(GSVA)

library(BiocParallel)

library(dplyr)

library(tibble)

library(ggplot2)

if ("collection" %in% names(formals(msigdbr::msigdbr))) {
    hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
} else {
    hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
}

dim(hallmark_df)

hallmark_gene_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

hallmark_gene_sets <- lapply(hallmark_gene_sets, unique)

length(hallmark_gene_sets)

names(hallmark_gene_sets)


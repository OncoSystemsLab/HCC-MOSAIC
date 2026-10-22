if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

check_geneset_coverage <- function(gene_sets, background_genes) {
    bg <- unique(background_genes)
    result <- lapply(names(gene_sets), function(gs_name) {
        genes <- unique(gene_sets[[gs_name]])
        present_genes <- intersect(genes, bg)
        missing_genes <- setdiff(genes, bg)
        data.frame(gene_set = gs_name, n_total = length(genes), n_present = length(present_genes), n_missing = length(missing_genes), 
            coverage = length(present_genes)/length(genes), missing_genes = paste(missing_genes, collapse = ";"), 
            stringsAsFactors = FALSE)
    })
    bind_rows(result)
}

hallmark_coverage <- check_geneset_coverage(hallmark_gene_sets, background_genes)

hallmark_coverage <- hallmark_coverage %>% arrange(coverage)

hallmark_coverage %>% select(gene_set, n_total, n_present, coverage, missing_genes)

hallmark_coverage %>% filter(gene_set %in% c("HALLMARK_BILE_ACID_METABOLISM", "HALLMARK_WNT_BETA_CATENIN_SIGNALING", 
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE"))

ici_coverage <- check_geneset_coverage(ici_gene_sets, background_genes)

ici_coverage <- ici_coverage %>% arrange(coverage)

ici_coverage

hallmark_low_coverage <- hallmark_coverage %>% filter(coverage < 0.8 | n_present < 10)

hallmark_low_coverage

trim_genesets_to_background <- function(gene_sets, background_genes) {
    lapply(gene_sets, function(gs) {
        intersect(unique(gs), background_genes)
    })
}

hallmark_gene_sets_bg <- trim_genesets_to_background(hallmark_gene_sets, background_genes)

ici_gene_sets_bg <- trim_genesets_to_background(ici_gene_sets, background_genes)

sapply(hallmark_gene_sets_bg, length)

sapply(ici_gene_sets_bg, length)

all_gene_sets <- c(hallmark_gene_sets_bg, ici_gene_sets_bg)

length(all_gene_sets)

anyDuplicated(names(all_gene_sets))


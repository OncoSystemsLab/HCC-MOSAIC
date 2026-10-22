if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

TCGA.CTNNB1_BA_gene_test <- bind_rows(lapply(BA_genes, test_BA_gene))

TCGA.CTNNB1_BA_gene_test <- TCGA.CTNNB1_BA_gene_test %>% mutate(FDR = p.adjust(P_value, method = "BH"), 
    direction = case_when(delta_mean > 0 ~ "Higher_in_CTNNB1_MUT", delta_mean < 0 ~ "Lower_in_CTNNB1_MUT", 
        TRUE ~ "No_change")) %>% arrange(FDR, P_value)

TCGA.BA_gene_effect <- bind_rows(lapply(BA_genes, function(gene) {
    x <- TCGA.BA_expr[gene, , drop = TRUE]
    mut <- TCGA.CTNNB1_status_vector == "MUT"
    data.frame(Gene = gene, standardized_difference = (mean(x[mut], na.rm = TRUE) - mean(x[!mut], na.rm = TRUE))/sd(x, 
        na.rm = TRUE))
}))

TCGA.CTNNB1_BA_gene_test <- TCGA.CTNNB1_BA_gene_test %>% left_join(TCGA.BA_gene_effect, by = "Gene")

TCGA.CTNNB1_BA_gene_test <- TCGA.CTNNB1_BA_gene_test %>% arrange(FDR, desc(abs(standardized_difference)))

TCGA.CTNNB1_BA_gene_test %>% head(30)

TCGA.CTNNB1_BA_significant <- TCGA.CTNNB1_BA_gene_test %>% filter(FDR < 0.05)

nrow(TCGA.CTNNB1_BA_significant)

TCGA.CTNNB1_BA_significant

TCGA.CTNNB1_BA_up <- TCGA.CTNNB1_BA_significant %>% filter(delta_mean > 0) %>% arrange(FDR, desc(standardized_difference))

TCGA.CTNNB1_BA_up

TCGA.CTNNB1_BA_down <- TCGA.CTNNB1_BA_significant %>% filter(delta_mean < 0) %>% arrange(FDR, standardized_difference)

TCGA.CTNNB1_BA_down

TCGA.BA_gene_direction_summary <- data.frame(total_BA_genes = length(BA_genes), significant = nrow(TCGA.CTNNB1_BA_significant), 
    higher_in_MUT = nrow(TCGA.CTNNB1_BA_up), lower_in_MUT = nrow(TCGA.CTNNB1_BA_down))

TCGA.BA_gene_direction_summary

BA_key_genes <- c("CYP7A1", "CYP8B1", "CYP27A1", "HSD3B7", "AKR1D1", "SLC27A5", "BAAT")

TCGA.CTNNB1_BA_gene_test %>% filter(Gene %in% BA_key_genes) %>% arrange(FDR)

TCGA.BA_volcano_data <- TCGA.CTNNB1_BA_gene_test %>% mutate(minus_log10_FDR = -log10(pmax(FDR, 1e-300)), 
    significant = FDR < 0.05, label = ifelse(Gene %in% BA_key_genes | (FDR < 0.01 & abs(standardized_difference) > 
        0.4), Gene, NA))

p_BA_volcano <- ggplot(TCGA.BA_volcano_data, aes(x = delta_mean, y = minus_log10_FDR)) + geom_point(aes(shape = significant), 
    size = 2.5, alpha = 0.8) + geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = -log10(0.05), 
    linetype = 2) + geom_text_repel(aes(label = label), na.rm = TRUE, size = 3.5, max.overlaps = Inf) + 
    theme_classic(base_size = 13) + labs(x = expression(Delta * " mean log"[2] * "(TPM+1), MUT - WT"), 
    y = expression(-log[10] * "(FDR)"), title = "CTNNB1-associated bile acid metabolic genes")

p_BA_volcano

TCGA.BAAT_data <- data.frame(Patient_ID = colnames(TCGA.expr), BAAT = TCGA.log2_tpm["BAAT", colnames(TCGA.expr)], 
    CTNNB1_status = ifelse(TCGA.mut["CTNNB1", colnames(TCGA.expr)] > 0, "MUT", "WT"))

wilcox.test(BAAT ~ CTNNB1_status, data = TCGA.BAAT_data, exact = FALSE)

TCGA.BAAT_data %>% group_by(CTNNB1_status) %>% summarise(n = n(), mean = mean(BAAT), median = median(BAAT), 
    SD = sd(BAAT))

cor.test(TCGA.BAAT_data$BAAT, TCGA.BA4_features$HALLMARK_BILE_ACID_METABOLISM[match(TCGA.BAAT_data$Patient_ID, 
    TCGA.BA4_features$Patient_ID)], method = "spearman")

cor.test(TCGA.BAAT_data$BAAT, TCGA.BA4_features$BA4_ICBScore[match(TCGA.BAAT_data$Patient_ID, TCGA.BA4_features$Patient_ID)], 
    method = "spearman")


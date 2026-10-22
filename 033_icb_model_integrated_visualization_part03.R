if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

TCGA.p_gene_BA4_effect <- TCGA.p_gene_BA4_effect + theme(axis.text.x = element_text(size = 10, color = "black"), 
    axis.text.y = element_text(size = 10, color = "black", face = "italic"), axis.title.x = element_text(size = 12), 
    axis.title.y = element_text(size = 12))

TCGA.effect_dir <- file.path("ICB.model", "TCGA_visualization")

dir.create(TCGA.effect_dir, recursive = TRUE, showWarnings = FALSE)

write.table(TCGA.gene_BA4_effects, file.path(TCGA.effect_dir, "TCGA_gene_BA4_effects.tsv"), sep = "\t", 
    quote = FALSE, row.names = FALSE, na = "NA")

write.table(TCGA.effect_plot_data, file.path(TCGA.effect_dir, "TCGA_gene_BA4_effect_plot_data.tsv"), 
    sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

ggsave(filename = file.path(TCGA.effect_dir, "TCGA_gene_BA4_effect_forest.pdf"), plot = TCGA.p_gene_BA4_effect, 
    width = 4.5, height = max(3.5, 0.28 * nrow(TCGA.effect_plot_data) + 2), limitsize = FALSE)


if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

ici_gene_sets <- list(ICI__GEP18_T_CELL_INFLAMED = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", 
    "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", 
    "STAT1", "TIGIT"), ICI__IFNG6 = c("CXCL9", "CXCL10", "IDO1", "IFNG", "HLA-DRA", "STAT1"), ICI__CYT2 = c("GZMA", 
    "PRF1"), ICI__HCC_INFLAMMATORY4 = c("CD274", "CD8A", "LAG3", "STAT1"), ICI__HCC_TEFF3 = c("CXCL9", 
    "PRF1", "GZMB"), ICI__HCC_ABRS10 = c("CXCR2P1", "ICOS", "TIMD4", "CTLA4", "PAX5", "KLRC3", "FCRL3", 
    "AIM2", "GBP5", "CCL4"), ICI__HCC_ISS10 = c("ITGAL", "IKZF3", "CD8A", "TOX", "GBP5", "PRF1", "CCL5", 
    "LCK", "SLAMF7", "CXCL9"))

sapply(ici_gene_sets, length)


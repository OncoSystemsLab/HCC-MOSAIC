if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

GSE202069.ann <- read.table("raw_data/GSE202069/Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", 
    quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

GSE202069.bulk <- read.table("raw_data/GSE202069/GSE202069_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, 
    sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

GSE202069.clin <- read_xlsx("raw_data/GSE202069/GSE202069_clin.xlsx", sheet = "Sheet1")

{
    GSE202069.annotation.length <- GSE202069.ann %>% select(GeneID, Symbol, Length) %>% filter(!is.na(GeneID), 
        GeneID != "", !is.na(Symbol), Symbol != "", !is.na(Length), Length > 0) %>% distinct(GeneID, 
        .keep_all = TRUE)
    GSE202069.count <- GSE202069.bulk %>% inner_join(GSE202069.annotation.length, by = "GeneID")
    sample_cols <- setdiff(colnames(GSE202069.count), c("GeneID", "Symbol", "Length"))
    GSE202069.rpk <- GSE202069.count
    GSE202069.rpk[sample_cols] <- sweep(as.matrix(GSE202069.count[sample_cols]), 1, GSE202069.count$Length/1000, 
        "/")
    GSE202069.tpm <- GSE202069.rpk
    GSE202069.tpm[sample_cols] <- sweep(as.matrix(GSE202069.rpk[sample_cols]), 2, colSums(as.matrix(GSE202069.rpk[sample_cols]), 
        na.rm = TRUE), "/") * 1e+06
    head(GSE202069.tpm)
    colSums(GSE202069.tpm[sample_cols])
    GSE202069.tpm.symbol <- GSE202069.tpm %>% select(Symbol, all_of(sample_cols)) %>% group_by(Symbol) %>% 
        summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")
    sum(duplicated(GSE202069.tpm.symbol$Symbol))
    colSums(GSE202069.tpm.symbol[-1])
    gse202069.tpm = GSE202069.tpm.symbol
    gse202069.tpm <- gse202069.tpm %>% select(Symbol, all_of(GSE202069.clin$sample))
}


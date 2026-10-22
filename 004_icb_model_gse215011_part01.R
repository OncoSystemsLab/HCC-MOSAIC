if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

GSE215011.ann <- read.table("raw_data/GSE215011/Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", 
    quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

GSE215011.bulk <- read.table("raw_data/GSE215011/GSE215011_raw_counts_GRCh38.p13_NCBI.tsv", header = TRUE, 
    sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

GSE215011.clin <- read_xlsx("raw_data/GSE215011/GSE215011_clin.xlsx", sheet = "Sheet1")

{
    GSE215011.annotation.length <- GSE215011.ann %>% select(GeneID, Symbol, Length) %>% filter(!is.na(GeneID), 
        GeneID != "", !is.na(Symbol), Symbol != "", !is.na(Length), Length > 0) %>% distinct(GeneID, 
        .keep_all = TRUE)
    GSE215011.count <- GSE215011.bulk %>% inner_join(GSE215011.annotation.length, by = "GeneID")
    sample_cols <- setdiff(colnames(GSE215011.count), c("GeneID", "Symbol", "Length"))
    GSE215011.rpk <- GSE215011.count
    GSE215011.rpk[sample_cols] <- sweep(as.matrix(GSE215011.count[sample_cols]), 1, GSE215011.count$Length/1000, 
        "/")
    GSE215011.tpm <- GSE215011.rpk
    GSE215011.tpm[sample_cols] <- sweep(as.matrix(GSE215011.rpk[sample_cols]), 2, colSums(as.matrix(GSE215011.rpk[sample_cols]), 
        na.rm = TRUE), "/") * 1e+06
    head(GSE215011.tpm)
    colSums(GSE215011.tpm[sample_cols])
    GSE215011.tpm.symbol <- GSE215011.tpm %>% select(Symbol, all_of(sample_cols)) %>% group_by(Symbol) %>% 
        summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")
    sum(duplicated(GSE215011.tpm.symbol$Symbol))
    colSums(GSE215011.tpm.symbol[-1])
    GSE215011.tpm = GSE215011.tpm.symbol
}


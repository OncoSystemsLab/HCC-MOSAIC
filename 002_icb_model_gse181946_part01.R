if (!exists("hcc_input_file", mode = "function")) source("000_project_config.R")

GSE181946.ann <- read.table("raw_data/GSE181946/Human.GRCh38.p13.annot.tsv", header = TRUE, sep = "\t", 
    quote = "", comment.char = "", stringsAsFactors = FALSE, check.names = FALSE)

GSE181946.bulk <- read_xlsx("raw_data/GSE181946/GSE181946.xlsx", sheet = "Sheet1")

GSE181946.clin <- read_xlsx("raw_data/GSE181946/GSE181946.xlsx", sheet = "clinical")

{
    GSE181946.annotation.length <- GSE181946.ann %>% select(EnsemblGeneID, Symbol, Length) %>% filter(!is.na(EnsemblGeneID), 
        EnsemblGeneID != "", !is.na(Symbol), Symbol != "", !is.na(Length), Length > 0) %>% distinct(EnsemblGeneID, 
        .keep_all = TRUE)
    GSE181946.count <- GSE181946.bulk %>% mutate(EnsemblGeneID = sub("\\..*$", "", Tag)) %>% inner_join(GSE181946.annotation.length, 
        by = "EnsemblGeneID")
    sample_cols <- setdiff(colnames(GSE181946.count), c("Tag", "EnsemblGeneID", "Symbol", "Length"))
    GSE181946.rpk <- GSE181946.count
    GSE181946.rpk[sample_cols] <- sweep(as.matrix(GSE181946.count[sample_cols]), 1, GSE181946.count$Length/1000, 
        "/")
    GSE181946.tpm <- GSE181946.rpk
    GSE181946.tpm[sample_cols] <- sweep(as.matrix(GSE181946.rpk[sample_cols]), 2, colSums(as.matrix(GSE181946.rpk[sample_cols]), 
        na.rm = TRUE), "/") * 1e+06
    head(GSE181946.tpm)
    colSums(GSE181946.tpm[sample_cols])
    GSE181946.tpm.symbol <- GSE181946.tpm %>% select(Symbol, all_of(sample_cols)) %>% group_by(Symbol) %>% 
        summarise(across(all_of(sample_cols), sum, na.rm = TRUE), .groups = "drop")
    sum(duplicated(GSE181946.tpm.symbol$Symbol))
    colSums(GSE181946.tpm.symbol[-1])
    gse181946.tpm = GSE181946.tpm.symbol
}


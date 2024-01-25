# fluoxetine signature analysis
library(signatureSearch)

# fluox_gct <- system.file("data", "lincs_fluoxetine_signatures.gct", package = "signatureSearch")
#
# h5file <- tempfile(fileext=".h5")
# gctx2h5("./data/lincs_fluoxetine_signatures.gct", cid=1,
#         new_cid=c('fluoxetine__trt_cp'),
#         h5file=h5file, overwrite=TRUE)
#
# fluox_sigs <- rhdf5::h5ls(h5file)
# data(h5file)

# disease signatures
precystic_sig <- read_csv("./res/deseq2_outputs/signature_annotation_p70.csv")
cyst_p21_sig <- read_csv("./res/deseq2_outputs/signature_annotation_p21.csv")
cyst_p28_sig <- read_csv("./res/deseq2_outputs/signature_annotation_p28.csv")

fluox_sig <- read_delim("./data/L1000_LINCS_DCIC_CPC009_HA1E_6H_J24_fluoxetine_10uM.tsv")

fluox_up <- fluox_sig %>%
  slice_max(order_by = `CD-coefficient`, n = 500)

fluox_down <- fluox_sig %>%
  slice_min(order_by = `CD-coefficient`, n = 500)

fluox_sig_genes <- c(fluox_up$symbol, fluox_down$symbol)

# function to find reverse gene relationships # can add arguments for gene column names and gene values
find_reverse_genes <- function(disease_sig, drug_sig){
  rev_genes <- disease_sig %>%
    inner_join(drug_sig, by = c("hgnc_symbol" = "symbol")) %>%
    select(hgnc_symbol, log2FoldChange, padj, CD = `CD-coefficient`) %>%
    filter(log2FoldChange > 0 & CD < 0 | log2FoldChange < 0 & CD > 0)
  return(rev_genes)
}

# filter precystic sig for fluox sig genes
precyst_fluox_reversion <- find_reverse_genes(precystic_sig, rbind(fluox_up, fluox_down))

cyst_p21_reversion <- find_reverse_genes(cyst_p21_sig, rbind(fluox_up, fluox_down))

cyst_p28_reversion <- find_reverse_genes(cyst_p28_sig, rbind(fluox_up, fluox_down))

# save results
write_csv(precyst_fluox_reversion, "./res/fluoxetine/precystic_fluoxetine_reversed_genes.csv")
write_csv(cyst_p21_reversion, "./res/fluoxetine/cysticp21_fluoxetine_reversed_genes.csv")
write_csv(cyst_p28_reversion, "./res/fluoxetine/cysticp28_fluoxetine_reversed_genes.csv")

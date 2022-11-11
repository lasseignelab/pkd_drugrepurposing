
#### Custom functions

#__salmon2dds__: Salmon merged gene counts from nf-core RNAseq pipeline (results/star_salmon/salmon.merged.gene_counts.rds) and SRA metadata as inputs, function data wrangles and runs DESeq2, resulting in a dds object with design by genotype. Just for dataset  GSE69556 (aka 56, p21), set argument data = not "regular", for extra characters in counts rownames just in that dataset

#Salmon merged gene counts to dds results
salmon2dds <- function(counts, metadata, dataset = "regular") {
  #count matrix in the sumarized experiment object
  cts <- counts@assays@data@listData[["counts"]]
  #round, salmon outputs have decimals
  cts <- round(cts)
  if (dataset == "dataset56"){
    colnames(cts) <- str_sub(colnames(cts), 1, -2) ##ALTERATION just for GSE69556
  }

  #remove version numbers from transcript id's
  rownames(cts) <- gsub("\\..*", "", rownames(cts))

  #remove all non-alphanumeric characters from Genotype
  metadata$Genotype <- str_replace_all(metadata$Genotype, "[^[:alnum:]]", "")

  #remove IFT88 for dataset39
  #if (dataset == "dataset39"){
  #metadata <- dplyr::filter(metadata, Genotype != "Ift88flflPkd2flflPax8rtTATetOcre")

  #select only non-IFT88 samples from counts
  #cts <- dplyr::select(cts, c(rownames(metadata)))
  #}

  #contrasts need to be factor levels
  metadata$Genotype <- as.factor(metadata$Genotype)

  #check ordering matches
  print("Checking order of count column names to metadata rownames... If false, fix and rerun!")
  print(colnames((cts)) == rownames(metadata))
  dds <- DESeqDataSetFromMatrix(cts, metadata, design = ~ Genotype)

  return(dds)
}




#__uplfc__ & __downlfc__: Formats deseq2 res() output to select for upregulated and downregulated genes by log2 fold change threshold > 2 and < -2 and p-adjusted < 0.05.
uplfc <- function(res, lfcutoff){
  #NAs removed first
  res <- res[!is.na(res$padj),]
  #p-adjusted value cutoff 0.05 and log2fc greater than 2.0 and less than -2.0
  upfc <- res[res$log2FoldChange > lfcutoff & res$padj < 0.05,]

  #pull just logfold change values from DESeqResults object and name the values by associated gene
  #up <- upfc$log2FoldChange
  #names(up) <- rownames(upfc)
  #up <- cbind(Mouse.gene.stable.ID = rownames(upfc), LFC = upfc$log2FoldChange, padj = upfc$padj)
  #  up <- upfc %>% as.data.frame() %>% dplyr::select(LFC = log2FoldChange, padj = upfc$padj) %>% tibble::rownames_to_column(var = "Mouse.gene.stable.ID")
  up <- upfc %>% tibble::rownames_to_column(var = "Mouse.gene.stable.ID")

  return(up)
}

downlfc <- function(res, lfcutoff){
  #NAs removed first
  res <- res[!is.na(res$padj),]
  #p-adjusted value cutoff 0.05 and log2fc greater than 2.0 and less than -2.0
  downfc <- res[res$log2FoldChange < lfcutoff & res$padj < 0.05,]

  #pull just logfold change values from DESeqResults object and name the values by associated gene
  #down <- downfc$log2FoldChange
  #names(down) <- rownames(downfc)
  #down <- cbind(Mouse.gene.stable.ID = rownames(downfc), LFC = down)
  down <- downfc %>% tibble::rownames_to_column(var = "Mouse.gene.stable.ID")

  return(down)
}


#__convertMouseGeneList__: Uses Ensembl human GRCh38.p13 and Ensembl mouse C57BL_6NJ_v1 genes and to find human orthologous genes (more details [here](https://docs.google.com/document/d/1jU3EOVaZXJMzxwH9SsssHoArz8pMdcfW8MSr8hbeE4Q/edit) ), and biomart to convert to human entrez
convertMouseGeneList <- function(lfc){
  #read in mose to human ens annotation
  ens_annot <- read.csv(here("data", "annot_ens_humanmouse.csv"))
  #conversion table
  ens_mouse2human <- merge(lfc, ens_annot, by = "Mouse.gene.stable.ID")
  #biomart to conver human ens to human entrez
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  genes <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","entrezgene_id",  "hgnc_symbol", "description"), values = ens_mouse2human$Gene.stable.ID, mart = human)

  #names(ens_mouse2human) <- c("mouse_ensembl_gene_id", "log2FoldChange", "ensembl_gene_id")
  conv_tbl <- merge(genes, ens_mouse2human, by.x = "ensembl_gene_id", by.y = "Gene.stable.ID")
  return(conv_tbl)
}


#__toplincs__: Uses uplfc and downlfc outputs to then find top 100 genes that are available in LINCS (best input for signature reversion)
eh <- ExperimentHub()
lincs <- eh[["EH3226"]]; lincs_expr <- eh[["EH3227"]]

lincs_genes <- rhdf5:::h5read(lincs, "rownames", drop=TRUE)

toplincs <- function(conv_tbl, direction){
  lincs_genes <- cbind(entrezgene_id = lincs_genes)
  conv_tbl <- merge(conv_tbl, lincs_genes, by = "entrezgene_id")
  #conv_tbl$log2FoldChange <- as.numeric(conv_tbl$log2FoldChange)
  if(direction == "up"){
    top100 <- slice_max(conv_tbl, n = 100, order_by = log2FoldChange)
  }
  else{
    top100 <- slice_min(conv_tbl, n = 100, order_by = log2FoldChange)
  }
  return(top100)
}

#__lincsGenes__: Wraps convertMouseGeneList to map orthologs, and toplincs to filter for LINCS-available genes
lincsGenes <- function(lfc, direction){
  conv_tbl <- convertMouseGeneList(lfc)
  top100 <- toplincs(conv_tbl, direction = direction)
  return(top100)
}



## Custom functions for FEA

#Enrichment function for up and down genes
enrich <- function(upgenes, downgenes, topreturned) {#, filepath) {
  up_fea <- gost(upgenes, multi_query = FALSE, correction_method = "bonferroni", evcodes = TRUE, sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP"))
  #remove abitrary pathways -- any pathways with > 1000 genes
  up_fea <- up_fea$result %>% dplyr::filter(., term_size < 1000 & term_size > 5) %>% dplyr::mutate(., .keep = "all", direction = "upregulated")
  up_mat <- up_fea  %>% dplyr::slice_max(., recall, n = topreturned)

  down_fea <- gost(downgenes, multi_query = FALSE, correction_method = "bonferroni", evcodes = TRUE, sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP"))
  #remove abitrary pathways -- any pathways with > 1000 genes
  down_fea <- down_fea$result %>% dplyr::filter(., term_size < 1000 & term_size > 5) %>% dplyr::mutate(., .keep = "all", direction = "downregulated")
  down_mat <- down_fea  %>% dplyr::slice_max(., recall, n = topreturned)
  combined_fea <- rbind(up_mat, down_mat)

  fea_res <- list(downregulated = down_fea, upregulated = up_fea, topcompared = combined_fea)

  return(fea_res)
}


#MOUSE GENES Enrichment function for up and down genes
#Runs an ordered query
fea <- function(genelist, sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP")){
  fea <- gost(genelist, organism = "mmusculus", ordered_query = FALSE, multi_query = FALSE, evcodes = TRUE, correction_method = "bonferroni", sources = sources)
  #remove abitrary pathways -- any pathways with > 1000 genes
  if(length(fea) > 0){
    fea <- fea$result %>% dplyr::filter(., term_size < 1000 & term_size > 10)
  }
  return(fea)
}

mouseenrich <- function(upgenes, downgenes, topreturned, sources = c("GO:BP", "GO:MF", "GO:CC", "REAC", "WP")) {

  up_fea <- fea(upgenes, sources = sources)
  down_fea <- fea(downgenes, sources = sources)

  if(length(up_fea) > 0){
    up_fea <- dplyr::mutate(up_fea, .keep = "all", direction = "upregulated")
    up_fea <- up_fea  %>% dplyr::arrange(., desc(recall)) %>% dplyr::slice_head(., n = topreturned)
  }
  if(length(down_fea) > 0){
    down_fea <- dplyr::mutate(down_fea, .keep = "all", direction = "downregulated")
    down_fea <- down_fea  %>% dplyr::arrange(., desc(recall)) %>% dplyr::slice_head(., n = topreturned)
  }

  combined_fea <- rbind(up_fea, down_fea)
  fea_res <- list(downregulated = down_fea, upregulated = up_fea, topcompared = combined_fea)

  return(fea_res)
}



#Custom background for FEA (for enrichment using mouse to human to entrez to lincs conversion -- background gene list is all possible measured genes that have gone through those conversions)

#Needs to only be possible genes -- all genes measured in salmon gene counts that map to human entrez
#mesgenes <- data@NAMES
#remove version numbers
#mesgenes <- gsub("\\..*", "", mesgenes)

#genebg <- convertMouseGeneList(mesgenes)
#write.csv2(genebg, file = "./data/background_genelist_genestoLINCSgenes.csv")

genebg <- read.csv(file = here("data", "background_genelist_genestoLINCSgenes.csv"))



#bubbleplot, up and down for x axis and terms on y axis

bubbleplot <- function(combined_fea){
  bubbleplot <- ggplot(combined_fea, aes(x=direction, y=term_name, size = recall, fill = p_value)) +
    geom_point(alpha=0.7, shape = 21) +
    scale_size(range = c(2, 10), name = "Recall") +
    scale_fill_distiller(palette = "Purples", direction = 1) +
    labs(x = "Differential Expression Direction", y = "Functional Enrichment Terms") +
    theme_minimal(base_size = 8) #+ labs(title = "Top Enriched Terms")
  #  ggsave(filepath, bubbleplot, width = 10, height = 7)
  return(bubbleplot)
}


#barplot, up and down merged, colored by p-value, terms on y axis and # genes on x axis

barplot <- function(combined_fea){
  barplot <- ggplot(combined_fea, aes(x=intersection_size, y=term_name,  fill = recall)) +
    geom_bar( stat = "identity") +
    scale_fill_continuous(type = "viridis") +
    labs(x = "Number of Genes Matched to Term", y = "Functional Terms") +
    theme_minimal(base_size = 8)
  return(barplot)
}



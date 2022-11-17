
#### Custom functions

#__salmon2dds__: Salmon merged gene counts from nf-core RNAseq pipeline (results/star_salmon/salmon.merged.gene_counts.rds) and SRA metadata as inputs, function data wrangles and runs DESeq2, resulting in a dds object with design by genotype. Just for dataset  GSE69556 (aka 56, p21), set argument data = not "regular", for extra characters in counts rownames just in that dataset

#Salmon merged gene counts to dds results
salmon2dds <- function(counts, metadata) {
  #round, salmon outputs have decimals
  cts <- round(counts)

  #remove version numbers from transcript id's
  rownames(cts) <- gsub("\\..*", "", rownames(cts))

  #remove all non-alphanumeric characters from Genotype
  metadata$Genotype <- str_replace_all(metadata$Genotype, "[^[:alnum:]]", "")

  #select only counts that match metadata
  cts <- dplyr::select(cts, c(rownames(metadata)))

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

convertMouseGeneList_2 <- function(mousegenes){
  #read in mouse to human ens annotation
  ens_annot <- read.csv(here("data", "annot_ens_humanmouse.csv"))
  #conversion table
  ens_mouse2human <- dplyr::filter(ens_annot, Mouse.gene.stable.ID %in% mousegenes)
  #biomart to conver human ens to human entrez
  human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  genes <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "entrezgene_id",  "hgnc_symbol", "description"), values = ens_mouse2human$Gene.stable.ID, mart = human)

  conv_tbl <- merge(genes, ens_mouse2human, by.x = "ensembl_gene_id", by.y = "Gene.stable.ID")
  return(conv_tbl)
}

##converting gene IDs for normalized (rlog) gene counts
rlog_wrangle <- function(rlog, geneconv){
  #retrieve counts data and ENMUSG #
  cts_df <- as.data.frame(assays((rlog)))
  #gene_id into a column
  cts_df <- tibble:::rownames_to_column(cts_df, var = "ensembl_gene_id")
  #Merge gene symbols into counts and gene ID
  conv_df <- merge(cts_df, geneconv, by = "ensembl_gene_id")
  #Filter out repeating data
  conv_df <- dplyr::distinct(conv_df, mgi_symbol, .keep_all = TRUE)
  #Convert back to rownames but with gene symbols
  conv_df <- tibble:::column_to_rownames(conv_df, "mgi_symbol")
  #Remove unnecessary data
  conv_df <- conv_df[,!(colnames(conv_df) %in% c("ensembl_gene_id", "group", "group_name", "entrezgene_id", "description"))]
  return(conv_df)
}

#Conversion table from mouse ensembl to human symbol and find targets that are expressed above a LFC cutoff
targets_in_degs <- function(conv_tbl, degs, targets, cutoff){
  #convert degs gene id's
  degs <- tibble:::rownames_to_column(degs, var = "Mouse.gene.stable.ID")
  conv_tbl <- merge(conv_tbl, degs, by = "Mouse.gene.stable.ID")
  print(str(conv_tbl))

  #remove duplicates due to multiple homolog mapping
  conv_degs <- conv_tbl %>% dplyr::arrange(order_by = desc(abs(log2FoldChange))) %>% dplyr::distinct(hgnc_symbol, .keep_all = TRUE)

  #combine degs with drug target data
  targetdegs <- merge(targets, conv_degs, by.x = "t_gn_sym", by.y = "hgnc_symbol")
  targetdegs <- tidyr::drop_na(targetdegs)

  targetdegs <- dplyr::filter(targetdegs, log2FoldChange > cutoff)
  print(unique(targetdegs$t_gn_sym))
  print(unique(targetdegs))
  print(unique(targetdegs$pert))

  targetdeglist <- list(degs = conv_degs, targetdegs = targetdegs)

  return(targetdeglist)
}


target_deg_plot <- function(targetdegs, cutoff, title){

  #targetdegs <- dplyr::distinct(targetdegs, t_gn_sym, .keep_all = TRUE)

  volcanoplot <- EnhancedVolcano(targetdegs$degs, lab = targetdegs$degs$hgnc_symbol,
  x = 'log2FoldChange', y = 'padj', selectLab = unique(targetdegs$targetdegs$t_gn_sym),
  pCutoff = 0.05, FCcutoff = cutoff, pointSize = 3.0, labSize = 5.0,
  col = c("#440154", "#21918c", "#fde725", "#5ec962"), title = title,
  labFace = 'bold', boxedLabels = TRUE,drawConnectors = TRUE,widthConnectors = 0.5,colConnectors = 'black',
  max.overlaps = 25)
  return(volcanoplot)
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


##converting gene IDs for normalized (rlog) gene counts
rlog_wrangle <- function(rlog){
  #retrieve counts data and ENMUSG #
  cts_df <- as.data.frame(assays((rlog)))
  #gene_id into a column
  cts_df <- tibble:::rownames_to_column(cts_df, var = "ensembl_gene_id")
  #Merge gene symbols into counts and gene ID
  conv_df <- merge(cts_df, gene_desc, by = "ensembl_gene_id")
  #Filter out repeating data
  conv_df <- dplyr::distinct(conv_df, mgi_symbol, .keep_all = TRUE)
  #Convert back to rownames but with gene symbols
  conv_df <- tibble:::column_to_rownames(conv_df, "mgi_symbol")
  #Remove unnecessary data
  conv_df <- conv_df[,!(colnames(conv_df) %in% c("ensembl_gene_id", "group", "group_name", "entrezgene_id", "description"))]
  return(conv_df)
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



#venn diagram for 3 data sets/signatures/etc
venn3 <- function(set1, set2, set3, categories, filename){
  myCol <- c("FDE725FF", "287D8EFF", "481567FF")

  venn.diagram(x = list(set1, set2, set3),
               category.names = categories,
               filename = filename, # imagetype = "png",
               output=TRUE,

               # Output features
               imagetype="png" ,
               height = 480 ,
               width = 570 ,
               resolution = 300,
               compression = "lzw",

               # Circles
               lwd = 2,
               col=c("#440154ff", '#21908dff', '#fde725ff'),
               #lty = 'blank',
               fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3), alpha('#fde725ff',0.3)),

               # Numbers
               cex = .6,
               fontface = "bold",
               fontfamily = "sans",

               # Set names
               cat.cex = 0.55,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-27, 27, 135),
               cat.dist = c(0.055, 0.055, 0.085),
               cat.fontfamily = "sans",
               rotation = 1 )
}


go_bp_sim <- function(fea_obj, orgdb, threshold){
  go_up <- dplyr::filter(fea_obj$upregulated, source == "GO:BP")
  go_down <- dplyr::filter(fea_obj$downregulated, source == "GO:BP")
  go_comb <- rbind(go_up, go_down)

  simMatrix <- calculateSimMatrix(go_comb$term_id, orgdb = orgdb, ont = "BP", method = "Wang")

  scores <- setNames(-log10(go_comb$p_value), go_comb$term_id)
  reducedTerms <- reduceSimMatrix(simMatrix, scores, threshold= threshold, orgdb=orgdb)

  sim_obj <- list(simMatrix = simMatrix, reducedTerms = reducedTerms)
  return(sim_obj)
}

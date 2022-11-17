
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

treemap2 <- function (dtf, index, vSize, vColor = NULL, stdErr = NULL, type = "index",
                      fun.aggregate = "sum", title = NA, title.legend = NA, algorithm = "pivotSize",
                      sortID = "-size", mirror.x = FALSE, mirror.y = FALSE, palette = NA,
                      palette.HCL.options = NULL, range = NA, mapping = NA, n = 7,
                      na.rm = TRUE, na.color = "#DDDDDD", na.text = "Missing",
                      fontsize.title = 14, fontsize.labels = 11, fontsize.legend = 12,
                      fontcolor.labels = NULL, fontface.labels = c("bold", rep("plain",
                                                                               length(index) - 1)), fontfamily.title = "sans", fontfamily.labels = "sans",
                      fontfamily.legend = "sans", border.col = "black", border.lwds = c(length(index) +
                                                                                          1, (length(index) - 1):1), lowerbound.cex.labels = 0.4,
                      inflate.labels = TRUE, bg.labels = NULL, force.print.labels = FALSE,
                      overlap.labels = 0.5, align.labels = c("center", "center"),
                      xmod.labels = 0, ymod.labels = 0, eval.labels = FALSE, position.legend = NULL,
                      reverse.legend = FALSE, format.legend = NULL, drop.unused.levels = TRUE,
                      aspRatio = NA, vp = NULL, draw = TRUE, ...)
{
  s <- NULL
  vColor.temp <- i <- w <- h <- NULL
  if (!exists("dtf"))
    stop("Dataframe <dtf> not defined")
  if (!exists("index"))
    stop("Attribute <index> not defined")
  if (!exists("vSize"))
    stop("Attribute <vSize> not defined")
  if (!inherits(dtf, "data.frame"))
    stop("Object <dtf> is not a data.frame")
  if (nrow(dtf) == 0)
    stop("data.frame doesn't have any rows")
  if (any(!index %in% names(dtf)))
    stop("<index> contains invalid column names")
  depth <- length(index)
  if (length(vSize) != 1)
    stop("vSize should be one column name")
  if (!vSize %in% names(dtf))
    stop("vSize is invalid column name")
  if (!is.numeric(dtf[[vSize]]))
    stop(paste("Column(s) in vSize not numeric", sep = ""))
  if (!is.null(stdErr)) {
    if (length(stdErr) != 1)
      stop("stdErr should be one column name")
    if (!stdErr %in% names(dtf))
      stop("stdErr is invalid column name")
    if (!is.numeric(dtf[[stdErr]]))
      stop(paste("Column(s) in stdErr not numeric", sep = ""))
  }
  else {
    stdErr <- vSize
  }
  vColorMplySplit <- function(vColor) {
    divided <- 0
    vColorMply <- unlist(strsplit(vColor, split = "*", fixed = TRUE))
    if (length(vColorMply) == 1) {
      vColorMply <- unlist(strsplit(vColor, split = "/",
                                    fixed = TRUE))
      if (length(vColorMply) == 1) {
        vColorMply <- c(vColorMply, 1)
      }
      else {
        vColorMply[2] <- (1/as.numeric(vColorMply[2]))
        divided <- 1
      }
    }
    return(c(vColorMply, divided))
  }
  if (type %in% c("index", "depth"))
    vColor <- NULL
  if (!is.null(vColor)) {
    if (length(vColor) != 1)
      stop("length of vColor should be one")
    vColor2 <- vColorMplySplit(vColor)
    vColorX <- as.numeric(vColor2[2])
    if (is.na(vColorX))
      stop("vColor: invalid scale factor")
    vColor <- vColor2[1]
    if (!(vColor %in% names(dtf)))
      stop("Invalid column name in vColor.")
    vColorDiv <- as.logical(as.numeric(vColor2[3]))
  }
  if (!type %in% c("value", "categorical", "comp", "dens",
                   "index", "depth", "color", "manual"))
    stop("Invalid type")
  if (!is.function(match.fun(fun.aggregate)))
    stop("fun.aggregate is not a function")
  if (type == "dens")
    fun.aggregate <- "weighted.mean"
  if (!is.na(title[1]) && length(title) != 1) {
    warning("Length of title should be 1")
    title <- NA
  }
  if (is.na(title[1])) {
    title <- vSize
  }
  if (!is.na(title.legend[1]) && length(title.legend) != 1) {
    warning("Length of title.legend should be 1")
    title.legend <- NA
  }
  formatColorTitle <- function(var, varX = NA, var2 = NA, var2X = NA,
                               div) {
    if (!is.na(var2)) {
      if (var2X != 1) {
        if (div)
          var2 <- paste(1/var2X, var2, sep = "*")
        else var <- paste(var2X, var, sep = "*")
      }
      var <- paste(var, "per", var2, sep = " ")
    }
    else {
      if (varX != 1) {
        if (div)
          var <- paste(var, 1/varX, sep = "/")
        else var <- paste(varX, var, sep = "*")
      }
    }
    var
  }
  if (is.na(title.legend[1])) {
    suppressWarnings({
      if (!is.null(vColor)) {
        if (type == "dens")
          title.legend <- formatColorTitle(var = vColor,
                                           var2 = vSize, var2X = vColorX, div = vColorDiv)
        else title.legend <- formatColorTitle(var = vColor,
                                              varX = vColorX, div = vColorDiv)
      }
      else title.legend <- ""
    })
  }
  if (!algorithm %in% c("pivotSize", "squarified"))
    stop("Invalid algorithm")
  if (length(sortID) != 1)
    stop("sortID should be of length one")
  ascending <- substr(sortID, 1, 1) != "-"
  if (!ascending)
    sortID <- substr(sortID, 2, nchar(sortID))
  if (sortID == "size")
    sortID <- vSize
  if (sortID == "color")
    sortID <- vColor
  if (sortID == "se")
    sortID <- stdErr
  if (!(sortID %in% names(dtf)))
    stop("Incorrect sortID")
  if (is.na(palette[1])) {
    if (type == "comp") {
      palette <- brewer.pal(11, "RdYlGn")
    }
    else if (type == "dens") {
      palette <- brewer.pal(9, "OrRd")
    }
    else if (type == "depth") {
      palette <- brewer.pal(8, "Set2")
    }
    else if (type == "index") {
      palette <- "HCL"
    }
    else if (type == "value") {
      palette <- brewer.pal(11, "RdYlGn")
    }
    else if (type == "categorical") {
      palette <- "HCL"
    }
    else if (type == "manual") {
      stop("For \"manual\" treemaps, a palette should be provided.")
    }
  }
  else {
    reverse <- (substr(palette[1], 1, 1) == "-")
    if (reverse)
      palette[1] <- substr(palette[1], 2, nchar(palette[1]))
    if ((length(palette) == 1) && (palette[1] %in% row.names(brewer.pal.info))) {
      palette <- brewer.pal(brewer.pal.info[palette, "maxcolors"],
                            palette)
      if (reverse)
        palette <- rev(palette)
    }
    else {
      if (palette[1] == "HCL" && !(type %in% c("depth",
                                               "index", "categorical"))) {
        stop("HCL palette only applicable for treemap types \"depth\", \"index\" and \"categorical\".")
      }
      if (palette[1] != "HCL" & inherits(try(col2rgb(palette),
                                             silent = TRUE), "try-error")) {
        stop("color palette is not correct")
      }
    }
  }
  palette.HCL.options <- treemap:::tmSetHCLoptions(palette.HCL.options)
  if (!all(is.na(range))) {
    if (length(range) != 2)
      stop("length range should be 2")
    if (!is.numeric(range))
      stop("range is not numeric")
  }
  else range <- c(NA, NA)
  if (!all(is.na(mapping))) {
    if (!length(mapping) %in% c(2, 3))
      stop("length range should be 2 or 3")
    if (!is.numeric(mapping))
      stop("range is not numeric")
    if (length(mapping) == 2) {
      mapping <- c(mapping[1], mean(mapping), mapping[2])
    }
  }
  else mapping <- c(NA, NA, NA)
  if (length(fontsize.title) != 1 || !is.numeric(fontsize.title))
    stop("Invalid fontsize.title")
  if (title == "")
    fontsize.title <- 0
  if (!is.numeric(fontsize.labels))
    stop("Invalid fontsize.labels")
  fontsize.labels <- rep(fontsize.labels, length.out = depth)
  cex_indices <- fontsize.labels/12
  if (length(fontsize.legend) != 1 || !is.numeric(fontsize.legend))
    stop("Invalid fontsize.legend")
  if (!missing(fontcolor.labels))
    if (length(fontcolor.labels) != depth)
      fontcolor.labels <- rep(fontcolor.labels, length.out = depth)
  if (length(fontface.labels) != depth)
    fontface.labels <- rep(fontface.labels, length.out = depth)
  if (length(border.col) != depth)
    border.col <- rep(border.col, length.out = depth)
  if (length(border.lwds) != depth)
    border.lwds <- rep(border.lwds, length.out = depth)
  if (length(lowerbound.cex.labels) != 1 || !is.numeric(lowerbound.cex.labels))
    stop("Invalid lowerbound.cex.labels")
  if (lowerbound.cex.labels < 0 || lowerbound.cex.labels >
      1)
    stop("lowerbound.cex.labels not between 0 and 1")
  if (length(inflate.labels) != 1 || class(inflate.labels) !=
      "logical")
    stop("Invalid inflate.labels")
  if (missing(bg.labels)) {
    bg.labels <- ifelse(type %in% c("value", "categorical"),
                        "#CCCCCCDC", 220)
  }
  else {
    if (length(bg.labels) != 1)
      stop("Invalid bg.labels")
    if (!is.numeric(bg.labels)) {
      if (class(try(col2rgb(bg.labels), silent = TRUE)) ==
          "try-error")
        stop("Invalid bg.labels")
    }
    else {
      if (bg.labels < 0 || bg.labels > 255)
        stop("bg.labels should be between 0 and 255")
    }
  }
  if (length(force.print.labels) != 1 || class(force.print.labels) !=
      "logical")
    stop("Invalid force.print.labels")
  if (length(overlap.labels) != 1 || !is.numeric(overlap.labels))
    stop("Invalid overlap.labels")
  if (overlap.labels < 0 || overlap.labels > 1)
    stop("overlap.labels should be between 0 and 1")
  if (!is.list(align.labels))
    align.labels <- list(align.labels)
  if (length(align.labels) != depth)
    align.labels <- rep(align.labels, length.out = depth)
  lapply(align.labels, function(al) if (!(al[1] %in% c("left",
                                                       "center", "centre", "right") && al[2] %in% c("top", "center",
                                                                                                    "centre", "bottom")))
    stop("incorrect align.labels"))
  if (is.list(xmod.labels))
    xmod.labels <- as.list(xmod.labels)
  if (is.list(ymod.labels))
    ymod.labels <- as.list(ymod.labels)
  if (length(xmod.labels) != depth)
    xmod.labels <- rep(xmod.labels, length.out = depth)
  if (length(ymod.labels) != depth)
    ymod.labels <- rep(ymod.labels, length.out = depth)
  if (missing(position.legend)) {
    position.legend <- switch(type, categorical = "right",
                              depth = "right", index = "none", color = "none",
                              "bottom")
  }
  else {
    if (!position.legend %in% c("right", "bottom", "none"))
      stop("Invalid position.legend")
  }
  if (length(drop.unused.levels) != 1 || class(drop.unused.levels) !=
      "logical")
    stop("Invalid drop.unused.levels")
  if (length(aspRatio) != 1 || (!is.na(aspRatio[1]) && !is.numeric(aspRatio)))
    stop("Invalid aspRatio")
  args <- list(...)
  args$na.rm <- na.rm
  if (inherits(dtf, c("tbl_df"))) {
    dtfDT <- data.table::as.data.table(data.frame(dtf[, c(index, vSize,
                                                          vColor, sortID, stdErr)]))
  }
  else if (data.table::is.data.table(dtf)) {
    dtfDT <- copy(dtf[, c(index, vSize, vColor, sortID, stdErr),
                      with = FALSE])
  }
  else {
    dtfDT <- data.table::as.data.table(dtf[, c(index, vSize, vColor,
                                               sortID, stdErr)])
  }
  if (is.null(vColor)) {
    vColor <- "vColor.temp"
    vColorX <- 1
    dtfDT[, `:=`(vColor.temp, 1)]
    data.table::setcolorder(dtfDT, c(1:(ncol(dtfDT) - 3), ncol(dtfDT),
                                     ncol(dtfDT) - 2, ncol(dtfDT) - 1))
  }
  indexList <- paste0("index", 1:depth)
  setnames(dtfDT, old = 1:ncol(dtfDT), new = c(indexList, "s",
                                               "c", "i", "se"))
  if (vColorX != 1)
    dtfDT[, `:=`(c, c/vColorX)]
  if (fun.aggregate == "weighted.mean") {
    if ("w" %in% names(args)) {
      if (is.character(args$w)) {
        dtfDT[, `:=`(w, eval(parse(text = args$w)))]
      }
      else {
        dtfDT[, `:=`(w, args$w)]
      }
    }
    else {
      dtfDT[, `:=`(w, s)]
    }
  }
  else {
    dtfDT[, `:=`(w, 1)]
  }
  for (d in 1:depth) {
    if (is.numeric(dtfDT[[d]])) {
      fact <- factor(dtfDT[[d]], levels = sort(unique(dtfDT[[d]])))
      dtfDT[, `:=`((d), fact)]
    }
    else if (!is.factor(dtfDT[[d]])) {
      fact <- factor(dtfDT[[d]])
      dtfDT[, `:=`((d), fact)]
    }
  }
  if (is.character(dtfDT[["c"]])) {
    dtfDT[, `:=`(c, factor(c))]
  }
  if (!is.numeric(dtfDT[["i"]])) {
    warning("sortID must be a numeric variable")
    dtfDT[, `:=`(i, integer(nrow(dtfDT)))]
  }
  if (!is.null(stdErr) && !is.numeric(dtfDT[["se"]])) {
    warning("stdErr must be a numeric variable")
    dtfDT[, `:=`("se", integer(dtfDT[["se"]]))]
  }
  setkeyv(dtfDT, indexList)
  datlist <- treemap:::tmAggregate(dtfDT, indexList, type, ascending,
                                   drop.unused.levels, fun.aggregate, args)
  catLabels <- switch(type, categorical = levels(datlist$c),
                      index = levels(datlist$index1), depth = index, standErr = datlist$se,
                      NA)
  if (!draw)
    position.legend <- "none"
  vps <- treemap:::tmGetViewports(vp, fontsize.title, fontsize.labels,
                                  fontsize.legend, position.legend, type, aspRatio, title.legend,
                                  catLabels)
  if (draw)
    treemap:::tmPrintTitles(vps, title, title.legend, position.legend,
                            fontfamily.title, fontfamily.legend)
  if (type == "color") {
    datlist$color <- as.character(datlist$c)
    datlist$colorvalue <- NA
  }
  else {
    attr(datlist, "range") <- 1:2
    datlist <- treemap:::tmColorsLegend(datlist, vps, position.legend,
                                        type, palette, range, mapping, indexNames = index,
                                        palette.HCL.options = palette.HCL.options, border.col,
                                        fontfamily.legend, n, na.color, na.text, format.legend,
                                        reverse.legend)
  }
  datlist <- treemap:::tmGenerateRect(datlist, vps, indexList, algorithm)
  if (mirror.x)
    datlist <- within(datlist, x0 <- 1 - x0 - w)
  if (mirror.y)
    datlist <- within(datlist, y0 <- 1 - y0 - h)
  if (draw) {
    treemap:::tmDrawRect(datlist, vps, indexList, lowerbound.cex.labels,
                         inflate.labels, bg.labels, force.print.labels, cex_indices,
                         overlap.labels, border.col, border.lwds, fontcolor.labels,
                         fontface.labels, fontfamily.labels, align.labels,
                         xmod.labels, ymod.labels, eval.labels)
  }
  upViewport(0 + !is.null(vp))
  tm <- datlist[, c(indexList, "s", "c", "se", "colorvalue",
                    "l", "x0", "y0", "w", "h", "color"), with = FALSE]
  if (type == "dens")
    tm[, `:=`(c, c * s)]
  data.table::setnames(tm, c(index, "vSize", "vColor", "stdErr", "vColorValue",
                             "level", "x0", "y0", "w", "h", "color"))
  tmSave <- list(tm = as.data.frame(tm), type = type, vSize = vSize,
                 vColor = ifelse(vColor == "vColor.temp", NA, vColor),
                 stdErr = stdErr, algorithm = algorithm, vpCoorX = vps$vpCoorX,
                 vpCoorY = vps$vpCoorY, aspRatio = vps$aspRatio, range = range,
                 mapping = mapping, draw = draw)
  invisible(tmSave)
}

#rrvgo altered treemap function for label background transparancy and viridis palette
treemapPlot2 <- function (reducedTerms, size = "score", title = "", ...)
{
  if (!all(sapply(c("treemap"), requireNamespace, quietly = TRUE))) {
    stop("Package treemap and/or its dependencies not available. ",
         "Consider installing it before using this function.",
         call. = FALSE)
  }
  treemap(reducedTerms, index = c("parentTerm", "term"),
          vSize = size, type = "index", title = title, palette = #metafolio::gg_color_hue(length(unique(reducedTerms$parent))),
            viridis(n = length(unique(reducedTerms$parent)), end = 0.9),
          fontcolor.labels = c("#FFFFFFDD", "#00000080"), bg.labels = 200,
          border.col = "#00000080", ...)
}

#rrvgo altered scatterplot function for max overlaps
scatterPlot2 <- function (simMatrix, reducedTerms, size = "score", addLabel = TRUE,
                          labelSize = 3)
{
  if (!all(sapply(c("ggplot2", "ggrepel"), requireNamespace,
                  quietly = TRUE))) {
    stop("Packages ggplot2, ggrepel and/or its dependencies not available. ",
         "Consider installing them before using this function.",
         call. = FALSE)
  }
  x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE,
                k = 2)
  df <- cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points),
                                                          reducedTerms$go), c("term", "parent", "parentTerm", "size")])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = V1, y = V2, color = parentTerm)) +
    ggplot2::geom_point(ggplot2::aes(size = size), alpha = 0.5) +
    ggplot2::scale_color_discrete(guide = "none") + ggplot2::scale_size_continuous(guide = "none",
                                                                                   range = c(0, 25)) + ggplot2::scale_x_continuous(name = "") +
    ggplot2::scale_y_continuous(name = "") + ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank())
  if (addLabel) {
    p + ggrepel::geom_label_repel(aes(label = parentTerm),
                                  data = subset(df, parent == rownames(df)), box.padding = grid::unit(0.5,
                                                                                                      "lines"), size = labelSize, max.overlaps = 20)
  }
  else {
    p
  }
}


#__ssRun__: Wrapper for qsiq(), signature query, using up and down genes, and gess_lincs(), the 'lincs' signatureSearch method for gene expression signature search, sorting the results by WTCS and calculating Tau values, and then filtering for inversely related signatures and for kidney cell lines (HA1E and NKDBA) only
ssRun <- function(upgenes, downgenes, cells){
  q <- qSig(query = list(upset=as.character(upgenes), downset=as.character(downgenes)), gess_method="LINCS", refdb= lincs)
  gess <- gess_lincs(q, sortby="WTCS", tau=TRUE, workers=1)
  gess_res <- result(gess)
  gess_res <- dplyr::filter(gess_res, WTCS < 0) #filter for drugs with only inverse signatures
  if(cells == "kidney"){
    gess_kidney <- gess_res %>% dplyr::filter(cell == "HA1E" | cell == "NKDBA")
    return(gess_kidney)
  }
  else{
    return(gess_res)
  }
}

#adaptation of Jen's function for checking FDA approval
#This function is check the FDA-approval of candidates. The approved list comes from Drugs@FDA
FDA_APPROVAL_CHECK<- function(drug_list, fda_approved){
  return_list <- c()
  for (i in 1:length(drug_list)){
    test <- agrep( drug_list[i], fda_approved$ActiveIngredient, ignore.case= TRUE, max.distance = 0.1)
    if (length(test)>0 ){
      return_list[i] <- TRUE
    }else{
      return_list[i] <- FALSE
    }

  }
  names(return_list)<- drug_list
  return(return_list)
}


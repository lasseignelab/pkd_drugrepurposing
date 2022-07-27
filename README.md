# pkd_drugrepurposing

PKD drug repurposing analysis. Drug candidate discovery and prioritization for PKD using publicly available PKD2 KO mouse data sets.
- GSE149739 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149739)
- GSE134719 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134719)
- GSE69556 (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69556)

## code/. 
### DESeq2. 
PKD drug repurposing analysis. DESeq2 analysis was run in separate markdowns for each dataset.
- deseq2_GSE149739.Rmd
- deseq2_GSE134719.Rmd
- deseq2_GSE69556.Rmd

### FEA and Comparisons. 
Functional enrichment analyses and data set comparisons.
- compare_datasets.Rmd

### Signature Reversion and Drug Annotations. 
- sigsearch_analyses.Rmd 

## data/. 
- Metadata for all 3 datasets (as .txt files)
- Processed STAR-salmon count data for all 3 datasets
- Ensembl human and mouse gene conversion annotation
- Background gene list (for FEA)

## res/. 
- deseq2_outputs (degs)
- sigsearch_outputs
- fea

### Versions  
Platform
```
x86_64-apple-darwin17.0 
R version 4.1.1 (2021-08-10)
```
Packages
```
##              biomaRt                rhdf5        ExperimentHub 
##             "2.50.3"             "2.38.1"              "2.2.1" 
##        AnnotationHub        BiocFileCache               dbplyr 
##              "3.2.2"              "2.2.1"              "2.1.1" 
##      signatureSearch  signatureSearchData                 Rcpp 
##              "1.8.2"              "1.8.4"            "1.0.8.3" 
##              ggrepel              ggplot2                dplyr 
##              "0.9.1"              "3.3.5"              "1.0.8" 
##             pheatmap               apeglm              stringr 
##             "1.0.12"             "1.16.0"              "1.4.0" 
##               DESeq2 SummarizedExperiment              Biobase 
##             "1.34.0"             "1.24.0"             "2.54.0" 
##       MatrixGenerics          matrixStats        GenomicRanges 
##              "1.6.0"             "0.61.0"             "1.46.1" 
##         GenomeInfoDb              IRanges            S4Vectors 
##             "1.30.1"             "2.28.0"             "0.32.4" 
##         BiocGenerics 
##             "0.40.0"
```

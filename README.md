# pkd_drugrepurposing
### Authors
- Elizabeth J. Wilk, Timothy C. Howton, Jennifer L. Fisher, Vishal Oza, Ryan Brownlee, Kasi McPherson, Hanna Cleary, Brittany N. Lasseigne

PKD drug repurposing analyses. Drug candidate discovery and prioritization for PKD using publicly available pre-cystic and cystic kidney PKD2 KO mouse data sets.
- Pre-Cystic P70: [GSE149739](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149739)
- Cystic P21: [GSE134719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134719)
- Cystic P28: [GSE69556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69556)

![alt text](https://github.com/lasseignelab/pkd_drugrepurposing/res/biorender/repo_workflow.png?raw=true)

## code/. 
### DESeq2. 
PKD drug repurposing analysis. DESeq2 analyses and visualizations, as well as drug target comparison (from pre-cystic p70 data ran through sigsearch_analyses.Rmd) to differential expression of each dataset. 
- deseq2_analyses.Rmd

### FEA and Comparisons. 
Data set comparisons and functional enrichment analyses with [gprofiler2](https://academic.oup.com/nar/article/47/W1/W191/5486750)
- compare_datasets.Rmd

### Signature Reversion and Drug Annotations. 
Analyses using [signatureSearch](https://github.com/girke-lab/signatureSearch) for signature reversion and the [Drug Repurposing Hub](https://www.nature.com/articles/nm.4306) and [onSides](https://github.com/tatonetti-lab/onsides) 
- sigsearch_analyses.Rmd 

## data/. 
- Metadata for all 3 datasets (as .txt files) (*_metadata.txt)
- Processed STAR-salmon count data for all 3 datasets (*_salmonmerged_gene_counts.rds)
- Ensembl human and mouse gene conversion annotation (annot_ens_humanmouse.csv)
- Background gene list (for FEA) 
- Side effects daata from onSides  
- Data from the Drug Repurposing Hub  

## res/. 
- deseq2_outputs (degs)
- sigsearch_outputs
- fea

## Acknowledgements

 - Michal Mrug
 - Brad Yoder

### Versions  

Dependencies
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

Platform
```
x86_64-apple-darwin17.0 
R version 4.1.1 (2021-08-10)
```

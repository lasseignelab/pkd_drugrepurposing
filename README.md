# Drug Repurposing Analyses in Polycystic Kidney Disease
### Authors
- Elizabeth J. Wilk, Timothy C. Howton, Jennifer L. Fisher, Vishal Oza, Ryan Brownlee, Kasi McPherson, Hanna Cleary, Brittany N. Lasseigne

## Overview

Autosomal dominant polycystic kidney disease (ADPKD) is characterized by renal cyst expansion and is primarily caused by variants in the PKD1 or PKD2 gene that encode the transmembrane proteins Polycystin-1 (PC1) and Polycystin-2 (PC2), respectively.[1](https://pubmed.ncbi.nlm.nih.gov/29326913/) The current only FDA approved drug for ADPKD is Tolvaptan, a vasopressin receptor 2 antagonist that cannot be used as long-term treatment due to liver toxicity side effects.[2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6873754/) Molecular signature reversion is a method of drug repurposing that could lead to new treatment options for ADPKD patients significantly sooner than traditional drug discovery. As gene expression profiles from ADPKD patients and preclinical models have significantly altered transcriptomic signatures, this approach presents an opportunity to compare ADPKD disease signatures to drug response signatures from cell lines treated with drugs and identify candidates that may reverse ADPKD-associated cellular phenotypes, ultimately slowing or reducing kidney cyst growth. In these analyses, we evaluated transcriptomic signatures to determine drug repurposing candidates for ADPKD using publicly available Pkd2 mouse data by detecting inversely related gene expression signatures from the Library of Integrated Network-Based Cellular Signatures (LINCS)[3](https://pubmed.ncbi.nlm.nih.gov/29195078/) database. Drug candidates were further prioritized based on their known mechanism of action (MOA), FDA status, targets, side effects,[4](https://github.com/tatonetti-lab/onsides) and functional enrichment analysis. 
  
 ### Data Sets
- Pre-Cystic P70: [GSE149739](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149739)
- Cystic P21: [GSE134719](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134719)
- Cystic P28: [GSE69556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE69556) 

 
## code/. 
![alt text](res/biorender/repo_workflow.png)
 
### DESeq2. 
PKD drug repurposing analysis. DESeq2 analyses and visualizations, as well as drug target comparison (from drug candidates found from pre-cystic p70 signature reversion in sigsearch_analyses.Rmd) to differentially expressed genes in each dataset. 
- deseq2_analyses.Rmd

### FEA and Comparisons. 
Data set comparisons and functional enrichment analyses with [gprofiler2](https://academic.oup.com/nar/article/47/W1/W191/5486750)
- compare_datasets.Rmd

### Signature Reversion and Drug Annotations. 
Analyses using [signatureSearch](https://github.com/girke-lab/signatureSearch) for signature reversion and the [Drug Repurposing Hub](https://www.nature.com/articles/nm.4306) and [onSides](https://github.com/tatonetti-lab/onsides) 
- sigsearch_analyses.Rmd 

## data/. 
- Metadata each dataset (*_metadata.txt)
- Processed STAR-salmon count data (*_salmonmerged_gene_counts.rds)
- Ensembl human and mouse gene conversion annotation (annot_ens_humanmouse.csv)
- Background gene list (for FEA) (background_genelist_genestoLINCSgenes.csv)
- Side effects data from onSides (onsides_v01)
  - Adverse reactions from "ADVERSE REACTIONS" section of most recent labels (adverse_reactions.csv)
  - Adverse reactions from "BOXED WARNINGS" section most recent labels (boxed_warnings.csv)
  - Most recent label for each drug/combination of drugs (latest_labels_bydrug.csv)
- Drug data incuding clinical trial status, MOA, original indication, etc. (Repurposing_Hub_export.txt)


## res/. 
- deseq2_outputs
- sigsearch_outputs
- FEA
- biorender (repo workflow image)

## Acknowledgements

 - Michal Mrug
 - Brad Yoder

## License

[MIT](https://choosealicense.com/licenses/mit/)

### Versions  

Dependencies
```
##  signatureSearchData      EnhancedVolcano                 here 
##             "1.10.0"             "1.14.0"              "1.0.1" 
##              viridis          viridisLite              biomaRt 
##              "0.6.2"              "0.4.0"             "2.52.0" 
##                rhdf5        ExperimentHub        AnnotationHub 
##             "2.40.0"              "2.4.0"              "3.4.0" 
##        BiocFileCache               dbplyr      signatureSearch 
##              "2.4.0"              "2.2.1"             "1.10.0" 
##                 Rcpp              ggrepel              ggplot2 
##              "1.0.9"              "0.9.1"              "3.3.6" 
##                dplyr             pheatmap               apeglm 
##              "1.0.9"             "1.0.12"             "1.18.0" 
##              stringr               DESeq2 SummarizedExperiment 
##              "1.4.0"             "1.36.0"             "1.26.1" 
##              Biobase       MatrixGenerics          matrixStats 
##             "2.56.0"              "1.8.1"             "0.62.0" 
##        GenomicRanges         GenomeInfoDb              IRanges 
##             "1.48.0"             "1.32.2"             "2.30.0" 
##            S4Vectors         BiocGenerics 
##             "0.34.0"             "0.42.0"

```

Platform
```
x86_64-apple-darwin17.0 
R version 4.1.1 (2021-08-10)
```

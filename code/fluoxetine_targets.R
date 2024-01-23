# script to pull fluoxetine drug targets and compare by tissue and sex specificity
library(ExperimentHub); library(rhdf5)
library(signatureSearch)
library(readxl)
library(readr)
library(tidyverse)
library(gprofiler2)
library(ComplexHeatmap)
# devtools::install_version("dbplyr", version = "2.3.4")
#
# # set cache location
# ExperimentHub(setExperimentHubOption("CACHE", "./bin/docker/"))
# # signaturesearch drug data
# eh <- ExperimentHub()
#
# query(eh, c("signatureSearchData", "lincs2"))
# lincs2_path <- eh[['EH7297']]
# rhdf5::h5ls(lincs2_path)

# read in gtex specificity data
t <- read_rds("./data/sexde.ez.fl_pc3.v_mle.mash_model.FILTERED2.rds")
# lfsr = local false sign rate, used with 0.05 cutoff in gtex paper
# PosteriorMean = Mash 'beta' used in gtex paper
gtex_sexde <- as.data.frame(gtex_sexde$result)
# remove enembl gene versions
rownames(gtex_sexde) <- gsub("\\.\\d+$", "", rownames(gtex_sexde))

gtex_readin <- function(gtex_acron, exppath = NULL, metapath = NULL){
  gtex <- read_excel(exppath, sheet = gtex_acron, col_names = TRUE)
  #gtex <- gtex %>% dplyr::mutate(tissue = tissue)
  meta <- read_excel(metapath, sheet = 2, col_names = TRUE)

  meta_sub <- meta %>% dplyr::filter(`Abbreviation Id` == gtex_acron)
  gtex <- gtex %>% dplyr::mutate(Tissue = meta_sub$Id) %>% dplyr::left_join( meta_sub, by = c("Tissue" = "Id"))

  return(gtex)
}

kidney <- gtex_readin("KDNCTX", "./data/gtex_sexspecificity_s2.xlsx", "./data/gtex_sexspecificity_s1.xlsx")

# read in previous drug - target data
drugtargets_p70 <- read_csv("res/sigsearch_outputs/drugtargets_p70.csv")

# filter for fluoxetine
fluox <- drugtargets_p70[,-1] %>%
  filter(pert == "fluoxetine")

fluox_targets <- fluox$t_gn_sym
fluox_ensmbl <- gprofiler2::gconvert(fluox_targets)

# filter GTEx sex-specific degs for fluox targets
gtex_fluox <- gtex_sexde[c(fluox_ensmbl$target),]

fluox_beta <- gtex_fluox %>%
  select(starts_with("PosteriorMean")) %>%
#  select(sum(is.na(everything)) != nrow(.))
  purrr::discard(~all(is.na(.))) %>%
  rename_all(~stringr::str_remove(., "PosteriorMean."))

fluox_sd <- gtex_fluox %>%
  select(starts_with("PosteriorSD"))

fluox_lfsr <- gtex_fluox %>%
  select(starts_with("lfsr")) %>%
  rename_all(~stringr::str_remove(., "lfsr.")) %>%
  select(contains(colnames(fluox_beta)))


png("./res/test.png", width = 25, height = 40, units = "cm", res = 300)

ComplexHeatmap::Heatmap(as.matrix(fluox_beta), cluster_rows = FALSE, cluster_columns = FALSE,
                        row_labels = fluox_ensmbl$input,
                        cell_fun = function(j,i,x,y,w,h,fill){
                          if(is.na(as.matrix(fluox_lfsr)[i,j])){
                            grid.text(" ", x, y)
                          } else if(as.matrix(fluox_lfsr)[i,j] < 0.05){
                            grid.text("*", x, y)
                          }
                          else{
                          grid.text(" ", x, y)
                        }
                          }
                        ) #annotation_colors = color_anno, #annotation_row = Target,
                         #annotation_col = colnames(target_tissue_exp),
#                         cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE,
                         #labels_col = stringr::str_trunc(colnames(target_tissue_exp_hp), 18), color = myColor, #color = brewer.pal(100, name = "PRGn"),
#                         fontsize = 20, fontface = "bold", fontsize_row = 16, fontsize_col = 16, annotation_col = tissue_annot,# viridis(n= 10000, alpha = 1, begin = 0, end = 1, option = "viridis"),
#                         border_color = "NA", main = "Sex-Biased DE of \n Drug Targets by Tissue", angle_col = "90" )#,
                         #breaks = myBreaks, heatmap_legend_param = list(title = "Expression", at = c(-1, 0, 1), labels = c("Male", "None", "Female")))


dev.off()


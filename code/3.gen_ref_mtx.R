library(Seurat)
library(readxl)
library(edgeR)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(matrixStats)
library(data.table)
library(SingleCellExperiment)
library(Matrix.utils)

####generate ref mtx
####get selected genes in each cluster
tissue <- "tissue"
deg <- read.csv(paste0(tissue, "/celltype/", tissue, "_50deg_rfe_roc_summary.txt"), sep="\t")

files <- list.files(path = paste0(tissue, "/celltype"), pattern = "top50_deg_ref.txt$", recursive = T, full.names = T)
df <- data.frame(tissue = character(),
                 cluster = character(),
                 gene = character(),
                 var_imp = numeric())
for (file in files){
  tdf <- read.csv(file, sep="\t")
  ####get the preferred num of deg from the roc results
  num <- filter(deg, cluster==unique(tdf$cluster)) %>% pull(selected)
  if (nrow(tdf)<num){
    print(paste0(unique(tdf$cluster), " varImp reported ", dim(tdf)[1], " < ", num, " genes"))
  } else {
    tdf_sum <- slice(tdf, 1: num) %>% mutate(cluster=unique(tdf$cluster))
    df <- bind_rows(df, tdf_sum)}}
print(paste0(tissue, " ref mtx has ", nrow(df), " genes; unique ", length(unique(df$gene)), " genes"))
write.table(df, file = paste0(tissue, "/celltype/", tissue, "_deg_rfe_varimp.txt"), sep = "\t", quote = F, row.names = F)

####calculate gene exp across cell types
deg <- read.table(paste0(tissue, "/celltype/", tissue, "_deg_rfe_varimp.txt"), sep = "\t", check.names = F, header = T)
count.norm <- read.table(paste0(tissue, "/", tissue, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", check.names = F)
genes <- unique(deg$gene)
clusters <- unique(deg$cluster)
ref <- data.frame(GeneSymbol = character())
for (type in clusters){
  idx <- which(tstrsplit(colnames(count.norm), "_")[[1]] == type) 
  exp <- count.norm[genes, idx] 
  exp <- exp %>% mutate(GeneSymbol=rownames(exp), !!type:=rowMeans(exp)) %>% select(GeneSymbol, !!type)
  ref <- full_join(ref, exp, by="GeneSymbol") %>% replace(is.na(.), 0)}
write.table(ref, file = paste0(, tissue, "/", tissue, "_ref_matrix.txt"), sep = "\t", quote = F, row.names = F)
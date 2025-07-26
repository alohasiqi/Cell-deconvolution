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

setwd("your_path/")

#### prepare your own single-cell seurat object
#### prep ref matrix
obj <- readRDS("your_path/your_obj.rds")

tissue <- "tissue"
celltype <- "cluster_annot"
sample_id <- "sample_name"

#### Create SCE 
counts <- obj@assays$RNA@counts 
metadata <- obj@meta.data
metadata$cluster_id <- factor(obj@meta.data[,celltype]) # Set up metadata as desired for aggregation and DE analysis
metadata$sample_id <- factor(metadata[,sample_id])

sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
sce <- sce[rowSums(counts(sce) > 0) >= 10, ] # Remove low expressed genes which have less than 10 cells with any counts

cluster_names <- levels(colData(sce)$cluster_id)
length(cluster_names) 
sample_names <- levels(colData(sce)$sample_id)
length(sample_names) 

# Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
groups <- colData(sce)[, c("cluster_id", "sample_id")]
# Aggregate across cluster-sample groups
# transposing row/columns to have cell_ids as row names matching those of groups
aggr_counts <- aggregate.Matrix(t(counts(sce)), 
                                groupings = groups, fun = "sum") 
# Transpose aggregated matrix to have genes as rows and samples as columns
aggr_counts <- t(aggr_counts)
print(paste0(tissue, " gene# is ", dim(aggr_counts)[1], "; ", tissue, " sample# is ", dim(aggr_counts)[2]))

df <- as.data.frame(as.matrix(aggr_counts))
dim(df)
write.table(df, file = paste0(tissue, "/celltype/", tissue, "_pseudobulk_raw_counts_cell_type.txt"), sep = "\t", quote = F, col.names=colnames(aggr_counts))

#clean up aggr_counts, filter genes
var_regex <- "^HLA-|^IG[HJKL]|^RNA|^MT|^RP" # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
var_idx <- which(rownames(aggr_counts) %in% grep(var_regex, rownames(aggr_counts), invert=T, value=T))
aggr_counts <- aggr_counts[var_idx,] 
dge_obj <- DGEList(counts=aggr_counts)

#remove genes with <20% sample counts
m <- sweep(dge_obj$counts, 2, 1e6 / dge_obj$samples$lib.size, `*`)

sample_counts <- dge_obj$samples %>% mutate(cell_type=tstrsplit(rownames(dge_obj$samples), "_")[[1]])
sample_sum <- sample_counts %>% dplyr::count(cell_type)
print(paste0(tissue, " ", sample_sum[which(sample_sum$n==min(sample_sum$n)), 1], " has the minimal number ", min(sample_sum$n), " samples"))

min_count_no <- 0.2*min(sample_sum$n)
ridx <- rowSums(m > 1) >= min_count_no
exp_fil <- table(ridx) 
names(exp_fil) <- c('removed','retained') 
print(paste0(tissue, ", removed ", exp_fil[1], ", keep ", exp_fil[2]))

dge_obj <- dge_obj[ridx,] 
dim(dge_obj$counts)

# TMM normilization
dge_obj <- calcNormFactors(dge_obj,method="TMM")

#design matrix
nf <- calcNormFactors(dge_obj$counts,method="TMM")

#estimate dispersion
dge_obj <- estimateCommonDisp(dge_obj,verbose=TRUE) #overdispersion,biological/sample-sample variability

count.norm <- edgeR::cpm(dge_obj, log = TRUE, prior.count=1)
write.table(count.norm, file = paste0(tissue,"/celltype/", tissue, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", quote = F)

print(paste0("finished ", tissue, " count mtx, start DEG"))

#run voom and return deg
clusters <- unique(tstrsplit(rownames(dge_obj$samples), "_")[[1]])

for (cluster in clusters){
  dge_obj$samples <- mutate(dge_obj$samples, sample=tstrsplit(rownames(dge_obj$samples), "_")[[2]], 
                            group=tstrsplit(rownames(dge_obj$samples), "_")[[1]], 
                            comp_group=ifelse(group==cluster, T, F))
  design <- model.matrix(~0+comp_group,data=dge_obj$samples)
  
  #deg using edgeR
  y <- voom(dge_obj$counts,design, plot=TRUE, lib.size=colSums(dge_obj$counts)*nf, save.plot=T)
  fit <- lmFit(y,design)
  contrast <- makeContrasts(comp_groupTRUE-comp_groupFALSE, levels=design)
  fit2 <- contrasts.fit(fit,contrast)
  fit2 <- eBayes(fit2)
  print(paste0("finished ", cluster, " deg analysis"))
  cluster_deg <- topTable(fit2,coef=1,adjust="BH",number=Inf,confint=TRUE, sort.by = "p")
  
  #save deg for each cluster
  write.table(cluster_deg, file = paste0(tissue, "/celltype/", tissue, "_pseudobulk_", cluster, "_deg.txt"), sep = "\t", quote = F)}

#select top 50 deg from each cluster for each tissue and save
df <- data.frame(cluster = character(),
                 gene = character(),
                 logFC = numeric(),
                 P.Value = numeric())

files <- list.files(path = paste0(tissue, "/celltype"), pattern = "*_deg.txt$", recursive = T, full.names = T)

print(paste0(tissue, ": ", length(files)))

for (file in files){
  all_deg <- read.table(file)
  type <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][3] #get the cluster name
  print(type)
  tdf <- all_deg %>% mutate(gene=rownames(all_deg), cluster=type) %>% filter(logFC>0) %>% 
    filter(gene %in% grep(var_regex, gene, invert=T, value=T)) %>% top_n(-50, P.Value) %>% 
    dplyr::select(cluster, gene, logFC, P.Value)
  df <- bind_rows(df, tdf)}

write.table(df, file = paste0(tissue, "/celltype/", tissue, "_pseudobulk_top50_clusters_deg.txt"), sep = "\t", quote = F, row.names = F)

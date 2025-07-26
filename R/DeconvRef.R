#' This function allows to build user specified CRM using their own scRNA-seq
#' It would return the CRM files and all the intermediate files.
#' Txt files Txt files of the intermediate results can be downloaded in the sub-directory.
#' @title Generate CRM
#' @importFrom Seurat
#' @importFrom tidyverse
#' @importFrom tibble
#' @importFrom stringr
#' @importFrom readxl
#' @importFrom readr
#' @importFrom dplyr
#' @importFrom SingleCellExperiment
#' @importFrom matrixStats
#' @importFrom Matrix.utils
#' @importFrom edgeR       
#' @importFrom e1071
#' @importFrom caret
#' @importFrom randomForest 
#' @importFrom PRROC
#' @importFrom Metrics 
#' @importFrom data.table  
#' @param analysis_name A character vector representing name of analyses
#' @param file A character of a single path to scRNA-seq in rds format
#' @param path A character of a single path to save all the results
#' @param celltype A character specify the layer to the scRNA-seq to read cell type information
#' @param sample_id A character specify the layer to the scRNA-seq to read sample name information
#' @param top_deg_num A numeric specify how many top DEGs will be selected
#' @examples
#' deconvref(analysis_name="test", 
#' file="data/test.nsclc.seur.meta.norm.qc.rds", 
#' path="/your_path/, 
#' celltype="cluster_annot", 
#' sample_id="sample_name", 
#' top_deg_num=50)
#' @export
deconvref <- function(analysis_name="test", 
          file="data/test.nsclc.seur.meta.norm.qc.rds", 
          path="/your_path/", 
          celltype="cluster_annot", 
          sample_id="sample_name", 
          top_deg_num=50){
  pseudobulk_deg(analysis_name, 
                 file, 
                 path, 
                 celltype, 
                 sample_id, 
                 top_deg_num)
  cell_type_markers(path,analysis_name, top_deg_num)
  gen_crm(path, analysis_name) 
}

#' This function allows to upload their own scRNA-seq, build the pseudobulk and find DEGs.
#' It would return the pseudobulk files and DEGs from each cluster.
#' Txt files Txt files of the intermediate results can be downloaded in the sub-directory.
#' @title Generate DEGs
#' @importFrom Seurat
#' @importFrom tidyverse
#' @importFrom dplyr
#' @importFrom readr
#' @importFrom stringr
#' @importFrom readxl
#' @importFrom SingleCellExperiment
#' @importFrom matrixStats
#' @importFrom Matrix.utils
#' @importFrom edgeR 
#' @param analysis_name A character vector representing name of analyses
#' @param file A character of a single path to scRNA-seq in rds format
#' @param path A character of a single path to save all the results
#' @param celltype A character specify the layer to the scRNA-seq to read cell type information
#' @param sample_id A character specify the layer to the scRNA-seq to read sample name information
#' @param top_deg_num A numeric specify how many top DEGs will be selected
#' @examples
#' pseudobulk_deg(analysis_name="test", 
#' file="data/test.nsclc.seur.meta.norm.qc.rds", 
#' path="/your_path/, 
#' celltype="cluster_annot", 
#' sample_id="sample_name", 
#' top_deg_num=50)
#' @export

#### prepare your own single-cell seurat object
#### Create SCE 
pseudobulk_deg <- function(analysis_name="test", 
                       file="data/test.nsclc.seur.meta.norm.qc.rds", 
                       path="/your_path/", 
                       celltype="cluster_annot", 
                       sample_id="sample_name", 
                       top_deg_num=50){
  
  dir.create(file.path(paste0(path, analysis_name, "/celltype/")), recursive = TRUE) #directory for cell type analyses
  dir.create(file.path(paste0(path, analysis_name, "/intermediates/")), recursive = TRUE) #directory for saving intermediate files
  
  obj <- readRDS(file)
  counts <- obj@assays$RNA@counts 
  metadata <- obj@meta.data
  metadata$cluster_id <- factor(obj@meta.data[,celltype]) # Set up metadata as desired for aggregation and DE analysis
  metadata$sample_id <- factor(metadata[,sample_id])
  
  sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
  sce <- sce[rowSums(counts(sce) > 0) >= 10, ] # Remove low expressed genes which have less than 10 cells with any counts
  
  cluster_names <- levels(colData(sce)$cluster_id)
  sample_names <- levels(colData(sce)$sample_id)
  # Subset metadata to include only the variables you want to aggregate across (here, we want to aggregate by sample and by cluster)
  groups <- colData(sce)[, c("cluster_id", "sample_id")]
  # Aggregate across cluster-sample groups
  # transposing row/columns to have cell_ids as row names matching those of groups
  aggr_counts <- Matrix.utils::aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  # Transpose aggregated matrix to have genes as rows and samples as columns
  aggr_counts <- t(aggr_counts)
  print(paste0(analysis_name, " includes ", dim(aggr_counts)[1], " genes and ", dim(aggr_counts)[2], " samples"))
  
  df <- as.data.frame(as.matrix(aggr_counts))
  write.table(df, file = paste0(analysis_name, "/intermediates/", analysis_name, "_pseudobulk_raw_counts_cell_type.txt"), sep = "\t", quote = F, col.names=colnames(aggr_counts))
  
  #clean up aggr_counts, filter genes
  var_regex <- "^HLA-|^IG[HJKL]|^RNA|^MT|^RP" # remove HLA, immunoglobulin, RNA, MT, and RP genes based on HUGO gene names
  var_idx <- which(rownames(aggr_counts) %in% grep(var_regex, rownames(aggr_counts), invert=T, value=T))
  aggr_counts <- aggr_counts[var_idx,] 
  dge_obj <- DGEList(counts=aggr_counts)
  
  #remove genes with <20% sample counts
  m <- sweep(dge_obj$counts, 2, 1e6 / dge_obj$samples$lib.size, `*`)
  
  sample_counts <- dge_obj$samples %>% mutate(cell_type=tstrsplit(rownames(dge_obj$samples), "_")[[1]])
  sample_sum <- sample_counts %>% dplyr::count(cell_type)
  print(paste0(analysis_name, " cell type ", sample_sum[which(sample_sum$n==min(sample_sum$n)), 1], " has the minimal number ", min(sample_sum$n), " samples"))
  
  min_count_no <- 0.2*min(sample_sum$n)
  ridx <- rowSums(m > 1) >= min_count_no
  exp_fil <- table(ridx) 
  names(exp_fil) <- c('removed','retained') 
  print(paste0("removed ", exp_fil[1], " and, keep ", exp_fil[2], " genes for QC"))
  
  dge_obj <- dge_obj[ridx,] 
  dim(dge_obj$counts)
  
  # TMM normilization
  dge_obj <- calcNormFactors(dge_obj,method="TMM")
  
  #design matrix
  nf <- calcNormFactors(dge_obj$counts,method="TMM")
  
  #estimate dispersion
  dge_obj <- estimateCommonDisp(dge_obj,verbose=TRUE) #overdispersion,biological/sample-sample variability
  
  count.norm <- edgeR::cpm(dge_obj, log = TRUE, prior.count=1)
  write.table(count.norm, file = paste0(path, analysis_name,"/intermediates/", analysis_name, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", quote = F)
  
  print(paste0("finished building ", analysis_name, " pseudobulk count matrix"))

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
    write.table(cluster_deg, file = paste0(path, analysis_name, "/celltype/", analysis_name, "_pseudobulk_", cluster, "_deg.txt"), sep = "\t", quote = F)}
  
  #select top deg from each cluster for each analysis and save
  df <- data.frame(cluster = character(),
                   gene = character(),
                   logFC = numeric(),
                   P.Value = numeric())
  
  files <- list.files(path = paste0(path, analysis_name, "/celltype"), pattern = "*_deg.txt$", recursive = T, full.names = T)
  for (file in files){
    all_deg <- read.table(file)
    type <- strsplit(strsplit(file, "/")[[1]][3], "_")[[1]][3] #get the cluster name
    tdf <- all_deg %>% mutate(gene=rownames(all_deg), cluster=type) %>% filter(logFC>0) %>% 
      filter(gene %in% grep(var_regex, gene, invert=T, value=T)) %>% top_n(-top_deg_num, P.Value) %>% 
      dplyr::select(cluster, gene, logFC, P.Value)
    df <- bind_rows(df, tdf)}
  write.table(df, file = paste0(path, analysis_name, "/intermediates/", analysis_name, "_pseudobulk_topdeg_clusters_deg.txt"), sep = "\t", quote = F, row.names = F)
}

#functions, prepare for different training datasets
ds_train <- function(ds, cluster, deg){
  cluster_idx <- which(tstrsplit(rownames(ds), "_")[[1]] == cluster)
  cluster_counts <- ds[cluster_idx, ] %>% dplyr::select(c(deg, sample)) %>% mutate(class="T") 
  rest_counts <- ds[-cluster_idx, ] %>% dplyr::select(c(deg, sample)) %>% mutate(class="F") 
  train <- rbind(cluster_counts, rest_counts)
  #levels(train$class) <- c("T", "F")
  train$class <- as.factor(train$class)
  return(train)
}

#rfe-rf model
twoClassSummaryCustom <- function (data, lev = NULL, model = NULL) {
  if (length(lev) > 2) {
    stop(paste("Your outcome has", length(lev), "levels. The twoClassSummary() function isn't appropriate."))
  }
  caret:::requireNamespaceQuietStop(c("pROC","PRROC"))
  
  if (!all(levels(data[, "pred"]) == lev)) {
    stop("levels of observed and predicted data do not match")
  }
  rocObject <- try(pROC::roc(data$obs, data[, lev[1]], direction = ">", quiet = TRUE), silent = TRUE)
  rocAUC <- if (inherits(rocObject, "try-error")) 
    NA
  else rocObject$auc
  prrocObject <- try(PRROC::pr.curve(scores.class0=data[, lev[1]],scores.class1=data[, lev[2]]), silent = TRUE)
  prAUC <- if (inherits(rocObject, "try-error")) 
    NA
  else prrocObject$auc.integral
  out <- c(rocAUC, prAUC)
  names(out) <- c("ROC", "PRROC")
}

#' This function will apply recursive feature elimination (RFE) with a random forest model to find representative marker genes.
#' It would return the marker genes for each cluster.
#' Txt files of the intermediate results can be downloaded in this sub-directory.
#' @title Generate markers from DEGs
#' @importFrom tidyverse
#' @importFrom tibble
#' @importFrom stringr
#' @importFrom readxl
#' @importFrom readr 
#' @importFrom dplyr
#' @importFrom e1071
#' @importFrom caret
#' @importFrom randomForest 
#' @importFrom PRROC
#' @importFrom Metrics 
#' @importFrom data.table 
#' @param path A character of a single path to save all the results
#' @param analysis_name A character vector representing name of analyses
#' @param top_deg_num A numeric specify how many top DEGs will be selected
#' @examples
#' cell_type_markers(path="/your_path/",analysis_name="test", 
#' top_deg_num=50)
#' @export
#' 
cell_type_markers <- function(path="/your_path/", analysis_name="test", top_deg_num=50){
  caret:::caretFuncs$summary <- twoClassSummaryCustom
  
  ctrl <- rfeControl(functions=caret:::caretFuncs, 
                     method = "repeatedcv",
                     repeats =2, number = 5,
                     returnResamp="all", verbose = TRUE)
  
  trainctrl <- trainControl(classProbs= TRUE,
                            summaryFunction = twoClassSummaryCustom)
  
  ####clean up the normalized count matrix
  count.norm <- read.table(paste0(path,analysis_name, "/intermediates/", analysis_name, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", check.names = F)
  ds <- as.data.frame(as.matrix(t(count.norm)))
  ds$sample <- rownames(ds)
  all_deg <- read.table(paste0(path,analysis_name, "/intermediates/", analysis_name, "_pseudobulk_topdeg_clusters_deg.txt"), header = TRUE, sep="\t")
  
  ####markers roc
  for (cell_type in unique(all_deg$cluster)){
    top_deg <- filter(all_deg, cluster==cell_type) %>% top_n(-top_deg_num,P.Value) %>% pull(gene)
    train <- ds_train(ds, cell_type, top_deg)
    
    ####only top markers, rfe-rf model
    results <- rfe(train[,-c(top_deg_num+1, top_deg_num+2)],train[,top_deg_num+2],size=c(5:(ncol(train)-2)),
                   rfeControl=ctrl,
                   method="rf",
                   metric = c("ROC"),
                   trControl = trainctrl,
                   allowParallel = T)
    
    print(paste0(cell_type, " RFE finished"))
    save(results,file=paste0(path, analysis_name, "/celltype/", analysis_name, "_pseudobulk_", cell_type, "_deg_rfe.RData"))
    ####save feature importance
    tdf <- varImp(results) %>% dplyr::mutate(cluster=cell_type) %>% dplyr::rename(var_imp=Overall)
    tdf$gene <- rownames(tdf)
    write.table(tdf, file = paste0(path, analysis_name, "/celltype/", analysis_name, "_pseudobulk_", cell_type, "_topdeg_deg_ref.txt"), sep = "\t", quote = F)}

  df <- data.frame(cluster = character(),
                   DEG = numeric(),
                   ROC = numeric())
  
  files <- list.files(path = paste0(path, analysis_name, "/celltype"), pattern = "_deg_rfe.RData$", recursive = T, full.names = T)
  for (file in files){
    load(file)
    cluster <- strsplit(file, "_")[[1]][3]
    tdf <- results$results %>% mutate(cluster=cluster) %>% dplyr::rename("DEG" = "Variables") %>% select(cluster, DEG, ROC)
    df <- bind_rows(df, tdf)}
  write.table(df, file = paste0(path, analysis_name, "/intermediates/", analysis_name, "_topdeg_rfe_roc.txt"), sep = "\t", quote = F, row.names = F)
  
  ####get the min DEG# of the max-ROC and the summary stats
  tdf_sum <- df %>% group_by(cluster, .add=T) %>% summarise(maxROC=max(ROC), minDEG=DEG[which.max(ROC)], 
                                                            ROC5=ROC[which(DEG==5)], ROC10=ROC[which(DEG==10)], 
                                                            selected=first(sort(DEG[ROC>.9])), ROC_selected=ROC[which(DEG==selected)]) 
  
  write.table(tdf_sum, file = paste0(path, analysis_name, "/intermediates/", analysis_name, "_topdeg_rfe_roc_summary.txt"), sep = "\t", quote = F, row.names = F)

  deg <- read.csv(paste0(path, analysis_name, "/intermediates/", analysis_name, "_topdeg_rfe_roc_summary.txt"), sep="\t")
  files <- list.files(path = paste0(path, analysis_name, "/celltype"), pattern = "topdeg_deg_ref.txt$", recursive = T, full.names = T)
  df <- data.frame(analysis_name = character(),
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
  write.table(df, file = paste0(path, analysis_name, "/intermediates/", analysis_name, "_deg_rfe_varimp.txt"), sep = "\t", quote = F, row.names = F)
}

#' This function will use the selected markers to build CRM
#' It would return the final CRM.
#' Txt files of the CRM can be downloaded in this sub-directory.
#' @title Generate DEGs
#' @importFrom data.table 
#' @importFrom dplyr
#' @param analysis_name A character vector representing name of analyses
#' @param path A character of a single path to save all the results
#' @examples
#' cell_type_markers(path="/your_path/", analysis_name="test")
#' @export
#' 
####calculate gene exp across cell types
gen_crm <- function(path="/your_path/", analysis_name="test"){
  deg <- read.table(paste0(path, analysis_name, "/intermediates/", analysis_name, "_deg_rfe_varimp.txt"), sep = "\t", check.names = F, header = T)
  count.norm <- read.table(paste0(path, analysis_name, "/intermediates/", analysis_name, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", check.names = F)
  genes <- unique(deg$gene)
  clusters <- unique(deg$cluster)
  ref <- data.frame(GeneSymbol = character())
  for (type in clusters){
    idx <- which(tstrsplit(colnames(count.norm), "_")[[1]] == type) 
    exp <- count.norm[genes, idx] 
    exp <- exp %>% mutate(GeneSymbol=rownames(exp), !!type:=rowMeans(exp)) %>% select(GeneSymbol, !!type)
    ref <- full_join(ref, exp, by="GeneSymbol") %>% replace(is.na(.), 0)}
  write.table(ref, file = paste0(path, analysis_name, "/", analysis_name, "_ref_matrix.txt"), sep = "\t", quote = F, row.names = F)}

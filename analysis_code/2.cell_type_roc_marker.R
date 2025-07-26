library(e1071)
library(tibble)
library(Metrics)
library(data.table)
library(dplyr)
library(caret)
library(PRROC)
library(randomForest)
library(foreach)
library(readr)
library(doParallel)

setwd("your_path")
tissue <- "tissue"
type <- "cell_type"

#functions, prepare for different training datasets
ds_train <- function(ds, cluster, deg){
  cluster_idx <- which(tstrsplit(rownames(ds), "_")[[1]] == cluster)
  cluster_counts <- ds[cluster_idx, ] %>% dplyr::select(c(deg, sample)) %>% mutate(class="T") 
  rest_counts <- ds[-cluster_idx, ] %>% dplyr::select(c(deg, sample)) %>% mutate(class="F") 
  train <- rbind(cluster_counts, rest_counts)
  #levels(train$class) <- c("T", "F")
  train$class <- as.factor(train$class)
  return(train)}

#rfe-rf model
twoClassSummaryCustom = function (data, lev = NULL, model = NULL) {
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
  out}

caretFuncs$summary <- twoClassSummaryCustom

ctrl <- rfeControl(functions=caretFuncs, 
                   method = "repeatedcv",
                   repeats =2, number = 5,
                   returnResamp="all", verbose = TRUE)

trainctrl <- trainControl(classProbs= TRUE,
                          summaryFunction = twoClassSummaryCustom)

####clean up the normalized count matrix
count.norm <- read.table(paste0(tissue, "_pseudobulk_norm_counts_cpm_cell_type.txt"), sep = "\t", check.names = F)
ds <- count.norm
ds <- as.data.frame(as.matrix(t(count.norm)))
ds$sample <- rownames(ds)
all_deg <- read.table(paste0(tissue, "_pseudobulk_top50_clusters_deg.txt"), header = TRUE, sep="\t")

####markers roc
set.seed(1)
top_deg <- filter(all_deg, cluster==type) %>% top_n(-50,P.Value) %>% pull(gene)
train <- ds_train(ds, type, top_deg)

####only top 50 markers, rfe-rf model
cl2 <- makeCluster(10)
registerDoParallel(cl2)
print(paste0(type, " rfe started"))
results <- rfe(train[,-c(51, 52)],train[,52],size=c(5:(ncol(train)-2)),
               rfeControl=ctrl,
               method="rf",
               metric = c("ROC"),
               trControl = trainctrl,
               allowParallel = T)
               
print(paste0(type, " rfe finished"))

####save feature importance
tdf <- varImp(results) %>% dplyr::mutate(cluster=type) %>% dplyr::rename(var_imp=Overall)
tdf$gene <- rownames(tdf)
write.table(tdf, file = paste0("celltype/", tissue, "_pseudobulk_", type, "_top50_deg_ref.txt"), sep = "\t", quote = F)

####get the summarized DEG-ROC df for each tissue
df <- data.frame(type = character(),
                 cluster = character(),
                 DEG = numeric(),
                 ROC = numeric())

files <- list.files(path = paste0(tissue, "/celltype"), pattern = "_deg_rfe.RData$", recursive = T, full.names = T)
for (file in files){
  load(file)
  cluster <- strsplit(file, "_")[[1]][3]
  tdf <- results$results %>% mutate(type=tissue, cluster=cluster) %>% dplyr::rename("DEG" = "Variables") %>% select(type, cluster, DEG, ROC)
  df <- bind_rows(df, tdf)}
write.table(df, file = paste0(tissue, "/celltype/", tissue, "_50deg_rfe_roc.txt"), sep = "\t", quote = F, row.names = F)

####get the min DEG# of the max-ROC and the summary stats
tdf_sum <- df %>% group_by(cluster, .add=T) %>% summarise(maxROC=max(ROC), minDEG=DEG[which.max(ROC)], 
                                                          ROC5=ROC[which(DEG==5)], ROC10=ROC[which(DEG==10)], 
                                                          selected=first(sort(DEG[ROC>.9])), ROC_selected=ROC[which(DEG==selected)]) %>% 
  mutate(tissue=unique(tdf$type))
write.table(tdf_sum, file = paste0(tissue, "/celltype/", tissue, "_50deg_rfe_roc_summary.txt"), sep = "\t", quote = F, row.names = F)

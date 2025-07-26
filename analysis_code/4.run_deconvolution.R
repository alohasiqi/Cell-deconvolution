bulk <- "your_bulk_transcriptomics"
library(dplyr)
library(purrr)
library(MAST)
library(readr)
library(foreach)
library(doParallel)
library(Matrix.utils)
library(pheatmap)
library(e1071)
library(som)
library(methods)
library(MASS)

source("/your_code_path/deconvolution_abbCell_rlm_mass.R")
source("/your_code_path/deconvolution_rlm_mass.R")

cl <- makeCluster("the_number_of_your_ref_mtx")
registerDoParallel(cl)

models <- c("your_ref_mtx1","your_ref_mtx2","your_ref_mtx3")

foreach (model = models, .packages = c("dplyr", "purrr", "abbCell")) %dopar% {
  dir.create(file.path("/your_path/", bulk), recursive = TRUE)
  mixDataDir <- paste0("/your_path_of_bulk_transcriptomics/", bulk, ".txt")
  outputDirectory <- paste0("/your_path/", bulk)
  file <- list.files(path = paste0("/your_path_of_ref_mtx/", model), pattern = "_ref_matrix.txt$", recursive = T, full.names = T)
  #run abbcell
  abbCell_rlm(mixDataDir=mixDataDir,cellMatrixDir=file,dataDescription=bulk,projectName=model,outputDirectory=outputDirectory,annotationLegend=TRUE)
}

library(Seurat)
library(readxl)
library(edgeR)
library(dplyr)
library(tidyverse)
library(stringr)
library(readr)
library(matrixStats)
library(data.table)
library(foreach)
library(doParallel)
library(SingleCellExperiment)
library(Matrix.utils)
library(RColorBrewer)

####get deconvolution results
ds <- "your_bulk_transcriptomics"

df <- data.frame(sampleID = character(),
                 RMSE = numeric(),
                 Cor = numeric(),
                 Pvalue = numeric(),
                 model = character())

files <- list.files(path = paste0("/your_path/", ds), pattern = "_abbCell.txt$", recursive = T, full.names = T)
for(file in files){
  model_name <- strsplit(strsplit(file, split = "CellFraction_")[[1]][2], split = "_abbCell")[[1]][1] #get the ref matrix name
  f <- read.table(file, sep = "\t", header = T)
  tdf <- as.data.frame(f) %>% select(sampleID, RMSE, Cor, Pvalue) %>% mutate(model=model_name)
  df <- bind_rows(df, tdf)
}
write.table(df, file = paste0("/your_path/summary/", ds, "_mtx_abbcell_corr_summary.txt"), sep = "\t", quote = F)


####plotting
df$model <- factor(df$model, levels = c("your_ref_mtx1","your_ref_mtx2","your_ref_mtx3"))
cor_boxplot <- ggplot(df, aes(x=model, y=Cor, color=model)) +
  scale_color_manual(values=c(brewer.pal(3, "Reds")[1], brewer.pal(5, "Blues")[2])) +
  labs(title=ds, y = "Goodness-of-fit", color="Reference matrices") +
  geom_boxplot(outlier.shape = NA) +
  ylim(.1, 1) +
  theme_minimal() +
  theme(title = element_text(size = 10), plot.title = element_text(hjust = 0.5), 
        strip.text = element_text(size = 12), legend.position = "right", 
        axis.text.x = element_text(angle = 45, hjust = 1))
print(cor_boxplot)
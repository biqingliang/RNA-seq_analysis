# Check if packages exist, if not, install packages
if (!requireNameSpace("BiocManager", quietly = TRUE)
  install.packages("BiocManager")

BiocManager::install(c("DESeq2","edgeR","limma","tximport")
install.packages(c("ggplot2","pheatmap","RColorBrewer","ggrepel","dplyr")

# Load packages
library(DESeq2)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(dplyr)


counts <- read.table("counts/count_matrix.txt",
header = TRUE,
row.names = 1,
sep = "\t")

colnames(counts) <- gsub("aligned_reads/", "", colnames(counts))
colnames(counts) <- gsub("_sorted.bam", "", colnames(counts))

sample_info <- data.frame(
sample = colnames(counts),
condition = c(rep("condition1", 4), rep("condition2", 4)),
replicate = rep(1:4, 2)
)
rownames(sample_info) <- sample_info$sample

all(rownames(sample_info) == colnames(counts))

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

# Importing Count Matrix                 
counts <- read.table("count/count_matrix.txt", header = TRUE, row.names =1. sep = "\t")

# Normalizing column names
colnames(counts) <- gsub("aligned_reads.", "", colnames(counts))
colnames(counts) <- gsub("_sorted.bam", "", colnames(counts))

# Create sample data
sample_info <- data.frame(
sample = colnames(counts),
condition = c(rep("condition1", 2), rep("condition2", 2)),
replicate = rep(1:2, 2)
)
rownames(sample_info) <- sample_info$sample

#!!! Verify wheter your data order matches
all(rownames(sample_info) == colnames(counts))

source("setup.R")

# Create DESeq2 object from the count matrix
dds <- DESeqDataSetFromMatrix(
countData = counts,
colData = sample_info,
design = ~ condition
)

# Remove genes with very low counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "condition2", "condition1")) # change the field for "conditon1" and "condition2" for your sample conditions
res_ordered <- res[order(res$padj),] # order by adjusted p-value
summary(res)
write.csv(as.data.frame(res_ordered),file = "results/DESeq2_results.csv")

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts,
file = "results/normalized_counts.csv")

# Define significance thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Filter significant genes
sig_genes <- subset(res_ordered,
padj < padj_cutoff & abs(log2FoldChange) > log2fc_cutoff)

# Separate up and down-regulated
sig_up <- subset(sig_genes, log2FoldChange > log2fc_cutoff)
sig_down <- subset(sig_genes, log2FoldChange < -log2fc_cutoff)
cat(sprintf("Upregulated genes: %d\n", nrow(sig_up)))
cat(sprintf("Downregulated genes: %d\n", nrow(sig_down)))

# Save significant genes
write.csv(as.data.frame(sig_genes),
file = "results/significant_genes.csv")
write.csv(as.data.frame(sig_up),
file = "results/upregulated_genes.csv")
write.csv(as.data.frame(sig_down),
file = "results/downregulated_genes.csv")

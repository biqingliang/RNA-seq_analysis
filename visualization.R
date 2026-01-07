source("setup.R")
source("DESeq2_analysis.R")
source("chromosome.R")

# 1. PCA plot to coursely visualize difference between replicates of two conditions
# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA plot with labels
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
png("figures/PCA_labeled.png",  width = 8, height = 6, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 3) +
  geom_text_repel() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot with Sample Labels")
dev.off()


# 2. Sample distance heatmap
# Calculate sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)

# Create annotation for heatmap
annotation_col <- data.frame(Condition = sample_info$condition)
rownames(annotation_col) <- colnames(sampleDistMatrix)

# Plot heatmap
png("figures/sample_distance_heatmap.png", width = 8, height = 7, units = "in", res = 300)
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  annotation_col = annotation_col,
  main = "Sample Distance Heatmap")
dev.off()


# 3. Dispersion Plot
png("figures/dispersion_plot.png", width = 8, height = 6, units = "in", res = 300)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

#-----Visualizing differential expression-----#

# 1. Volcano Plots
# Prepare data for volcano plot
volcano_data <- as.data.frame(res_ordered)
volcano_data$gene <- rownames(volcano_data)
volcano_data$significant <- ifelse(
  volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 1,
  ifelse(volcano_data$log2FoldChange > 1, "Upregulated", "Downregulated"),
  "Not Significant"
)

# Create volcano plot
png("figures/volcano_plot.png", width = 10, height = 8, , units = "in", res = 300)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                    "Downregulated" = "blue",
                    "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(title = "Volcano Plot of Differential Expression",
  x = "Log2 Fold Change",
  y = "-Log10 Adjusted P-value",
  color = "Regulation")
dev.off()

# Volcano plot with top gene labels
top_genes <- volcano_data %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  arrange(padj) %>%
  head(20)

png("figures/volcano_plot_labeled.png", width = 10, height = 8, units = "in", res = 300)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(padj)),
      color = "black", size = 2) +
  geom_text_repel(data = top_genes,
          aes(label = gene),
          color = "black",
          size = 3,
          max.overlaps = 20) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "Volcano Plot with Top 20 Genes Labeled",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
dev.off()

# MA plot 
png("figures/MA_plot.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res, ylim = c(-5, 5), main = "MA Plot", alpha = 0.05)
dev.off()


# 2. Gene filtering and Heatmap
# Select 50 most differentially expressed genes
top_genes_ids <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]

# Get normalized counts for top genes
top_genes_counts <- normalized_counts[top_genes_ids, ]

# Scale by row (z-score)
top_genes_scaled <- t(scale(t(top_genes_counts)))

# Create annotation for samples
annotation_col <- data.frame(
  Condition = sample_info$condition
)
rownames(annotation_col) <- colnames(top_genes_scaled)

# Plot heatmap
png("figures/heatmap_screened_genes.png", width = 8, height = 12, units = "in", res = 300)
pheatmap(top_genes_scaled,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6,
  main = "Differentially Expressed Genes Heatmap")
dev.off()


# 3. Volcano plot of each chromosome
png("figures/volcano_plot_chr1.png", width = 10, height = 8, units = "in", res = 300)
chr1_data <- volcano_data[rownames(chr1_genes), ]
ggplot(chr1_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "Volcano Plot - Chromosome 1",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
dev.off()


png("figures/volcano_plot_chr2.png", width = 10, height = 8, units = "in", res = 300)
chr2_data <- volcano_data[rownames(chr2_genes), ]
ggplot(chr2_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_bw() +
  labs(title = "Volcano Plot - Chromosome 2",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")
dev.off()

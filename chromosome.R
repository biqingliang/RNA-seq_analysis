source("setup.R")
source("DESeq2_analysis.R")

# Load annotation data
annotation <- read.table("reference/gene_chromosome_map.txt", header = TRUE, sep = "\t")
# Merge annotation data with results
res_with_chr <- merge(as.data.frame(res_ordered),
annotation,
by.x = "row.names",
by.y = "gene_id",
all.x = TRUE)
rownames(res_with_chr) <- res_with_chr$Row.names
res_with_chr$Row.names <- NULL

# Assign results to each chromosomes
# make sure to change "chr1" and "chr2" field to your specific chromosome names, to see the names, simply type (grep "^>" *.fna) in your reference folder
chr1_genes <- res_with_chr[res_with_chr$chromosome == "chr1", ] 
chr2_genes <- res_with_chr[res_with_chr$chromosome == "chr2", ]

# Statistics for differential regulation between chromosomes (change "chr1" and "chr2" field as well)
chr_comparison <- data.frame(
chromosome = c("chr1", "chr2"),
n_total = c(nrow(chr1_genes), nrow(chr2_genes)),
n_sig = c(
sum(chr1_genes$padj < 0.05 & abs(chr1_genes$log2FoldChange) > 1, na.rm = TRUE),
sum(chr2_genes$padj < 0.05 & abs(chr2_genes$log2FoldChange) > 1, na.rm = TRUE)
),
mean_log2fc = c(
mean(chr1_genes$log2FoldChange, na.rm = TRUE),
mean(chr2_genes$log2FoldChange, na.rm = TRUE)
)
)
chr_comparison$percent_sig <- (chr_comparison$n_sig / chr_comparison$n_total) * 100
print(chr_comparison)
write.csv(chr_comparison, "results/chromosome_comparison.csv")


# MA plots [log2 fold change (y) vs. average log2 expression (x) for each gene] in each chromosome 
dir.create("figures", showWarnings = FALSE)
png("figures/MA_plot_by_chromosome.png", width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))

# Chromosome 1
plotMA(res_ordered[rownames(chr1_genes), ], ylim = c(-5, 5), main = "Chromosome 1", alpha = 0.05)

# Chromosome 2
plotMA(res_ordered[rownames(chr2_genes), ], ylim = c(-5, 5), main = "Chromosome 2", alpha = 0.05)
dev.off()

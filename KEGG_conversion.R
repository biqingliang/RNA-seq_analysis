# Load required libraries
library(tidyverse)

# Function to extract attribute values from GTF
extract_gtf_attribute <- function(attribute_string, attribute_name) {
  pattern <- paste0(attribute_name, ' "([^"]+)"')
  matches <- str_match(attribute_string, pattern)
  return(matches[, 2])
}

# Read GTF file
cat("Reading GTF file...\n")
gtf <- read_tsv("reference/reference.gtf", 
                comment = "#", 
                col_names = c("seqname", "source", "feature", "start", "end", 
                              "score", "strand", "frame", "attribute"),
                col_types = cols(.default = "c"))

# Filter for gene features only (to avoid duplicates)
gtf_genes <- gtf %>%
  filter(feature == "gene")

# Extract locus_tag and old_locus_tag from attributes
cat("Extracting locus tags...\n")
locus_mapping <- gtf_genes %>%
  mutate(
    locus_tag = extract_gtf_attribute(attribute, "locus_tag"),
    old_locus_tag = extract_gtf_attribute(attribute, "old_locus_tag")
  ) %>%
  select(locus_tag, old_locus_tag) %>%
  filter(!is.na(locus_tag))  # Remove any rows without locus_tag

# Display some examples of the mapping
cat("\nExample locus tag mappings:\n")
print(head(locus_mapping, 10))

# Count how many genes have old_locus_tag
n_with_old <- sum(!is.na(locus_mapping$old_locus_tag))
n_total <- nrow(locus_mapping)
cat(sprintf("\nFound %d genes with old_locus_tag out of %d total genes (%.1f%%)\n", 
            n_with_old, n_total, 100*n_with_old/n_total))

# Read DESeq2 results
cat("\nReading DESeq2 results...\n")
deseq_results <- read_csv("RNAseq_results/significant_genes.csv")

# Check the first column name and rename if needed
if (names(deseq_results)[1] == "" || names(deseq_results)[1] == "...1") {
  names(deseq_results)[1] <- "locus_tag"
}

cat(sprintf("DESeq2 results contain %d genes\n", nrow(deseq_results)))

# Merge the datasets
cat("\nMerging datasets...\n")
merged_results <- deseq_results %>%
  left_join(locus_mapping, by = "locus_tag") %>%
  select(locus_tag, old_locus_tag, everything())  # Reorder columns

# Check how many significant genes have old_locus_tag
n_sig_with_old <- sum(!is.na(merged_results$old_locus_tag))
cat(sprintf("\n%d out of %d significant genes have old_locus_tag (%.1f%%)\n", 
            n_sig_with_old, nrow(merged_results), 100*n_sig_with_old/nrow(merged_results)))

# Display first few rows
cat("\nFirst few rows of merged results:\n")
print(head(merged_results, 10))

# Save the merged results
output_file <- "significant_genes_old_locus.csv"
cat(sprintf("\nSaving merged results to: %s\n", output_file))
write_csv(merged_results, output_file)

# Optional: Create a summary of genes without old_locus_tag
genes_without_old <- merged_results %>%
  filter(is.na(old_locus_tag)) %>%
  select(locus_tag, baseMean, log2FoldChange, padj)

if (nrow(genes_without_old) > 0) {
  cat(sprintf("\nWarning: %d genes do not have old_locus_tag mappings:\n", nrow(genes_without_old)))
  print(head(genes_without_old, 5))
  
  # Save this list for reference
  write_csv(genes_without_old, "genes_without_old_locus.csv")
}

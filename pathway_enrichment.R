source("setup.R")
source("DESeq2_analysis.R")
source("chromosome.R")
source("visualization.R")

library(clusterProfile)
library(ggplot2)
library(dplyr)
library(pathview)

# Prepare dataset for KEGG analysis
data <- read.csv("*.csv") # this should be your file of significant genes with your gene locus id (old gene locus format for KEGG)
sig_gene_df <- subset(data, padj < 0.05 & abs(log2FoldChange) > 1)
genes_to_test <- sig_genes_df$old_locus_tag
all_fc <- sig_gene_df$log2FoldChange
names(all_fc) <- genes_to_test

# Perform KEGG analysis
kegg_results <- enrichKEGG(gene = genes_to_test,
                           organism = "cvi", # KEGG code for Chromobacterium violaceum
                           keyType = "kegg",
                           pvalueCutoff = 0.05)

# Function for CNET plot
directional_cnet <- function(kegg_res, fc_vector, title_text, filename) {
  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    return(message(paste("No significant pathway found for: ", title_text)))
    }
  num_pathways <- 5
  p <- cnetplot(kegg_res, foldChange = fc_vector, showCategory = num_pathways, nodel_label = "none", cex_category = 3.5, cex_gene = 0.8)
  p_final <- p + geom_node_text(aes(label = ifelse(seq_along(name) > num_pathways, name, "")),
                                size = 2.5, color = "grey30", repel = TRUE, max.overlap = 10)
               + geom_node_text(aes(label = ifelse(seq_along(name) <= num_pathways, name, "")),
                                size =5, fontface = "bold", color = "black", repel = TRUE, box.padding = 1.5, point.padding = 0.5)
               + ggtitle(title_text) 
              + theme(plot.title = element_text(hjust = 0.5, size = 16))
  ggsave(filename, plot = P_final, width =15, height = 10. dpi = 300)
  return(p_final)
 } 

if(!is.null(kegg_results)) {
  dot_all_gene <- dotplot(kegg_results, showCategory = 20) + ggtitle("Regulated KEGG pathways")
  ggsave("dotplot.png, plot = dot_all_gene, width = 8, height = 10)
  plot_directional_cnet(kegg_results, all_fc, "Regulated Gene-pathway Network", "cnet_regulated.png")

# Prepare data for regulated gene visualization in pathway
pw_data <- sig_genes_df$log2FoldChange
names(pw_data) <- sig_genes_df$olf_locus_tag

# Pathway gene regulation output
pw.out <- pathview(gene.data = pw_data, 
         pathway.id = ***, # insert your pathway of interest here
         species = "cvi", gene.idtype ="KEGG", keggnative = TRUE, 
         limit = list(gene = 2, cpd = 1), low = "lightblue", mid = "gray", high = "red")
          
                    



# RNA-seq analysis for Chromobacterium subtsugae
## 1 -  Requirement
### Files: 
 - Raw read files (fastq or fastq.gz files)
 - Reference genomes (*.fna and *.gtf or *.gtt files)
### Softwares:
  - Quality Control: FastQC and MultiQC
  - Trimming & filtering: Fastp
  - Reads alignment: HISAT2
  - Quantification: featureCounts
  - Differential analysis and plotting: R (DESeq2, edgeR, ggplot2, dplyr etc...)

## 2 - Preprocessing and Alignment 
### 2.1 Environment Setup
```bash
# Create a new environment and download everythinh you need from bioconda

conda create -n rnaseq -c bioconda -c conda-forge \
    python=3.10 \
    fastqc \
    multiqc \
    trimmomatic \
    hisat2 \
    samtools \
    subread
    
conda activate rnaseq
```
### 2.2 Review raw reads: quality control
This yields -
* Per base sequence quality (Phread > 30)
* Adapter content
* GC / sequence
* Duplication level

```bash
# Create QC directories (make sure to put all your raw reads in a raw_read directory beforehand)
mkdir -p fastqc_raw multiqc_raw

# FastQC on your raw reads
fastqc -t 8 -o fastqc_raw/ raw_reads/*.fastq.gz

# Compiling results with MultiQC
multiqc fastqc_raw/ -o multiqc_raw/
```

### 2.3 Trim and filter reads
Use fastp for sequence trimming: 
- you should have R1 and R2 reads for each sample
- and duplicates for each condition
- sample1-4 are specific names (e.g. WT-1, WT-2, LasR-1, LasR-2), only specify names in these fields
  
```bash
for sample in sample1 sample2 sample 3 sample 4; do 
fastp -i raw_reads/${sample}_R1.fastq.gz \
      -I raw_reads/${sample}_R2.fastq.gz \
      -o trimmed_reads/${sample}_R1_trimmed.fastq.gz \
      -O trimmed_reads/${sample}_R2_trimmed.fastq.gz \
      -h trummed_reads/${sample}_report.html \
      -j trummed_reads/${sample}_report.json \
      - thread 8 \
      - qualified_quality_phread_20 \
      - length_required 36
done

# Perform a post-trimming QC
mkdir -p fastqc_trimmed multiqc_trimmed
fastqc -t 8 -o fastq_trimmed/trimmed_reads/*_paired.fastq.gz
multiqc fastqc_trimmed/ -o multiqc_trimmed/
```
### 2.4 Prepare the reference genome
Download the .gtf/.gtt and .fna files for Chromobacterium subtsugae (https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/676/875/GCF_001676875.1_ASM167687v1/)
```bash
# Make reference directory, move reference files in, and check
mkdir -p reference
cd reference

# Check the number of chromosomes if you have more than one
grep "^>" *.fna

# Build the genome index for HISAT2
mkdir -p hisat2_index
hisat2-build -p 8 \
  reference/*.fna \
  hisat2_index/genome_index
```

### 2.5 Use HISAT2 to align reads
```bash
mkdir -p aligned_reads

# Read alignment
for sample in sample1 sample 2 sample 3 sample 4; do
hisat2 -p 8 \
        -dta \
        -x hisat2_index/genome_index \
        -1 trimmed_reads/${sample}_R1_paired.fastq.gz \
        -2 trimmed_reads/${sample}_R2_paired.fastq.gz \
        -S aligned_reads/${sample}.sam \
        2> aligned_reads/${sample}_alignment_stats.txt


# SAM to BAM (save space)
samtools view -@ 8 -bS aligned_reads/${sample}.sam|\
samtools sort -@ 8 -o aligned_reads/${sample}_sorted.bam

samtools index aligned_reads/${sample}_sorted.bam

rm aligned_reads/${sample}.sam

done

#------------------------------------------------------
# Post-alignment statistics
mkdir -p alignment_stats

for bam in aligned_reads/*_sorted.bam; do
  samtools flagstat $bam > alignment_stats/$(basename $bam.bam)_flagstat.txt
done
```
### 2.6 Read Quantification
```bash
mkdir -p counts
featureCounts \
  -T 8
  -p \
  -s 2 \
  -t gene \
  -g gene_id \
  -a reference/*.gtf \
  -o counts/raw_counts.txt \
  aligned_reads/*_sorted.bam

# Extract count matrix
cut -f1,7- counts/raw_counts.txt > counts/count_matrix.txt
```

## 3 - Differential Expression Analysis
### 3.1 - R packages installation and loading 

```R
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

```

### 3.2 - Import Count Matrix 
```R
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
```

### 3.3 - Fold Change Analysis, Chromosomes Assignment, and Visualization
- DESEq2_analysis.R : DESeq2 for fold change analysis
- chromosome.R : assign results to each chromosome (reference gene files needed)
- visualization.R : generates PCA, MA, Dispersion, and Volcano plots, as well as heatmaps for visualization. 


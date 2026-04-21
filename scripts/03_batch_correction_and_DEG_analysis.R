## This script merges two preprocessed HCC transcriptomic datasets (downloaded and processed in 01_HCCdata_preparation.R), 
## performs batch correction to remove inter-dataset technical variation, 
## and conducts differential expression analysis between tumor and non-tumor samples to identify DEGs and generate PCA, heatmap, and volcano plot visualizations.

# Code block 1: Batch correction

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R packages
library(GEOquery)
library(data.table)
library(dplyr)
library(tidyverse)
# install.packages("ggplot2")
library(ggplot2)
## BiocManager::install("limma")
library(limma)
## Load plotting-related packages
library(patchwork)
library(reshape2)
library(tools)
library(ggpubr)
library(RColorBrewer)

## Set working directory
getwd()
### Be sure to set your own working directory
setwd("~/Desktop/YUE_KAUST/tutorial/AI_omic/Hua_course/Hua_AI2025.09/q6：批次矫正和差异分析--复现文章Fig3/数据/")

exp_GSE25097 = read.csv("GSE25097_exp.csv", row.names = 1)
phe_GSE25097 = read.csv("GSE25097_pdata.csv", row.names = 1)
exp_GSE36376 = read.csv("GSE36376/GSE36376_exp.csv", row.names = 1)
phe_GSE36376 = read.csv("GSE36376/GSE36376_pdata.csv", row.names = 1)

#################################################################################
#### Data normalization

# Check whether log2 transformation is needed
exp = exp_GSE25097
quantiles = quantile(exp, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
needLog = (quantiles[5] > 100) || ((quantiles[6] - quantiles[1]) > 50 && quantiles[2] > 0)
if (needLog) {
  exp[exp < 0] = 0  # Set values smaller than 0 to 0
  exp = log2(exp + 1)  # Perform log2 transformation
}

# Normalize the data
normalizedData = normalizeBetweenArrays(exp)
exp_GSE25097 = normalizedData

exp = exp_GSE36376
quantiles = quantile(exp, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE)
needLog = (quantiles[5] > 100) || ((quantiles[6] - quantiles[1]) > 50 && quantiles[2] > 0)
if (needLog) {
  exp[exp < 0] = 0  # Set values smaller than 0 to 0
  exp = log2(exp + 1)  # Perform log2 transformation
}

# Normalize the data
normalizedData = normalizeBetweenArrays(exp)
exp_GSE36376 = normalizedData

#############################################################
## Merge the two datasets

gene_inter = intersect(rownames(exp_GSE25097), rownames(exp_GSE36376))
exp_GSE25097 = as.data.frame(exp_GSE25097)
exp_GSE36376 = as.data.frame(exp_GSE36376)
merge_data = merge(x = exp_GSE25097[gene_inter, ], y = exp_GSE36376[gene_inter, ], by = 0)
rownames(merge_data) = merge_data$Row.names
merge_data = merge_data[, -1]

##################################################################
## Remove batch effects

# BiocManager::install("sva")
library(sva)

# Perform batch effect correction on the merged data
?ComBat
batchType = c(rep("GSE25097", 557), rep("GSE36376", 433))
batchType <- as.factor(batchType)
merge_data1 = na.omit(merge_data)
correctedData = ComBat(merge_data1, batchType, par.prior = TRUE)

getwd()
setwd("~/Desktop/YUE_KAUST/tutorial/AI_omic/Hua_course/Hua_AI2025.09/q6：批次矫正和差异分析--复现文章Fig3/数据/")
write.csv(correctedData, file = "correctedData.csv")

huage_pca_plot <- function(mat, batch_info, title, file = NULL) {
  data_pca <- prcomp(t(mat), scale. = TRUE)
  pc_var <- data_pca$sdev^2 / sum(data_pca$sdev^2)  # Proportion of explained variance
  
  pca_df <- data.frame(
    Sample = rownames(data_pca$x),
    PC1 = data_pca$x[, 1],
    PC2 = data_pca$x[, 2],
    Batch = factor(batch_info, levels = unique(batch_info))
  )
  
  n_batch <- length(unique(batch_info))
  mycol <- c("#e1615e", "#26a3a3")
  if (n_batch > 8) mycol <- colorRampPalette(brewer.pal(9, "Set1"))(n_batch)
  
  pc1_lab <- paste0("PC1 (", round(pc_var[1] * 100, 1), "%)")
  pc2_lab <- paste0("PC2 (", round(pc_var[2] * 100, 1), "%)")
  
  p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 2, alpha = 0.90) +
    stat_ellipse(aes(fill = Batch), geom = "polygon", alpha = 0.14, linetype = 2, show.legend = FALSE) +
    scale_color_manual(values = mycol) +
    scale_fill_manual(values = mycol) +
    labs(title = title, x = pc1_lab, y = pc2_lab, color = "Batch") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_line(colour = "gray90", linetype = 2),
      legend.title = element_text(face = "bold"),
      legend.background = element_rect(colour = "black", fill = NA, size = 0.25),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(face = "bold")
    )
  
  if (!is.null(file)) {
    ggsave(file, p, width = 8, height = 6)
  }
  
  return(p)
}

# Export PCA plots as PDF
# Prepare data for PCA analysis
mode(correctedData) <- "numeric"  # Make sure it is a numeric matrix

huage_pca_plot(
  merge_data1,
  batchType,
  "PCA_Before Batch Effect",
  file.path("PCA_Before_Batch.pdf")
)

huage_pca_plot(
  correctedData,
  batchType,
  "PCA_After Batch Effect",
  file.path("PCA_After_Batch.pdf")
)


# Code block 2: Differential expression analysis and visualization

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R packages
library(GEOquery)
library(data.table)
library(dplyr)
library(tidyverse)
# install.packages("ggplot2")
library(ggplot2)
## BiocManager::install("limma")
library(limma)
## Load plotting-related packages
library(patchwork)
library(reshape2)
library(tools)
library(ggpubr)
library(RColorBrewer)

## Set working directory
getwd()
### Be sure to set your own working directory
setwd("~/Desktop/YUE_KAUST/tutorial/AI_omic/Hua_course/Hua_AI2025.09/q6：批次矫正和差异分析--复现文章Fig3/数据/")

phe_GSE25097 = read.csv("GSE25097_pdata.csv", row.names = 1)
table(phe_GSE25097$source_name_ch1)
s1_n = filter(phe_GSE25097, source_name_ch1 == "non_tumor")
s1_t = filter(phe_GSE25097, source_name_ch1 == "tumor")

phe_GSE36376 = read.csv("GSE36376/GSE36376_pdata.csv", row.names = 1)
table(phe_GSE36376$tissue.ch1)
s2_n = filter(phe_GSE36376, tissue.ch1 == "adjacent non-tumor liver")
s2_t = filter(phe_GSE36376, tissue.ch1 == "liver tumor")

correctedData = read.csv("correctedData.csv", row.names = 1)

### Reorder samples so that all adjacent non-tumor samples come first
### and all tumor samples come afterward
length(c(rownames(s1_n), rownames(s2_n)))
length(c(rownames(s1_t), rownames(s2_t)))

col_index = c(rownames(s1_n), rownames(s2_n), rownames(s1_t), rownames(s2_t))
correctedData = correctedData[, col_index]

############################################################################
############################################################################
############################################################################

## Create group information:
## the first 436 samples are controls, and the next 508 are tumor samples
length(c(rownames(s1_n), rownames(s2_n)))
length(c(rownames(s1_t), rownames(s2_t)))

group_list = c(rep("Nor", 436), rep("Tumor", 508))

## Convert to factor
## You do not need to fully understand this now;
## it is used to construct the design matrix for differential analysis
## Explicitly enforce the group order
group_list <- factor(group_list, levels = c("Nor", "Tumor"))

## Construct the design matrix for differential expression analysis
design = model.matrix(~ group_list)
colnames(design) <- c("Nor", "Tumor")
rownames(design) <- colnames(correctedData)

## lmFit(): build a linear model
fit = lmFit(correctedData, design)  # Fit one linear model for each gene

## eBayes(): apply empirical Bayes moderation to stabilize standard errors
fit2 <- eBayes(fit)  # Empirical Bayes correction makes the statistics more stable

## topTable(): return a ranked list of genes that are most likely differentially expressed
allDiff = topTable(fit2, coef = 2, adjust.method = "fdr", number = Inf)
write.csv(allDiff, file = "allDiff.csv")

############################################################################
############################################################################
############################################################################

### Differential gene expression analysis
# Transcriptomic data were analyzed using the limma package.
# Differentially expressed genes (DEGs) were identified with thresholds of
# FDR-adjusted P < 0.05 and |log2FC| > 0.585 (1.5-fold change).
# Results were visualized using ggplot2.

select.log2FC <- abs(allDiff$logFC) > 0.585
table(select.log2FC)

select.qval <- allDiff$adj.P.Val < 0.05
table(select.qval)

select.vec = select.log2FC & select.qval
table(select.vec)

deg = allDiff[select.vec, ]
write.csv(deg, file = "deg.csv")

############################################################################
############################################################################
############################################################################

### Plot heatmap

ordered_DEGs <- allDiff[order(as.numeric(as.vector(allDiff$logFC))), ]
ordered_gene_names <- rownames(ordered_DEGs)
total_DEG_count <- length(ordered_gene_names)

### The paper displays the top 25 upregulated and top 25 downregulated genes
display_genes_num = 25

selected_gene_set <- ordered_gene_names[c(
  1:display_genes_num,
  (total_DEG_count - display_genes_num + 1):total_DEG_count
)]

heatmap_exp <- correctedData[selected_gene_set, ]
col_index = c(rownames(s1_n), rownames(s2_n), rownames(s1_t), rownames(s2_t))

## Draw heatmap
length(c(rownames(s1_n), rownames(s2_n)))
length(c(rownames(s1_t), rownames(s2_t)))

annotation_col1 = data.frame(
  Type = c(rep("Nor", 436), rep("Tumor", 508)),
  Project = c(
    rep("GSE25097", 243),
    rep("GSE36376", 193),
    rep("GSE25097", 268),
    rep("GSE36376", 240)
  )
)

rownames(annotation_col1) = colnames(heatmap_exp)

# install.packages("pheatmap")
library(pheatmap)

pheatmap::pheatmap(
  heatmap_exp,                 # Data used for the heatmap
  cluster_rows = TRUE,         # Cluster rows
  cluster_cols = FALSE,        # Do not cluster columns, so sample separation can be seen directly
  annotation_col = annotation_col1,
  show_colnames = FALSE,
  scale = "row",               # Row-wise standardization
  color = colorRampPalette(c("blue", "white", "red"))(100)
)

############################################################################
############################################################################
############################################################################

### Plot volcano plot

library(ggpubr)
library(ggplot2)
library(ggrepel)

data <- allDiff
data$significant = "stable"
data$significant[data$logFC >= 0.585 & data$adj.P.Val < 0.05] = "up"

# data$significant[ ] = "up"
x = data$logFC >= 1.5 & data$adj.P.Val < 0.05
table(x)
x[1:20]

data$significant[data$logFC <= -0.585 & data$adj.P.Val < 0.05] = "down"
table(data$significant)

log10(0.001)
log10(0.00001)
# ggplot(data, aes(x = logFC, y = -1 * log10(adj.P.Val)))
log10(0.0001)
log10(0.05)
colnames(data)

ggplot(data, aes(x = logFC, y = -1 * log10(adj.P.Val))) +
  xlim(-4, 3) +
  ylim(0, 310) +
  geom_point(aes(color = significant), size = 0.8) +
  theme_classic() +
  scale_color_manual(values = c("green", "grey", "red")) +
  geom_hline(yintercept = 1.3, linetype = 4, size = 0.8) +
  geom_vline(xintercept = c(-0.585, 0.585), linetype = 4, size = 0.8) +
  theme(title = element_text(size = 12), text = element_text(size = 12)) +
  labs(x = "log2 fold change", y = "-log10(p_value)") +
  theme_bw()

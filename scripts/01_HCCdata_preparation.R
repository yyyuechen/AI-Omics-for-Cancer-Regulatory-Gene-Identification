## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R packages
# BiocManager::install('seandavi/GEOquery')
library(GEOquery)
library(dplyr)
library(tidyverse)

##download datasets from https://www.ncbi.nlm.nih.gov/geo/
# GSE25097
# GSE36376

## Set working directory
getwd()

### Make sure to set your own working directory
setwd("~/path/to/your/directory/data/")

# ***Extract expression matrix***

gset = getGEO(filename = 'GSE25097_series_matrix.txt.gz', destdir = ".", getGPL = FALSE)

save(gset, file = "GSE25097.rdata")

### Open the working directory and you will see that GSE25097.rdata has been created
## Save the downloaded data as an rdata file; it can also be saved as an rds file
########### At this point, you can clear the environment and reload the saved data using load()
load("GSE25097.rdata")

############ After loading, the gset object will appear in the environment again

exp.s4 = gset[[1]]
exp.s4 = gset

# ***0 3 In-depth interpretation of the S4 object***

## This code extracts the expression matrix from an ExpressionSet-like object
## and converts it into a data frame
exp = exp.s4@assayData[["exprs"]]
df.exp = as.data.frame(exp)
class(exp)
class(df.exp)

## Save phenotype/group information
pdata = exp.s4@phenoData@data
write.csv(pdata, file = "GSE25097_pdata.csv")
pdata1 = read.csv("GSE25097_pdata.csv")

# ***0 4 Convert probe IDs to gene symbols***

# install.packages("data.table")
library(data.table)

anno = fread("GPL10687.txt", sep = "\t", header = TRUE, data.table = FALSE) # Microarray platform annotation file

# Extract key information from the annotation file
colnames(anno)
gpl_gene = anno[, c(1, 3)]
# gpl_gene <- anno[, c("ID", "Gene Symbol")]
colnames(gpl_gene)

# Merge the expression matrix with the platform annotation information

### exp is a matrix, but merge() requires data frames
################################################ Demonstration of merge() usage
### Perform probe annotation conversion
View(df.exp)
?merge
exp.anno = merge(x = gpl_gene, y = df.exp, by.x = 1, by.y = 0)
x = exp.anno$GeneSymbol
x[1:50]

### At this point, both the gene symbol annotation file and the expression matrix are ready
############################################################################
### Organize the expression matrix
exp1 = exp.anno
colnames(exp1)
rownames(exp1) = exp1$GeneSymbol

### The error occurs because there are duplicated gene names:
### multiple probes may correspond to the same gene
exp2.1 = distinct(exp1, GeneSymbol, .keep_all = FALSE)
exp2 = distinct(exp1, GeneSymbol, .keep_all = TRUE)

### At this point, the row count changes from 37582 to 18077
rownames(exp2) = exp2$GeneSymbol
View(exp2)

### Remove the first and second columns
exp3 = exp2[, -c(1, 2)]
write.csv(exp3, file = "GSE25097_exp.csv")

### do the same for GSE36376

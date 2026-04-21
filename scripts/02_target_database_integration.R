## download data about AFB1 from following 2 datasets:
# ChEMBL: targets that have already been reported in the literature or existing databases
# SwissTargetPrediction: predicts potential targets based on molecular similarity
# PharmMapper: predicts potential targets based on 3D pharmacophore matching

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R package
library(data.table)

setwd("~/path/to/your/directory/data/")

# Read input data
chembl_gene <- fread("chembl/DOWNLOAD_gene.tsv", sep = "\t", header = TRUE)
ref <- fread("chembl/chembl_ref.tsv", sep = "\t", header = TRUE)

x1 = strsplit(chembl_gene$`UniProt Accessions`, split = "|", fixed = TRUE)

## https://github.com/Aleksobrad/single-cell-rcc-pipeline/blob/master/singlecell_gex_viper_analysis.R
## Code adapted from lines 941 to 946
# Convert ChEMBL target UniProt IDs into gene names, then export the final gene list
f1 = x1[[1]][1]

for (i in 2:72) {
  f2 = x1[[i]][1]
  f1 = c(f1, f2)
}

f1

# This code keeps only the first UniProt ID from each row
# If a ChEMBL entry contains multiple accessions, all IDs after the first one are discarded

chembl_gene$gene.all = f1
# Add the extracted UniProt ID back into the table
# Here, gene.all actually stores only the first UniProt accession from each row

inters = intersect(f1, ref$Entry)
# Find the overlap between ChEMBL UniProt IDs and the reference table

result <- merge(chembl_gene, ref,
                by.x = "gene.all", by.y = "Entry",
                all.x = TRUE, sort = FALSE)

gene_last = result$`Gene Names`
gene_last
gene_last = unique(gene_last)
gene_last
gene_last = as.character(na.omit(gene_last))
gene_last
write.csv(gene_last, file = "chembl_gene_last.csv")

## "DDR1 /// MIR4640" is a string
### String splitting / string processing
### There are dedicated tutorials on string handling

x3 = strsplit("RPL21 /// RPL21P28 /// SNORA27 /// SNORD102", split = " /// ", fixed = TRUE)

for (i in 1:10) {
  k1 = i + 5
  k2 = k1 + 10
  k3 = k2 - 3
  print(k3)
}

i = 1
k3 = 13

i = 2
k3 = 14

...
i = 10
k3 = 22

################################################################################################
################################################################################################
################################################################################################
## Code block 2: combine gene sets from three databases

## Set system messages to English
Sys.setenv(LANGUAGE = "en")

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R package
library(data.table)

setwd("~/Desktop/YUE_KAUST/tutorial/AI_omic/Hua_course/Hua_AI2025.09/q5.数据库整理--文章Fig2复现/数据/")

# Read ChEMBL gene list
ChEMBL_gene = read.csv("chembl/chembl_gene_last.csv")
ChEMBL_gene = ChEMBL_gene[, 2]

# Remove leading and trailing spaces
ChEMBL_gene <- trimws(ChEMBL_gene)

# Remove empty strings
ChEMBL_gene <- ChEMBL_gene[ChEMBL_gene != ""]

# Remove duplicated genes
ChEMBL_gene <- unique(ChEMBL_gene)

################################################
# Read PharmMapper gene list
PharmMapper_gene <- fread("pharmmapper/pharmmapper_gene.txt", sep = "\t", header = FALSE)
PharmMapper_gene = PharmMapper_gene$V1
PharmMapper_gene <- trimws(PharmMapper_gene)
PharmMapper_gene <- PharmMapper_gene[PharmMapper_gene != ""]
PharmMapper_gene <- unique(PharmMapper_gene)

################################################
# Read SwissTargetPrediction gene list
Swiss_gene <- fread("Swiss/swiss_gen.txt", sep = "\t", header = FALSE)
Swiss_gene = Swiss_gene$V1
Swiss_gene <- trimws(Swiss_gene)
Swiss_gene <- Swiss_gene[Swiss_gene != ""]
Swiss_gene <- unique(Swiss_gene)

################################################
################################################
### Combine genes from the three databases into a list
geneList = list(ChEMBL = ChEMBL_gene,
                PharmMapper = PharmMapper_gene,
                Swiss = Swiss_gene)

# Draw a Venn diagram
# install.packages("ggvenn")
library(ggvenn)

ggvenn(geneList, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#e28764", "#728fc6", "#f3b670")[seq_along(geneList)],
       set_name_color = c("#e28764", "#728fc6", "#f3b670")[seq_along(geneList)],
       set_name_size = 4, text_size = 3.5)

# Calculate the union of all gene names from the three databases
unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)
unionCount

# Save the result
write.csv(unionGenes, file = "unionGenes.csv")

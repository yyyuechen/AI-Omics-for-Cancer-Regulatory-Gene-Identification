## This script integrates HCC-related candidate genes and AFB1-related target genes, 
## constructs a PPI network for the selected gene set, and performs GO/KEGG enrichment analysis to characterize the functional pathways associated with the candidate genes.
## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required package
library(data.table)

setwd("~/path/to/your/directory/data/")

# Read HCC candidate genes
HCC_gene = read.csv("all.unionGenes.csv")
# all.unionGenes.csv = DEG genes + brown-module genes from WGCNA
HCC_gene = HCC_gene[,2]

# Remove duplicated genes
HCC_gene <- unique(HCC_gene)

################################################
# Read AFB1-related genes
AFB1_gene <- read.csv("~/path/to/your/directory/data/unionGenes.csv")
AFB1_gene = AFB1_gene[,2]  # Gene set collected from the three target-prediction databases
AFB1_gene <- unique(AFB1_gene)

################################################
################################################
### Organize the two gene sets into geneList
geneList = list(AFB1_gene = AFB1_gene,
                HCC_gene = HCC_gene)

# Draw Venn diagram
# install.packages("ggvenn")
library(ggvenn)

ggvenn(geneList, show_percentage = TRUE,
       stroke_color = "white", stroke_size = 0.5,
       fill_color = c("#e28764", "#728fc6", "#f3b670")[seq_along(geneList)],
       set_name_color = c("#e28764", "#728fc6", "#f3b670")[seq_along(geneList)],
       set_name_size = 4, text_size = 3.5)

last = c("RND3", "AURKA", "PCK1", "BCAT2", "UCK2", "CCNB1")

intersect(HCC_gene, last)
intersect(last, AFB1_gene)

# These six genes are the final core genes reported in the paper.
# They are included here to show that even though the current Venn analysis
# does not produce exactly the same gene set as the paper,
# the final key genes were not missed.

# Check the overlap between HCC-related and AFB1-related genes
intersect(HCC_gene, AFB1_gene)

# Calculate the union of all gene names across the input gene sets
unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)
unionCount

# Save the result

setwd("~/path/to/your/directory/data/")

wenzhang.gene = read.table("genes_list.txt")
# genes_list.txt contains the genes shown in the published PPI figure.
# Because the original paper did not provide the code,
# the exact same gene list cannot be fully reproduced from the previous step.
# Therefore, to keep the downstream machine-learning results consistent with the paper,
# the genes appearing in the published PPI figure are directly included here.

wenzhang.gene = wenzhang.gene$V1

all.gene = union(wenzhang.gene, unique(c(last, intersect(HCC_gene, AFB1_gene))))
write.csv(all.gene, file = "all.genes.csv")


# PPI network

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
## Load plotting packages
library(patchwork)
library(reshape2)
library(tools)
library(ggpubr)
library(RColorBrewer)

## Set working directory
getwd()

### Be sure to set your own working directory
setwd("~/path/to/your/directory/data/")

deg = read.csv("~/path/to/your/directory/data/allDiff.csv", row.names = 1)
gene = read.csv("all.genes.csv", header = TRUE)

deg = deg[gene$x, ]
deg$group = ifelse(deg$logFC > 0, "up", "down")

#####################################################################################

## https://string-db.org/cgi/network?pollingId=bUpew6ZPM4MT&sessionId=bq3ONpHVnJ16

library(ggraph)
library(igraph)

nodes = data.frame(V1 = rownames(deg), Name = rownames(deg), group = deg$group)
links <- read.table("string_interactions_short.tsv", header = FALSE, sep = "\t")

nodes
links
net <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)

while (!is.null(dev.list())) dev.off()

par(mar = c(1, 1, 1, 1))
plot(net,
     vertex.size = 8,
     vertex.label.cex = 0.7,
     edge.arrow.size = 0.3)

?ggraph
ggraph(net, layout = "linear", circular = TRUE) +
  geom_edge_link(color = "blue") +
  geom_node_text(aes(label = Name)) +
  geom_node_point(size = 10, color = "red", alpha = 0.5) +
  theme_void()

## layout = circlepack
ggraph(net, layout = "treemap") +
  geom_edge_link(color = "blue") +
  geom_node_text(aes(label = Name)) +
  geom_node_point(size = 10, color = "red", alpha = 0.5) +
  theme_void()

ggraph(net, layout = "kk") +
  geom_edge_link(color = "blue") +
  geom_node_text(aes(label = Name)) +
  geom_node_point(size = 10, color = "red", alpha = 0.5) +
  theme_void()

# After comparing several layouts, layout = "kk" looks better.

## Next, color the nodes according to upregulated and downregulated groups
nodes
links
nodes$Name <- nodes$V1
nodes$Group <- nodes$group
net <- graph_from_data_frame(d = links, vertices = nodes, directed = TRUE)

ggraph(net, layout = "kk") +
  geom_edge_link() +
  geom_node_point(size = 10, aes(color = Group)) +
  geom_node_text(aes(label = Name)) +
  theme_void()

# Draw the network
ggraph(net, layout = "circle") +
  geom_edge_link(color = "gray") +
  geom_node_point(size = 10, aes(color = Group)) +
  geom_node_text(aes(label = Name), repel = TRUE, size = 3) +
  scale_color_manual(values = c("green", "red")) +
  theme_void()

# Assume that links and nodes data frames are already available
# links: from, to
# nodes: Name, Group

# Build the network object
net <- graph_from_data_frame(d = links, vertices = nodes, directed = FALSE)

# Group information
V(net)$Group <- nodes$Group

# Create a two-ring layout
layout_circle2 <- function(graph) {
  vnames <- V(graph)$name
  groups <- V(graph)$Group
  
  # Set radii: outer ring = 1.0, inner ring = 0.5
  radius_inner <- 0.5
  radius_outer <- 1.0
  
  # Get groups
  group_levels <- unique(groups)
  
  # Nodes in each group
  group1_nodes <- vnames[groups == group_levels[1]]
  group2_nodes <- vnames[groups == group_levels[2]]
  
  n1 <- length(group1_nodes)
  n2 <- length(group2_nodes)
  
  # Coordinates (convert polar coordinates to Cartesian coordinates)
  angle1 <- seq(0, 2*pi, length.out = n1 + 1)[-(n1 + 1)]
  angle2 <- seq(0, 2*pi, length.out = n2 + 1)[-(n2 + 1)]
  
  layout <- matrix(NA, nrow = length(vnames), ncol = 2)
  rownames(layout) <- vnames
  
  layout[group1_nodes, 1] <- radius_outer * cos(angle1)
  layout[group1_nodes, 2] <- radius_outer * sin(angle1)
  
  layout[group2_nodes, 1] <- radius_inner * cos(angle2)
  layout[group2_nodes, 2] <- radius_inner * sin(angle2)
  
  return(layout)
}

# Get manual layout
my_layout <- layout_circle2(net)

# Plot the network
ggraph(net, layout = "manual", x = my_layout[,1], y = my_layout[,2]) +
  geom_edge_link(color = "gray", alpha = 0.6) +
  geom_node_point(aes(color = Group), size = 10) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  scale_color_manual(values = c("green", "red")) +
  theme_void()

ggraph(net, layout = "manual", x = my_layout[,1], y = my_layout[,2]) +
  geom_edge_link(color = "gray", alpha = 0.6) +
  geom_node_point(aes(color = Group), size = 10) +
  geom_node_text(
    aes(label = name),
    repel = FALSE,
    size = 2,
    vjust = 0.5,
    hjust = 0.5,
    color = "black"
  ) +
  scale_color_manual(values = c("green", "red")) +
  theme_void()

ggraph(net, layout = "manual", x = my_layout[,1], y = my_layout[,2]) +
  geom_edge_link(color = "gray", alpha = 0.6) +
  geom_node_point(aes(color = Group), size = 10) +
  geom_node_text(
    aes(label = name),
    repel = FALSE,
    size = 4,
    vjust = 0.5,
    hjust = 0.5,
    color = "black"
  ) +
  scale_color_manual(values = c("green", "red")) +
  theme_void()

ggraph(net, layout = "manual", x = my_layout[,1], y = my_layout[,2]) +
  geom_edge_link(color = "black", alpha = 0.6) +
  geom_node_point(aes(color = Group), size = 10) +
  geom_node_text(
    aes(label = name),
    repel = FALSE,
    size = 4,
    vjust = 0.5,
    hjust = 0.5,
    color = "black"
  ) +
  scale_color_manual(values = c("green", "red")) +
  theme_void()

## Ask GPT to help center the text labels
ggraph(net, layout = "manual", x = my_layout[,1], y = my_layout[,2]) +
  geom_edge_link(color = "gray", alpha = 0.6) +
  geom_node_point(aes(color = Group), size = 10) +
  geom_node_text(aes(label = name), size = 3, color = "black", hjust = 0.5, vjust = 0.5) +
  scale_color_manual(values = c("down" = "green", "up" = "red")) +
  theme_void()


# Enrichment analysis

## Set system messages to English
Sys.setenv(LANGUAGE = "en")

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required R packages
library(GEOquery)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

## Set working directory
getwd()

### Be sure to set your own working directory
setwd("~/path/to/your/directory/data/")

gene = read.csv("~/path/to/your/directory/data/all.genes.csv", row.names = 1)

degs.list = gene$x
degs.list

erich.go.BP = enrichGO(gene = degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
dotplot(erich.go.BP, showCategory = 15)
?dotplot
barplot(erich.go.BP, showCategory = 10)
erich.go.BP = erich.go.BP@result
write.table(erich.go.BP, "erich.go.BP.txt", sep = "\t", col.names = NA)
write.csv(erich.go.BP, "erich.go.BP.csv")

erich.go.CC = enrichGO(gene = degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

p1 = dotplot(erich.go.CC, showCategory = 10)
p2 = barplot(erich.go.CC, showCategory = 8)
p1 + p2

erich.go.MF = enrichGO(gene = degs.list,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
p3 = dotplot(erich.go.MF, showCategory = 10)
p4 = barplot(erich.go.MF, showCategory = 8)
p3 + p4

####################################################

erich.go.all = enrichGO(gene = degs.list,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        pvalueCutoff = 0.05,
                        qvalueCutoff = 0.05,
                        ont = "all")

## ont = "all": calculate GO enrichment for BP, CC, and MF at the same time

# Convert results into a data frame
erich.go.all.df <- as.data.frame(erich.go.all@result)

barplot(erich.go.all, drop = TRUE, showCategory = 10,
        split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = "free")

barplot(erich.go.all, drop = TRUE, showCategory = 10, split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = "free") +
  theme(
    axis.text.y = element_text(size = 8, angle = 0, hjust = 1),
    text = element_text(lineheight = 0.8)
  )

####################################################
####################################################
####################################################

kegg = read.gmt("c2.cp.kegg_medicus.v2025.1.Hs.symbols.gmt")
length(unique(kegg$term))

kegg_res = enricher(gene = degs.list, TERM2GENE = kegg,
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)
dotplot(kegg_res, showCategory = 20)

dotplot(kegg_res, showCategory = 20) +
  theme(
    axis.text.y = element_text(size = 8, angle = 0, hjust = 1),
    text = element_text(lineheight = 0.8)
  )

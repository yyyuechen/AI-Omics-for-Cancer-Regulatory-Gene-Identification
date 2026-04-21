# 04_WGCNA_module_analysis_and_DEG_integration.R

# This script performs WGCNA on the batch-corrected HCC expression matrix to identify co-expression modules, 
# evaluate module-trait relationships, extract module genes, and integrate key modules with DEG results to prioritize tumor-associated candidate genes.

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required packages
library(WGCNA)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(grid)
library(reshape2)
library(RColorBrewer)

## Allow multi-threading in WGCNA if supported
allowWGCNAThreads()

## Set working directories
getwd()
setwd("~/path/to/your/directory/data/")

###############################################################################
## 1. Read batch-corrected expression matrix and sample metadata
###############################################################################

correctedData <- read.csv("/path/to/your/directory/data/correctedData.csv", row.names = 1)

phe_GSE25097 <- read.csv("/path/to/your/directory/data/GSE25097_pdata.csv", row.names = 1)
phe_GSE36376 <- read.csv("/path/to/your/directory/data/GSE36376/GSE36376_pdata.csv", row.names = 1)

## Identify normal and tumor samples in each dataset
s1_n <- filter(phe_GSE25097, source_name_ch1 == "non_tumor")
s1_t <- filter(phe_GSE25097, source_name_ch1 == "tumor")

s2_n <- filter(phe_GSE36376, tissue.ch1 == "adjacent non-tumor liver")
s2_t <- filter(phe_GSE36376, tissue.ch1 == "liver tumor")

## Reorder samples: all normal samples first, then all tumor samples
col_index <- c(rownames(s1_n), rownames(s2_n), rownames(s1_t), rownames(s2_t))
correctedData1 <- correctedData[, col_index]

## Create trait matrix
group_list <- c(rep("Nor", length(c(rownames(s1_n), rownames(s2_n)))),
                rep("Tumor", length(c(rownames(s1_t), rownames(s2_t)))))

design <- model.matrix(~0 + group_list)
colnames(design) <- c("Nor", "Tumor")
rownames(design) <- colnames(correctedData1)

###############################################################################
## 2. Filter genes and prepare expression matrix for WGCNA
###############################################################################

## Filter out low-variance genes
filteredData <- correctedData1[apply(correctedData1, 1, sd, na.rm = TRUE) > 0.5, ]

## Transpose matrix for WGCNA:
## rows = samples, columns = genes
exprMatrix <- t(filteredData)

## Check for problematic samples or genes
gsg <- goodSamplesGenes(exprMatrix, verbose = 3)
if (!gsg$allOK) {
  exprMatrix <- exprMatrix[gsg$goodSamples, gsg$goodGenes]
}

###############################################################################
## 3. Sample clustering for outlier detection
###############################################################################

sampleTree <- hclust(dist(exprMatrix), method = "average")

pdf("Sample_Clustering.pdf", width = 8, height = 6)
par(cex = 0.8)
par(mar = c(4, 4, 2, 1))
plot(sampleTree,
     main = "Sample Clustering for Outlier Detection",
     sub = "",
     xlab = "",
     cex.lab = 1.2,
     cex.axis = 1.0,
     cex.main = 1.5)
abline(h = 40, col = "red", lwd = 2)
dev.off()

## Plot dendrogram with trait heatmap
traitColors <- numbers2colors(design, signed = FALSE)
pdf("Sample_Dendrogram_TraitHeatmap.pdf", width = 8, height = 6)
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = colnames(design),
                    main = "Sample Dendrogram and Trait Heatmap",
                    dendroLabels = FALSE)
dev.off()

###############################################################################
## 4. Pick soft-thresholding power
###############################################################################

powers <- 1:20
sft <- pickSoftThreshold(exprMatrix, powerVector = powers, verbose = 5)

pdf("SoftThreshold_Selection.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     type = "n",
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale independence")
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers,
     col = "red")
abline(h = 0.90, col = "red")

plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     type = "n",
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = "Mean connectivity")
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     col = "red")
dev.off()

softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 6
}

###############################################################################
## 5. Construct network and detect modules
###############################################################################

adjacency <- adjacency(exprMatrix, power = softPower)
tom_sim <- TOMsimilarity(adjacency)
tom_dis <- 1 - tom_sim

geneTree <- hclust(as.dist(tom_dis), method = "average")

pdf("Gene_Clustering_TOM_Dissimilarity.pdf", width = 8, height = 6)
plot(geneTree,
     xlab = "",
     sub = "",
     main = "Gene Clustering on TOM-based Dissimilarity",
     labels = FALSE,
     hang = 0.04)
dev.off()

## Dynamic tree cut
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = tom_dis,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)

dynamicColors <- labels2colors(dynamicMods)

pdf("Dynamic_Tree_Cut_Modules.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05,
                    main = "Gene Dendrogram and Module Colors")
dev.off()

## TOM heatmap
plot_sim <- -(1 - tom_sim)
diag(plot_sim) <- NA
pdf("Network_Heatmap_TOM.pdf", width = 7, height = 7)
TOMplot(plot_sim, geneTree, dynamicColors,
        main = "Network Heatmap Plot")
dev.off()

###############################################################################
## 6. Calculate module eigengenes and merge similar modules
###############################################################################

MEList <- moduleEigengenes(exprMatrix, colors = dynamicColors)
MEs <- MEList$eigengenes

ME_cor <- cor(MEs)

pdf("Module_Eigengene_Clustering.pdf", width = 8, height = 6)
METree <- hclust(as.dist(1 - ME_cor), method = "average")
plot(METree, main = "Clustering of Module Eigengenes", xlab = "", sub = "")
abline(h = 0.4, col = "blue")
dev.off()

merge_module <- mergeCloseModules(exprMatrix,
                                  dynamicColors,
                                  cutHeight = 0.4,
                                  verbose = 3)

mergedColors <- merge_module$colors
mergedMEs <- merge_module$newMEs

pdf("Merged_Modules.pdf", width = 8, height = 6)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged Modules"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
dev.off()

## Update final module assignments
moduleColors <- mergedColors
moduleEigengenes <- mergedMEs

###############################################################################
## 7. Module–trait relationship analysis
###############################################################################

moduleTraitCor <- cor(moduleEigengenes, design, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(exprMatrix))

write.table(moduleTraitCor, "moduleTraitCor.txt", sep = "\t", col.names = NA, quote = FALSE)
write.table(moduleTraitPvalue, "moduleTraitPvalue.txt", sep = "\t", col.names = NA, quote = FALSE)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

pdf("Module_Trait_Relationships.pdf", width = 8, height = 8)
par(mar = c(5, 10, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(design),
               yLabels = colnames(moduleEigengenes),
               ySymbols = colnames(moduleEigengenes),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               cex.text = 0.7,
               zlim = c(-1, 1),
               main = "Module-Trait Relationships")
dev.off()

## Enhanced ggplot heatmap
corDF <- melt(moduleTraitCor)
pvalDF <- melt(moduleTraitPvalue)
heatmapDF <- merge(corDF, pvalDF, by = c("Var1", "Var2"))
colnames(heatmapDF) <- c("Module", "Trait", "Correlation", "Pvalue")

p_trait <- ggplot(heatmapDF, aes(x = Trait, y = Module, fill = Correlation)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#0D63A5", high = "#cb4936", mid = "white",
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = sprintf("%.2f\n(%.1e)", Correlation, Pvalue)), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Module-Trait Relationships")

ggsave("Module_Trait_Relationships_ggplot.pdf", p_trait, width = 8, height = 6)

###############################################################################
## 8. Gene-level statistics within modules
###############################################################################

geneModuleMembership <- as.data.frame(cor(exprMatrix, moduleEigengenes, use = "p"))
geneMMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(exprMatrix)))

geneTraitSignificance <- as.data.frame(cor(exprMatrix, design, use = "p"))
geneGSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(exprMatrix)))

geneNames <- colnames(exprMatrix)
geneInfo <- data.frame(Gene = geneNames, Module = moduleColors)

for (mod in colnames(moduleEigengenes)) {
  geneInfo[, paste0("MM_", mod)] <- geneModuleMembership[, mod]
  geneInfo[, paste0("p.MM_", mod)] <- geneMMPvalue[, mod]
}

for (trait in colnames(design)) {
  geneInfo[, paste0("GS_", trait)] <- geneTraitSignificance[, trait]
  geneInfo[, paste0("p.GS_", trait)] <- geneGSPvalue[, trait]
}

geneInfo <- geneInfo[order(geneInfo$Module), ]
write.table(geneInfo, file = "GeneInfo_Modules.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

###############################################################################
## 9. Export genes for each module
###############################################################################

uniqueModules <- unique(moduleColors)

for (mod in uniqueModules) {
  moduleGenes <- geneNames[moduleColors == mod]
  write.table(moduleGenes,
              file = paste0("ModuleGenes_", mod, ".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)
}

moduleSizes <- table(moduleColors)
moduleSizeDF <- as.data.frame(moduleSizes)
colnames(moduleSizeDF) <- c("Module", "GeneCount")

p_size <- ggplot(moduleSizeDF, aes(x = Module, y = GeneCount, fill = Module)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Gene Counts per Module", x = "Module", y = "Gene Count")

ggsave("Module_Gene_Counts.pdf", p_size, width = 8, height = 5)

###############################################################################
## 10. MM vs GS plots for all modules
###############################################################################

targetTrait <- "Tumor"

for (mod in uniqueModules) {
  pdf(file = paste0("MM_vs_GS_", mod, ".pdf"), width = 6, height = 6)

  inModule <- (moduleColors == mod)
  mmColName <- paste0("ME", mod)
  gsColName <- targetTrait

  if (!(mmColName %in% colnames(geneModuleMembership))) {
    dev.off()
    next
  }

  MM <- as.numeric(geneModuleMembership[inModule, mmColName])
  GS <- as.numeric(geneTraitSignificance[inModule, gsColName])

  if (sum(!is.na(MM) & !is.na(GS)) >= 2) {
    corTest <- cor.test(MM, GS)
    corVal <- corTest$estimate
    pVal <- corTest$p.value

    plot(MM, GS,
         xlab = paste("Module Membership in", mod, "module"),
         ylab = paste("Gene Significance for", targetTrait),
         main = paste0("Module: ", mod,
                       "\ncor = ", signif(corVal, 3),
                       ", p = ", format(pVal, scientific = TRUE, digits = 2)),
         pch = 21,
         bg = adjustcolor(mod, alpha.f = 0.6),
         col = "black",
         cex = 1.5)

    abline(lm(GS ~ MM), col = "blue", lwd = 2, lty = 2)
  }

  dev.off()
}

###############################################################################
## 11. Module-level summary plots
###############################################################################

df <- data.frame(
  Module = moduleColors,
  MM = sapply(seq_along(moduleColors), function(i) {
    geneModuleMembership[i, paste0("ME", moduleColors[i])]
  }),
  GS = geneTraitSignificance[, "Tumor"]
)

p_box <- ggplot(df, aes(x = Module, y = MM, fill = Module)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Boxplot of Module Membership across Modules",
       y = "Module Membership") +
  theme(legend.position = "none")

ggsave("Module_Membership_Boxplot.pdf", p_box, width = 8, height = 5)

p_violin <- ggplot(df, aes(x = Module, y = GS, fill = Module)) +
  geom_violin(alpha = 0.7) +
  theme_minimal() +
  labs(title = "Violin Plot of Gene Significance across Modules",
       y = "Gene Significance") +
  theme(legend.position = "none")

ggsave("Gene_Significance_Violin.pdf", p_violin, width = 8, height = 5)

###############################################################################
## 12. Correlation among modules
###############################################################################

moduleCor <- cor(moduleEigengenes, use = "p")
moduleCor_melt <- melt(moduleCor)

p_modcor <- ggplot(moduleCor_melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), name = "Correlation") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  theme_minimal() +
  ggtitle("Module-to-Module Correlation Heatmap") +
  xlab("Module Eigengenes") +
  ylab("Module Eigengenes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Module_to_Module_Correlation.pdf", p_modcor, width = 8, height = 6)

###############################################################################
## 13. TOM plot for the selected target module
###############################################################################

## The brown module is manually selected here based on its relevance
## to the tumor trait and/or the published study, rather than being
## automatically chosen by the script.
targetModule <- "brown"

inModule <- (moduleColors == targetModule)
modGenes <- geneNames[inModule]
modTOM <- tom_sim[inModule, inModule]

geneDendroModule <- hclust(as.dist(1 - modTOM), method = "average")

pdf(paste0("TOMplot_", targetModule, "_module.pdf"), width = 7, height = 7)
TOMplot(modTOM, geneDendroModule,
        main = paste("TOM Plot for", targetModule, "Module"))
dev.off()

###############################################################################
## 14. Heatmap + eigengene barplot for each module
###############################################################################

## Use the module name itself as the plotting color whenever possible
for (targetModule in uniqueModules) {
  mod_genes <- colnames(exprMatrix)[moduleColors == targetModule]

  if (length(mod_genes) < 2) next

  sample_list <- rownames(exprMatrix)
  me_col <- paste0("ME", targetModule)

  if (!(me_col %in% colnames(moduleEigengenes))) next

  mod_expr <- exprMatrix[sample_list, mod_genes, drop = FALSE]
  eigengene <- moduleEigengenes[sample_list, me_col, drop = TRUE]

  bar_df <- data.frame(
    Sample = factor(sample_list, levels = sample_list),
    Eigengene = eigengene
  )

  heatmap_grob <- grid.grabExpr({
    pheatmap(
      t(mod_expr),
      cluster_cols = FALSE,
      cluster_rows = TRUE,
      show_rownames = FALSE,
      show_colnames = FALSE,
      scale = "row",
      color = colorRampPalette(c("skyblue", "white", "orange"))(100),
      border_color = NA,
      legend = FALSE,
      fontsize = 10
    )
  })

  mod_color <- tryCatch(targetModule, error = function(e) "#888888")

  bar_p <- ggplot(bar_df, aes(x = Sample, y = Eigengene)) +
    geom_bar(stat = "identity", width = 1, fill = mod_color) +
    labs(y = "Eigengene", x = "Sample") +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      axis.title.y = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      plot.margin = margin(0, 5, 0, 5),
      axis.line = element_line(linewidth = 0.6)
    )

  pdf_file <- paste0("Module_", targetModule, "_heatmap_eigengene.pdf")
  pdf(pdf_file, width = 9, height = 7)

  grid.newpage()
  pushViewport(viewport(layout = grid.layout(
    2, 1,
    heights = unit(c(0.58, 0.38), "npc"),
    widths = unit(1, "npc")
  )))

  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.draw(heatmap_grob)
  grid.text(targetModule, x = 0.5, y = unit(0.98, "npc"),
            gp = gpar(fontface = "bold", fontsize = 20, col = mod_color))
  popViewport()

  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(bar_p, newpage = FALSE)
  popViewport()

  dev.off()
}

###############################################################################
## 15. Integrate DEG results with the selected WGCNA module
###############################################################################

setwd("~/path/to/your/directory/data/")

deg_gene <- read.csv("/path/to/your/directory/data/deg.csv")
deg_gene <- unique(deg_gene[, 1])

wgcna_gene <- fread("/path/to/your/directory/data/ModuleGenes_brown.txt", sep = "\t", header = FALSE)$V1
wgcna_gene <- unique(wgcna_gene)

geneList <- list(DEG = deg_gene, BrownModule = wgcna_gene)

library(ggvenn)
ggvenn(geneList,
       show_percentage = TRUE,
       stroke_color = "white",
       stroke_size = 0.5,
       fill_color = c("#e28764", "#728fc6"),
       set_name_color = c("#e28764", "#728fc6"),
       set_name_size = 4,
       text_size = 3.5)

# Calculate the union of gene names across all input gene sets
unionGenes <- Reduce(union, geneList)
unionCount <- length(unionGenes)
unionCount

# Save the union gene list
write.csv(unionGenes, file = "all.unionGenes.csv")

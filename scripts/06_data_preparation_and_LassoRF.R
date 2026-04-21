## This script prepares multiple HCC transcriptomic datasets, 
## converts probe-level data to gene-level expression matrices, merges and batch-corrects them across cohorts, 
## selects previously prioritized candidate genes, and constructs training and external test datasets for machine learning. 
## It then applies Lasso feature selection followed by random forest modeling to predict tumor status, 
## outputs sample-level risk scores and class labels, evaluates model performance by cohort using AUC, 
## and identifies the top genes contributing to the final model.

###############################################
# 06_multicohort_preparation_and_LassoRF_modeling.R
###############################################

###############################################
# Code block 1: Download, preprocess, and organize multiple datasets
###############################################

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required packages
library(GEOquery)
library(dplyr)
library(tidyverse)
library(data.table)
library(sva)

## Set working directory
getwd()

###############################################
# Part 1. Download and preprocess GSE14811
###############################################

### Be sure to set your own working directory
setwd("~/path/to/your/directory/GSE14811/")

### Read GEO series matrix
## If direct GEO download fails, the local file can be read directly by using filename=
gset = getGEO(filename = "GSE14811_series_matrix.txt.gz", destdir = ".", getGPL = FALSE)

#### Extract the expression matrix
exp = gset@assayData[["exprs"]]
df.exp = as.data.frame(exp)

## Extract phenotype / group information
pdata = gset@phenoData@data
table(pdata$source_name_ch2)
table(pdata$`tissue:ch2`)
## "liver" = adjacent non-tumor
## "hepatocellular carcinoma" = tumor

write.csv(pdata, file = "GSE14811.csv")

### Prepare platform annotation for probe-to-gene conversion
anno = fread("GPL8177.txt", sep = "\t", header = TRUE, data.table = FALSE)

colnames(anno)
gpl_gene = anno[, c(1, 5)]
colnames(gpl_gene)[2] = "GeneSymbol"

## Merge platform annotation with the expression matrix
exp.anno = merge(x = gpl_gene, y = df.exp, by.x = 1, by.y = 0)

### Convert probes to gene symbols
exp1 = exp.anno
rownames(exp1) = exp1$GeneSymbol

## Duplicate gene names occur because multiple probes can map to the same gene
exp2 = distinct(exp1, GeneSymbol, .keep_all = TRUE)

rownames(exp2) = exp2$GeneSymbol

## Remove the first two annotation columns
exp3 = exp2[, -c(1, 2)]

write.csv(exp3, file = "GSE14811_exp.csv")
GSE14811_exp = exp3


###############################################
# Part 2. Download and preprocess GSE54236
###############################################

setwd("~/path/to/your/directory/GSE54236/")

### Read GEO series matrix
gset = getGEO(filename = "GSE54236_series_matrix.txt.gz", destdir = ".", getGPL = FALSE)

#### Extract expression matrix
exp = gset@assayData[["exprs"]]
df.exp = as.data.frame(exp)

## Extract phenotype / group information
pdata = gset@phenoData@data
table(pdata$source_name_ch2)
table(pdata$`tissue type:ch1`)
## "Biopsy of cirrhotic non-malignant tissue" = adjacent non-tumor
## "Biopsy of tumor tissue" = tumor

write.csv(pdata, file = "GSE54236.csv")
GSE54236_pdata = pdata

### Prepare platform annotation
anno = fread("GPL6480-9577.txt", sep = "\t", header = TRUE, data.table = FALSE)

colnames(anno)
gpl_gene = anno[, c(1, 7)]
colnames(gpl_gene)[2] = "GeneSymbol"

## Merge platform annotation with the expression matrix
exp.anno = merge(x = gpl_gene, y = df.exp, by.x = 1, by.y = 0)

### Convert probes to gene symbols
exp1 = exp.anno
rownames(exp1) = exp1$GeneSymbol

## Multiple probes may correspond to the same gene
exp2 = distinct(exp1, GeneSymbol, .keep_all = TRUE)

rownames(exp2) = exp2$GeneSymbol

## Remove the first two annotation columns
exp3 = exp2[, -c(1, 2)]

write.csv(exp3, file = "GSE54236_exp.csv")
GSE54236_exp = exp3


###############################################
# Part 3. Download and preprocess GSE76427
###############################################

setwd("~/path/to/your/directory/GSE76427/")

### Read GEO series matrix
gset = getGEO(filename = "GSE76427_series_matrix.txt.gz", destdir = ".", getGPL = FALSE)

#### Extract expression matrix
exp = gset@assayData[["exprs"]]
df.exp = as.data.frame(exp)

## Extract phenotype / group information
pdata = gset@phenoData@data
table(pdata$characteristics_ch1.2)
## "adjacent non-tumor liver tissue" = adjacent non-tumor
## "primary hepatocellular carcinoma tumor" = tumor

write.csv(pdata, file = "GSE76427.csv")
GSE76427_pdata = pdata

### Prepare platform annotation
anno = fread("GPL10558-50081.txt", sep = "\t", header = TRUE, data.table = FALSE)

colnames(anno)
gpl_gene = anno[, c(1, 13)]
colnames(gpl_gene)[2] = "GeneSymbol"

## Merge platform annotation with the expression matrix
exp.anno = merge(x = gpl_gene, y = df.exp, by.x = 1, by.y = 0)

### Convert probes to gene symbols
exp1 = exp.anno
rownames(exp1) = exp1$GeneSymbol

## Multiple probes may correspond to the same gene
exp2 = distinct(exp1, GeneSymbol, .keep_all = TRUE)

rownames(exp2) = exp2$GeneSymbol

## Remove the first two annotation columns
exp3 = exp2[, -c(1, 2)]

write.csv(exp3, file = "GSE76427_exp.csv")
GSE76427_exp = exp3


###############################################
# Part 4. Merge datasets and perform batch correction
###############################################

setwd("~/path/to/your/directory/")

### Read the previously merged and batch-corrected dataset
## This dataset corresponds to the earlier combined training cohorts:
## GSE25097 + GSE36376
exp_1 = read.csv(
  "~/path/to/your/directory/correctedData.csv",
  row.names = 1
)

## Take the intersection of genes shared by all datasets
gene_collections = list(
  GSE14811_name = rownames(GSE14811_exp),
  GSE54236_name = rownames(GSE54236_exp),
  GSE76427_name = rownames(GSE76427_exp),
  exp_1 = rownames(exp_1)
)

common_genes <- Reduce(intersect, gene_collections)

## Combine all datasets using only shared genes
exp_all = cbind(
  GSE14811_exp[common_genes, ],
  GSE54236_exp[common_genes, ],
  GSE76427_exp[common_genes, ],
  exp_1[common_genes, ]
)

## Define batch labels for ComBat
batch_vector = c(
  rep("data1", length(colnames(GSE14811_exp))),
  rep("data2", length(colnames(GSE54236_exp))),
  rep("data3", length(colnames(GSE76427_exp))),
  rep("data4", length(colnames(exp_1)))
)

## Perform batch correction
corrected_mat <- ComBat(dat = exp_all, batch = batch_vector, par.prior = TRUE)


###############################################
# Part 5. Restrict to preselected candidate genes
###############################################

## Read the final candidate gene set selected from previous analyses
## Instead of using all genes directly for machine learning,
## only previously prioritized candidate genes are retained here
selected_gene = read.csv(
  "~/path/to/your/directory/all.genes.csv"
)

## all.genes.csv comes from the previous integration step
selected_gene = selected_gene$x

## Keep only the selected genes that are present in the corrected matrix
selected_gene_s = intersect(selected_gene, rownames(corrected_mat))

selected_data <- corrected_mat[selected_gene_s, ]
selected_data <- t(selected_data)   # rows = samples, columns = genes


###############################################
# Part 6. Split into training and test sets
###############################################

### Use the earlier combined dataset (exp_1) as the training set
### Use the remaining three GEO cohorts as external test sets
train_exp = as.data.frame(selected_data[colnames(exp_1), ])
test_exp  = as.data.frame(selected_data[setdiff(rownames(selected_data), colnames(exp_1)), ])


###############################################
# Part 7. Add class labels to the training set
###############################################

## Training cohorts: GSE25097 + GSE36376
## Assign labels: tumor = 1, adjacent normal = 0

phe_GSE25097 = read.csv(
  "~/path/to/your/directory/GSE25097_pdata.csv",
  row.names = 1
)
table(phe_GSE25097$source_name_ch1)
s1_n = filter(phe_GSE25097, source_name_ch1 == "non_tumor")
s1_t = filter(phe_GSE25097, source_name_ch1 == "tumor")

phe_GSE36376 = read.csv(
  "~/path/to/your/directory/GSE36376/GSE36376_pdata.csv",
  row.names = 1
)
table(phe_GSE36376$tissue.ch1)
s2_n = filter(phe_GSE36376, tissue.ch1 == "adjacent non-tumor liver")
s2_t = filter(phe_GSE36376, tissue.ch1 == "liver tumor")

exp1_col_index = c(
  rownames(s1_n), rownames(s2_n),
  rownames(s1_t), rownames(s2_t)
)

train_exp = train_exp[exp1_col_index, ]
train_exp$Type = c(
  rep(0, length(c(rownames(s1_n), rownames(s2_n)))),
  rep(1, length(c(rownames(s1_t), rownames(s2_t))))
)


###############################################
# Part 8. Add class labels to the test set
###############################################

## External test cohorts: GSE14811, GSE54236, GSE76427

# GSE14811
GSE14811_pdata = read.csv(
  "~/path/to/your/directory/GSE14811/GSE14811.csv",
  row.names = 1
)
table(GSE14811_pdata$tissue.ch2)
GSE14811_n = filter(GSE14811_pdata, tissue.ch2 == "liver")
GSE14811_t = filter(GSE14811_pdata, tissue.ch2 == "hepatocellular carcinoma")

# GSE54236
GSE54236_pdata = read.csv(
  "~/path/to/your/directory/GSE54236/GSE54236.csv",
  row.names = 1
)
table(GSE54236_pdata$tissue.type.ch1)
GSE54236_n = filter(GSE54236_pdata, tissue.type.ch1 == "Biopsy of cirrhotic non-malignant  tissue")
GSE54236_t = filter(GSE54236_pdata, tissue.type.ch1 == "Biopsy of tumor tissue")

# GSE76427
GSE76427_pdata = read.csv(
  "~/path/to/your/directory/GSE76427/GSE76427.csv",
  row.names = 1
)
table(GSE76427_pdata$characteristics_ch1.2)
GSE76427_n = filter(GSE76427_pdata, characteristics_ch1.2 == "tissue: adjacent non-tumor liver tissue")
GSE76427_t = filter(GSE76427_pdata, characteristics_ch1.2 == "tissue: primary hepatocellular carcinoma tumor")

test_col_index = c(
  rownames(GSE14811_n), rownames(GSE54236_n), rownames(GSE76427_n),
  rownames(GSE14811_t), rownames(GSE54236_t), rownames(GSE76427_t)
)

test_exp = test_exp[test_col_index, ]
test_exp$Type = c(
  rep(0, length(c(rownames(GSE14811_n), rownames(GSE54236_n), rownames(GSE76427_n)))),
  rep(1, length(c(rownames(GSE14811_t), rownames(GSE54236_t), rownames(GSE76427_t))))
)

## Save training and test sets
setwd("~/path/to/your/directory/")
write.csv(test_exp,  file = "test.csv")
write.csv(train_exp, file = "train.csv")


###############################################
# Code block 2: Full machine-learning workflow
###############################################

## Reset environment
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

## Load required packages
library(dplyr)
library(tidyverse)
library(data.table)
library(glmnet)
library(randomForestSRC)
library(pROC)

## Set working directory
getwd()
setwd("~/path/to/your/directory/")


###############################################
# Part 9. Read training and test files
###############################################

Train_data <- read.table("train.csv", header = TRUE, sep = ",",
                         check.names = FALSE, row.names = 1,
                         stringsAsFactors = FALSE)

Train_expr = Train_data[, 1:(ncol(Train_data) - 1), drop = FALSE]
Train_class = Train_data[, ncol(Train_data), drop = FALSE]

Test_data <- read.table("test.csv", header = TRUE, sep = ",",
                        check.names = FALSE, row.names = 1,
                        stringsAsFactors = FALSE)

Test_expr = Test_data[, 1:(ncol(Test_data) - 1), drop = FALSE]
Test_class = Test_data[, ncol(Test_data), drop = FALSE]


###############################################
# Part 10. Add cohort information to the test set
###############################################

## Because the test set contains three cohorts,
## cohort information is added for later cohort-wise evaluation

# GSE14811
GSE14811_pdata = read.csv(
  "~/path/to/your/directory/GSE14811/GSE14811.csv"
)
GSE14811_pdata$Cohort = "GSE14811"
GSE14811_pdata = GSE14811_pdata[, c(1, 41)]

# GSE54236
GSE54236_pdata = read.csv(
  "~/path/to/your/directory/GSE54236/GSE54236.csv"
)
GSE54236_pdata$Cohort = "GSE54236"
GSE54236_pdata = GSE54236_pdata[, c(1, 47)]

# GSE76427
GSE76427_pdata = read.csv(
  "~/path/to/your/directory/GSE76427/GSE76427.csv"
)
GSE76427_pdata$Cohort = "GSE76427"
GSE76427_pdata = GSE76427_pdata[, c(1, 54)]

## Combine phenotype tables from all external cohorts
pdata_all = rbind(GSE14811_pdata, GSE54236_pdata, GSE76427_pdata)
table(pdata_all$Cohort)

## Merge cohort information into Test_class
Test_class = merge(Test_class, pdata_all, by.x = 0, by.y = 1)
rownames(Test_class) = Test_class$Row.names
Test_class = Test_class[, -1]
Test_class = Test_class[, c("Cohort", "Type")]

## Each test sample now has:
## Type   = tumor or adjacent normal
## Cohort = source dataset


###############################################
# Part 11. Define scaling functions
###############################################

## Standardization helper
standarize.fun <- function(indata, centerFlag, scaleFlag) {
  scale(indata, center = centerFlag, scale = scaleFlag)
}

## Standardize either:
## - the whole dataset together, or
## - each cohort separately
scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL) {
  samplename = rownames(data)

  if (is.null(cohort)) {
    data <- list(data)
    names(data) = "training"
  } else {
    data <- split(as.data.frame(data), cohort)
  }

  if (is.null(centerFlags)) {
    centerFlags = FALSE
    message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags) == 1) {
    centerFlags = rep(centerFlags, length(data))
    message("Set centerFlags for all cohorts as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))) {
    names(centerFlags) <- names(data)
    message("Match centerFlags with cohorts by order")
  }

  if (is.null(scaleFlags)) {
    scaleFlags = FALSE
    message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags) == 1) {
    scaleFlags = rep(scaleFlags, length(data))
    message("Set scaleFlags for all cohorts as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))) {
    names(scaleFlags) <- names(data)
    message("Match scaleFlags with cohorts by order")
  }

  centerFlags <- centerFlags[names(data)]
  scaleFlags  <- scaleFlags[names(data)]

  outdata <- mapply(
    standarize.fun,
    indata = data,
    centerFlag = centerFlags,
    scaleFlag = scaleFlags,
    SIMPLIFY = FALSE
  )

  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}


###############################################
# Part 12. Harmonize features and standardize data
###############################################

## Keep only genes shared by training and test sets
comgene <- intersect(colnames(Train_expr), colnames(Test_expr))

## The model must use exactly the same features in training and test data
Train_expr <- as.matrix(Train_expr[, comgene])
Test_expr  <- as.matrix(Test_expr[, comgene])

## Standardize the training set using z-score normalization
Train_set = scaleData(data = Train_expr, centerFlags = TRUE, scaleFlags = TRUE)

## Add a small amount of random noise to the training set
## This is mainly for tutorial demonstration and reproducibility
set.seed(123)
noise <- matrix(
  rnorm(n = nrow(Train_set) * ncol(Train_set), mean = 0, sd = 0.01),
  nrow = nrow(Train_set),
  ncol = ncol(Train_set)
)
Train_set <- Train_set + noise

## Standardize the test set by cohort
Test_set = scaleData(data = Test_expr,
                     cohort = Test_class$Cohort,
                     centerFlags = TRUE,
                     scaleFlags = TRUE)

## Add a small amount of random noise to the test set
## Strictly speaking, adding noise to the test set is not standard practice
noise_test <- matrix(
  rnorm(n = nrow(Test_set) * ncol(Test_set), mean = 0, sd = 0.01),
  nrow = nrow(Test_set),
  ncol = ncol(Test_set)
)
Test_set <- Test_set + noise_test


###############################################
# Part 13. Lasso feature selection
###############################################

## y must be numeric 0/1
y_train <- Train_class$Type

## Lasso (alpha = 1) with 10-fold cross-validation
cvfit <- cv.glmnet(
  x = as.matrix(Train_set),
  y = y_train,
  family = "binomial",
  alpha = 1,
  nfolds = 10
)

## Extract coefficients at lambda.min
coef_mat <- coef(cvfit, s = "lambda.min")
coef_mat <- as.matrix(coef_mat)

## Keep non-zero coefficients
lasso_features <- rownames(coef_mat)[which(as.numeric(coef_mat) != 0)]

## Remove intercept
lasso_features <- setdiff(lasso_features, c("(Intercept)", "Intercept"))

cat("Number of features selected by Lasso:", length(lasso_features), "\n")

## Subset the data using selected features
Train_sel <- as.data.frame(Train_set[, lasso_features, drop = FALSE])
Test_sel  <- as.data.frame(Test_set[,  lasso_features, drop = FALSE])


###############################################
# Part 14. Random Forest model training
###############################################

rf_train_df <- cbind(Train_sel, Type = factor(Train_class$Type))

set.seed(123)
rf_fit <- rfsrc(
  Type ~ .,
  data       = rf_train_df,
  ntree      = 1000,
  nodesize   = 5,
  importance = TRUE,
  proximity  = FALSE,
  forest     = TRUE
)


###############################################
# Part 15. Predict risk scores (probabilities)
###############################################

## Combine Train + Test so that predictions can later be split by cohort
all_sel <- rbind(Train_sel, Test_sel)
rf_pred <- predict(rf_fit, newdata = all_sel, na.action = "na.impute")

## Robustly extract the probability of class 1 (tumor)
if (is.matrix(rf_pred$predicted)) {
  rf_prob <- if ("1" %in% colnames(rf_pred$predicted)) {
    rf_pred$predicted[, "1"]
  } else {
    rf_pred$predicted[, 2]
  }
} else {
  rf_prob <- as.numeric(rf_pred$predicted)
}

## Split probabilities back into train and test
n_tr <- nrow(Train_sel)
rf_prob_train <- rf_prob[seq_len(n_tr)]
rf_prob_test  <- rf_prob[-seq_len(n_tr)]


###############################################
# Part 16. Predict class labels
###############################################

rf_cls_train <- ifelse(rf_prob_train > 0.5, "1", "0")
rf_cls_test  <- ifelse(rf_prob_test  > 0.5, "1", "0")

## Save risk and class matrices
risk_mat <- data.frame(
  id = c(rownames(Train_sel), rownames(Test_sel)),
  RF_Prob = rf_prob,
  stringsAsFactors = FALSE
)

class_mat <- data.frame(
  id = c(rownames(Train_sel), rownames(Test_sel)),
  RF_Class = c(rf_cls_train, rf_cls_test),
  stringsAsFactors = FALSE
)

write.table(risk_mat,  "LassoRF.riskMatrix.txt",  sep = "\t", row.names = FALSE, quote = FALSE)
write.table(class_mat, "LassoRF.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)


###############################################
# Part 17. Calculate cohort-wise AUC
###############################################

## Mark the training set as cohort "Train"
lab_df <- rbind(
  data.frame(Cohort = "Train", Type = Train_class$Type, stringsAsFactors = FALSE),
  data.frame(Cohort = Test_class$Cohort, Type = Test_class$Type, stringsAsFactors = FALSE)
)

table(lab_df$Cohort)

## Convert labels to 0/1 if needed
y_all <- lab_df$Type
if (!all(y_all %in% c(0, 1))) {
  y_all <- as.numeric(factor(y_all)) - 1
}

## Calculate ROC AUC for each cohort separately
auc_by_cohort <- tapply(seq_along(y_all), lab_df$Cohort, function(idx) {
  yy <- y_all[idx]
  pp <- rf_prob[idx]

  ## AUC requires both classes to be present
  if (length(unique(yy[!is.na(yy)])) < 2) return(NA_real_)

  as.numeric(auc(roc(response = yy, predictor = pp, quiet = TRUE)))
})

auc_df <- data.frame(
  Cohort = names(auc_by_cohort),
  AUC = as.numeric(auc_by_cohort),
  row.names = NULL
)

print(auc_df)
write.table(auc_df, "LassoRF.AUC.byCohort.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

## Optional overall AUC across all samples
overall_auc <- as.numeric(auc(roc(response = y_all, predictor = rf_prob, quiet = TRUE)))
cat("Overall AUC (Lasso + RF):", sprintf("%.4f", overall_auc), "\n")


###############################################
# Part 18. Extract RF variable importance and top 10 genes
###############################################

if (is.null(rf_fit$importance)) {
  stop("rf_fit$importance is empty. Please ensure importance = TRUE was set in rfsrc().")
}

vimp_vec <- rf_fit$importance

if (is.matrix(vimp_vec) || is.data.frame(vimp_vec)) {
  if ("all" %in% colnames(vimp_vec)) {
    vimp_vec <- vimp_vec[, "all"]
  } else {
    vimp_vec <- rowMeans(as.matrix(vimp_vec))
  }
}

## Keep only features actually used in the model
vimp_vec <- vimp_vec[names(vimp_vec) %in% colnames(Train_sel)]

## Rank features by importance
ord <- order(vimp_vec, decreasing = TRUE, na.last = NA)
topN <- min(10, length(ord))
top10_genes <- names(vimp_vec)[ord][seq_len(topN)]
top10_imp   <- vimp_vec[top10_genes]

cat("Top 10 important genes:\n")
print(top10_genes)
print(top10_imp)


###############################################
# Notes
###############################################

## Adding noise to the test set is not considered standard rigorous practice.
## This is mainly included for tutorial demonstration.
## In formal analysis, it is generally better to:
## - consider augmentation or resampling in the training set if needed
## - keep the test set as unchanged as possible

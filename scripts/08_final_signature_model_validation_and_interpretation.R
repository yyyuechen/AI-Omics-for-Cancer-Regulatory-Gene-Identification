# This script performs final validation and interpretation of a 6-gene HCC signature using logistic regression, 
# ROC/DCA/calibration analyses, multi-model comparison, and SHAP-based explainability.

# In 07_machine_learning_model_comparison.R, the workflow focused on:
#
# - variable-selection methods such as Lasso / Enet / Stepglm
# - model methods such as RF / SVM / GBM
# - comparing many model combinations by AUC
# - identifying the best-performing strategy
# - retaining the final core genes
#
# So that step answered:
#
# Which genes should be kept?
# Which modeling strategy performs best?
#
# That was the feature-selection stage.
#
# What does 08_final_signature_model_validation_and_interpretation.R do?
#
# This step is more like:
#
# Now that the final 6 genes have been determined,
# treat them as the final signature and build a final diagnostic model
# for validation and interpretation.
#
# In other words, this step is no longer about selecting genes.
# Instead, it asks:
#
# - Do the 6 genes individually have diagnostic power?
# - Do they perform better when combined?
# - How can the final 6-gene model be interpreted?
# - Which genes contribute most to the final model?
#
# So this step answers:
#
# How well does the final 6-gene diagnostic model perform,
# and how should its performance be interpreted?
#
# This is the finalization + validation + interpretation stage.

## Set system messages to English
Sys.setenv(LANGUAGE = "en")

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required packages
library(GEOquery)
library(dplyr)
library(tidyverse)
library(rms)
library(regplot)
library(rmda)
library(pROC)
library(tools)
library(ggplot2)

## Set working directory
getwd()

# Because the preprocessing was not exactly the same as in the paper,
# the key genes selected by our own model are not identical to those reported in the paper:
# ALDH1B1  glmBoost+Stepglm[both]
# AURKA    glmBoost+Stepglm[both]
# BCAT2    glmBoost+Stepglm[both]
# GSTP1    glmBoost+Stepglm[both]
# PCK1     glmBoost+Stepglm[both]
# RND3     glmBoost+Stepglm[both]
# CCNB2    glmBoost+Stepglm[both]
# CDA      glmBoost+Stepglm[both]
# CTSB     glmBoost+Stepglm[both]

############################################################################
############################################################################
############################################################################
### 2. Read data and preprocess
setwd("~/path/to/your/directory/")

test_exp = read.csv("test.csv", row.names = 1)
colnames(test_exp)

## Read the final selected 6 genes
selected_gene_s = c("RND3", "PCK1", "AURKA", "CCNB1", "BCAT2", "UCK2", "Type")
test_selected_data <- test_exp[, selected_gene_s]

train_exp = read.csv("train.csv", row.names = 1)
colnames(train_exp)

## Read the final selected 6 genes
selected_gene_s = c("RND3", "PCK1", "AURKA", "CCNB1", "BCAT2", "UCK2", "Type")
train_selected_data <- train_exp[, selected_gene_s]

## Combine training and test data
df_expr <- rbind(train_selected_data, test_selected_data)
colnames(df_expr)[7] = "GroupType"

setwd("~/path/to/your/directory/")
write.csv(df_expr, file = "df_expr.csv")

############################################################################
############################################################################
############################################################################
## 3. Prepare model input
ddinfo <- datadist(df_expr)
options(datadist = "ddinfo")

############################################################################
############################################################################
############################################################################
## 4. Logistic regression modeling and nomogram
model_vars <- setdiff(colnames(df_expr), "GroupType")
reg_formula <- as.formula(paste("GroupType ~", paste(model_vars, collapse = " + ")))

# Use the 6 genes RND3, PCK1, AURKA, CCNB1, BCAT2, and UCK2
# to predict GroupType

# Fit the model and draw a nomogram
# rms::lrm() is used because the outcome is binary
lrm_fit <- rms::lrm(reg_formula, data = df_expr, x = TRUE, y = TRUE)

nomo_obj <- rms::nomogram(
  lrm_fit,
  fun = plogis,
  fun.at = c(0.05, 0.2, 0.5, 0.8, 0.95),
  lp = FALSE,
  funlabel = "Disease Risk"
)

plot(nomo_obj)

############################################################################
############################################################################
############################################################################
## 5. Calibration curve
# Use bootstrap resampling (1000 times) for internal validation
# to reduce overfitting from evaluating only on the training data
calibrate_obj <- calibrate(lrm_fit, method = "boot", B = 1000)

plot(calibrate_obj, xlab = "Predicted", ylab = "Observed", sub = FALSE)

# Ideal: ideal calibration line
# Apparent: apparent performance
# Bias-corrected: bootstrap-corrected performance
#
# If the bias-corrected curve is close to the ideal line,
# the predicted probabilities are relatively reliable.

############################################################################
############################################################################
############################################################################
## 6. DCA (Decision Curve Analysis)
# Explicitly convert the label to 0/1,
# which is the most stable input format for rmda::decision_curve()
# ROC evaluates statistical discrimination,
# whereas DCA evaluates net clinical benefit
df_expr$GroupType <- as.numeric(factor(df_expr$GroupType)) - 1
set.seed(123456)

dca_obj <- decision_curve(
  formula = reg_formula,
  data = df_expr,
  thresholds = seq(0, 1, by = 0.01),
  family = binomial(link = "logit"),
  bootstraps = 100
)

plot_decision_curve(
  dca_obj,
  xlab = "Threshold Probability",
  col = "orange",
  confidence.intervals = TRUE,
  standardize = TRUE,
  cost.benefit.axis = TRUE,
  legend.position = "bottomleft"
)

############################################################################
############################################################################
############################################################################
## 7. Export tables

## Output the cleaned expression matrix
setwd("~/path/to/your/directory/")
write.csv(df_expr, file = "Filtered_Expression_Matrix.csv", quote = FALSE)

## DCA-derived data
Decision_Curve_Data_df = as.data.frame(dca_obj$derived.data)
write.csv(Decision_Curve_Data_df, file = "Decision_Curve_Data.csv")

## Predicted probability for each sample
pred_probs <- predict(lrm_fit, newdata = df_expr, type = "fitted")
df_out <- data.frame(
  Sample = rownames(df_expr),
  GroupType = df_expr$GroupType,
  Predicted_Prob = pred_probs
)
write.csv(df_out, file = "Sample_Predicted_Probabilities.csv")

############################################################################
############################################################################
############################################################################
## 8. Export coefficient table based on the logistic regression model
# Generate coefficients (beta), OR, 95% CI, and significance table

coef_df <- as.data.frame(summary(lrm_fit))
if ("Effect" %in% names(coef_df) & "S.E." %in% names(coef_df)) {
  coef_df$Effect <- as.numeric(as.character(coef_df$Effect))
  coef_df$`S.E.` <- as.numeric(as.character(coef_df$`S.E.`))
  coef_df$OR <- exp(coef_df$Effect)
  coef_df$OR_low <- exp(coef_df$Effect - 1.96 * coef_df$`S.E.`)
  coef_df$OR_high <- exp(coef_df$Effect + 1.96 * coef_df$`S.E.`)

  write.csv(coef_df, file = "Model_Coefficients.csv", row.names = FALSE)

  if ("Lower 0.95" %in% names(coef_df) & "Upper 0.95" %in% names(coef_df)) {
    write.csv(coef_df, file = "Model_Coefficients_With_CI.csv", row.names = FALSE)
  }

  or_selected <- coef_df[!is.na(coef_df$OR) & coef_df$OR > 1, ]
  write.csv(or_selected, file = "Model_OR_GT1.csv", row.names = FALSE)

  if ("P" %in% names(coef_df)) {
    sig_coef <- coef_df[!is.na(coef_df$P) & coef_df$P < 0.05, ]
    write.csv(sig_coef, file = "Model_Significant_Coefficients.csv", row.names = FALSE)
  }
}

############################################################################
############################################################################
############################################################################
## 9. Enhanced nomogram using regplot ("virtual patient")

median_point <- sapply(df_expr[, model_vars, drop = FALSE], median)
median_df <- as.data.frame(t(median_point))
median_df$GroupType <- 1
rownames(median_df) <- "AllMedian"

regplot(
  lrm_fit,
  showP = TRUE,
  rank = "sd",
  distribution = TRUE,
  observation = median_df,
  title = "Prediction Nomogram"
)

############################################################################
############################################################################
############################################################################
## 10. Combined ROC curve for the logistic regression model

roc_y <- as.numeric(df_expr$GroupType)
roc_pred <- pred_probs
roc_obj <- roc(roc_y, roc_pred, levels = c(0, 1), direction = "<")

# Define a set of fixed distinguishable colors
my_cols <- c("red", "gold", "forestgreen", "deepskyblue",
             "blue", "purple")

while (length(my_cols) < length(model_vars) + 1) {
  my_cols <- c(my_cols, rainbow(length(model_vars) + 1 - length(my_cols)))
}

# Calculate AUC
auc_val <- as.numeric(auc(roc_obj))

# Calculate 95% confidence interval of AUC
auc_ci <- ci(roc_obj)

# Visualization
par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)

plot(1 - roc_obj$specificities, roc_obj$sensitivities, type = "l",
     col = my_cols[1], lwd = 4, lty = 1,
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = "Model ROC Curve",
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45,
     xlim = c(0, 1), ylim = c(0, 1))

abline(0, 1, lty = 2, col = "gray70", lwd = 2)

legend(
  "bottomright",
  legend = sprintf("AUC = %.3f (95%% CI: %.3f, %.3f)", auc_val, auc_ci[1], auc_ci[2]),
  col = my_cols[1],
  lwd = 4,
  lty = 1,
  bty = "n",
  cex = 1.2
)

############################################################################
############################################################################
############################################################################
## 11. Multi-color ROC curves for individual genes

stopifnot(length(my_cols) >= length(model_vars) + 1)

expression_threshold <- median(
  as.matrix(df_expr[, model_vars, drop = FALSE]),
  na.rm = TRUE
)

roc_list <- list()
auc_list <- c()
leglab <- c()

par(mar = c(5, 6, 4, 2) + 0.1, cex = 1.3)
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = expression("1 - Specificity"), ylab = "Sensitivity",
     main = "Gene ROC Curves",
     cex.lab = 1.4, cex.axis = 1.15, cex.main = 1.45)
abline(0, 1, lty = 2, col = "gray70", lwd = 2)

lab <- df_expr$GroupType

for (i in seq_along(model_vars)) {
  g <- model_vars[i]
  gene_exp <- df_expr[[g]]

  keep <- !(is.na(gene_exp) | is.na(lab))
  gene_exp_i <- gene_exp[keep]
  lab_i <- lab[keep]

  if (length(gene_exp_i) < 5 || length(unique(lab_i)) < 2) {
    message(sprintf("Skip %s: not enough valid samples or label has a single class.", g))
    next
  }

  mexp <- mean(gene_exp_i, na.rm = TRUE)
  if (is.na(mexp) || is.na(expression_threshold)) {
    message(sprintf("Skip %s: mean or threshold is NA.", g))
    next
  }

  roc_y <- if (mexp > expression_threshold) 1 - as.numeric(lab_i) else as.numeric(lab_i)

  cur_roc <- pROC::roc(response = roc_y, predictor = gene_exp_i, quiet = TRUE, na.rm = TRUE)

  if (is.null(cur_roc)) {
    message(sprintf("Skip %s: ROC failed.", g))
    next
  }

  cur_auc <- as.numeric(pROC::auc(cur_roc))

  if (!is.na(cur_auc) && cur_auc < 0.5) {
    cur_roc <- pROC::roc(response = roc_y, predictor = -gene_exp_i, quiet = TRUE, na.rm = TRUE)
    cur_auc <- as.numeric(pROC::auc(cur_roc))
  }

  lines(1 - cur_roc$specificities, cur_roc$sensitivities,
        col = my_cols[i + 1], lwd = 3)

  roc_list[[g]] <- cur_roc
  auc_list <- c(auc_list, cur_auc)
  leglab <- c(leglab, sprintf("%s   AUC=%.3f", g, cur_auc))
}

legend(
  "bottomright",
  legend = leglab,
  col = my_cols[seq_along(leglab) + 1],
  lwd = 3,
  cex = 0.9,
  bty = "n"
)

gene_auc_df <- data.frame(Gene = model_vars, AUC = auc_list)
write.csv(gene_auc_df, file = "IndividualGenes_AUC.csv", row.names = FALSE)

############################################################
# Code block 2: Fig5DEFG
############################################################

## Set system messages to English
Sys.setenv(LANGUAGE = "en")

## Prevent automatic conversion of strings to factors
options(stringsAsFactors = FALSE)

## Clear the workspace
rm(list = ls())

## Load required packages
library(caret)
library(DALEX)
library(ggplot2)
library(randomForest)
library(kernlab)
library(pROC)
library(shapviz)
library(xgboost)
library(klaR)
library(RColorBrewer)
library(pheatmap)

setwd("~/path/to/your/directory/")
set.seed(12345678)

# Create output folder
if (!dir.exists("results")) dir.create("results")

# ===================== Step 3: Read data and preprocess =====================

test_exp = read.csv("~/path/to/your/directory/test.csv", row.names = 1)
colnames(test_exp)

## Read the final selected 6 genes
selected_gene_s = c("RND3", "PCK1", "AURKA", "CCNB1", "BCAT2", "UCK2", "Type")
test_selected_data <- test_exp[, selected_gene_s]

train_exp = read.csv("~/path/to/your/directory/train.csv", row.names = 1)
colnames(train_exp)

## Read the final selected 6 genes
selected_gene_s = c("RND3", "PCK1", "AURKA", "CCNB1", "BCAT2", "UCK2", "Type")
train_selected_data <- train_exp[, selected_gene_s]

## Combine training and test data
df_expr <- rbind(train_selected_data, test_selected_data)
colnames(df_expr)[7] = "Group"
table(df_expr$Group)

df_expr$Group = ifelse(df_expr$Group == 0, "Control", "Treatment")
table(df_expr$Group)

df_expr$Group <- as.factor(df_expr$Group)
expr_data = na.omit(df_expr)

# ===================== Step 4: Split into training and test sets =====================

if (length(unique(expr_data$Group)) != 2) stop("The number of groups is not 2, so binary classification cannot be performed.")

split_idx <- createDataPartition(y = expr_data$Group, p = 0.7, list = FALSE)
train_set <- expr_data[split_idx, ]
test_set <- expr_data[-split_idx, ]
test_labels <- test_set$Group
test_features <- test_set[, -ncol(test_set)]

# ===================== Step 5: Model training and evaluation =====================

ml_methods <- data.frame(
  ModelName = c("RF", "SVM", "XGB", "GBM", "KNN"),
  MethodID = c("rf", "svmRadial", "xgbTree", "gbm", "knn")
)

color_palette = c('#cb4936', '#5ac1b3', '#1616c8',
                  '#A12568', '#F499C1', '#F7C394', '#B2A157', '#ade87c',
                  '#000000', '#03C4A1', '#bb4316', '#7382BC', '#F0E442',
                  '#3B185F', '#0d6c0d', '#FEC260', '#FD7014', '#a67063',
                  '#1B9B8A', '#D0EBE7', '#713045', '#F6E0EA', '#AD6D28',
                  '#EAB67D', '#5ac1b3', '#EE4590')

auc_results <- c()
model_auc_map <- list()
roc_colors <- color_palette[1:nrow(ml_methods)]
roc_ci_list <- list()

pdf_width <- 6
pdf_height <- 6
pdf(file = "model_roc_curve.pdf", width = 6, height = 6)
par(mar = c(5, 5, 4, 2) + 0.1)

for (i in 1:nrow(ml_methods)) {
  cv_folds <- 5
  mdl_name <- ml_methods$ModelName[i]
  mdl_id <- ml_methods$MethodID[i]

  if (nrow(train_set) < 10) stop("Too few samples in the training set.")

  if (mdl_id == "svmRadial") {
    mdl <- train(
      Group ~ ., data = train_set, method = mdl_id, prob.model = TRUE,
      trControl = trainControl(method = "repeatedcv", number = cv_folds, savePredictions = TRUE)
    )
  } else {
    mdl <- train(
      Group ~ ., data = train_set, method = mdl_id,
      trControl = trainControl(method = "repeatedcv", number = cv_folds, savePredictions = TRUE)
    )
  }

  pred_prob <- predict(mdl, newdata = test_features, type = "prob")
  if (!"Treatment" %in% colnames(pred_prob)) stop("The probability output does not contain a 'Treatment' column.")

  roc_obj <- roc(as.numeric(test_labels) - 1, as.numeric(pred_prob[, "Treatment"]))
  auc_val <- as.numeric(roc_obj$auc)
  auc_ci <- ci.auc(roc_obj, conf.level = 0.95)

  roc_ci_list[[mdl_id]] <- auc_ci
  auc_results <- c(
    auc_results,
    sprintf("%s: %.03f [%.03f, %.03f]", mdl_name, auc_val, auc_ci[1], auc_ci[3])
  )
  model_auc_map[[mdl_id]] <- auc_val

  if (i == 1) {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "",
         col = roc_colors[i], lwd = 3)
  } else {
    plot(roc_obj, print.auc = FALSE, legacy.axes = TRUE, main = "",
         col = roc_colors[i], lwd = 3, add = TRUE)
  }
}

legend("bottomright", auc_results, col = roc_colors, lwd = 3, bty = "n", cex = 0.9)
dev.off()

# ===================== Step 6: Best-model selection and SHAP calculation =====================

auc_vec <- unlist(model_auc_map)
if (length(auc_vec) == 0) stop("AUC results are empty.")
best_method <- names(which.max(auc_vec))
cat("The model with the highest AUC is:", best_method, "\n")

## 0) Preliminary checks: ensure the response variable and features are valid
stopifnot(ncol(train_set) >= 2)

# Response variable: Group
# Two classes: Treatment / Control
train_set$Group <- factor(train_set$Group)
if (!all(c("Treatment", "Control") %in% levels(train_set$Group))) {
  stop("The levels of Group must include 'Treatment' and 'Control'. Please adjust according to your data.")
}
train_set$Group <- relevel(train_set$Group, ref = "Treatment")

# Use numeric features only (required by SVM-RBF), and ensure no NA/Inf
X_train <- train_set[, -ncol(train_set), drop = FALSE]
if (!all(vapply(X_train, is.numeric, logical(1)))) {
  stop("svmRadial requires all features to be numeric. Please convert character/factor variables before fitting.")
}
stopifnot(!anyNA(X_train))
stopifnot(all(is.finite(as.matrix(X_train))))

## 1) Training: enable class probabilities and let ksvm use probability modeling
library(caret)
set.seed(12345)

ctrl <- trainControl(
  method = "repeatedcv",
  number = cv_folds,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = TRUE
)

final_model <- train(
  Group ~ ., data = data.frame(X_train, Group = train_set$Group),
  method = best_method,
  preProcess = c("center", "scale"),
  tuneLength = 10,
  metric = "ROC",
  trControl = ctrl,
  prob.model = TRUE
)

## 2) Check probability prediction: NA is not allowed
probs <- predict(final_model, newdata = X_train, type = "prob")
if (!"Treatment" %in% colnames(probs)) {
  stop("The probability output does not contain a 'Treatment' column. Please check the naming of Group levels.")
}
if (anyNA(probs[, "Treatment"])) {
  stop("Probability prediction contains NA, so SHAP cannot be calculated. Please check the data or use another model such as xgbTree or rf.")
}

## 3) kernelshap: provide a small background sample bg_X for better stability and speed
library(kernelshap)
set.seed(12345)

bg_id <- sample(seq_len(nrow(X_train)), size = min(50, nrow(X_train)))
bg_X <- X_train[bg_id, , drop = FALSE]

pred_fun <- function(model, newdata) {
  p <- predict(model, newdata = newdata, type = "prob")[, "Treatment"]
  if (anyNA(p)) stop("pred_fun returned NA, so kernelshap cannot continue.")
  as.numeric(p)
}

shap_fit <- kernelshap(
  object = final_model,
  X = X_train,
  bg_X = bg_X,
  pred_fun = pred_fun
)

stopifnot(!is.null(shap_fit$S))
stopifnot(!anyNA(shap_fit$S))
stopifnot(is.finite(shap_fit$baseline))

# ===================== Step 7: SHAP visualization =====================

shap_obj <- shapviz(
  shap_fit,
  X_pred = train_set[, -ncol(train_set)],
  X = train_set[, -ncol(train_set)],
  interactions = TRUE
)

shap_importance <- sort(colMeans(abs(shap_obj$S)), decreasing = TRUE)
top_features <- names(shap_importance)

custom_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right"
  )

options(shapviz.colors = color_palette)

# SHAP importance bar plot
print(sv_importance(shap_obj, kind = "bar", show_numbers = TRUE) + custom_theme)
dev.off()

# SHAP importance beeswarm plot
print(sv_importance(shap_obj, kind = "bee", show_numbers = TRUE) + custom_theme)

# Dependence plots for all top features
print(sv_dependence(shap_obj, v = top_features) + custom_theme)

# Force plot for sample 17
print(sv_force(shap_obj, row_id = 17) + custom_theme)
rownames(train_set)[17]

## 06_data_preparation_and_LassoRF.R

## asks:

## Can Lasso + RF classify tumor vs normal?
## What is the AUC of this one pipeline?
## What are the top features in this one model?
## 07_machine_learning_model_comparison.R

## asks:

## Among many ML pipelines, which performs best?
## Which feature-selection strategy is most effective?
## Which classifier generalizes best across cohorts?
## Which genes are consistently selected across models?
## Can a simpler final logistic model be built on top of selected features?

## This script performs multi-model benchmarking and feature-stability analysis on the HCC machine-learning dataset, 
## comparing different feature-selection/classification pipelines and identifying robust predictive genes across cohorts.

###############################################
# 07_multimodel_benchmarking_and_feature_stability.R
###############################################

Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list = ls())

###############################################
# Load required packages
###############################################

library(dplyr)
library(tidyverse)
library(data.table)
library(tgp)
library(snowfall)        # parallel computing support
library(seqinr)          # sequence analysis utilities
library(plyr)            # data processing tools
library(randomForestSRC) # random forest
library(glmnet)          # generalized linear models with regularization
library(plsRglm)         # partial least squares regression for GLM
library(gbm)             # gradient boosting machine
library(mboost)          # boosting models
library(e1071)           # SVM and related algorithms
library(BART)            # Bayesian additive regression trees
library(MASS)            # statistical methods
library(xgboost)         # extreme gradient boosting
library(caret)           # model training and testing
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)
library(circlize)

###############################################
# Set working directory
###############################################

getwd()
setwd("~//path/to/your/directory/")

############################################################
## Define model-specific functions
############################################################

### Model 1: Lasso regression
RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

### Model 2: Ridge regression
RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

### Model 3: Elastic Net
RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = cv.glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    nfolds = 10
  )
  
  fit = glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    lambda = cv.fit$lambda.min
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 4: Stepwise GLM
RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- step(
    glm(
      formula = Train_label[[classVar]] ~ .,
      family = "binomial",
      data = as.data.frame(Train_set)
    ),
    direction = direction,
    trace = 0
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 5: PLSRGLM
RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  cv.plsRglm.res = cv.plsRglm(
    formula = Train_label[[classVar]] ~ .,
    data = as.data.frame(Train_set),
    nt = 10,
    verbose = FALSE
  )
  
  fit <- plsRglm(
    Train_label[[classVar]],
    as.data.frame(Train_set),
    modele = "pls-glm-logistic",
    verbose = FALSE,
    sparse = TRUE
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 6: GLMBoost
RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  fit <- glmboost(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    family = Binomial()
  )
  
  cvm <- cvrisk(
    fit,
    papply = lapply,
    folds = cv(model.weights(fit), type = "kfold")
  )
  
  fit <- glmboost(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    family = Binomial(),
    control = boost_control(mstop = max(mstop(cvm), 40))
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 7: Linear Discriminant Analysis
RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit = train(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    method = "lda",
    trControl = trainControl(method = "cv")
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 8: Naive Bayes
RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  fit <- naiveBayes(
    eval(parse(text = paste(classVar, "~."))),
    data = data
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 9: Random Forest
RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- rfsrc(
    formula = formula(paste0(classVar, "~.")),
    data = cbind(Train_set, Train_label[classVar]),
    ntree = 1000,
    nodesize = rf_nodesize,
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 10: GBM
RunGBM <- function(Train_set, Train_label, mode, classVar){
  fit <- gbm(
    formula = Train_label[[classVar]] ~ .,
    data = as.data.frame(Train_set),
    distribution = "bernoulli",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = 6
  )
  
  best <- which.min(fit$cv.error)
  
  fit <- gbm(
    formula = Train_label[[classVar]] ~ .,
    data = as.data.frame(Train_set),
    distribution = "bernoulli",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    n.cores = 8
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 11: XGBoost
RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  indexes = createFolds(Train_label[[classVar]], k = 5, list = TRUE)
  
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgb.DMatrix(data = Train_set[-pt, ], label = Train_label[-pt, ])
    dtest  = xgb.DMatrix(data = Train_set[pt, ],  label = Train_label[pt, ])
    watchlist <- list(train = dtrain, test = dtest)
    
    bst <- xgb.train(
      data = dtrain,
      max.depth = 2,
      eta = 1,
      nthread = 2,
      nrounds = 10,
      watchlist = watchlist,
      objective = "binary:logistic",
      verbose = FALSE
    )
    
    which.min(bst$evaluation_log$test_logloss)
  }))
  
  nround <- as.numeric(names(which.max(table(CV))))
  
  fit <- xgboost(
    data = Train_set,
    label = Train_label[[classVar]],
    max.depth = 2,
    eta = 1,
    nthread = 2,
    nrounds = nround,
    objective = "binary:logistic",
    verbose = FALSE
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

### Model 12: SVM
RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit = svm(
    formula = eval(parse(text = paste(classVar, "~."))),
    data = data,
    probability = TRUE
  )
  
  fit$subFeature = colnames(Train_set)
  
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

############################################################
## General helper functions
############################################################

RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  method = gsub(" ", "", method)
  method_name = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  
  method_param = switch(
    EXPR = method_name,
    "Enet" = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, "; using ", ncol(Train_set), " variables")
  
  args = list(
    "Train_set" = Train_set,
    "Train_label" = Train_label,
    "mode" = mode,
    "classVar" = classVar
  )
  
  args = c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name), args = args)
  
  if (mode == "Variable"){
    message(length(obj), " variables retained;\n")
  } else {
    message("\n")
  }
  
  return(obj)
}

standarize.fun <- function(indata, centerFlag, scaleFlag){
  scale(indata, center = centerFlag, scale = scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  
  if (is.null(cohort)){
    data <- list(data)
    names(data) = "training"
  } else {
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)){
    centerFlags = FALSE
    message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags) == 1){
    centerFlags = rep(centerFlags, length(data))
    message("Set centerFlags for all cohorts as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data)
    message("Match centerFlags with cohorts by order\n")
  }
  
  if (is.null(scaleFlags)){
    scaleFlags = FALSE
    message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags) == 1){
    scaleFlags = rep(scaleFlags, length(data))
    message("Set scaleFlags for all cohorts as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data)
    message("Match scaleFlags with cohorts by order\n")
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

ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
    "glm" = names(coef(fit)),
    "svm.formula" = fit$subFeature,
    "train" = fit$coefnames,
    "glmboost" = names(coef(fit)[abs(coef(fit)) > 0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
    "rfsrc" = names(which(fit$importance[, 1] > 0.01)),
    "gbm" = rownames(summary.gbm(fit, plotit = FALSE))[summary.gbm(fit, plotit = FALSE)$rel.inf > 0],
    "xgb.Booster" = fit$subFeature,
    "naiveBayes" = fit$subFeature
  ))
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = predict(fit, type = "response", as.matrix(new_data)),
    "glm" = predict(fit, type = "response", as.data.frame(new_data)),
    "svm.formula" = {
      prob <- attr(predict(fit, as.data.frame(new_data), probability = TRUE), "probabilities")
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    },
    "train" = {
      prob <- predict(fit, new_data, type = "prob")
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    },
    "glmboost" = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc" = {
      pred <- predict(fit, as.data.frame(new_data))$predicted
      if (is.matrix(pred)) pred[, "1"] else pred
    },
    "gbm" = predict(fit, type = "response", as.data.frame(new_data)),
    "xgb.Booster" = predict(fit, as.matrix(new_data)),
    "naiveBayes" = {
      prob <- predict(object = fit, type = "raw", newdata = new_data)
      if ("1" %in% colnames(prob)) prob[, "1"] else prob[, 1]
    }
  ))
  
  RS <- as.numeric(RS)
  
  if (length(RS) != nrow(new_data)) {
    warning(paste0("Prediction length mismatch in model [", class(fit)[1], "] → ",
                   "RS length: ", length(RS), " vs samples: ", nrow(new_data)))
    RS <- rep(NA, nrow(new_data))
  }
  
  names(RS) <- rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet" = predict(fit, type = "class", as.matrix(new_data)),
    "glm" = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "svm.formula" = predict(fit, as.data.frame(new_data), decision.values = TRUE),
    "train" = predict(fit, new_data, type = "raw"),
    "glmboost" = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "rfsrc" = predict(fit, as.data.frame(new_data))$class,
    "gbm" = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "xgb.Booster" = ifelse(predict(fit, as.matrix(new_data)) > 0.5, "1", "0"),
    "naiveBayes" = predict(object = fit, type = "class", newdata = new_data)
  ))
  
  label <- as.character(label)
  
  if (length(label) != nrow(new_data)) {
    warning(paste0("PredictClass length mismatch in model [", class(fit)[1], "] → ",
                   "predicted: ", length(label), " vs actual: ", nrow(new_data)))
    label <- rep(NA, nrow(new_data))
  }
  
  names(label) <- rownames(new_data)
  return(label)
}

RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){
  
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column."))
  }
  
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(
      Train_set[, fit$subFeature],
      Test_set[, fit$subFeature]
    )
    
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    
    Test_label <- rbind.data.frame(
      Train_label[, c(cohortVar, classVar)],
      Test_label[, c(cohortVar, classVar)]
    )
    
    Test_label[, 1] <- factor(
      Test_label[, 1],
      levels = c(unique(Train_label[, cohortVar]),
                 setdiff(unique(Test_label[, cohortVar]), unique(Train_label[, cohortVar])))
    )
  } else {
    new_data <- Test_set[, fit$subFeature]
  }
  
  RS <- suppressWarnings(as.numeric(CalPredictScore(fit = fit, new_data = new_data)))
  
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[, cohortVar])
  
  unlist(lapply(Predict.out, function(data){
    if (!is.numeric(data$RS) || all(is.na(data$RS))) {
      warning("AUC skipped: prediction RS invalid.")
      return(NA)
    }
    if (length(unique(na.omit(data[[classVar]]))) < 2) {
      warning("AUC skipped: only one class in this cohort.")
      return(NA)
    }
    tryCatch({
      auc(suppressMessages(roc(response = data[[classVar]], predictor = as.numeric(data$RS))))
    }, error = function(e){
      warning(paste0("AUC failed: ", e$message))
      return(NA)
    })
  }))
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows){
  
  col_ha = columnAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    show_annotation_name = FALSE
  )
  
  row_ha = rowAnnotation(
    bar = anno_barplot(
      avg_Cindex,
      bar_width = 0.8,
      border = FALSE,
      gp = gpar(fill = barCol, col = NA),
      add_numbers = TRUE,
      numbers_offset = unit(-10, "mm"),
      axis_param = list("labels_rot" = 0),
      numbers_gp = gpar(fontsize = 9, col = "white"),
      width = unit(3, "cm")
    ),
    show_annotation_name = FALSE
  )
  
  Heatmap(
    as.matrix(Cindex_mat),
    name = "AUC",
    right_annotation = row_ha,
    top_annotation = col_ha,
    col = c("#4195C1", "#FFFFFF", "#FFBC90"),
    rect_gp = gpar(col = "black", lwd = 1),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
    height = unit(cellheight * nrow(Cindex_mat), "cm"),
    column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(
        label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
        x, y, gp = gpar(fontsize = 10)
      )
    }
  )
}

quiet <- function(..., messages = FALSE, cat = FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

#########################################################################
# Read training and test data
#########################################################################

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

#########################################################################
# Add cohort information to the external test set
#########################################################################

GSE14811_pdata = read.csv("~/path/to/your/directory/GSE14811/GSE14811.csv")
GSE14811_pdata$Cohort = "GSE14811"
GSE14811_pdata = GSE14811_pdata[, c(1, 41)]

GSE54236_pdata = read.csv("~/path/to/your/directory/GSE54236/GSE54236.csv")
GSE54236_pdata$Cohort = "GSE54236"
GSE54236_pdata = GSE54236_pdata[, c(1, 47)]

GSE76427_pdata = read.csv("~/path/to/your/directory/GSE76427/GSE76427.csv")
GSE76427_pdata$Cohort = "GSE76427"
GSE76427_pdata = GSE76427_pdata[, c(1, 54)]

pdata_all = rbind(GSE14811_pdata, GSE54236_pdata, GSE76427_pdata)

Test_class = merge(Test_class, pdata_all, by.x = 0, by.y = 1)
rownames(Test_class) = Test_class$Row.names
Test_class = Test_class[, -1]
Test_class = Test_class[, c("Cohort", "Type")]

#########################################################################
# Keep shared genes and standardize data
#########################################################################

comgene <- intersect(colnames(Train_expr), colnames(Test_expr))

Train_expr <- as.matrix(Train_expr[, comgene])
Test_expr  <- as.matrix(Test_expr[, comgene])

Train_set = scaleData(data = Train_expr, centerFlags = TRUE, scaleFlags = TRUE)

set.seed(123)
noise <- matrix(
  rnorm(n = nrow(Train_set) * ncol(Train_set), mean = 0, sd = 0.01),
  nrow = nrow(Train_set),
  ncol = ncol(Train_set)
)
Train_set <- Train_set + noise

Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = TRUE, scaleFlags = TRUE)

noise_test <- matrix(
  rnorm(n = nrow(Test_set) * ncol(Test_set), mean = 0, sd = 0.01),
  nrow = nrow(Test_set),
  ncol = ncol(Test_set)
)
Test_set <- Test_set + noise_test

names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = TRUE, scaleFlags = TRUE)

#########################################################################
# Read method list
#########################################################################

methodRT <- read.table("refer.txt", header = TRUE, sep = "\t", check.names = FALSE)
methods = methodRT$Model
methods <- gsub("-| ", "", methods)

#########################################################################
# Prepare model parameters
#########################################################################

classVar = "Type"
Variable = colnames(Train_set)

preTrain.method = strsplit(methods, "\\+")
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method))

###################### Variable selection stage ######################

preTrain.var <- list()
set.seed(seed = 123)

for (method in preTrain.method){
  preTrain.var[[method]] = RunML(
    method = method,
    Train_set = Train_set,
    Train_label = Train_class,
    mode = "Variable",
    classVar = classVar
  )
}

preTrain.var[["simple"]] <- colnames(Train_set)

###################### Model fitting stage ######################

model <- list()
set.seed(seed = 123)

Train_set_bk <- Train_set
min.selected.var = 2

for (method in methods) {
  cat(match(method, methods), ":", method, "\n")
  
  parts <- strsplit(method, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  
  vars <- preTrain.var[[parts[1]]]
  
  if (length(vars) <= min.selected.var) {
    message("SKIP ", parts[1], " → only ", length(vars), " variables\n")
    next
  }
  
  ts <- Train_set_bk[, vars, drop = FALSE]
  
  fit <- RunML(
    method = parts[2],
    Train_set = ts,
    Train_label = Train_class,
    mode = "Model",
    classVar = classVar
  )
  
  if (length(ExtractVar(fit)) <= min.selected.var) {
    message("DROP ", method, " → only ", length(ExtractVar(fit)), " vars after modelling\n")
  } else {
    model[[method]] <- fit
  }
}

Train_set <- Train_set_bk

saveRDS(model, "model.MLmodel.rds")

#########################################################################
# Optional final logistic model
#########################################################################

FinalModel <- c("panML", "multiLogistic")[2]

if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){
    tmp <- glm(
      formula = Train_class[[classVar]] ~ .,
      family = "binomial",
      data = as.data.frame(Train_set[, ExtractVar(fit)])
    )
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}

saveRDS(logisticmodel, "model.logisticmodel.rds")

#########################################################################
# Compute prediction scores for all models
#########################################################################

model <- readRDS("model.MLmodel.rds")
methodsValid <- names(model)

all(model[[method]]$subFeature %in% colnames(Train_set))
all(model[[method]]$subFeature %in% colnames(Test_set))

fit <- model[[1]]
nd <- rbind(Train_set[, fit$subFeature], Test_set[, fit$subFeature])
RS <- CalPredictScore(fit, new_data = nd)

length(RS)
length(rownames(nd))

for (method in methodsValid) {
  cat("Checking:", method, "\n")
  fit <- model[[method]]
  nd <- rbind(
    Train_set[, fit$subFeature, drop = FALSE],
    Test_set[, fit$subFeature, drop = FALSE]
  )
  
  tryCatch({
    RS <- CalPredictScore(fit = fit, new_data = nd)
    if (length(RS) != nrow(nd)) {
      cat("Length mismatch →", "RS:", length(RS), "vs", "samples:", nrow(nd), "\n")
    } else {
      cat("Passed:", method, "\n")
    }
  }, error = function(e){
    cat("ERROR in", method, "→", e$message, "\n")
  })
}

RS_list <- list()
for (method in methodsValid){
  RS_list[[method]] <- CalPredictScore(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}

riskTab = as.data.frame(t(do.call(rbind, RS_list)))
riskTab = cbind(id = row.names(riskTab), riskTab)
write.table(riskTab, "model.riskMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#########################################################################
# Predict class labels for all models
#########################################################################

Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(
    fit = model[[method]],
    new_data = rbind.data.frame(Train_set, Test_set)
  )
}

Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
classTab = cbind(id = row.names(Class_mat), Class_mat)
write.table(classTab, "model.classMatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#########################################################################
# Extract selected features from each model
#########################################################################

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})

fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file = "model.genes.txt", sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

#########################################################################
# Calculate AUC for each model across cohorts
#########################################################################

table(Test_class$Cohort, Test_class$Type)

AUC_list <- list()
for (method in methodsValid){
  AUC_list[[method]] <- RunEval(
    fit = model[[method]],
    Test_set = Test_set,
    Test_label = Test_class,
    Train_set = Train_set,
    Train_label = Train_class,
    Train_name = "Train",
    cohortVar = "Cohort",
    classVar = classVar
  )
}

AUC_mat <- do.call(rbind, AUC_list)
aucTab = cbind(Method = row.names(AUC_mat), AUC_mat)
write.table(aucTab, "model.AUCmatrix.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#########################################################################
# Draw AUC heatmap
#########################################################################

AUC_mat <- read.table("model.AUCmatrix.txt", header = TRUE, sep = "\t",
                      check.names = FALSE, row.names = 1, stringsAsFactors = FALSE)

avg_AUC <- apply(AUC_mat, 1, mean)
avg_AUC <- sort(avg_AUC, decreasing = TRUE)
AUC_mat <- AUC_mat[names(avg_AUC), ]

fea_sel <- fea_list[[rownames(AUC_mat)[1]]]
avg_AUC <- as.numeric(format(avg_AUC, digits = 3, nsmall = 3))

CohortCol <- brewer.pal(n = ncol(AUC_mat), name = "Paired")
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1
cellheight = 0.5

hm <- SimpleHeatmap(
  Cindex_mat = AUC_mat,
  avg_Cindex = avg_AUC,
  CohortCol = CohortCol,
  barCol = "steelblue",
  cellwidth = cellwidth,
  cellheight = cellheight,
  cluster_columns = FALSE,
  cluster_rows = FALSE
)

pdf(file = "model.AUCheatmap.pdf",
    width = cellwidth * ncol(AUC_mat) + 6,
    height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm, heatmap_legend_side = "right", annotation_legend_side = "right")
dev.off()

#########################################################################
# Summarize the top 10 repeatedly selected genes and build importance matrix
#########################################################################

all_features <- unlist(fea_list, use.names = FALSE)
top_genes_tab <- sort(table(all_features), decreasing = TRUE)
top_genes <- names(top_genes_tab)[1:10]

importance_mat <- matrix(0, nrow = length(methodsValid), ncol = length(top_genes))
rownames(importance_mat) <- methodsValid
colnames(importance_mat) <- top_genes

for (m in methodsValid) {
  fit <- model[[m]]
  cls <- class(fit)[1]
  
  imp <- rep(0, length(top_genes))
  names(imp) <- top_genes
  
  if (cls == "lognet") {
    coefs <- coef(fit)
    coef_val <- as.numeric(coefs)
    names(coef_val) <- rownames(coefs)
    imp_val <- coef_val[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
    
  } else if (cls == "glm") {
    coefs <- coef(fit)
    coef_val <- as.numeric(coefs)
    names(coef_val) <- names(coefs)
    imp_val <- coef_val[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
    
  } else if (cls == "glmboost") {
    coefs <- coef(fit)
    imp_val <- as.numeric(coefs[top_genes])
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
    
  } else if (cls == "rfsrc") {
    if (!is.null(fit$importance)) {
      use_genes <- intersect(top_genes, rownames(fit$importance))
      if (length(use_genes) > 0) {
        imp_val <- fit$importance[use_genes, 1]
        imp[use_genes] <- imp_val
      }
    }
    
  } else if (cls == "gbm") {
    gbm_imp <- summary.gbm(fit, plotit = FALSE)
    vals <- gbm_imp$rel.inf
    names(vals) <- gbm_imp$var
    imp_val <- vals[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
    
  } else if (cls == "xgb.Booster") {
    imp_tab <- xgb.importance(feature_names = fit$subFeature, model = fit)
    vals <- imp_tab$Gain
    names(vals) <- imp_tab$Feature
    imp_val <- vals[top_genes]
    imp[!is.na(imp_val)] <- imp_val[!is.na(imp_val)]
  }
  
  importance_mat[m, ] <- imp
}

importance_mat[is.na(importance_mat)] <- 0

## you may ask:

## Why is my model performance on the test set different from the paper?

# 1) The preprocessing is not exactly the same as in the paper
# a. Was log2 transformation applied?
#    For example, the paper may have used: log2(x + 1)
# b. Is the normalization method consistent?
#    The normalization method I used is defined in this function: scaleData
# c. Was batch effect correction performed?
# d. How were missing values, duplicated genes, and low-expression genes handled?

# 2) Gene IDs / probes / features do not match
# a. If the same gene appears in multiple rows, how should they be merged?
# b. After mapping probes to gene symbols, should values be averaged or should the maximum be used?
# c. Differences in letter case
# d. Some gene names may use different annotation versions
# e. The number of intersecting genes between the training and test sets is different from the paper

# How can I compare different preprocessing strategies to see which one performs better?
# For each normalization strategy:
# perform 5-fold or 10-fold cross-validation on the training set
# calculate validation AUC for each fold
# then compare:
# the mean AUC
# the standard deviation of AUC

# The strategy with higher average AUC and smaller fluctuation is more worth keeping.
# Then use the new model to run prediction on the test set.

#!/usr/bin/env Rscript

# Classify resistant versus susceptible individuals.

suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("kernlab"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("randomForest"))
suppressPackageStartupMessages(library("C50"))
suppressPackageStartupMessages(library("penalizedLDA"))
# Loaded by caret when running certain algorithms, recommended to load before
# dplyr.
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("limma"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

set.seed(12345)

# Input
v <- readRDS(file.path(data_dir, "results-limma-voom.rds"))
anno <- read.delim(file.path(data_dir, "experiment-info-filtered.txt"),
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(colnames(v$E) == rownames(anno))
results <- readRDS(file.path(data_dir, "results-limma-stats.rds"))

# Functions

calc_f1 <- function(recall, precision) {
  # Calculate F1 score.
  #
  # recall - sensitivity, TP / (TP + FN)
  # precision - specificity, TN / (TN + FP)
  #
  # http://stackoverflow.com/a/36843900/2483477
  #
  f1 <- (2 * precision * recall) / (precision + recall)
  return(f1)
}

extract_kappa <- function(results, bestTune) {
  # Extract the Kappa statistic from the model with the best results
  #
  # results - output from caret::train
  # bestTune - output from caret::train
  #
  results_final <- results
  for (parameter in names(bestTune)) {
    param_filter <- results_final[, parameter] == bestTune[1, parameter]
    results_final <- results_final[param_filter, ]
  }
  stopifnot(nrow(results_final) == 1)
  return(results_final$Kappa)
}

# Algorithms to try
algos <- c("svmLinear", "glmnet", "rf")
# algos <- c("svmLinear", "glmnet", "C5.0", "rf", "PenalizedLDA")
# getModelInfo(algos[3])[[algos[3]]]
# svalue cutoffs
sval_cut <- seq(0.05, 0.25, by = 0.05)

# Use leave-one-out-cross-validation
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePred = "final")

predictions <- vector(length = length(algos), mode = "list")
names(predictions) <- algos

for (alg in algos) {
  checkInstall(getModelInfo(alg)$library)
  predictions[[alg]] <- vector(length = length(sval_cut), mode = "list")
  names(predictions[[alg]]) <- paste0("s", sval_cut)
  for (cutoff in sval_cut) {
    cutoff_name <- paste0("s", cutoff)
    genes_for_classifer <- v$E[results[["diff_before"]]$svalue < cutoff, ]
    genes_for_classifer <- t(genes_for_classifer)
    genes_for_classifer <- as.data.frame(genes_for_classifer)
    genes_for_classifer$status <- anno$status
    genes_for_classifer <- genes_for_classifer[anno$treatment == "none", ]
    fit <- train(status ~ ., data = genes_for_classifer, method = alg,
                 trControl = ctrl, metric = "Kappa")
    fit$num_genes <- ncol(genes_for_classifer) - 1
    data <- fit$pred[, "pred"]
    reference <- fit$pred[, "obs"]
    fit$recall <- sensitivity(data, reference, positive = "tb")
    fit$precision <- specificity(data, reference, negative = "contact")
    fit$f1 <- calc_f1(fit$recall, fit$precision)
    fit$kappa <- extract_kappa(fit$results, fit$bestTune)
    fit$pred$id <- rownames(genes_for_classifer)[fit$pred$rowIndex]
    predictions[[alg]][[cutoff_name]] <- fit
  }
}

# Save results
saveRDS(predictions, file.path(data_dir, "classifier-predictions.rds"))

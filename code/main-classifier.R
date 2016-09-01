#!/usr/bin/env Rscript

# Classify resistant versus susceptible individuals.

suppressPackageStartupMessages(library("caret"))
suppressPackageStartupMessages(library("kernlab"))
suppressPackageStartupMessages(library("glmnet"))
suppressPackageStartupMessages(library("randomForest"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

set.seed(12345)

# Input ------------------------------------------------------------------------

regressed <- read.delim(file.path(data_dir, "combined-regressed.txt"),
                row.names = 1, check.names = FALSE)
anno_comb <- read.delim(file.path(data_dir, "combined-annotation.txt"),
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(colnames(regressed) == rownames(anno_comb))
results <- readRDS(file.path(data_dir, "combined-limma.rds"))
stopifnot(rownames(regressed) == rownames(results[["suscept_ni"]]))

# Split into training and test sets. Only include noninfected samples.
training <- regressed[, anno_comb$study == "current" &
                        anno_comb$treatment == "noninf"]
anno_train <- anno_comb[anno_comb$study == "current" &
                          anno_comb$treatment == "noninf", ]
training <- data.frame(status = anno_train$status, t(training))

testing <- regressed[, anno_comb$study == "lbb2012" &
                       anno_comb$treatment == "noninf"]
anno_test <- anno_comb[anno_comb$study == "lbb2012" &
                       anno_comb$treatment == "noninf", ]
testing <- data.frame(t(testing))

# All potential genes for classifier
genes_all <- rownames(results[["status_ni"]])

# Scoring functions --------------------------------------------------------------------

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

# Model ------------------------------------------------------------------------

algos <- c("svmLinear", "glmnet", "rf")
# getModelInfo(algos[3])[[algos[3]]]

# qvalue cutoffs
qval_cut <- seq(0.05, 0.25, by = 0.05)

# Use leave-one-out-cross-validation
ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePred = "final")

predictions <- vector(length = length(algos), mode = "list")
names(predictions) <- algos
predictions_lbb <- vector(length = length(algos), mode = "list")
names(predictions_lbb) <- algos

for (alg in algos) {
  checkInstall(getModelInfo(alg)$library)
  predictions[[alg]] <- vector(length = length(qval_cut), mode = "list")
  names(predictions[[alg]]) <- sprintf("q%.2f", qval_cut)
  predictions_lbb[[alg]] <- vector(length = length(qval_cut), mode = "list")
  names(predictions_lbb[[alg]]) <- sprintf("q%.2f", qval_cut)
  for (cutoff in qval_cut) {
    cutoff_name <- sprintf("q%.2f", cutoff)
    genes_for_classifier <- genes_all[results[["status_ni"]]$qvalue < cutoff]
    training_sub <- training[, c("status", genes_for_classifier)]
    fit <- train(status ~ ., data = training_sub, method = alg,
                 trControl = ctrl, metric = "Kappa")
    fit$num_genes <- length(genes_for_classifier)
    data <- fit$pred[, "pred"]
    reference <- fit$pred[, "obs"]
    fit$recall <- sensitivity(data, reference, positive = "suscep")
    fit$precision <- specificity(data, reference, negative = "resist")
    fit$f1 <- calc_f1(fit$recall, fit$precision)
    fit$kappa <- extract_kappa(fit$results, fit$bestTune)
    fit$pred$id <- rownames(training_sub)[fit$pred$rowIndex]
    fit$separation <- mean(fit$pred[fit$pred[, "obs"] == "resist", "resist"]) -
                      mean(fit$pred[fit$pred[, "obs"] == "suscep", "resist"])
    predictions[[alg]][[cutoff_name]] <- fit
    # Predict in Barreiro et al., 2012
    testing_sub <- testing[, genes_for_classifier]
    stopifnot(ncol(testing_sub) == ncol(training_sub) - 1)
    result_test <- predict(fit, testing_sub, type = "prob")
    rownames(result_test) <- rownames(testing_sub)
    predictions_lbb[[alg]][[cutoff_name]] <- result_test
  }
}

# Save results -----------------------------------------------------------------

saveRDS(predictions, file.path(data_dir, "classifier-predictions.rds"))
saveRDS(predictions_lbb, file.path(data_dir, "classifier-predictions-lbb.rds"))

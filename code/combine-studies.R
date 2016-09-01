#!/usr/bin/env Rscript

# Combine expression data from the current study with Barreiro et al., 2012.

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("ashr"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

# Input ------------------------------------------------------------------------

# Current study
counts <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                     check.names = FALSE, row.names = 1)
anno <- read.delim(file.path(data_dir, "experiment-info-filtered.txt"),
                   stringsAsFactors = FALSE, row.names = 1)
anno$study <- "current"
stopifnot(colnames(counts) == rownames(anno))

# TMM-normalized log2 counts per million. The default prior counts is 0.25.
norm_factors <- calcNormFactors(counts)
counts_cpm <- cpm(counts, lib.size = colSums(counts) * norm_factors,
                  log = TRUE)

# lbb2012
load(file.path(data_dir, "Exp_final_Batch_corrected.Rdata"))
lbb2012 <- Exp_final_Batch_corrected[, -1:-2]
rownames(lbb2012) <- Exp_final_Batch_corrected$Ensembl_ID
id <- sub("Mtb", "", colnames(lbb2012))
id <- sub("_", "", id)
id <- sub("\\.neg", "", id)
id <- sub("\\.pos", "", id)
id <- sub("\\.TB", "", id)
id <- sub("\\.1", "", id)
id <- sub("\\.2", "", id)
id <- as.numeric(id)
stopifnot(id[1:65] == id[66:(65*2)])
anno_lbb <- data.frame(individual = sprintf("u%02d", id),
                       status = "unknow", study= "lbb2012",
                       treatment = ifelse(grepl("neg", colnames(lbb2012)),
                                              "noninf", "infect"))
colnames(lbb2012) <- paste(anno_lbb$individual,
                           anno_lbb$status,
                           anno_lbb$treatment, sep = "-")
rownames(anno_lbb) <- colnames(lbb2012)
stopifnot(ncol(lbb2012) == 65 * 2, table(id) == 2)

# Intersect gene lists ---------------------------------------------------------

common_genes <- intersect(rownames(counts_cpm), rownames(lbb2012))
# length(common_genes)
counts_cpm_sub <- counts_cpm[rownames(counts_cpm) %in% common_genes, ]
counts_cpm_sub <- counts_cpm_sub[order(rownames(counts_cpm_sub)), ]
lbb2012_sub <- lbb2012[rownames(lbb2012) %in% common_genes, ]
lbb2012_sub <- lbb2012_sub[order(rownames(lbb2012_sub)), ]
stopifnot(nrow(counts_cpm_sub) == length(common_genes),
          rownames(counts_cpm_sub) == rownames(lbb2012_sub))

# Combine the two data sets
exp_combined <- cbind(counts_cpm_sub, lbb2012_sub)
anno_comb <- rbind(anno[, c("individual", "status", "treatment", "study")],
                   anno_lbb)
stopifnot(colnames(exp_combined) == rownames(anno_comb))

# Normalize gene expression distributions---------------------------------------

# Normalize to standard normal: N(0, 1) with qqnorm
normalize <- function(x) {qqnorm(x, plot.it = FALSE)$x}

exp_combined_norm <- apply(exp_combined, 2, normalize)
rownames(exp_combined_norm) <- rownames(exp_combined)

# Regress PC1 ------------------------------------------------------------------

pca <- prcomp(t(exp_combined_norm), scale. = TRUE)
variances <- pca$sdev^2
explained <- variances / sum(variances)
pca_data <- cbind(anno_comb, pca$x[, 1:10])
# plot(pca, main = "Variance per PC")

regressed <- removeBatchEffect(exp_combined_norm, covariates = pca_data$PC1)

# Perform PCA on regressed data
pca2 <- prcomp(t(regressed), scale. = TRUE)
variances2 <- pca2$sdev^2
explained2 <- variances2 / sum(variances2)
pca_data2 <- cbind(anno_comb, pca2$x[, 1:10])
# plot(pca2, main = "Variance per PC")

# Redo limma analysis with adjusted, non-count data ----------------------------

# Extract the adjusted data for the current study
current_adj <- regressed[, anno_comb$study == "current"]
stopifnot(colnames(current_adj) == colnames(counts_cpm),
          colnames(current_adj) == rownames(anno))

# Design matrix:
# Use the standard treatment-contrasts parametrization. See Ch. 9 of limma
# User's Guide.
anno$status <- factor(anno$status, levels = c("resist", "suscep"))
anno$treatment <- factor(anno$treatment, levels = c("noninf", "infect"))
design <- model.matrix(~treatment + status + status:treatment,
                       data = anno)
colnames(design)[1] <- "Intercept"
colnames(design) <- gsub("treatment", "", colnames(design))
colnames(design) <- gsub("status", "", colnames(design))
colnames(design) <- gsub(":", ".", colnames(design))
# colSums(design)

# Fit model:
# Model individual as a random effect.
corfit <- duplicateCorrelation(current_adj, design, block = anno$individual)
# corfit$consensus
fit <- lmFit(current_adj, design, block = anno$individual,
             correlation = corfit$consensus)

# Fit contrasts
contrasts_mat <- makeContrasts(
  # Difference between susceptible and resistant individuals in the noninfected state
  status_ni = suscep,
  # Difference between susceptible and resistant individuals in the infected state
  status_ii = suscep + infect.suscep,
  # Effect of treatment in resistant individuals
  treat_resist = infect,
  # Effect of treatment in susceptible individuals
  treat_suscep = infect + infect.suscep,
  # Difference in effect of treatment between susceptible and resistant
  # individuals
  interact = infect.suscep,
  levels = design)
fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2)

# Get model results:
get_results <- function(x, number = nrow(x$coefficients), sort.by = "none",
                        ...) {
  # x - object MArrayLM from eBayes output
  # ... - additional arguments passed to topTable
  stopifnot(class(x) == "MArrayLM")
  results <- topTable(x, number = number, sort.by = sort.by, ...)
  return(results)
}

run_ash <- function(x, coef) {
  # Perform multiple testing correction with adaptive shrinkage (ASH)
  #
  # x - object MArrayLM from eBayes output
  # coef - coefficient tested by eBayes
  stopifnot(class(x) == "MArrayLM", coef %in% colnames(x$coefficients))
  result <- ash(betahat = x$coefficients[, coef],
                sebetahat = x$stdev.unscaled[, coef] * sqrt(x$s2.post),
                df = x$df.total)
  return(result)
}

tests <- colnames(fit2$coefficients)
results <- vector(length = length(tests), mode = "list")
names(results) <- tests

for (test in tests) {
  # Extract limma results
  results[[test]] <- get_results(fit2, coef = test)
  # Add mutliple testing correction with ASH
  ash_log <- capture.output(output_ash <- run_ash(fit2, coef = test))
  results[[test]] <- cbind(results[[test]], lfsr = output_ash$lfsr,
                           lfdr = output_ash$lfdr, qvalue = output_ash$qvalue,
                           svalue = output_ash$svalue)
}

# Save for creating figures ----------------------------------------------------

# Combined data
write.table(exp_combined, file.path(data_dir, "combined-raw.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(anno_comb, file.path(data_dir, "combined-annotation.txt"),
            quote = FALSE, sep = "\t", col.names = NA)

# Normalized data
write.table(exp_combined_norm, file.path(data_dir, "combined-normalized.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(pca_data, file.path(data_dir, "combined-pca.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
saveRDS(explained, file.path(data_dir, "combined-pca-explained.txt"))

# Regressed data
write.table(regressed, file.path(data_dir, "combined-regressed.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(pca_data2, file.path(data_dir, "combined-pca-regressed.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
saveRDS(explained2, file.path(data_dir, "combined-pca-explained-regressed.txt"))

# limma results
saveRDS(results, file.path(data_dir, "combined-limma.rds"))

# ## ----partition-training--------------------------------------------------
# training <- regressed[, anno$study == "current"]
#
# ## ----loocv---------------------------------------------------------------
# ctrl <- trainControl(method = "LOOCV", classProbs = TRUE, savePred = "final")
#
# ## ----genes-for-classifier------------------------------------------------
# g <- rownames(training)[results[["diff_before"]]$svalue < .05]
# length(g)
# genes_for_classifer <- training[g, ]
# genes_for_classifer <- t(genes_for_classifer)
# genes_for_classifer <- as.data.frame(genes_for_classifer)
# genes_for_classifer$status <- anno$status[anno$study == "current"]
# genes_for_classifer <- genes_for_classifer[anno$treatment[anno$study == "current"] == "none", ]
#
# ## ----glmnet--------------------------------------------------------------
# fit <- train(status ~ ., data = genes_for_classifer, method = "glmnet",
#              trControl = ctrl, metric = "Kappa")
# fit
#
# ## ----vis-----------------------------------------------------------------
# pred_training <- fit$pred
# pred_training <- pred_training[order(pred_training$contact), ]
# plot(pred_training$contact, col = ifelse(pred_training$obs == "tb",
#                                          "red", "black"),
#      ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of current study")
# abline(h = 0.8, col = "blue")
#
# ## ----training-barplot----------------------------------------------------
# pred_training$classification <- NA
# for (i in 1:nrow(pred_training)) {
#   if (pred_training$contact[i] >= 0.9) {
#     pred_training$classification[i] <- "resistant"
#   } else if (pred_training$contact[i] <= 0.1) {
#     pred_training$classification[i] <- "susceptible"
#   } else {
#     pred_training$classification[i] <- "ambiguous"
#   }
# }
# pred_training$classification <- factor(pred_training$classification,
#                                        levels = c("resistant", "ambiguous", "susceptible"))
# barplot(table(pred_training$classification), xlab = "Classification",
#         ylab = "Number of individuals",
#         main = "Interpretation of classifier results\nof current study")
#
# ## ----test----------------------------------------------------------------
# testing <- regressed[g, anno$study == "lbb2012"]
# testing <- t(testing)
# testing <- testing[anno$treatment[anno$study == "lbb2012"] == "none", ]
# fit_lbb <- predict(fit, testing, type = "prob")
# fit_lbb
#
# ## ----vis-lbb-------------------------------------------------------------
# pred_testing <- fit_lbb
# pred_testing <- pred_testing[order(pred_testing$contact), ]
# plot(pred_testing$contact, ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of Barreiro et al., 2012")
# abline(h = 0.8, col = "blue")
#
# ## ----testing-barplot-----------------------------------------------------
# pred_testing$classification <- NA
# for (i in 1:nrow(pred_testing)) {
#   if (pred_testing$contact[i] >= 0.9) {
#     pred_testing$classification[i] <- "resistant"
#   } else if (pred_testing$contact[i] <= 0.1) {
#     pred_testing$classification[i] <- "susceptible"
#   } else {
#     pred_testing$classification[i] <- "ambiguous"
#   }
# }
# pred_testing$classification <- factor(pred_testing$classification,
#                                        levels = c("resistant", "ambiguous", "susceptible"))
# barplot(table(pred_testing$classification), xlab = "Classification",
#         ylab = "Number of individuals",
#         main = "Interpretation of classifier results\nof Barreiro et al., 2012")
#
# ## ----compare-exp-level---------------------------------------------------
# training_medians <- apply(genes_for_classifer[-ncol(genes_for_classifer)],
#                           2, median)
# testing_medians <- apply(testing, 2, median)
# plot(training_medians, testing_medians,
#      xlab = "Current study", ylab = "Barreiro et al., 2012",
#      main = "Comparing median expression levels\nof genes used in classifer")
# abline(0, 1, col = "red")
#
# ## ----hclust--------------------------------------------------------------
# clust_current <- genes_for_classifer[, -ncol(genes_for_classifer)]
# labels <- ifelse(grepl("tb", genes_for_classifer$status),
#                  "susceptible", "resistant")
# labels[labels == "susceptible"] <- paste(labels[labels == "susceptible"],
#                                          1:sum(labels == "susceptible"),
#                                          sep = ".")
# labels[labels == "resistant"] <- paste(labels[labels == "resistant"],
#                                          1:sum(labels == "resistant"),
#                                          sep = ".")
# rownames(clust_current) <- labels
# d1 <- dist(clust_current)
# h1 <- hclust(d1)
# ggdendrogram(h1, rotate = TRUE)
# clust_lbb <- testing
# labels_lbb <- paste0("Luis.", 1:nrow(clust_lbb))
# rownames(clust_lbb) <- labels_lbb
# clust_all <- rbind(clust_current, clust_lbb)
# d2 <- dist(clust_all)
# h2 <- hclust(d2)
# ggdendrogram(h2, rotate = TRUE)
#
# ## ----multi-eps-----------------------------------------------------------
# postscript("../figure/classifier.eps", width = 14, height = 7)
# par(mfrow = c(1, 2))
# plot(pred_training$contact, col = ifelse(pred_training$obs == "tb",
#                                          "red", "black"),
#      ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of current study")
# abline(h = 0.8, col = "blue")
# plot(pred_testing$contact, ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of Barreiro et al., 2012")
# abline(h = 0.8, col = "blue")
# invisible(dev.off())
#
# ## ----multi-png-----------------------------------------------------------
# png("../figure/classifier.png", width = 14, height = 7,
#     units = "in", res = 72)
# par(mfrow = c(1, 2))
# plot(pred_training$contact, col = ifelse(pred_training$obs == "tb",
#                                          "red", "black"),
#      ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of current study")
# abline(h = 0.8, col = "blue")
# plot(pred_testing$contact, ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of Barreiro et al., 2012")
# abline(h = 0.8, col = "blue")
# invisible(dev.off())
#
# ## ----multi-pdf-----------------------------------------------------------
# pdf("../figure/classifier.pdf", width = 14, height = 7)
# par(mfrow = c(1, 2))
# plot(pred_training$contact, col = ifelse(pred_training$obs == "tb",
#                                          "red", "black"),
#      ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of current study")
# abline(h = 0.8, col = "blue")
# plot(pred_testing$contact, ylim = c(0, 1), xlab = "Individual",
#      ylab = "Predicted probability of TB resistance",
#      main = "Classifier results of Barreiro et al., 2012")
# abline(h = 0.8, col = "blue")
# invisible(dev.off())


#!/usr/bin/env Rscript

# Combine expression data from the current study with Thuong et al., 2008.

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

# thuong2008
thuong2008 <- read.delim(file.path(data_dir, "thuong2008-exp.txt"),
                         check.names = FALSE, row.names = 1)
anno_thuong <- read.delim(file.path(data_dir, "thuong2008-anno.txt"),
                          stringsAsFactors = FALSE, row.names = 1)
anno_thuong$study <- "thuong2008"
stopifnot(colnames(thuong2008) == rownames(anno_thuong))

# Intersect gene lists ---------------------------------------------------------

common_genes <- intersect(rownames(counts_cpm), rownames(thuong2008))
# length(common_genes)
counts_cpm_sub <- counts_cpm[rownames(counts_cpm) %in% common_genes, ]
counts_cpm_sub <- counts_cpm_sub[order(rownames(counts_cpm_sub)), ]
thuong2008_sub <- thuong2008[rownames(thuong2008) %in% common_genes, ]
thuong2008_sub <- thuong2008_sub[order(rownames(thuong2008_sub)), ]
stopifnot(nrow(counts_cpm_sub) == length(common_genes),
          rownames(counts_cpm_sub) == rownames(thuong2008_sub))

# Combine the two data sets
exp_combined <- cbind(counts_cpm_sub, thuong2008_sub)
anno_cols <- c("individual", "status", "treatment", "study")
anno_comb <- rbind(anno[, anno_cols], anno_thuong[, anno_cols])
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
# plot(pca$x[, 1], pca$x[, 2],
#      col = ifelse(anno_comb$study == "current", "red", "blue"),
#      pch = ifelse(anno_comb$treatment == "noninf", 1, 19))

regressed <- removeBatchEffect(exp_combined_norm, covariates = pca_data$PC1)

# Perform PCA on regressed data
pca2 <- prcomp(t(regressed), scale. = TRUE)
variances2 <- pca2$sdev^2
explained2 <- variances2 / sum(variances2)
pca_data2 <- cbind(anno_comb, pca2$x[, 1:10])
# plot(pca2, main = "Variance per PC")
# plot(pca2$x[, 1], pca2$x[, 2],
#      col = ifelse(anno_comb$study == "current", "red", "blue"),
#      pch = ifelse(anno_comb$treatment == "noninf", 1, 19))

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
write.table(exp_combined, file.path(data_dir, "combined-raw-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(anno_comb, file.path(data_dir, "combined-annotation-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)

# Normalized data
write.table(exp_combined_norm, file.path(data_dir, "combined-normalized-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(pca_data, file.path(data_dir, "combined-pca-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
saveRDS(explained, file.path(data_dir, "combined-pca-explained-thuong.txt"))

# Regressed data
write.table(regressed, file.path(data_dir, "combined-regressed-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(pca_data2, file.path(data_dir, "combined-pca-regressed-thuong.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
saveRDS(explained2, file.path(data_dir, "combined-pca-explained-regressed-thuong.txt"))

# limma results
saveRDS(results, file.path(data_dir, "combined-limma-thuong.rds"))

#!/usr/bin/env Rscript

# Perform differential expression analysis with limma.

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("ashr"))

# Input ------------------------------------------------------------------------

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

counts <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                     check.names = FALSE, row.names = 1)
anno <- read.delim(file.path(data_dir, "experiment-info-filtered.txt"),
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(colnames(counts) == rownames(anno))

# Design matrix ----------------------------------------------------------------

# Use the standard treatment-contrasts parametrization. See Ch. 9 of limma
# User's Guide.
anno$status <- factor(anno$status, levels = c("contact", "tb"))
anno$treatment <- factor(anno$treatment, levels = c("none", "infected"))
design <- model.matrix(~treatment + status + status:treatment,
                       data = anno)
colnames(design)[1] <- "Intercept"
colnames(design) <- gsub("treatment", "", colnames(design))
colnames(design) <- gsub("status", "", colnames(design))
colnames(design) <- gsub(":", ".", colnames(design))
# colSums(design)

# Fit model --------------------------------------------------------------------

# Model individual as a random effect.
# Recommended to run both voom and duplicateCorrelation twice.
# https://support.bioconductor.org/p/59700/#67620

y <- DGEList(counts)
y <- calcNormFactors(y)
v1 <- voom(y, design)
corfit1 <- duplicateCorrelation(v1, design, block = anno$individual)
# corfit1$consensus
v2 <- voom(y, design, block = anno$individual, correlation = corfit1$consensus)
corfit2 <- duplicateCorrelation(v2, design, block = anno$individual)
# corfit2$consensus
fit <- lmFit(v2, design, block = anno$individual,
             correlation = corfit2$consensus)

# Fit contrasts ----------------------------------------------------------------

contrasts_mat <- makeContrasts(
  # Difference between susceptible and resistant individuals before treatment
  diff_before = tb,
  # Difference between susceptible and resistant individuals after treatment
  diff_after = tb + infected.tb,
  # Effect of treatment in resistant individuals
  treat_resist = infected,
  # Effect of treatment in susceptible individuals
  treat_suscept = infected + infected.tb,
  # Difference in effect of treatment between susceptible and resistant
  # individuals
  diff_treat = infected.tb,
  levels = design)
fit2 <- contrasts.fit(fit, contrasts_mat)
fit2 <- eBayes(fit2)

# Get model results

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
  output_ash <- run_ash(fit2, coef = test)
  results[[test]] <- cbind(results[[test]], lfsr = output_ash$lfsr,
                           lfdr = output_ash$lfdr, qvalue = output_ash$qvalue,
                           svalue = output_ash$svalue)
}

# Save results from analysis with limma and ash.
saveRDS(v2, file.path(data_dir, "results-limma-voom.rds"))
saveRDS(fit2, file.path(data_dir, "results-limma-fit.rds"))
saveRDS(results, file.path(data_dir, "results-limma-stats.rds"))

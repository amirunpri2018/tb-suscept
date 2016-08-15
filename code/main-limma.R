#!/usr/bin/env Rscript

# Perform differential expression analysis with limma.

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("ashr"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

# Input
counts <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                     check.names = FALSE, row.names = 1)
anno <- read.delim(file.path(data_dir, "experiment-info-filtered.txt"),
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(colnames(counts) == rownames(anno))

# Design matrix Use the group-means parametrization: The mean of each group is
# estimated, and the differences are tested with contasts. See Ch. 9 of limma
# User's Guide.
combined <- paste(anno$status, anno$treatment, sep = ".")
combined <- factor(combined, levels = c("contact.none", "tb.none",
                                        "contact.infected", "tb.infected"))
design <- model.matrix(~0 + combined)
colnames(design) <- gsub("combined", "", colnames(design))
# colnames(design)

# Fit model
y <- DGEList(counts)
y <- calcNormFactors(y)
v <- voom(y, design)
# Model individual as a random effect
corfit <- duplicateCorrelation(v, design, block = anno$individual)
# corfit$consensus
v <- voom(y, design, block = anno$individual, correlation = corfit$consensus)
fit <- lmFit(v, design, block = anno$individual, correlation = corfit$consensus)

# Fit contrasts
contrasts_mat <- makeContrasts(
  # Difference between susceptible and resistant individual before treatment
  diff_before = tb.none - contact.none,
  # Difference between susceptible and resistant individual after treatment
  diff_after = tb.infected - contact.infected,
  # Effect of treatment in resistant individuals
  treat_resist = contact.infected - contact.none,
  # Effect of treatment in susceptible individuals
  treat_suscept = tb.infected - tb.none,
  # Difference in effect of treatment between susceptible and resistant
  # individuals
  diff_treat = (tb.infected - tb.none) - (contact.infected - contact.none),
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
saveRDS(v, file.path(data_dir, "results-limma-voom.rds"))
saveRDS(fit2, file.path(data_dir, "results-limma-fit.rds"))
saveRDS(results, file.path(data_dir, "results-limma-stats.rds"))

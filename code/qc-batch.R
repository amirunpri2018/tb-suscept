#!/usr/bin/env Rscript

# Check for technical batch effects.

suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("tidyr"))

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

# PCA --------------------------------------------------------------------------

# Calculate counts per million.
counts_cpm <- cpm(counts, log = TRUE)

# Perform PCA
pca <- prcomp(t(counts_cpm), scale. = TRUE)
variances <- pca$sdev^2
explained <- variances / sum(variances)

# Calculate the relationship between each recorded covariate and the top 6 PCs.
p_comps <- 1:6
pc_cov_cor <- matrix(nrow = ncol(anno), ncol = length(p_comps),
                     dimnames = list(colnames(anno), colnames(pca$x)[p_comps]))
for (pc in p_comps) {
  for (covariate in 1:ncol(anno)) {
    lm_result <- lm(pca$x[, pc] ~ anno[, covariate])
    r2 <- summary(lm_result)$r.squared
    pc_cov_cor[covariate, pc] <- r2
  }
}

# Convert from wide to long format
pc_cov_cor_2 <- as.data.frame(pc_cov_cor)
pc_cov_cor_2$covariate <- rownames(pc_cov_cor)
pc_cov_cor_2 <- gather(pc_cov_cor_2, key = "pc", value = "cor", -covariate)

# Save PCA results
pca_data <- cbind(anno, pca$x[, p_comps])
write.table(pca_data, file = file.path(data_dir, "results-pca.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
saveRDS(explained, file = file.path(data_dir, "results-pca-explained.rds"))
write.table(pc_cov_cor_2, file = file.path(data_dir, "results-pca-covariates.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

#!/usr/bin/env Rscript

# Identify outliers based on PCA.

suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

# Input
counts <- read.delim(file.path(data_dir, "counts.txt"),
                     check.names = FALSE, row.names = 1)
info <- read.delim(file.path(data_dir, "experiment-info.txt"), stringsAsFactors = FALSE)
rownames(info) <- info$id
info <- info %>% select(-id, -tube_ind)
stopifnot(colnames(counts) == rownames(info))

# Calculate log counts per million
counts_cpm <- cpm(counts, log = TRUE)

# Perform PCA
pca <- prcomp(t(counts_cpm), scale. = TRUE)
d <- cbind(info, pca$x)

# Function for identifying outliers
identify_outliers <- function(pc, multiplier) {
  # Returns character vector of outliers
  #
  # pc - numeric vector of PC values with sample names
  # multiplier - The number of standard deviations from the mean to be labeled
  #              an outlier.
  stopifnot(is.numeric(pc), !is.null(names(pc)),
            is.numeric(multiplier), length(multiplier) == 1)
  mu <- mean(pc)
  sd <- sd(pc)
  cutoff_lower <- mu - sd * multiplier
  cutoff_upper <- mu + sd * multiplier
  outliers <- names(pc[pc < cutoff_lower | pc > cutoff_upper])
#   if (length(outliers) == 0) {
#     outliers <- ""
#   }
  return(outliers)
}

out_list <- vector(length = 6, mode = "list")
# Because the first PC separates the samples by treatment, perform the outlier
# detection separately for each treatment for this PC.
for(treatment in unique(d$treatment)) {
  pc <- d[d$treatment == treatment, "PC1"]
  names(pc) <- rownames(d)[d$treatment == treatment]
  out_list[[1]] <- c(out_list[[1]], identify_outliers(pc, 2))
}

# Search PCs 2-6 for outliers
for (i in 2:6) {
  pc <- d[, paste0("PC", i)]
  names(pc) <- rownames(d)
  out_list[[i]] <- identify_outliers(pc, 2)
}

outliers <- unique(unlist(out_list))

# Remove outliers
outlier_index <- rownames(info) %in% outliers
info <- info[!outlier_index, ]
counts <- counts[, !outlier_index]
stopifnot(sum(outlier_index) == length(outliers), nrow(info) == ncol(counts))
# table(info$status, info$treatment)

# Save filtered data
write.table(counts, file = file.path(data_dir, "counts-filtered.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(info, file = file.path(data_dir, "experiment-info-filtered.txt"),
            quote = FALSE, sep = "\t", col.names = NA)

# Save outliers per PC
col_pc <- numeric()
col_sample <- character()
for (i in 1:length(out_list)) {
  for (o in out_list[[i]]) {
    col_pc <- c(col_pc, i)
    col_sample <- c(col_sample, o)
  }
}
out_df <- data.frame(pc = col_pc, sample = col_sample)
write.table(out_df, file = file.path(data_dir, "outliers.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

#!/usr/bin/env Rscript

# Filter lowly expressed genes

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

# Input
data_raw <- fread(file.path(data_dir, "subread-counts-per-sample.txt"),
                  data.table = FALSE)
# data_raw[1:5, 1:5]
info <- read.delim(file.path(data_dir, "experiment-info.txt"),
                   stringsAsFactors = FALSE)

# Sort
data_raw <- data_raw %>% arrange(desc(status), individual, desc(treatment))
stopifnot(data_raw$inindividual == info$individual,
          data_raw$status == info$status,
          data_raw$treatment == info$treatment)

# Transpose counts
counts_raw <- data_raw %>% select(-(individual:treatment)) %>% t()
colnames(counts_raw) <- info$id
stopifnot(ncol(counts_raw) == 50)

# Calculate log counts per million
cpm_raw <- cpm(counts_raw, log = TRUE)
# Save for creating figure
saveRDS(cpm_raw, file.path(data_dir, "cpm-all.rds"))
# Calculate the median log2 cpm per gene
cpm_raw_median <- apply(cpm_raw, 1, median)

# Filter genes
cutoff <- 0
counts <- as.data.frame(counts_raw) %>% filter(cpm_raw_median > cutoff)
rownames(counts) <- rownames(counts_raw)[cpm_raw_median > cutoff]

# Save counts
write.table(counts, file.path(data_dir, "counts.txt"), quote = FALSE, sep = "\t",
            col.names = NA)

#!/usr/bin/env Rscript

#  Calculate the mean number of reads per sample per processing stage.

suppressPackageStartupMessages(library("dplyr"))

if(interactive()) {
  data_dir <- "../data"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
}
stopifnot(dir.exists(data_dir))

# Input ------------------------------------------------------------------------

total <- read.delim(file.path(data_dir, "total-counts.txt"),
                    stringsAsFactors = FALSE)
total <- total %>%
  mutate(stage = factor(stage, levels = c("raw", "mapped to genome", "mapped to exons")),
         counts_mil <- counts / 10^6)

# Calculate total number of reads per sample -----------------------------------

total_per_sample <- total %>%
  group_by(stage, id, status, treatment) %>%
  summarize(counts = sum(counts),
            counts_mil = counts / 10^6)

# Calculate mean number of reads per stage -------------------------------------

mean_per_stage <- total_per_sample %>%
  group_by(stage) %>%
  summarize(mean_counts = mean(counts),
            sem_counts = sd(counts) / sqrt(n()),
            mean_counts_mil = mean(counts_mil),
            sem_counts_mil = sd(counts_mil) / sqrt(n()))

# Save results -----------------------------------------------------------------

write.table(mean_per_stage, file.path(data_dir, "mean-counts-per-stage.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

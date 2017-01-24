#!/usr/bin/env Rscript

# Updating the naming of the fastq files such that they match the new scheme.
#
# Currently symlinks are created in a new directory
#
# Examples:
#
# From: 01-contact-infected-H7LVWBBXX-l2.fastq.gz
# To: r01-resist-infect-H7LVWBBXX-l2.fastq.gz
#
# From: 01-tb-infected-H7LVWBBXX-l2.fastq.gz
# To: s01-suscep-infect-H7LVWBBXX-l2.fastq.gz
#
# From: 01-tb-none-H7LVWBBXX-l2.fastq.gz
# To: s01-suscep-noninf-H7LVWBBXX-l2.fastq.gz
#
# https://github.com/jdblischak/tb-suscept/commit/f0857f4833e57047ffbe5c2270876df60dd75fa3

# Run from the code subdirectory
stopifnot(grepl("code$", getwd()))

library("dplyr")
library("stringr")

dir_orig <- "/project/gilad/jdblischak/tb-suscept/fastq"
dir_new <- "/project/gilad/jdblischak/tb-suscept/fastq-symlink"
stopifnot(dir.exists(dir_orig))
if (!dir.exists(dir_new)) {
  dir.create(dir_new)
}

info <- read.delim("../data/experiment-info.txt", stringsAsFactors = FALSE)

fname_orig <- list.files(path = dir_orig,
                         pattern = "^[^undetermined].*fastq.gz$")
# 40 samples were sequenced on 3 lanes and 10 samples were sequenced on 2 lanes
num_expected <- 40 * 3 + 10 * 2
stopifnot(length(fname_orig) == num_expected)

# Parse original filenames
d_orig <- str_replace(fname_orig, pattern = ".fastq.gz$", replacement = "") %>%
  str_split(pattern = "-", simplify = TRUE) %>%
  as.data.frame()
stopifnot(is.data.frame(d_orig), !is.na(d_orig))
colnames(d_orig) <- c("individual", "status", "treatment", "flow_cell", "lane")
d_orig <- mutate(d_orig, id = paste(individual, status, treatment, sep = "-"),
                 fname = paste0(id, "-", flow_cell, "-", lane, ".fastq.gz"))
stopifnot(d_orig$fname == fname_orig)

# Convert to new filenames
d_new <- d_orig %>%
  mutate(individual = paste0(ifelse(status == "tb", "s", "r"),
                             individual),
         status = ifelse(status == "tb", "suscep", "resist"),
         treatment = ifelse(treatment == "none", "noninf", "infect"),
         id = paste(individual, status, treatment, sep = "-"),
         fname = paste0(id, "-", flow_cell, "-", lane, ".fastq.gz"))
stopifnot(is.data.frame(d_new), ncol(d_new) == ncol(d_orig),
          d_new$id %in% info$id)

# Create symlinks
for (i in 1:nrow(d_new)) {
  cat(sprintf("%d:\norig: %s\nnew: %s\n\n", i, d_orig$fname[i], d_new$fname[i]))
  fromfile <- file.path(dir_orig, d_orig$fname[i])
  stopifnot(file.exists(fromfile))
  tofile <- file.path(dir_new, d_new$fname[i])
  file.symlink(from = fromfile,
               to = tofile)
}

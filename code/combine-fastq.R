#!/usr/bin/env Rscript

# Combine the fastq files per sample to upload to GEO.
#
# Submitted with:
#
# sbatch combine-fastq.R --time=1-12:00:00 --mem=8g
#
# Requests 8g of RAM for a potential max of 36 hours (the max for sandby, see
# `rcchelp qos`)

# Run from the code subdirectory
stopifnot(grepl("code$", getwd()))

library("dplyr")
library("stringr")

# Directory which contains the fastq files per lane
dir_per_lane <- "/project/gilad/jdblischak/tb-suscept/fastq-symlink"
# Directory which contains the fastq files per sample
dir_per_sample <- "/scratch/midway/jdblischak/tb-suscept-geo"
stopifnot(dir.exists(dir_per_lane))
if (!dir.exists(dir_per_sample)) {
  dir.create(dir_per_sample)
}

info <- read.delim("../data/experiment-info.txt", stringsAsFactors = FALSE)

for (i in 1:nrow(info)) {
  id <- info$id[i]
  print(id)
  files_per_lane <- list.files(path = dir_per_lane,
                               pattern = paste0("^", id))
  stopifnot(length(files_per_lane) %in% 2:3)
  target_file <- file.path(dir_per_sample, paste0(id, ".fastq.gz"))
  if (file.exists(target_file)) {
    warning(sprintf("File %s already exists. Skipping.", target_file))
    next()
  }
  cmd <- sprintf("zcat %s* | gzip > %s",
                 file.path(dir_per_lane, id), target_file)
  print(cmd)
  system(cmd)
}

files_per_sample <- list.files(path = dir_per_sample,
                             pattern = ".fastq.gz$",
                             full.names = TRUE)

if (length(files_per_sample) != nrow(info)) {
  warning("Missing samples")
}

for (f in files_per_sample) {
  size <- file.size(f)
  if (is.na(size)) {
    warning("Something wrong with file: ", f)
  }
  if (file.size(f) == 0) {
    warning("File is empty: ", f)
  }
}

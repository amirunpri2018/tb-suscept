#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("edgeR"))

if(interactive()) {
  data_dir <- "../data"
  fig_dir <- "../figure"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
  fig_dir <- args[2]
}
stopifnot(dir.exists(data_dir))
stopifnot(dir.exists(fig_dir))

# Figure settings --------------------------------------------------------------

# Dimensions in inches
w <- 7
h <- 7

theme_set(theme_bw(base_size = 12))
theme_update(panel.grid.minor.x = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.grid.major.x = element_blank(),
             panel.grid.major.y = element_blank())


# Filter genes -----------------------------------------------------------------
# See code/qc-genes.R

plot_densities <- function(x, include_median = TRUE, cutoff = NULL, ...) {
  # Plot the density of each sample (column)
  #
  # x - matrix of log2 cpm values
  #
  x_density <- apply(x, 2, density)
  plot(x_density[[1]], lty = "dotted", col = gray(0.05, alpha = 0.5),
       xlab = "log2 cpm", ...)
  suppress_output <- lapply(x_density, lines, lty = "dotted",
                            col = gray(0.05, alpha = 0.5))
  subtitle <- ""
  if (include_median) {
    x_median <- apply(x, 1, median)
    x_median_density <- density(x_median)
    lines(x_median_density, col = "red", lwd = 3)
    subtitle <- paste(subtitle, "red = median")
  }
  cutoff <- 0
  if (!is.null(cutoff)) {
    abline(v = cutoff, col = "blue")
    subtitle <- paste(subtitle, "blue = cutoff")
  }
  title(sub = subtitle)
}

# log2 cpm - all genes, all samples
cpm_raw <- readRDS(file.path(data_dir, "cpm-all.rds"))
pdf(file.path(fig_dir, "density-all-genes-all-samples.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_densities(cpm_raw, include_median = TRUE, cutoff = 0,
               main = "All genes, all samples",
               xlim = c(-11, 16), ylim = c(0, 0.25))
dev.off()

# log2 cpm - filtered genes, all samples
counts <- read.delim(file.path(data_dir, "counts.txt"),
                     check.names = FALSE, row.names = 1)
counts_cpm <- cpm(counts, log = TRUE)

pdf(file.path(fig_dir, "density-filt-genes-all-samples.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_densities(counts_cpm, include_median = TRUE, cutoff = 0,
               main = "Filtered genes, all samples",
               xlim = c(-11, 16), ylim = c(0, 0.25))
dev.off()

# log2 cpm - filtered genes, filtered samples
counts_filt <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                          check.names = FALSE, row.names = 1)
counts_filt_cpm <- cpm(counts_filt, log = TRUE)
pdf(file.path(fig_dir, "density-filt-genes-filt-samples.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_densities(counts_filt_cpm, include_median = TRUE, cutoff = 0,
               main = "Filtered genes, filtered samples",
               xlim = c(-11, 16), ylim = c(0, 0.25))
dev.off()

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("VennDiagram"))

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
invisible(dev.off())

# log2 cpm - filtered genes, all samples
counts <- read.delim(file.path(data_dir, "counts.txt"),
                     check.names = FALSE, row.names = 1)
counts_cpm <- cpm(counts, log = TRUE)

pdf(file.path(fig_dir, "density-filt-genes-all-samples.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_densities(counts_cpm, include_median = TRUE, cutoff = 0,
               main = "Filtered genes, all samples",
               xlim = c(-11, 16), ylim = c(0, 0.25))
invisible(dev.off())

# log2 cpm - filtered genes, filtered samples
counts_filt <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                          check.names = FALSE, row.names = 1)
counts_filt_cpm <- cpm(counts_filt, log = TRUE)
pdf(file.path(fig_dir, "density-filt-genes-filt-samples.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_densities(counts_filt_cpm, include_median = TRUE, cutoff = 0,
               main = "Filtered genes, filtered samples",
               xlim = c(-11, 16), ylim = c(0, 0.25))
invisible(dev.off())

# Differential expression with limma -------------------------------------------
# see code/main-limma.R

plot_ma <- function(x, qval) {
  # Create MA plot.
  #
  # x - data frame with topTable and ASH output
  #     (columns logFC, AveExpr, and qvalue)
  #
  stopifnot(is.data.frame(x), c("logFC", "AveExpr", "qvalue") %in% colnames(x),
            is.numeric(qval), qval <= 1, qval >= 0)
  x$highlight <- ifelse(x$qvalue < qval, "darkred", "gray75")
  ggplot(x, aes(x = AveExpr, y = logFC, color = highlight)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = "Average expression level", y = "Log fold change") +
  scale_color_identity()
#   scale_color_gradient(low = "red", high = "white", limits = c(0, 0.25))
}

plot_pval_hist <- function(x) {
  # Create histogram of p-values.
  #
  # x - data frame with topTable output (column P.Value)
  #
  stopifnot(is.data.frame(x), c("P.Value") %in% colnames(x))
  ggplot(x, aes(x = P.Value)) +
    geom_histogram(binwidth = 0.01) +
    labs(x = "p-value", y = "Number of genes")
}

results <- readRDS(file.path(data_dir, "results-limma-stats.rds"))

# Difference between susceptible and resistant individuals before treatment
pdf(file.path(fig_dir, "ma-diff-before.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_ma(results[["diff_before"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = "Difference between susceptible and resistant individuals before treatment")
invisible(dev.off())
pdf(file.path(fig_dir, "pval-diff-before.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["diff_before"]]) +
  labs(title = "Difference between susceptible and resistant individuals before treatment")
invisible(dev.off())

# Difference between susceptible and resistant individuals after treatment
pdf(file.path(fig_dir, "ma-diff-after.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_ma(results[["diff_after"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = "Difference between susceptible and resistant individuals after treatment")
invisible(dev.off())
pdf(file.path(fig_dir, "pva-diff-after.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["diff_after"]]) +
  labs(title = "Difference between susceptible and resistant individuals after treatment")
invisible(dev.off())

# Effect of treatment in resistant individuals
pdf(file.path(fig_dir, "ma-treat-resist.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_ma(results[["treat_resist"]], qval = 0.1) +
  labs(title = "Effect of treatment in resistant individuals")
invisible(dev.off())
pdf(file.path(fig_dir, "pval-treat-resist.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["treat_resist"]]) +
  labs(title = "Effect of treatment in resistant individuals")
invisible(dev.off())

# Effect of treatment in susceptible individuals
pdf(file.path(fig_dir, "ma-treat-suscept.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_ma(results[["treat_suscept"]], qval = 0.1) +
  labs(title = "Effect of treatment in susceptible individuals")
invisible(dev.off())
pdf(file.path(fig_dir, "pval-treat-suscept.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["treat_suscept"]]) +
  labs(title = "Effect of treatment in susceptible individuals")
invisible(dev.off())

# Difference in effect of treatment between susceptible and resistant
# individuals
pdf(file.path(fig_dir, "ma-diff-treat.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_ma(results[["diff_treat"]], qval = 0.1) +
  labs(title = "Difference in effect of treatment between susceptible and resistant individuals")
invisible(dev.off())
pdf(file.path(fig_dir, "pval-diff-treat.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["diff_treat"]]) +
  labs(title = "Difference in effect of treatment between susceptible and resistant individuals")
invisible(dev.off())

# Venn diagram of DE gene overlap
# Not saving because they aren't very informative
ash_diff_before <- results[["diff_before"]]$svalue < 0.1
ash_diff_after <- results[["diff_after"]]$svalue < 0.1
ash_treat_resist <- results[["treat_resist"]]$svalue < 0.1
ash_treat_suscept <- results[["treat_suscept"]]$svalue < 0.1
ash_diff_treat <- results[["diff_treat"]]$svalue < 0.1
# grid.newpage()
venn_status <- draw.triple.venn(
  area1 = sum(ash_diff_before),
  area2 = sum(ash_diff_after),
  area3 = sum(ash_diff_treat),
  n12 = sum(ash_diff_before & ash_diff_after),
  n23 = sum(ash_diff_after & ash_diff_treat),
  n13 = sum(ash_diff_before & ash_diff_treat),
  n123 = sum(ash_diff_before & ash_diff_after & ash_diff_treat),
  category = c("before", "after", "response"),
  #   fill = c("darkolivegreen3", "cadetblue3", "darkorchid"),
  alpha = c(0.5, 0.5, 0.5), col = "black", cex = 2, cat.cex = 1.5,
  euler.d = FALSE, scaled = FALSE, ind = FALSE)
# grid.draw(venn_status)

# grid.newpage()
venn_treat <- draw.triple.venn(
  area1 = sum(ash_diff_before),
  area2 = sum(ash_treat_resist),
  area3 = sum(ash_treat_suscept),
  n12 = sum(ash_diff_before & ash_treat_resist),
  n23 = sum(ash_treat_resist & ash_treat_suscept),
  n13 = sum(ash_diff_before & ash_treat_suscept),
  n123 = sum(ash_diff_before & ash_treat_resist & ash_treat_suscept),
  category = c("before", "Treatment in resist", "Treatment in suscept"),
  #   fill = c("darkolivegreen3", "cadetblue3", "darkorchid"),
  alpha = c(0.5, 0.5, 0.5), col = "black", cex = 2, cat.cex = 1.5,
  euler.d = FALSE, scaled = FALSE, ind = FALSE)
# grid.draw(venn_treat)

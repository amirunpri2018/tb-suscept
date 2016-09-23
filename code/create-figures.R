#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("flux")) # for auc
suppressPackageStartupMessages(library("cowplot"))

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

# theme_set(theme_classic(base_size = 12))
# theme_update(panel.grid.minor.x = element_blank(),
#              panel.grid.minor.y = element_blank(),
#              panel.grid.major.x = element_blank(),
#              panel.grid.major.y = element_blank())

# Color scheme:
color_scheme <- brewer.pal(n = 5, name = "Set1")
# non-infected is purple
col_noninf <- color_scheme[4]
# infected is green
col_infect <- color_scheme[3]
# susceptible is blue
col_suscep <- color_scheme[2]
# resistant is orange
col_resist <- color_scheme[5]

my_ggsave <- function(fname, dims = c(1, 1)) {
  # Save the last plot using these defaults.
  #
  # fname - character string of base filename (no path)
  # dims - multiplication factors of default width and height.
  #        numeric vector of length 2: (width, height)
  ggsave(fname, path = fig_dir,
         width = w * dims[1], height = h * dims[2], units = "in",
         dpi = 300)
}

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
# log2 cpm - filtered genes, all samples
counts <- read.delim(file.path(data_dir, "counts.txt"),
                     check.names = FALSE, row.names = 1)
counts_cpm <- cpm(counts, log = TRUE)
# log2 cpm - filtered genes, filtered samples
counts_filt <- read.delim(file.path(data_dir, "counts-filtered.txt"),
                          check.names = FALSE, row.names = 1)
counts_filt_cpm <- cpm(counts_filt, log = TRUE)

for (device in c("png", "pdf")) {
  if (device == "png") {
    png(file.path(fig_dir, "gene-exp-distribution.png"),
        width = 3 * w, height = h, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "gene-exp-distribution.pdf"),
        width = 3 * w, height = h, useDingbats = FALSE)
  }
  par(mfrow = c(1, 3), cex = 1.5)
  plot_densities(cpm_raw, include_median = TRUE, cutoff = 0,
                 main = "All genes, all samples",
                 xlim = c(-11, 16), ylim = c(0, 0.25))
  plot_densities(counts_cpm, include_median = TRUE, cutoff = 0,
                 main = "Filtered genes, all samples",
                 xlim = c(-11, 16), ylim = c(0, 0.25))
  plot_densities(counts_filt_cpm, include_median = TRUE, cutoff = 0,
                 main = "Filtered genes, filtered samples",
                 xlim = c(-11, 16), ylim = c(0, 0.25))
  invisible(dev.off())
}

# Remove outliers ----------------------------------------------------------------
# see code/qc-outliers.R

pca_outliers <- read.delim(file.path(data_dir, "pca-outliers.txt"),
                           stringsAsFactors = FALSE, row.names = 1)
explained_outliers <- readRDS(file.path(data_dir, "explained-outliers.rds"))
outliers <- read.delim(file.path(data_dir, "outliers.txt"),
                       stringsAsFactors = FALSE)
outliers_simple <- outliers %>%
  group_by(sample) %>%
  summarize(PC = paste(sort(pc), collapse = ", ")) %>%
  arrange(PC) %>%
  rename(Sample = sample)

# PCA

# PC1 versus PC2.
pc1v2_out <- ggplot(pca_outliers, aes(x = PC1, y = PC2, color = treatment)) +
  geom_text(aes(label = individual)) +
  scale_color_manual(values = c(col_infect, col_noninf)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained_outliers[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained_outliers[2] * 100, 2)))
# PC3 versus PC4.
pc3v4_out <- ggplot(pca_outliers, aes(x = PC3, y = PC4, color = treatment)) +
  geom_text(aes(label = individual)) +
  scale_color_manual(values = c(col_infect, col_noninf)) +
  labs(x = sprintf("PC%d (%.2f%%)", 3, round(explained_outliers[3] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 4, round(explained_outliers[4] * 100, 2)))
# PC5 versus PC6.
pc5v6_out <- ggplot(pca_outliers, aes(x = PC5, y = PC6, color = treatment)) +
  geom_text(aes(label = individual)) +
  scale_color_manual(values = c(col_infect, col_noninf)) +
  labs(x = sprintf("PC%d (%.2f%%)", 5, round(explained_outliers[5] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 6, round(explained_outliers[6] * 100, 2)))
# Table of outliers
outliers_grob <- tableGrob(outliers_simple,
                           rows = rep("", nrow(outliers_simple)),
                           theme = ttheme_default(base_size = 16))

multi_outliers <- plot_grid(pc1v2_out + theme(legend.position = "none"),
                            pc3v4_out + theme(legend.position = "none"),
                            pc5v6_out + theme(legend.position = "none"),
                            outliers_grob,
                            labels = letters[1:4])

my_ggsave("outliers.pdf", dims = c(2, 2))
my_ggsave("outliers.png", dims = c(2, 2))

# Hierarchical clustering

# Data already loaded above:
# counts_cpm - Filtered genes, all samples
# counts_filt_cpm - Filtered genes, filtered samples
stopifnot(ncol(counts_cpm) == 50,
          ncol(counts_filt_cpm) == 50 - nrow(outliers_simple))

for (device in c("png", "pdf")) {
  if (device == "png") {
    png(file.path(fig_dir, "heatmap-all-samples.png"),
        width = w, height = h, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "heatmap-all-samples.pdf"),
        width = w, height = h, useDingbats = FALSE)
  }
  heatmap.2(cor(counts_cpm), trace = "none", margins = c(7, 7),
            main = "All samples")
  invisible(dev.off())
}

for (device in c("png", "pdf")) {
  if (device == "png") {
    png(file.path(fig_dir, "heatmap-no-outliers.png"),
        width = w, height = h, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "heatmap-no-outliers.pdf"),
        width = w, height = h, useDingbats = FALSE)
  }
  heatmap.2(cor(counts_filt_cpm), trace = "none", margins = c(7, 7),
            main = "After removing outliers")
  invisible(dev.off())
}

# Batch effects ----------------------------------------------------------------
# see code/qc-batch.R

pca_data <- read.delim(file.path(data_dir, "results-pca.txt"),
                       stringsAsFactors = FALSE, row.names = 1)
pca_data$status <- factor(pca_data$status, levels = c("resist", "suscep"),
                          labels = c("resistant", "susceptible"))
pca_data$treatment <- factor(pca_data$treatment, levels = c("noninf", "infect"),
                             labels = c("non-infected", "infected"))
explained <- readRDS(file.path(data_dir, "results-pca-explained.rds"))

# PC1 versus PC2.
pc1v2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_text(aes(label = individual)) +
  scale_color_manual(values = c(col_infect, col_noninf)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained[2] * 100, 2)))

# PC3 versus PC4.
pc3v4 <- ggplot(pca_data, aes(x = PC3, y = PC4, color = infection)) +
  geom_text(aes(label = individual)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  labs(x = sprintf("PC%d (%.2f%%)", 3, round(explained[3] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 4, round(explained[4] * 100, 2)))

# PC5 versus PC6.
pc5v6 <- ggplot(pca_data, aes(x = PC5, y = PC6, color = infection)) +
  geom_text(aes(label = individual)) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  labs(x = sprintf("PC%d (%.2f%%)", 5, round(explained[5] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 6, round(explained[6] * 100, 2)))

d_heatmap <- read.delim(file.path(data_dir, "results-pca-covariates.txt"),
                        stringsAsFactors = FALSE)
d_heatmap$covariate <- factor(d_heatmap$covariate,
                             levels = c("rin", "master_mix",
                                        "extraction", "arrival", "infection",
                                        "individual", "status", "treatment"),
                             labels = c("RNA quality", "Library prep batch",
                                        "RNA extraction batch", "Arrival batch",
                                        "Infection batch", "Individual",
                                        "Susceptibility status", "Treatment"))
pca_heat <- ggplot(d_heatmap, aes(x = pc, y = covariate)) +
  geom_tile(aes(fill = cor), colour = "white") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), name = "r") +
  labs(x = "Principal Component", y = "",
       title = "Correlations with experimental covariates")

p_batch_1 <- plot_grid(pc1v2 + theme(legend.position = "none"),
                       pc3v4 + theme(legend.position = "none"),
                       pc5v6 + theme(legend.position = "none"),
                       pca_heat,
                       labels = letters[1:4])

my_ggsave("batch-pca.pdf", dims = c(2, 2))
my_ggsave("batch-pca.png", dims = c(2, 2))

# Further investigate correlation of infection batch with PC3 and PC5

pc3_infection <- ggplot(pca_data, aes(x = infection, y = PC3, fill = treatment)) +
  geom_boxplot() +
  scale_fill_manual(name = "Treatment", values = c(col_noninf, col_infect)) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "PC3 varies by date of infection",
       x = "Date of infection experiment")

pc3_status <- ggplot(pca_data, aes(x = reorder(individual, PC3), y = PC3, color = status)) +
  geom_point() +
  scale_color_manual(name = "Status", values = c(col_resist, col_suscep)) +
  facet_wrap(~treatment) +
  labs(title = "PC3 does not correlate with susceptibility status",
       x = "Individual") +
  theme(axis.text.x = element_text(angle = 90, size = rel(0.5)))

pc5_infection <- pc3_infection %+% aes(y = PC5) +
  labs(title = "PC5 varies by date of infection")

pc5_status <- pc3_status %+% aes(x = reorder(individual, PC5), y = PC5) +
  labs(title = "PC5 does not correlate with susceptibility status",
       x = "Individual")

p_batch_2 <- plot_grid(pc3_infection + theme(legend.position = "bottom"),
                       pc5_infection  + theme(legend.position = "bottom"),
                       pc3_status + theme(legend.position = "bottom"),
                       pc5_status + theme(legend.position = "bottom"),
                       labels = letters[1:4])

my_ggsave("batch-infection.pdf", dims = c(2, 2))
my_ggsave("batch-infection.png", dims = c(2, 2))

# Differential expression with limma -------------------------------------------
# see code/main-limma.R

plot_ma <- function(x, qval) {
  # Create MA plot.
  #
  # x - data frame with topTable and ASH output
  #     (columns logFC, AveExpr, and qvalue)
  # qval - qvalue cutoff for calling a gene DE
  #
  stopifnot(is.data.frame(x), c("logFC", "AveExpr", "qvalue") %in% colnames(x),
            is.numeric(qval), qval <= 1, qval >= 0)
  x$highlight <- ifelse(x$qvalue < qval, "darkred", "gray75")
  x$highlight <- factor(x$highlight, levels = c("darkred", "gray75"))
  ggplot(x, aes(x = AveExpr, y = logFC, color = highlight, shape = highlight)) +
    geom_point() +
    labs(x = "Average expression level", y = "Log fold change") +
    scale_color_identity(drop = FALSE) +
    scale_shape_manual(values = c(16, 1), drop = FALSE) +
    theme(legend.position = "none")
#   scale_color_gradient(low = "red", high = "white", limits = c(0, 0.25))
}

plot_volcano <- function(x, qval) {
  # Create volcano plot.
  #
  # x - data frame with topTable and ASH output
  #     (columns logFC, P.Value, and qvalue)
  # qval - qvalue cutoff for calling a gene DE
  #
  stopifnot(is.data.frame(x), c("logFC", "P.Value", "qvalue") %in% colnames(x),
            is.numeric(qval), qval <= 1, qval >= 0)
  x$highlight <- ifelse(x$qvalue < qval, "darkred", "gray75")
  x$highlight <- factor(x$highlight, levels = c("darkred", "gray75"))
  ggplot(x, aes(x = logFC, y = -log10(P.Value), color = highlight)) +
    geom_point(shape = 1) +
    labs(x = "Log fold change",
         y = expression(-log[10] * " p-value")) +
    scale_color_identity(drop = FALSE) +
    theme(legend.position = "none")
}

plot_pval_hist <- function(x, qval) {
  # Create histogram of p-values.
  #
  # x - data frame with topTable and ash output (columns P.Value and qvalue)
  # qval - qvalue cutoff for calling a gene DE
  #
  stopifnot(is.data.frame(x), c("P.Value", "qvalue") %in% colnames(x))
  x$highlight <- ifelse(x$qvalue < qval, "darkred", "gray75")
  x$highlight <- factor(x$highlight, levels = c("darkred", "gray75"))
  ggplot(x, aes(x = P.Value, fill = highlight)) +
    geom_histogram(position = "stack", binwidth = 0.01) +
    scale_fill_identity(drop = FALSE) +
    labs(x = "p-value", y = "Number of genes")
}

results <- readRDS(file.path(data_dir, "results-limma-stats.rds"))

# Difference between susceptible and resistant individuals in the non-infected state
title_status_ni <- "Difference between susceptible and resistant\nindividuals in the non-infected state"
ma_status_ni <- plot_ma(results[["status_ni"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = title_status_ni)
volc_status_ni <- plot_volcano(results[["status_ni"]], qval = 0.1) +
  xlim(-3.5, 3.5) + ylim(0, 6.5) +
  labs(title = title_status_ni)
pval_status_ni <- plot_pval_hist(results[["status_ni"]], qval = 0.1) +
  ylim(0, 325) +
  labs(title = title_status_ni)

# Difference between susceptible and resistant individuals in the infected state
title_status_ii <- "Difference between susceptible and resistant\nindividuals in the infected state"
ma_status_ii <- plot_ma(results[["status_ii"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = title_status_ii)
volc_status_ii <- plot_volcano(results[["status_ii"]], qval = 0.1) +
  xlim(-3.5, 3.5) + ylim(0, 6.5) +
  labs(title = title_status_ii)
pval_status_ii <- plot_pval_hist(results[["status_ii"]], qval = 0.1) +
  ylim(0, 325) +
  labs(title = title_status_ii)

# Effect of treatment in resistant individuals
title_treat_resist <- "Effect of treatment in resistant individuals"
ma_treat_resist <- plot_ma(results[["treat_resist"]], qval = 0.1) +
  ylim(-8.5, 8.5) +
  labs(title = title_treat_resist)
volc_treat_resist <- plot_volcano(results[["treat_resist"]], qval = 0.1) +
  xlim(-8.5, 8.5) + ylim(0, 45) +
  labs(title = title_treat_resist)
pval_treat_resist <- plot_pval_hist(results[["treat_resist"]], qval = 0.1) +
  ylim(0, 9000) +
  labs(title = title_treat_resist)

# Effect of treatment in susceptible individuals
title_treat_suscep <- "Effect of treatment in susceptible individuals"
ma_treat_suscep <- plot_ma(results[["treat_suscep"]], qval = 0.1) +
  ylim(-8.5, 8.5) +
  labs(title = title_treat_suscep)
volc_treat_suscep <- plot_volcano(results[["treat_suscep"]], qval = 0.1) +
  xlim(-8.5, 8.5) + ylim(0, 45) +
  labs(title = title_treat_suscep)
pval_treat_suscep <- plot_pval_hist(results[["treat_suscep"]], qval = 0.1) +
  ylim(0, 9000) +
  labs(title = title_treat_suscep)

# Difference in effect of treatment between susceptible and resistant
# individuals
title_interact <- "Difference in effect of treatment between\nsusceptible and resistant individuals"
ma_interact <- plot_ma(results[["interact"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = title_interact)
volc_interact <- plot_volcano(results[["interact"]], qval = 0.1) +
  xlim(-3.5, 3.5) + ylim(0, 6.5) +
  labs(title = title_interact)
pval_interact <- plot_pval_hist(results[["interact"]], qval = 0.1) +
  ylim(0, 325) +
  labs(title = title_interact)

# Multipanel main figure
limma_main <- plot_grid(pval_status_ni,
                     pval_status_ii,
                     volc_status_ni + labs(title = ""),
                     volc_status_ii + labs(title = ""),
                     labels = letters[1:4])
my_ggsave("limma.eps", dims = c(2, 2))
my_ggsave("limma.pdf", dims = c(2, 2))
my_ggsave("limma.png", dims = c(2, 2))

# Multipanel supplemental figure
limma_supp <- plot_grid(pval_treat_resist,
                        pval_treat_suscep,
                        volc_treat_resist + labs(title = ""),
                        volc_treat_suscep + labs(title = ""),
                        labels = letters[1:4])
my_ggsave("limma-supp.pdf", dims = c(2, 2))
my_ggsave("limma-supp.png", dims = c(2, 2))

# Venn diagram of DE gene overlap
# Not saving because they aren't very informative
ash_status_ni <- results[["status_ni"]]$qvalue < 0.1
ash_status_ii <- results[["status_ii"]]$qvalue < 0.1
ash_treat_resist <- results[["treat_resist"]]$qvalue < 0.1
ash_treat_suscep <- results[["treat_suscep"]]$qvalue < 0.1
ash_interact <- results[["interact"]]$qvalue < 0.1
# grid.newpage()
venn_status <- draw.triple.venn(
  area1 = sum(ash_status_ni),
  area2 = sum(ash_status_ii),
  area3 = sum(ash_interact),
  n12 = sum(ash_status_ni & ash_status_ii),
  n23 = sum(ash_status_ii & ash_interact),
  n13 = sum(ash_status_ni & ash_interact),
  n123 = sum(ash_status_ni & ash_status_ii & ash_interact),
  category = c("Difference in non-infected state",
               "Difference in infected state",
               "Difference in response to infection"),
  #   fill = c("darkolivegreen3", "cadetblue3", "darkorchid"),
  alpha = c(0.5, 0.5, 0.5), col = "black", cex = 2, cat.cex = 1.5,
  euler.d = FALSE, scaled = FALSE, ind = FALSE)
# grid.draw(venn_status)

# grid.newpage()
venn_treat <- draw.triple.venn(
  area1 = sum(ash_status_ni),
  area2 = sum(ash_treat_resist),
  area3 = sum(ash_treat_suscep),
  n12 = sum(ash_status_ni & ash_treat_resist),
  n23 = sum(ash_treat_resist & ash_treat_suscep),
  n13 = sum(ash_status_ni & ash_treat_suscep),
  n123 = sum(ash_status_ni & ash_treat_resist & ash_treat_suscep),
  category = c("Difference in non-infected state",
               "Effect of treatment in resistant individuals",
               "Effect of treatment in susceptible individuals"),
  #   fill = c("darkolivegreen3", "cadetblue3", "darkorchid"),
  alpha = c(0.5, 0.5, 0.5), col = "black", cex = 2, cat.cex = 1.5,
  euler.d = FALSE, scaled = FALSE, ind = FALSE)
# grid.draw(venn_treat)

# GWAS enrichment --------------------------------------------------------------
# see code/main-gwas.R

# Main figure (gwas) will be:
# (a) enrichment plot for Gambia GWAS and status_ni
# (b) boxplot of permuted auc's to assess significance of enrichment for Gambia
#     GWAS and all DE tests
#
# Supp figure (gwas-supp) will be:
# (a) enrichment plot for Ghana GWAS and status_ni
# (b) boxplot of permuted auc's to assess significance of enrichment for Ghana
#     GWAS and all DE tests

# Input detailed data for enrichment of status_ni and Gambia GWAS
enrich_gambia_status_ni <- read.delim(file.path(data_dir,
                                                "gambia-status_ni-enrichment.txt"))
intervals_gambia_status_ni <- read.delim(file.path(data_dir,
                                                   "gambia-status_ni-intervals.txt"))
sizes_gambia_status_ni <- read.delim(file.path(data_dir,
                                               "gambia-status_ni-sizes.txt"))

# Input detailed data for enrichment of status_ni and Gambia GWAS
enrich_ghana_status_ni <- read.delim(file.path(data_dir,
                                               "ghana-status_ni-enrichment.txt"))
intervals_ghana_status_ni <- read.delim(file.path(data_dir,
                                                  "ghana-status_ni-intervals.txt"))
sizes_ghana_status_ni <- read.delim(file.path(data_dir,
                                              "ghana-status_ni-sizes.txt"))

# Input enrichment data only for all the comparisons
gwas_list <- list()
gwas_signif <- numeric()
for (gwas in c("gambia", "ghana")) {
  for (test in c("status_ni", "status_ii", "treat_resist", "treat_suscep")) {
    fname_base <- file.path(data_dir, paste(gwas, test, sep = "-"))
    print(fname_base)
    if (gwas == "gambia") {
      enrich_result <- read.delim(file.path(data_dir,
                                            paste(gwas, test, "enrichment.txt",
                                                  sep = "-")))
    } else if (gwas == "ghana") {
      enrich_result <- read.delim(file.path(data_dir,
                                            paste(gwas, test, "enrichment.txt",
                                                  sep = "-")))
    }
    # http://artax.karlin.mff.cuni.cz/r-help/library/flux/html/auc.html
    record_auc <- apply(enrich_result, 2, function(r) auc(x = 1:length(r), y = r))
    gwas_list[[gwas]][[test]] <- record_auc
    gwas_signif[paste(gwas, test, sep = "-")] <-
      sum(record_auc[-1] > record_auc[1]) / length(record_auc[-1])
  }
}
gwas_list[["gambia"]] <- data.frame(status_ni = gwas_list[["gambia"]][["status_ni"]],
                                    status_ii = gwas_list[["gambia"]][["status_ii"]],
                                    treat_resist = gwas_list[["gambia"]][["treat_resist"]],
                                    treat_suscep = gwas_list[["gambia"]][["treat_suscep"]])
gwas_list[["ghana"]] <- data.frame(status_ni = gwas_list[["ghana"]][["status_ni"]],
                                   status_ii = gwas_list[["ghana"]][["status_ii"]],
                                   treat_resist = gwas_list[["ghana"]][["treat_resist"]],
                                   treat_suscep = gwas_list[["ghana"]][["treat_suscep"]])

# Function to create enrichment plot
#
# Input is output from enrich_full.
#
# enrichment - fold enrichment of a particular metric at increasing stringency
#              of another metric
# intervals - the cutoff value for the metric
# sizes - the number of genes considered at each interval
# ... - arguments passed to `plot`
enrich_plot <- function(enrichment, intervals, sizes, ...) {
  plot(enrichment[, 1], type = "l", col = "red",
       ylim = c(min(enrichment, na.rm = TRUE),
                max(enrichment, na.rm = TRUE)),
       ylab = "Fold enrichment of GWAS P < 0.05",
       xlab = "Number of genes (effect size cutoff)",
       xaxt = "n", ...)
  n <- nrow(enrichment)
  spacing <- seq(1, n, by = 50)
  axis(side = 1, at = (1:n)[spacing], padj = 0.5,
       labels = paste0(sizes[, 1], "\n(", round(intervals[, 1], digits = 2),
                       ")")[spacing])
  apply(enrichment[, -1], 2, lines, col = "grey75")
  lines(enrichment[, 1], col = "red")
  abline(h = 1, col = "blue", lty = 2)
}

# Function to create boxplot of fold enrichments
#
# x - data frame with 4 columns of enrichment values for each of the 4 DE tests
# ... - arguments passed to boxplot
enrich_boxplot <- function(x, ...) {
  stopifnot(colnames(x) == c("status_ni", "status_ii", "treat_resist", "treat_suscep"))
  boxplot(x[-1, ], ylab = "Area under fold enrichment curve", xaxt = "n",
          ylim = c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)), ...)
  axis(side = 1, at = 1:4, tck = 0, labels = c("suscep v resist\n(noninf)",
                                               "suscep v resist\n(infect)",
                                               "infect v. noninf\n(resist)",
                                               "infect v. noninf\n(suscep)"))
  points(1:4, x[1, ], col = "red", pch = 19)
}

# Plot titles shared across main and supp figures
enrich_title_raw <- "Enrichment of GWAS signal in genes differentially expressed between susceptible and resistant individuals in the non-infected state"
enrich_title <- paste(strwrap(enrich_title_raw, width = 60), collapse = "\n")
enrich_title_box_raw <- "Significant enrichment of GWAS signal in each differential expression test"
enrich_title_box <- paste(strwrap(enrich_title_box_raw, width = 60),
                          collapse = "\n")

# Main figure using Gambia data
for (device in c("png", "pdf", "eps")) {
  if (device == "png") {
    png(file.path(fig_dir, "gwas.png"),
        width = 2 * w, height = h, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "gwas.pdf"),
        width = 2 * w, height = h, useDingbats = FALSE)
  } else if (device == "eps") {
    postscript(file.path(fig_dir, "gwas.eps"),
               width = 2 * w, height = h)
    par(cex.axis = 0.7)
  }
  par(mfrow = c(1, 2), cex = 1)
  # Make enrichment plot for Gambia GWAS and status_ni
  enrich_plot(enrich_gambia_status_ni,
              intervals_gambia_status_ni,
              sizes_gambia_status_ni,
              main = enrich_title)
  mtext("a", side = 3, line = 1, adj = 0, outer = FALSE, font = 2)
  enrich_boxplot(gwas_list[["gambia"]],
                 main = enrich_title_box)
  mtext("b", side = 3, line = 1, adj = 0, outer = FALSE, font = 2)
  invisible(dev.off())
}

# Supp figure using Ghana data
for (device in c("png", "pdf")) {
  if (device == "png") {
    png(file.path(fig_dir, "gwas-supp.png"),
        width = 2 * w, height = h, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "gwas-supp.pdf"),
        width = 2 * w, height = h, useDingbats = FALSE)
  }
  par(mfrow = c(1, 2), cex = 1)
  # Make enrichment plot for Ghana GWAS and status_ni
  enrich_plot(enrich_ghana_status_ni,
              intervals_ghana_status_ni,
              sizes_ghana_status_ni,
              main = enrich_title)
  mtext("a", side = 3, line = 1, adj = 0, outer = FALSE, font = 2)
  enrich_boxplot(gwas_list[["ghana"]],
                 main = enrich_title_box)
  mtext("b", side = 3, line = 1, adj = 0, outer = FALSE, font = 2)
  invisible(dev.off())
}

# Combine studies --------------------------------------------------------------
# see code/combine-studies.R

# Combined data
exp_combined <- read.delim(file.path(data_dir, "combined-raw.txt"),
                           row.names = 1, check.names = FALSE)
anno_comb <- read.delim(file.path(data_dir, "combined-annotation.txt"),
                        stringsAsFactors = FALSE, row.names = 1)

# Normalized data
exp_combined_norm <- read.delim(file.path(data_dir, "combined-normalized.txt"),
                                row.names = 1, check.names = FALSE)
combined_pca <- read.delim(file.path(data_dir, "combined-pca.txt"),
                           row.names = 1, stringsAsFactors = FALSE)
combined_explained <- readRDS(file.path(data_dir, "combined-pca-explained.txt"))

# Regressed data
regressed <- read.delim(file.path(data_dir, "combined-regressed.txt"),
                        row.names = 1, check.names = FALSE)
combined_pca_regr <- read.delim(file.path(data_dir, "combined-pca-regressed.txt"),
                                row.names = 1, check.names = FALSE)
combined_explained_regr <- readRDS(file.path(data_dir, "combined-pca-explained-regressed.txt"))

# Plots of distributions
dist_current_pre <- apply(exp_combined[, anno_comb$study == "current"], 1, median)
dist_lbb2012_pre <- apply(exp_combined[, anno_comb$study == "lbb2012"], 1, median)
dist_current_norm <- apply(exp_combined_norm[, anno_comb$study == "current"],
                           1, median)
dist_lbb2012_norm <- apply(exp_combined_norm[, anno_comb$study == "lbb2012"],
                           1, median)
dist_data <- gather(data.frame(dist_current_pre, dist_lbb2012_pre,
                        dist_current_norm, dist_lbb2012_norm),
                    key = "dist", value = "gene_exp")
dist_data <- separate(dist_data, col = dist,
                      into = c("delete", "study", "normalized")) %>% select(-delete)
dist_data$study <- factor(dist_data$study, levels = c("current", "lbb2012"),
                          labels = c("Current", "Barreiro et al., 2012"))
dist_data$normalized <- factor(dist_data$normalized, levels = c("pre", "norm"),
                               labels = c("Original", "Normalized N(0, 1)"))

ggplot(dist_data, aes(x = gene_exp, color = study)) +
  geom_density() +
  facet_wrap(~normalized, scales = "free") +
  labs(x = "Gene expression level", y = "Density",
       title = "Normalization to compare gene expression\nmeasured with different technologies") +
  theme(legend.position = "bottom")

my_ggsave("combined-distributions.pdf", dims = c(2, 1))
my_ggsave("combined-distributions.png", dims = c(2, 1))

# PCA

p_combined_pca <- ggplot(combined_pca, aes(x = PC1, y = PC2,
                                           color = study, size = treatment)) +
  geom_point(shape = 1) +
  scale_color_discrete(name = "Study", breaks = c("current", "lbb2012"),
                       labels = c("Current", "Barreiro et al., 2012")) +
  scale_size_discrete(name = "Treatment", breaks = c("noninf", "infect"),
                     labels = c("non-infected", "infected")) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(combined_explained[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(combined_explained[2] * 100, 2)),
       title = "PCA before batch correction")

p_combined_pca_regr <- p_combined_pca %+% combined_pca_regr +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(combined_explained_regr[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(combined_explained_regr[2] * 100, 2)),
       title = "PCA after batch correction")

combined_pca_multi <- plot_grid(p_combined_pca + theme(legend.position = c(.5, .5)),
                                p_combined_pca_regr + theme(legend.position = "none"),
                                labels = letters[1:2])

my_ggsave("combined-pca.pdf", dims = c(2, 1))
my_ggsave("combined-pca.png", dims = c(2, 1))

# Classifier -------------------------------------------------------------------
# see code/main-classifier.R

predictions <- readRDS(file.path(data_dir, "classifier-predictions.rds"))

# Plot the metrics of model performance

# Convert nested list of lists into data frame for plotting
extract_train_metrics <- function(x) {
  # x - a list with compenents kappa, recall, precision, f1
  return(c(kappa = x$kappa,
           recall = x$recall,
           precision = x$precision,
           f1 = x$f1,
           separation = x$separation,
           num_genes = x$num_genes))
}
# Step 1: Convert to list of data frames
predict_tmp <- lapply(predictions, function(x) {
    ldply(x, extract_train_metrics, .id = "qvalue")
  })
# Step 2: Convert to data frame
predict_df <- ldply(predict_tmp, .id = "method")

predict_df$method <- factor(predict_df$method,
                            levels = c("glmnet", "svmLinear", "rf"),
                            labels = c("Elastic Net", "Support Vector Machine",
                                       "Random Forest"))
metrics <- c("kappa", "recall", "precision", "f1", "separation")
.simpleCap <- function(x) {
  # from ?chartr
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
for (m in metrics) {
  fname <- paste0("classifier-", m, ".pdf")
  pdf(file.path(fig_dir, fname),
      width = w, height = h, useDingbats = FALSE)
  p <- ggplot(predict_df, aes_string(x = "num_genes", y = m, color = "method")) +
    geom_point() +
    facet_wrap(~method) +
    labs(x = "Number of genes",
         y = .simpleCap(m)) +
    scale_x_continuous(breaks = unique(predict_df$num_genes)) +
    ylim(0, 1) +
    theme_linedraw() +
    theme(legend.position = "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          axis.text.x = element_text(angle = 300, vjust = 0.5))
  print(p)
  invisible(dev.off())
}

# Save a nice version of separation comparison
class_compare <- ggplot(predict_df,
                        aes(x = num_genes, y = separation, color = method)) +
  geom_point() +
  facet_wrap(~method) +
  labs(x = "Number of genes",
       y = "Difference in mean probability of TB resistance\nbetween resistant and susceptible individuals") +
  scale_x_continuous(breaks = unique(predict_df$num_genes)) +
  ylim(0, 1) +
  theme_linedraw(base_size = 22) +
  theme(legend.position = "none",
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_text(angle = 300, vjust = 0.5))

my_ggsave("classifier-compare.pdf", dims = c(3, 1))
my_ggsave("classifier-compare.png", dims = c(3, 1))

# Plot the class probabilites for the different models

# Convert nested list of lists into data frame for plotting
extract_class_probs <- function(x) {
  # x - a list with component pred
  return(x$pred)
}
# Step 1: Convert to list of data frames
class_prob_tmp <- lapply(predictions, function(x) {
  ldply(x, extract_class_probs, .id = "qvalue")
})
# Step 2: Convert to data frame
class_prob_df <- ldply(class_prob_tmp, .id = "method")

class_prob_df$method <- factor(class_prob_df$method,
                            levels = c("glmnet", "svmLinear", "rf"),
                            labels = c("Elastic Net", "Support Vector Machine",
                                       "Random Forest"))
class_prob_df$Observed <- factor(class_prob_df$obs,
                                 levels = c("suscep", "resist"),
                                 labels = c("Susceptible", "Resistant"))
for (method in levels(class_prob_df$method)) {
  for (qval in unique(class_prob_df$qvalue)) {
    fname <- paste0("class-prob-", gsub(" ", "-", tolower(method)),
                    "-", qval, ".pdf")
    pdf(file.path(fig_dir, fname),
        width = w, height = h, useDingbats = FALSE)
    d_tmp <- class_prob_df[class_prob_df$method == method &
                             class_prob_df$qvalue == qval, ]
    # Order the samples by their classified probably of being resistant
    d_tmp$id <- factor(d_tmp$id, levels = d_tmp$id[order(d_tmp$resist)])
    d_tmp <- d_tmp[order(d_tmp$resist), ]
    d_tmp$text_label <- substr(as.character(d_tmp$id), 1, 3)
    d_tmp$text_label[11:nrow(d_tmp)] <- ""
    p <- ggplot(d_tmp, aes(x = id, y = resist)) +
      geom_point(aes(color = Observed)) +
      geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
                size = rel(2)) +
      labs(x = "Individual",
           y = "Assigned probability of being TB resistant",
           title = paste(method, qval)) +
      scale_x_discrete(labels = NULL) +
      ylim(0, 1) +
      theme(legend.position = c(0.75, 0.5),
            axis.ticks.x = element_blank())
    print(p)
    invisible(dev.off())
  }
}

# Plot results of predictions in Barreiro et al., 2012

predictions_lbb <- readRDS(file.path(data_dir, "classifier-predictions-lbb.rds"))

predict_lbb_tmp <- lapply(predictions_lbb, function(x) {
  ldply(x, .id = "qvalue")
})

# Step 1: Convert to list of data frames
predict_lbb_tmp <- lapply(predictions_lbb, function(x) {
  ldply(x, .id = "qvalue")
})
# Step 2: Convert to data frame
predict_lbb_df <- ldply(predict_lbb_tmp, .id = "method")

predict_lbb_df$method <- factor(predict_lbb_df$method,
                               levels = c("glmnet", "svmLinear", "rf"),
                               labels = c("Elastic Net", "Support Vector Machine",
                                          "Random Forest"))
predict_lbb_df$id <- rownames(predictions_lbb[[3]][[1]])
for (method in levels(predict_lbb_df$method)) {
  for (qval in unique(predict_lbb_df$qvalue)) {
    fname <- paste0("class-prob-lbb-", gsub(" ", "-", tolower(method)),
                    "-", qval, ".pdf")
    pdf(file.path(fig_dir, fname),
        width = w, height = h, useDingbats = FALSE)
    d_tmp <- predict_lbb_df[predict_lbb_df$method == method &
                              predict_lbb_df$qvalue == qval, ]
    # Order the samples by their classified probably of being resistant
    d_tmp$id <- factor(d_tmp$id, levels = d_tmp$id[order(d_tmp$resist)])
    p <- ggplot(d_tmp, aes(x = id, y = resist)) +
      geom_point() +
      labs(x = "Individual",
           y = "Assigned probability of being TB resistant",
           title = paste(method, qval)) +
      scale_x_discrete(labels = NULL) +
      ylim(0, 1) +
      theme(legend.position = c(0.75, 0.5),
            axis.ticks.x = element_blank())
    print(p)
    invisible(dev.off())
  }
}

# Make nicer plots for paper
qval <- "q0.05"
class_prob_df$text_label <- substr(as.character(class_prob_df$id), 1, 3)

# Elastic Net
en <- class_prob_df[class_prob_df$method == "Elastic Net" &
                         class_prob_df$qvalue == qval, ]

p_base <- ggplot(en, aes(x = reorder(id, resist), y = resist)) +
  geom_point() +
  geom_hline(yintercept = 0.75, color = "red", linetype = "dashed") +
  labs(x = "Individual",
       y = "Assigned probability of being TB resistant") +
  scale_x_discrete(labels = NULL) +
  scale_y_continuous(limits = c(0, 1.01), breaks = seq(0, 1, by = 0.1)) +
  theme(axis.ticks.x = element_blank())
p_en <- p_base +
  geom_point(aes(color = Observed)) +
  scale_color_manual(values = c(col_suscep, col_resist)) +
  geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
            size = rel(2)) +
  theme(legend.position = c(0.75, 0.5)) +
  labs(title = "Training classifier on current data")

en_lbb <- predict_lbb_df[predict_lbb_df$method == "Elastic Net" &
                           predict_lbb_df$qvalue == qval, ]

p_en_lbb <- p_base %+% en_lbb +
  labs(title = "Testing classifier on Barreiro et al. data")

multi_en <- plot_grid(p_en,
                      p_en_lbb,
                      labels = letters[1:2])

my_ggsave("classifier-en.pdf", dims = c(2, 1))
my_ggsave("classifier-en.png", dims = c(2, 1))

# Support Vector Machine
svm <- class_prob_df[class_prob_df$method == "Support Vector Machine" &
                      class_prob_df$qvalue == qval, ]

p_svm <- p_base %+% svm +
  geom_point(aes(color = Observed)) +
  scale_color_manual(values = c(col_suscep, col_resist)) +
  geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
            size = rel(2)) +
  theme(legend.position = c(0.75, 0.5)) +
  labs(title = "Training classifier on current data")

svm_lbb <- predict_lbb_df[predict_lbb_df$method == "Support Vector Machine" &
                          predict_lbb_df$qvalue == qval, ]

p_svm_lbb <- p_base %+% svm_lbb +
  labs(title = "Testing classifier on Barreiro et al. data")

multi_svm <- plot_grid(p_svm,
                      p_svm_lbb,
                      labels = letters[1:2])

my_ggsave("classifier-svm.eps", dims = c(2, 1))
my_ggsave("classifier-svm.pdf", dims = c(2, 1))
my_ggsave("classifier-svm.png", dims = c(2, 1))

# Random Forest
rf <- class_prob_df[class_prob_df$method == "Random Forest" &
                       class_prob_df$qvalue == qval, ]

p_rf <- p_base %+% rf +
  geom_point(aes(color = Observed)) +
  scale_color_manual(values = c(col_suscep, col_resist)) +
  geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
            size = rel(2)) +
  theme(legend.position = c(0.75, 0.5)) +
  labs(title = "Training classifier on current data")

rf_lbb <- predict_lbb_df[predict_lbb_df$method == "Random Forest" &
                            predict_lbb_df$qvalue == qval, ]

p_rf_lbb <- p_base %+% rf_lbb +
  labs(title = "Testing classifier on Barreiro et al. data")

multi_rf <- plot_grid(p_rf,
                       p_rf_lbb,
                       labels = letters[1:2])

my_ggsave("classifier-rf.pdf", dims = c(2, 1))
my_ggsave("classifier-rf.png", dims = c(2, 1))

# Compare expression of classifier genes between the two studies

limma_comb <- readRDS(file.path(data_dir, "combined-limma.rds"))
genes_q05 <- rownames(limma_comb[["status_ni"]][limma_comb[["status_ni"]]$qvalue < 0.05, ])
stopifnot(length(genes_q05) == predict_df$num_genes[predict_df$qvalue == "q0.05"][1])
# regressed and anno_comb were loaded above
training_medians <- apply(regressed[genes_q05,
                                    anno_comb$treatment == "noninf" &
                                    anno_comb$study == "current"],
                          1, median)
testing_medians <- apply(regressed[genes_q05,
                                   anno_comb$treatment == "noninf" &
                                   anno_comb$study == "lbb2012"],
                         1, median)
p_class_exp <- ggplot(data.frame(training_medians, testing_medians),
                      aes(x = training_medians, y = testing_medians)) +
  geom_point(shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Current study", y = "Barreiro et al., 2012",
       title = "Comparing median expression levels\nof genes used in classifer")

my_ggsave("classifier-exp.pdf")
my_ggsave("classifier-exp.png")

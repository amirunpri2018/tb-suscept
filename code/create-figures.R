#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("gplots"))
suppressPackageStartupMessages(library("VennDiagram"))
suppressPackageStartupMessages(library("ggbeeswarm"))
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
outliers_simple

# PCA

# PC1 versus PC2.
pc1v2_out <- ggplot(pca_outliers, aes(x = PC1, y = PC2, color = treatment)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained_outliers[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained_outliers[2] * 100, 2)))
# PC3 versus PC4.
pc3v4_out <- ggplot(pca_outliers, aes(x = PC3, y = PC4, color = treatment)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 3, round(explained_outliers[3] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 4, round(explained_outliers[4] * 100, 2)))
# PC5 versus PC6.
pc5v6_out <- ggplot(pca_outliers, aes(x = PC5, y = PC6, color = treatment)) +
  geom_text(aes(label = individual)) +
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
                             labels = c("noninfected", "infected"))
explained <- readRDS(file.path(data_dir, "results-pca-explained.rds"))

# PC1 versus PC2.
pc1v2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained[2] * 100, 2)))

# PC3 versus PC4.
pc3v4 <- ggplot(pca_data, aes(x = PC3, y = PC4, color = infection)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 3, round(explained[3] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 4, round(explained[4] * 100, 2)))

# PC5 versus PC6.
pc5v6 <- ggplot(pca_data, aes(x = PC5, y = PC6, color = infection)) +
  geom_text(aes(label = individual)) +
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
  scale_fill_discrete(name = "Treatment") +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "PC3 varies by date of infection",
       x = "Date of infection experiment")

pc3_status <- ggplot(pca_data, aes(x = reorder(individual, PC3), y = PC3, color = status)) +
  geom_point() +
  scale_color_discrete(name = "Status") +
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

# Difference between susceptible and resistant individuals in the noninfected state
title_status_ni <- "Difference between susceptible and resistant\nindividuals in the noninfected state"
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
  category = c("Difference in noninfected state",
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
  category = c("Difference in noninfected state",
               "Effect of treatment in resistant individuals",
               "Effect of treatment in susceptible individuals"),
  #   fill = c("darkolivegreen3", "cadetblue3", "darkorchid"),
  alpha = c(0.5, 0.5, 0.5), col = "black", cex = 2, cat.cex = 1.5,
  euler.d = FALSE, scaled = FALSE, ind = FALSE)
# grid.draw(venn_treat)

# GWAS enrichment --------------------------------------------------------------
# see code/main-gwas.R

gwas_results <- read.delim(file.path(data_dir, "results-gwas.txt"),
                           stringsAsFactors = FALSE)
gwas_lm <- read.delim(file.path(data_dir, "results-gwas-lm.txt"),
                      stringsAsFactors = FALSE)

# Gambia
gwas_scatter_gambia <- ggplot(gwas_results,
                              aes(x = status_ni, y = gwas_p_gambia)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm", col = "red", se = FALSE) +
  labs(x = "|logFC| between susceptible and resistant\nindividuals in the noninfected state",
       y = "GWAS p-value",
       title = "GWAS p-value vs. DE effect size")

# Remove interaction term b/c basically the same as noninfected state
gwas_lm <- gwas_lm[gwas_lm$test != "interact", ]
# Remove mean_expression_level b/c no longer clear this is actually an
# informative comparison.
gwas_lm <- gwas_lm[gwas_lm$test != "mean_expression_level", ]
gwas_lm$test <- factor(gwas_lm$test,
                       levels = c("treat_suscep",
                                  "treat_resist",
                                  "status_ni",
                                  "status_ii"),
                       labels = c("|logFC| between\ninfected and noninfected states\nin susceptible individuals",
                                  "|logFC| between\ninfected and noninfected states\nin resistant individuals",
                                  "|logFC| between\nsusceptible and resistant\nindividuals in the noninfected state",
                                  "|logFC| between\nsusceptible and resistant\nindividuals in the infected state"))
gwas_lm$population <- factor(gwas_lm$population, levels = c("gambia", "ghana"),
                             labels = c("Gambia", "Ghana"))
gwas_slopes <- ggplot(gwas_lm, aes(x = test, y = slope, color = population)) +
  geom_point(position = position_dodge(width = .5)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(ymin = slope - 2 * slope_se, ymax = slope + 2 * slope_se),
                 position = position_dodge(width = .5)) +
  scale_color_discrete(name = "Population") +
  labs(x = "", y = "Slope of best fit line (+/- 2 x standard error)",
       title = "Relationship between GWAS p-value\nand |logFC| in DCs") +
  coord_flip() +
  theme(legend.position = "bottom")

gwas_multi_gambia <- plot_grid(gwas_scatter_gambia,
                        gwas_slopes,
                        labels = letters[1:2])

my_ggsave("gwas.eps", dims = c(2, 1))
my_ggsave("gwas.pdf", dims = c(2, 1))
my_ggsave("gwas.png", dims = c(2, 1))

# Correlation with number of SNPs per gene
# cor(gwas_results$n_snps, gwas_results$gwas_p_gambia)
gwas_n_snps_gambia <- ggplot(gwas_results, aes(x = n_snps, y = gwas_p_gambia)) +
  geom_point(shape = 1) +
  geom_smooth(method = "lm",col = "red", se = FALSE) +
  labs(x = "Number of SNPs nearby gene",
       y = "Minimum GWAS p-value",
       title = "Relationship between number of tested SNPs near gene\nand the minimum GWAS p-value of these SNPs - Gambia")
gwas_n_snps_ghana <- gwas_n_snps_gambia %+% aes(y = gwas_p_ghana) +
  labs(y = "Minimum GWAS p-value",
       title = "Relationship between number of tested SNPs near gene\nand the minimum GWAS p-value of these SNPs - Ghana")

n_snps_multi <- plot_grid(gwas_n_snps_gambia,
                          gwas_n_snps_ghana,
                          labels = letters[1:2])

my_ggsave("gwas-n-snps.pdf", dims = c(2, 1))
my_ggsave("gwas-n-snps.png", dims = c(2, 1))

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
metrics <- c("kappa", "recall", "precision", "f1")
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
    scale_x_continuous(breaks = predict_df$num_genes) +
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
    d_tmp$text_label <- substr(as.character(d_tmp$id), 1, 5)
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

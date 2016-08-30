#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("edgeR"))
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

# Batch effects ----------------------------------------------------------------
# see code/qc-batch.R

pca_data <- read.delim(file.path(data_dir, "results-pca.txt"),
                       stringsAsFactors = FALSE, row.names = 1)
pca_data$status <- factor(pca_data$status, levels = c("contact", "tb"),
                          labels = c("resistant", "susceptible"))
pca_data$treatment <- factor(pca_data$treatment, levels = c("none", "infected"))
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
                       labels = LETTERS[1:4])

postscript(file.path(fig_dir, "batch-pca.eps"),
           width = 2 * w, height = 2 * h)
p_batch_1
invisible(dev.off())

png(file.path(fig_dir, "batch-pca.png"),
           width = 2 * w, height = 2 * h, units = "in", res = 72)
p_batch_1
invisible(dev.off())

# Further investigate correlation of infection batch with PC3 and PC5

pc3_infection <- ggplot(pca_data, aes(x = infection, y = PC3, fill = treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "PC3 varies by date of infection")

pc3_status <- ggplot(pca_data, aes(x = reorder(individual, PC3), y = PC3, color = status)) +
  geom_point() +
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
                       labels = LETTERS[1:4])

postscript(file.path(fig_dir, "batch-infection.eps"),
           width = 2 * w, height = 2 * h)
p_batch_2
invisible(dev.off())

png(file.path(fig_dir, "batch-infection.png"),
    width = 2 * w, height = 2 * h, units = "in", res = 72)
p_batch_2
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
  geom_point(size = 0.5) +
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
ma_diff_before <- plot_ma(results[["diff_before"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = "Difference between susceptible and resistant\nindividuals before treatment")
ma_diff_before
invisible(dev.off())
pdf(file.path(fig_dir, "pval-diff-before.pdf"),
    width = w, height = h, useDingbats = FALSE)
pval_diff_before <- plot_pval_hist(results[["diff_before"]]) +
  ylim(0, 325) +
  labs(title = "Difference between susceptible and resistant\nindividuals before treatment")
pval_diff_before
invisible(dev.off())

# Difference between susceptible and resistant individuals after treatment
pdf(file.path(fig_dir, "ma-diff-after.pdf"),
    width = w, height = h, useDingbats = FALSE)
ma_diff_after <- plot_ma(results[["diff_after"]], qval = 0.1) +
  ylim(-3.5, 3.5) +
  labs(title = "Difference between susceptible and resistant\nindividuals after treatment")
ma_diff_after
invisible(dev.off())
pdf(file.path(fig_dir, "pval-diff-after.pdf"),
    width = w, height = h, useDingbats = FALSE)
pval_diff_after <- plot_pval_hist(results[["diff_after"]]) +
  ylim(0, 325) +
  labs(title = "Difference between susceptible and resistant\nindividuals after treatment")
pval_diff_after
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
  labs(title = "Difference in effect of treatment between\nsusceptible and resistant individuals")
invisible(dev.off())
pdf(file.path(fig_dir, "pval-diff-treat.pdf"),
    width = w, height = h, useDingbats = FALSE)
plot_pval_hist(results[["diff_treat"]]) +
  labs(title = "Difference in effect of treatment between\nsusceptible and resistant individuals")
invisible(dev.off())

# Multipanel figure
p_limma <- plot_grid(pval_diff_before,
                     pval_diff_after,
                     ma_diff_before + labs(title = ""),
                     ma_diff_after + labs(title = ""),
                     labels = letters[1:4])

postscript(file.path(fig_dir, "limma.eps"),
           width = 2 * w, height = 2 * h)
p_limma
invisible(dev.off())

png(file.path(fig_dir, "limma.png"),
    width = 2 * w, height = 2 * h, units = "in", res = 72)
p_limma
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

# GWAS enrichment --------------------------------------------------------------
# see code/main-gwas.R

gwas_results <- read.delim(file.path(data_dir, "results-gwas.txt"),
                           stringsAsFactors = FALSE)
gwas_lm <- read.delim(file.path(data_dir, "results-gwas-lm.txt"),
                      stringsAsFactors = FALSE)

gwas_scatter <- ggplot(gwas_results, aes(x = diff_before, y = gwas_p_gambia)) +
  geom_point(alpha = 0.25) +
  geom_smooth(method = "lm", col = "red") +
  labs(x = "|logFC| between susceptible and resistant\nindividuals in the noninfected state",
       y = "GWAS p-value",
       title = "GWAS p-value vs. DE effect size")

# Remove interaction term b/c basically the same as noninfected state
gwas_lm <- gwas_lm[gwas_lm$test != "diff_treat", ]
gwas_lm$test <- factor(gwas_lm$test,
                       levels = c("mean_expression_level",
                                  "treat_suscept",
                                  "treat_resist",
                                  "diff_before",
                                  "diff_after"),
                       labels = c("Mean expression level",
                                  "|logFC| between\ninfected and noninfected states\nin susceptible individuals",
                                  "|logFC| between\ninfected and noninfected states\nin resistant individuals",
                                  "|logFC| between\nsusceptible and resistant\nindividuals in the noninfected state",
                                  "|logFC| between\nsusceptible and resistant\nindividuals in the infected state"))
gwas_slopes_gambia <- ggplot(gwas_lm[gwas_lm$population == "gambia", ],
                             aes(x = test, y = slope)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_linerange(aes(ymin = slope - slope_se, ymax = slope + slope_se)) +
  labs(x = "", y = "Slope of best fit line (+/- standard error)",
       title = "Relationship between GWAS p-value\nand |logFC| in DCs") +
  coord_flip()

gwas_multi <- plot_grid(gwas_scatter, gwas_slopes_gambia, labels = letters[1:2])

pdf(file.path(fig_dir, "gwas.pdf"),
           width = 2 * w, height = h)
gwas_multi
invisible(dev.off())

png(file.path(fig_dir, "gwas.png"),
    width = 2 * w, height = h, units = "in", res = 72)
gwas_multi
invisible(dev.off())

# Correlation with number of SNPs per gene
# cor(gwas_results$n_snps, gwas_results$gwas_p_gambia)
# plot(gwas_results$n_snps, gwas_results$gwas_p_gambia,
#      xlab = "Number of SNPs nearby gene",
#      ylab = "Minimum GWAS p-value",
#      main = "Relationship between number of tested SNPs near gene\nand the minimum GWAS p-value of these SNPs")


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
    ldply(x, extract_train_metrics, .id = "svalue")
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
  ldply(x, extract_class_probs, .id = "svalue")
})
# Step 2: Convert to data frame
class_prob_df <- ldply(class_prob_tmp, .id = "method")

class_prob_df$method <- factor(class_prob_df$method,
                            levels = c("glmnet", "svmLinear", "rf"),
                            labels = c("Elastic Net", "Support Vector Machine",
                                       "Random Forest"))
class_prob_df$Observed <- factor(class_prob_df$obs,
                                 levels = c("tb", "contact"),
                                 labels = c("Susceptible", "Resistant"))
for (method in levels(class_prob_df$method)) {
  for (sval in unique(class_prob_df$svalue)) {
    fname <- paste0("class-prob-", gsub(" ", "-", tolower(method)),
                    "-", sval, ".pdf")
    pdf(file.path(fig_dir, fname),
        width = w, height = h, useDingbats = FALSE)
    d_tmp <- class_prob_df[class_prob_df$method == method &
                             class_prob_df$svalue == sval, ]
    # Order the samples by their classified probably of being resistant
    d_tmp$id <- factor(d_tmp$id, levels = d_tmp$id[order(d_tmp$contact)])
    d_tmp <- d_tmp[order(d_tmp$contact), ]
    d_tmp$text_label <- substr(as.character(d_tmp$id), 1, 5)
    d_tmp$text_label[11:nrow(d_tmp)] <- ""
    p <- ggplot(d_tmp, aes(x = id, y = contact)) +
      geom_point(aes(color = Observed)) +
      geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
                size = rel(2)) +
      labs(x = "Individual",
           y = "Assigned probability of being TB resistant",
           title = paste(method, sval)) +
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
  ldply(x, .id = "svalue")
})

# Step 1: Convert to list of data frames
predict_lbb_tmp <- lapply(predictions_lbb, function(x) {
  ldply(x, .id = "svalue")
})
# Step 2: Convert to data frame
predict_lbb_df <- ldply(predict_lbb_tmp, .id = "method")

predict_lbb_df$method <- factor(predict_lbb_df$method,
                               levels = c("glmnet", "svmLinear", "rf"),
                               labels = c("Elastic Net", "Support Vector Machine",
                                          "Random Forest"))
predict_lbb_df$id <- rownames(predictions_lbb[[3]][[1]])
for (method in levels(predict_lbb_df$method)) {
  for (sval in unique(predict_lbb_df$svalue)) {
    fname <- paste0("class-prob-lbb-", gsub(" ", "-", tolower(method)),
                    "-", sval, ".pdf")
    pdf(file.path(fig_dir, fname),
        width = w, height = h, useDingbats = FALSE)
    d_tmp <- predict_lbb_df[predict_lbb_df$method == method &
                              predict_lbb_df$svalue == sval, ]
    # Order the samples by their classified probably of being resistant
    d_tmp$id <- factor(d_tmp$id, levels = d_tmp$id[order(d_tmp$contact)])
    p <- ggplot(d_tmp, aes(x = id, y = contact)) +
      geom_point() +
      labs(x = "Individual",
           y = "Assigned probability of being TB resistant",
           title = paste(method, sval)) +
      scale_x_discrete(labels = NULL) +
      ylim(0, 1) +
      theme(legend.position = c(0.75, 0.5),
            axis.ticks.x = element_blank())
    print(p)
    invisible(dev.off())
  }
}

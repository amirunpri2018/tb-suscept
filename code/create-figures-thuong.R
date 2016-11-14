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

# Combine studies --------------------------------------------------------------
# see code/combine-studies.R

# Combined data
exp_combined <- read.delim(file.path(data_dir, "combined-raw-thuong.txt"),
                           row.names = 1, check.names = FALSE)
anno_comb <- read.delim(file.path(data_dir, "combined-annotation-thuong.txt"),
                        stringsAsFactors = FALSE, row.names = 1)

# Normalized data
exp_combined_norm <- read.delim(file.path(data_dir, "combined-normalized-thuong.txt"),
                                row.names = 1, check.names = FALSE)
combined_pca <- read.delim(file.path(data_dir, "combined-pca-thuong.txt"),
                           row.names = 1, stringsAsFactors = FALSE)
combined_explained <- readRDS(file.path(data_dir, "combined-pca-explained-thuong.txt"))

# Regressed data
regressed <- read.delim(file.path(data_dir, "combined-regressed-thuong.txt"),
                        row.names = 1, check.names = FALSE)
combined_pca_regr <- read.delim(file.path(data_dir, "combined-pca-regressed-thuong.txt"),
                                row.names = 1, check.names = FALSE)
combined_explained_regr <- readRDS(file.path(data_dir, "combined-pca-explained-regressed-thuong.txt"))

# Plots of distributions
dist_current_pre <- apply(exp_combined[, anno_comb$study == "current"], 1, median)
dist_thuong_pre <- apply(exp_combined[, anno_comb$study == "thuong2008"], 1, median)
dist_current_norm <- apply(exp_combined_norm[, anno_comb$study == "current"],
                           1, median)
dist_thuong_norm <- apply(exp_combined_norm[, anno_comb$study == "thuong2008"],
                           1, median)
dist_data <- gather(data.frame(dist_current_pre, dist_thuong_pre,
                        dist_current_norm, dist_thuong_norm),
                    key = "dist", value = "gene_exp")
dist_data <- separate(dist_data, col = dist,
                      into = c("delete", "study", "normalized")) %>% select(-delete)
dist_data$study <- factor(dist_data$study, levels = c("current", "thuong"),
                          labels = c("Current", "Thuong et al., 2008"))
dist_data$normalized <- factor(dist_data$normalized, levels = c("pre", "norm"),
                               labels = c("Original", "Normalized N(0, 1)"))

ggplot(dist_data, aes(x = gene_exp, color = study)) +
  geom_density() +
  facet_wrap(~normalized, scales = "free") +
  labs(x = "Gene expression level", y = "Density",
       title = "Normalization to compare gene expression\nmeasured with different technologies") +
  theme(legend.position = "bottom")

# my_ggsave("combined-distributions.pdf", dims = c(2, 1))
# my_ggsave("combined-distributions.png", dims = c(2, 1))

# PCA

p_combined_pca <- ggplot(combined_pca, aes(x = PC1, y = PC2,
                                           color = study, size = treatment)) +
  geom_point(shape = 1) +
  scale_color_discrete(name = "Study", breaks = c("current", "thuong2008"),
                       labels = c("Current", "Thuong et al., 2008")) +
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

# my_ggsave("combined-pca.pdf", dims = c(2, 1))
# my_ggsave("combined-pca.png", dims = c(2, 1))

# Classifier -------------------------------------------------------------------
# see code/main-classifier.R

predictions <- readRDS(file.path(data_dir, "classifier-predictions-train-thuong.rds"))

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
#   pdf(file.path(fig_dir, fname),
#       width = w, height = h, useDingbats = FALSE)
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
  # invisible(dev.off())
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

# my_ggsave("classifier-compare.pdf", dims = c(3, 1))
# my_ggsave("classifier-compare.png", dims = c(3, 1))

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
#     pdf(file.path(fig_dir, fname),
#         width = w, height = h, useDingbats = FALSE)
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
    # invisible(dev.off())
  }
}

# Plot results of predictions in Thuong et al., 2008

predictions_thuong <- readRDS(file.path(data_dir, "classifier-predictions-test-thuong.rds"))

predict_thuong_tmp <- lapply(predictions_thuong, function(x) {
  ldply(x, .id = "qvalue")
})

# Step 1: Convert to list of data frames
predict_thuong_tmp <- lapply(predictions_thuong, function(x) {
  ldply(x, .id = "qvalue")
})
# Step 2: Convert to data frame
predict_thuong_df <- ldply(predict_thuong_tmp, .id = "method")

predict_thuong_df$method <- factor(predict_thuong_df$method,
                               levels = c("glmnet", "svmLinear", "rf"),
                               labels = c("Elastic Net", "Support Vector Machine",
                                          "Random Forest"))
predict_thuong_df$id <- rownames(predictions_thuong[[3]][[1]])
predict_thuong_df <- separate(predict_thuong_df, col = id, remove = FALSE,
                              into = c("individual", "status", "treatment"))
predict_thuong_df$Observed <- factor(predict_thuong_df$status,
                                levels = c("suscep", "resist"),
                                labels = c("Susceptible", "Resistant"))
for (method in levels(predict_thuong_df$method)) {
  for (qval in unique(predict_thuong_df$qvalue)) {
    fname <- paste0("class-prob-thuong-", gsub(" ", "-", tolower(method)),
                    "-", qval, ".pdf")
#     pdf(file.path(fig_dir, fname),
#         width = w, height = h, useDingbats = FALSE)
    d_tmp <- predict_thuong_df[predict_thuong_df$method == method &
                              predict_thuong_df$qvalue == qval, ]
    # Order the samples by their classified probably of being resistant
    d_tmp$id <- factor(d_tmp$id, levels = d_tmp$id[order(d_tmp$resist)])
    p <- ggplot(d_tmp, aes(x = id, y = resist, color = status)) +
      geom_text(aes(label = individual)) +
      labs(x = "Individual",
           y = "Assigned probability of being TB resistant",
           title = paste(method, qval)) +
      scale_x_discrete(labels = NULL) +
      ylim(0, 1) +
      theme(legend.position = c(0.75, 0.5),
            axis.ticks.x = element_blank())
    print(p)
    # invisible(dev.off())
  }
}

# Make nicer plots for paper
qval <- "q0.05"
class_prob_df$text_label <- substr(as.character(class_prob_df$id), 1, 3)
predict_thuong_df$text_label <- predict_thuong_df$individual

# Elastic Net
en <- class_prob_df[class_prob_df$method == "Elastic Net" &
                         class_prob_df$qvalue == qval, ]

p_base <- ggplot(en, aes(x = reorder(id, 1 - resist), y = 1 - resist)) +
  geom_point() +
  geom_hline(yintercept = 0.25, color = "red", linetype = "dashed") +
  labs(x = "Individual",
       y = "Assigned probability of being TB susceptible") +
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

en_thuong <- predict_thuong_df[predict_thuong_df$method == "Elastic Net" &
                           predict_thuong_df$qvalue == qval, ]

p_en_thuong <- p_en %+% en_thuong +
  labs(title = "Testing classifier on Thuong et al. data")

multi_en <- plot_grid(p_en,
                      p_en_thuong,
                      labels = letters[1:2])

# my_ggsave("classifier-en.pdf", dims = c(2, 1))
my_ggsave("classifier-en-thuong.png", dims = c(2, 1))

# Support Vector Machine
svm <- class_prob_df[class_prob_df$method == "Support Vector Machine" &
                      class_prob_df$qvalue == qval, ]

p_svm <- p_base %+% svm +
  geom_point(aes(color = Observed)) +
  scale_color_manual(values = c(col_suscep, col_resist)) +
#   geom_text(aes(label = text_label), nudge_x = -0.1, nudge_y = 0.01,
#             size = rel(2)) +
  theme(legend.position = c(0.25, 0.75)) +
  labs(title = "Training classifier on current data")

svm_thuong <- predict_thuong_df[predict_thuong_df$method == "Support Vector Machine" &
                          predict_thuong_df$qvalue == qval, ]

p_svm_thuong <- p_svm %+% svm_thuong +
  labs(title = "Testing classifier on Thuong et al. data")

multi_svm <- plot_grid(p_svm,
                      p_svm_thuong,
                      labels = letters[1:2])

# my_ggsave("classifier-svm.eps", dims = c(2, 1))
# my_ggsave("classifier-svm.pdf", dims = c(2, 1))
my_ggsave("classifier-svm-thuong.png", dims = c(2, 1))

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

rf_thuong <- predict_thuong_df[predict_thuong_df$method == "Random Forest" &
                            predict_thuong_df$qvalue == qval, ]

p_rf_thuong <- p_rf %+% rf_thuong +
  labs(title = "Testing classifier on Thuong et al. data")

multi_rf <- plot_grid(p_rf,
                       p_rf_thuong,
                       labels = letters[1:2])

# my_ggsave("classifier-rf.pdf", dims = c(2, 1))
my_ggsave("classifier-rf-thuong.png", dims = c(2, 1))

# Compare expression of classifier genes between the two studies

limma_comb <- readRDS(file.path(data_dir, "combined-limma-thuong.rds"))
genes_q05 <- rownames(limma_comb[["status_ni"]][limma_comb[["status_ni"]]$qvalue < 0.05, ])
stopifnot(length(genes_q05) == predict_df$num_genes[predict_df$qvalue == "q0.05"][1])
# regressed and anno_comb were loaded above
training_medians <- apply(regressed[genes_q05,
                                    anno_comb$treatment == "noninf" &
                                    anno_comb$study == "current"],
                          1, median)
testing_medians <- apply(regressed[genes_q05,
                                   anno_comb$treatment == "noninf" &
                                   anno_comb$study == "thuong2008"],
                         1, median)
p_class_exp <- ggplot(data.frame(training_medians, testing_medians),
                      aes(x = training_medians, y = testing_medians)) +
  geom_point(shape = 1) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(x = "Current study", y = "Thuong et al., 2008",
       title = "Comparing median expression levels\nof genes used in classifer")

# my_ggsave("classifier-exp.pdf")
my_ggsave("classifier-exp-thuong.png")

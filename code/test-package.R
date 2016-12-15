pkg_relative_path <- "../pkg/tbsuscept/"
devtools::document(pkg_relative_path)
suppressPackageStartupMessages(library("data.table"))
library("stringr")

if(interactive()) {
  data_dir <- "../data"
  fig_dir <- "../figure"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
  fig_dir <- args[2]
}
stopifnot(dir.exists(data_dir))

# Input ------------------------------------------------------------------------

fit <- readRDS(file.path(data_dir, "results-limma-fit.rds"))

gwas_thye_ghana <- fread(file.path(data_dir, "OUT_PLINK_Ghana.txt"),
                         data.table = FALSE, verbose = FALSE)

gwas_thye_gambia <- fread(file.path(data_dir, "OUT_PLINK_Gambia.txt"),
                          data.table = FALSE, verbose = FALSE)

# Missing allele frequencies in column 4 are coded as '.'
gwas_height <- fread(file.path(data_dir, "GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt"),
                         data.table = FALSE, verbose = FALSE, na.strings = ".")

gwas_russia <- fread(file.path(data_dir, "Export_TB-GWAS.txt"), sep = ",",
                     data.table = FALSE, verbose = FALSE)

gwas_sobota <- fread(file.path(data_dir, "sobota2016.txt"),
                     data.table = FALSE, verbose = FALSE)

gene_names <- rownames(fit$coefficients)

tss_all_fname <- file.path(data_dir, "tss-all.rds")

snp_coords_fname_ghana <- file.path(data_dir, "/snp-coords-ghana.rds")

snp_coords_fname_gambia <- file.path(data_dir, "/snp-coords-gambia.rds")

snp_coords_fname_height <- file.path(data_dir, "/snp-coords-height.rds")

snp_coords_fname_russia <- file.path(data_dir, "/snp-coords-russia.rds")

snp_coords_fname_sobota <- file.path(data_dir, "/snp-coords-sobota.rds")

# Ghana ------------------------------------------------------------------------

if (!file.exists(file.path(data_dir, "gwas-ghana.rds"))) {
  result_ghana <- list()
  for (i in 1:4) {
    test_name <- colnames(fit$coefficients)[i]
    result_ghana[[test_name]] <- run_gwas_enrich(gene_names = gene_names,
                                                 tss_all_fname = tss_all_fname,
                                                 window_size = 50000,
                                                 rsid = gwas_thye_ghana$SNP,
                                                 snp_coords_fname = snp_coords_fname_ghana,
                                                 pval = gwas_thye_ghana$PVAL,
                                                 effect_size = fit$coefficients[, i])
    plot(result_ghana[[test_name]], main = sprintf("Ghana: %s", test_name))
  }
  hist(gwas_thye_ghana$PVAL)
  saveRDS(result_ghana, file = file.path(data_dir, "gwas-ghana.rds"))
} else {
  result_ghana <- readRDS(file.path(data_dir, "gwas-ghana.rds"))
}

# The Gambia -------------------------------------------------------------------

if (!file.exists(file.path(data_dir, "gwas-gambia.rds"))) {
  result_gambia <- list()
  for (i in 1:4) {
    test_name <- colnames(fit$coefficients)[i]
    result_gambia[[test_name]] <- run_gwas_enrich(gene_names = gene_names,
                                                  tss_all_fname = tss_all_fname,
                                                  window_size = 50000,
                                                  rsid = gwas_thye_gambia$SNP,
                                                  snp_coords_fname = snp_coords_fname_gambia,
                                                  pval = gwas_thye_gambia$PVAL,
                                                  effect_size = fit$coefficients[, i])
    plot(result_gambia[[test_name]], main = sprintf("The Gambia: %s", test_name))
  }
  hist(gwas_thye_gambia$PVAL)
  saveRDS(result_gambia, file = file.path(data_dir, "gwas-gambia.rds"))
} else {
  result_gambia <- readRDS(file.path(data_dir, "gwas-gambia.rds"))
}

# Height -----------------------------------------------------------------------

if (!file.exists(file.path(data_dir, "gwas-height.rds"))) {
  result_height <- list()
  for (i in 1:4) {
    test_name <- colnames(fit$coefficients)[i]
    result_height[[test_name]] <- run_gwas_enrich(gene_names = gene_names,
                                                  tss_all_fname = tss_all_fname,
                                                  window_size = 50000,
                                                  rsid = gwas_height$MarkerName,
                                                  snp_coords_fname = snp_coords_fname_height,
                                                  pval = gwas_height$p,
                                                  effect_size = fit$coefficients[, i])
    
    plot(result_height[[test_name]], main = sprintf("Height from Lango Allen et al., 2010: %s",
                                                    test_name))
  }
  hist(gwas_height$p)
  saveRDS(result_height, file = file.path(data_dir, "gwas-height.rds"))
} else {
  result_height <- readRDS(file.path(data_dir, "gwas-height.rds"))
}

# Russia -----------------------------------------------------------------------

if (!file.exists(file.path(data_dir, "gwas-russia.rds"))) {
  result_russia <- list()
  for (i in 1:4) {
    test_name <- colnames(fit$coefficients)[i]
    result_russia[[test_name]] <- run_gwas_enrich(gene_names = gene_names,
                                                  tss_all_fname = tss_all_fname,
                                                  window_size = 50000,
                                                  rsid = gwas_russia$SNP[grepl("rs", gwas_russia$SNP)],
                                                  snp_coords_fname = snp_coords_fname_russia,
                                                  pval = gwas_russia$P[grepl("rs", gwas_russia$SNP)],
                                                  effect_size = fit$coefficients[, i])
    
    plot(result_russia[[test_name]], main = sprintf("Russia: %s", test_name))
  }
  hist(gwas_russia$P[grepl("rs", gwas_russia$SNP)])
  saveRDS(result_russia, file = file.path(data_dir, "gwas-russia.rds"))
} else {
  result_russia <- readRDS(file.path(data_dir, "gwas-russia.rds"))
}

# Sobota -----------------------------------------------------------------------

if (!file.exists(file.path(data_dir, "gwas-sobota.rds"))) {
  gwas_sobota_filtered <- gwas_sobota[grepl("rs", gwas_sobota$SNP) & !is.na(gwas_sobota$P), ]
  result_sobota <- list()
  for (i in 1:4) {
    test_name <- colnames(fit$coefficients)[i]
    result_sobota[[test_name]] <- run_gwas_enrich(gene_names = gene_names,
                                                  tss_all_fname = tss_all_fname,
                                                  window_size = 50000,
                                                  rsid = gwas_sobota_filtered$SNP,
                                                  snp_coords_fname = snp_coords_fname_sobota,
                                                  pval = gwas_sobota_filtered$P,
                                                  effect_size = fit$coefficients[, i])
    
    plot(result_sobota[[test_name]], main = sprintf("Sobota: %s", test_name))
  }
  hist(gwas_sobota_filtered$P)
  saveRDS(result_sobota, file = file.path(data_dir, "gwas-sobota.rds"))
} else {
  result_sobota <- readRDS(file.path(data_dir, "gwas-sobota.rds"))
}

# Exploratory Figures ----------------------------------------------------------

pdf("gwas-figures.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))
plot_titles <- c("Difference between susceptible and resistant\nindividuals in the non-infected state",
                 "Difference between susceptible and resistant\nindividuals in the infected state",
                 "Effect of treatment in resistant individuals",
                 "Effect of treatment in susceptible individuals")
for (i in 1:4) {
  boxplot_enrich(auc = c(result_russia[[i]]$auc$main,
                         result_gambia[[i]]$auc$main,
                         result_ghana[[i]]$auc$main,
                         result_sobota[[i]]$auc$main,
                         result_height[[i]]$auc$main),
                 permutations = cbind(result_russia[[i]]$auc$permutations,
                                      result_gambia[[i]]$auc$permutations,
                                      result_ghana[[i]]$auc$permutations,
                                      result_sobota[[i]]$auc$permutations,
                                      result_height[[i]]$auc$permutations),
                 main = plot_titles[i],
                 xaxt = "n")
  axis(1, at = 1:5, labels = FALSE)
  lab <- c("TB Russia", "TB Gambia", "TB Ghana", "TB Uganda",
           "Height\nEurope")
  text(x = 1:5, y = par()$usr[3] - 0.1 * (par()$usr[4] - par()$usr[3]),
       labels = lab, srt = 45, adj = 1, xpd = TRUE)
}
for (i in 1:4) {
  test_name <- colnames(fit$coefficients)[i]
  plot(result_russia[[test_name]], main = sprintf("TB Russia: %s",
                                                 test_name))
}
for (i in 1:4) {
  test_name <- colnames(fit$coefficients)[i]
  plot(result_gambia[[test_name]], main = sprintf("TB Gambia: %s",
                                                 test_name))
}
for (i in 1:4) {
  test_name <- colnames(fit$coefficients)[i]
  plot(result_ghana[[test_name]], main = sprintf("TB Ghana: %s",
                                                  test_name))
}
for (i in 1:4) {
  test_name <- colnames(fit$coefficients)[i]
  plot(result_sobota[[test_name]], main = sprintf("TB Uganda/Tanzania: %s",
                                                 test_name))
}
for (i in 1:4) {
  test_name <- colnames(fit$coefficients)[i]
  plot(result_height[[test_name]], main = sprintf("Height: %s",
                                                  test_name))
}
dev.off()

# Final Figures ----------------------------------------------------------------

add_fig_label <- function(label, line = 3) {
  # Adds label to top left of plot
  mtext(text = label, side = 2, line = line, font = 2, las = 1,
        at = par("usr")[4] + par("usr")[4]*.1, cex = 1.5)
}

test_name <- "status_ni"
plot_titles <- c("Difference between susceptible and resistant individuals in the non-infected state",
                 "Difference between susceptible and resistant individuals in the infected state",
                 "Effect of treatment in resistant individuals",
                 "Effect of treatment in susceptible individuals")
plot_titles <- str_wrap(plot_titles, width = 30)
gwas_labs <- c("TB Russia", "TB Gambia", "TB Ghana", "TB Uganda",
               "Height\nEurope")
for (device in c("png", "pdf")) {
  if (device == "png") {
    png(file.path(fig_dir, "gwas-final.png"),
        width = 14, height = 14, units = "in", res = 300)
  } else if (device == "pdf") {
    pdf(file.path(fig_dir, "gwas-final.pdf"),
        width = 14, height = 14, useDingbats = FALSE)
  }

m1 <- rbind(c(0.0, 0.5, 0.5, 1.0),
            c(0.0, 0.5, 0.0, 0.5),
            c(0.5, 1.0, 0.0, 1.0))
m2 <- rbind(c(0.0, 0.5, 0.5, 1.0),
            c(0.5, 1.0, 0.5, 1.0),
            c(0.0, 0.5, 0.0, 0.5),
            c(0.5, 1.0, 0.0, 0.5))
split.screen(m1)
screen(1)
# par(cex.axis = 0.5, cex.main = 0.25)
par(font.main = 1, cex.main = 1, oma = c(0, 0, 0, 0))
par(mar = c(5.1, 4.1, 4.1, 3.1))
title_russia <- "Genes with increasing effect size between susceptible and resistant individuals in the non-infected state are enriched for nearby SNPs with marginally significant P from a GWAS of TB susceptibility in Russia"
plot(result_russia[[test_name]], main = str_wrap(title_russia, width = 80),
     las = 1, ymin = 0, ymax = 1.7)
# mtext(text="a", side=2, line = 3, font = 2, las = 1,
#       at = par("usr")[4] + par("usr")[4]*.1, cex = 1.5)
add_fig_label("a")
screen(2)
par(mar = c(5.1, 4.1, 4.1, 3.1))
title_height <- "Genes with increasing effect size between susceptible and resistant individuals in the non-infected state are NOT enriched for nearby SNPs with marginally significant P from a GWAS of height in Europe"
plot(result_height[[test_name]], main = str_wrap(title_height, width = 80),
     las = 1, ymin = 0, ymax = 1.7)
add_fig_label("b")
split.screen(m2, screen = 3)
screen(4)
par(mar = c(2, 4, 5, 0))
boxplot_enrich(auc = c(result_russia[[1]]$auc$main,
                       result_gambia[[1]]$auc$main,
                       result_ghana[[1]]$auc$main,
                       result_sobota[[1]]$auc$main,
                       result_height[[1]]$auc$main),
               permutations = cbind(result_russia[[1]]$auc$permutations,
                                    result_gambia[[1]]$auc$permutations,
                                    result_ghana[[1]]$auc$permutations,
                                    result_sobota[[1]]$auc$permutations,
                                    result_height[[1]]$auc$permutations),
               ymin = -30, ymax = 55,
               ylab = "Area under the curve (backgroud subtracted)",
               xaxt =  "n", las = 1)
title(plot_titles[1], line = 1)
add_fig_label("c", line = 0.5)
screen(5)
par(mar = c(2, 1, 5, 3))
boxplot_enrich(auc = c(result_russia[[2]]$auc$main,
                       result_gambia[[2]]$auc$main,
                       result_ghana[[2]]$auc$main,
                       result_sobota[[2]]$auc$main,
                       result_height[[2]]$auc$main),
               permutations = cbind(result_russia[[2]]$auc$permutations,
                                    result_gambia[[2]]$auc$permutations,
                                    result_ghana[[2]]$auc$permutations,
                                    result_sobota[[2]]$auc$permutations,
                                    result_height[[2]]$auc$permutations),
               ymin = -30, ymax = 55,
               xaxt =  "n", yaxt = "n")
add_fig_label("d", line = 0.5)
title(plot_titles[2], line = 1)
screen(6)
par(mar = c(5, 4, 2, 0))
boxplot_enrich(auc = c(result_russia[[3]]$auc$main,
                       result_gambia[[3]]$auc$main,
                       result_ghana[[3]]$auc$main,
                       result_sobota[[3]]$auc$main,
                       result_height[[3]]$auc$main),
               permutations = cbind(result_russia[[3]]$auc$permutations,
                                    result_gambia[[3]]$auc$permutations,
                                    result_ghana[[3]]$auc$permutations,
                                    result_sobota[[3]]$auc$permutations,
                                    result_height[[3]]$auc$permutations),
               ylab = "Area under the curve (backgroud subtracted)",
               xaxt =  "n", las = 1, ymin = -30, ymax = 55)
add_fig_label("e", line = 0.5)
title(plot_titles[3])
# Axis rotation from:
# https://www.r-bloggers.com/rotated-axis-labels-in-r-plots/
# https://cran.r-project.org/doc/FAQ/R-FAQ.html#How-can-I-create-rotated-axis-labels_003f
axis(1, at = 1:5, labels = FALSE)
# cat(sprintf("usr:\nx-axis %.2f-%.2f\ny-axis %.2f-%.2f\n",
#             par()$usr[1], par()$usr[2], par()$usr[3], par()$usr[4]))
anno_position <- par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3])
# print(sprintf("annotation at y = %.2f", anno_position))
text(x = 1:5, y = anno_position,
     labels = gwas_labs, srt = 45, adj = 1, xpd = TRUE)
screen(7)
par(mar = c(5, 1, 2, 3))
boxplot_enrich(auc = c(result_russia[[4]]$auc$main,
                       result_gambia[[4]]$auc$main,
                       result_ghana[[4]]$auc$main,
                       result_sobota[[4]]$auc$main,
                       result_height[[4]]$auc$main),
               permutations = cbind(result_russia[[4]]$auc$permutations,
                                    result_gambia[[4]]$auc$permutations,
                                    result_ghana[[4]]$auc$permutations,
                                    result_sobota[[4]]$auc$permutations,
                                    result_height[[4]]$auc$permutations),
               xaxt =  "n", yaxt = "n", ymin = -30, ymax = 55)
add_fig_label("f", line = 0.5)
title(plot_titles[4])
axis(1, at = 1:5, labels = FALSE)
text(x = 1:5, y = par()$usr[3] - 0.05 * (par()$usr[4] - par()$usr[3]),
     labels = gwas_labs, srt = 45, adj = 1, xpd = TRUE)
close.screen(all.screens = TRUE)
invisible(dev.off())
}

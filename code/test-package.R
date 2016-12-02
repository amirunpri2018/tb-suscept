pkg_relative_path <- "../pkg/tbsuscept/"
devtools::document(pkg_relative_path)
suppressPackageStartupMessages(library("data.table"))

if(interactive()) {
  data_dir <- "../data"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
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

gene_names <- rownames(fit$coefficients)

tss_all_fname <- file.path(data_dir, "tss-all.rds")

snp_coords_fname_ghana <- file.path(data_dir, "/snp-coords-ghana.rds")

snp_coords_fname_gambia <- file.path(data_dir, "/snp-coords-gambia.rds")

snp_coords_fname_height <- file.path(data_dir, "/snp-coords-height.rds")

snp_coords_fname_russia <- file.path(data_dir, "/snp-coords-russia.rds")

# Ghana ------------------------------------------------------------------------

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

# The Gambia -------------------------------------------------------------------

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

# Height -----------------------------------------------------------------------

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

# Russia -----------------------------------------------------------------------

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

# Figures ----------------------------------------------------------------------

boxplot_enrich(auc = c(result_russia[[1]]$auc$main,
                       result_gambia[[1]]$auc$main,
                       result_ghana[[1]]$auc$main,
                       result_height[[1]]$auc$main),
               permutations = cbind(result_russia[[1]]$auc$permutations,
                                    result_gambia[[1]]$auc$permutations,
                                    result_ghana[[1]]$auc$permutations,
                                    result_height[[1]]$auc$permutations))


boxplot_enrich(auc = c(result_russia[[2]]$auc$main,
                       result_gambia[[2]]$auc$main,
                       result_ghana[[2]]$auc$main,
                       result_height[[2]]$auc$main),
               permutations = cbind(result_russia[[2]]$auc$permutations,
                                    result_gambia[[2]]$auc$permutations,
                                    result_ghana[[2]]$auc$permutations,
                                    result_height[[2]]$auc$permutations))

boxplot_enrich(auc = c(result_russia[[1]]$auc$main,
                       result_russia[[2]]$auc$main,
                       result_russia[[3]]$auc$main,
                       result_russia[[4]]$auc$main),
               permutations = cbind(result_russia[[1]]$auc$permutations,
                                    result_russia[[2]]$auc$permutations,
                                    result_russia[[3]]$auc$permutations,
                                    result_russia[[4]]$auc$permutations))

boxplot_enrich(auc = c(result_height[[1]]$auc$main,
                       result_height[[2]]$auc$main,
                       result_height[[3]]$auc$main,
                       result_height[[4]]$auc$main),
               permutations = cbind(result_height[[1]]$auc$permutations,
                                    result_height[[2]]$auc$permutations,
                                    result_height[[3]]$auc$permutations,
                                    result_height[[4]]$auc$permutations))

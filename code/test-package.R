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

result_ghana <- run_gwas_enrich(gene_names = gene_names,
                                tss_all_fname = tss_all_fname,
                                window_size = 50000,
                                rsid = gwas_thye_ghana$SNP,
                                snp_coords_fname = snp_coords_fname_ghana,
                                pval = gwas_thye_ghana$PVAL,
                                effect_size = fit$coefficients[, 1])

plot(result_ghana, main = "Ghana")

# The Gambia -------------------------------------------------------------------

result_gambia <- run_gwas_enrich(gene_names = gene_names,
                                tss_all_fname = tss_all_fname,
                                window_size = 50000,
                                rsid = gwas_thye_gambia$SNP,
                                snp_coords_fname = snp_coords_fname_gambia,
                                pval = gwas_thye_gambia$PVAL,
                                effect_size = fit$coefficients[, 1])

plot(result_gambia, main = "The Gambia")

# Height -----------------------------------------------------------------------

result_height <- run_gwas_enrich(gene_names = gene_names,
                                 tss_all_fname = tss_all_fname,
                                 window_size = 50000,
                                 rsid = gwas_height$MarkerName,
                                 snp_coords_fname = snp_coords_fname_height,
                                 pval = gwas_height$p,
                                 effect_size = fit$coefficients[, 1])

plot(result_height, main = "Height from Lango Allen et al., 2010")

# Russia -----------------------------------------------------------------------

result_russia <- run_gwas_enrich(gene_names = gene_names,
                                 tss_all_fname = tss_all_fname,
                                 window_size = 50000,
                                 rsid = gwas_russia$SNP[grepl("rs", gwas_russia$SNP)],
                                 snp_coords_fname = snp_coords_fname_russia,
                                 pval = gwas_russia$P[grepl("rs", gwas_russia$SNP)],
                                 effect_size = fit$coefficients[, 1])

plot(result_russia, main = "Russia")

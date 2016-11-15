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

gene_names <- rownames(fit$coefficients)

tss_all_fname <- file.path(data_dir, "tss-all.rds")

snp_coords_fname_ghana <- file.path(data_dir, "/snp-coords-ghana.rds")

# Package ----------------------------------------------------------------------

result <- run_gwas_enrich(gene_names = gene_names,
                          tss_all_fname = tss_all_fname,
                          window_size = 50000,
                          rsid = gwas_thye_ghana$SNP,
                          snp_coords_fname = snp_coords_fname_ghana)


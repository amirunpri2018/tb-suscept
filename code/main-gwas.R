#!/usr/bin/env Rscript

# Test for enrichment of GWAS signal in differentially expressed genes.

suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("SNPlocs.Hsapiens.dbSNP144.GRCh38"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("GenomicRanges"))

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

# Assign SNPs to genes ---------------------------------------------------------

# for each gene:
# 1. Obtain the window of the TSS +/- 50 kb
# 2. Identify all the GWAS SNPs in this window
# 3. Assign the SNP with the lowest p-value

gene_names <- rownames(fit$coefficients)

# 1. Obtain the window of the TSS +/- 50 kb

# Ensembl 83, Dec 2015, grch38.p5, hg38
# ensembl <- useMart(host = "dec2015.archive.ensembl.org",
#                    biomart = "ENSEMBL_MART_ENSEMBL",
#                    dataset = "hsapiens_gene_ensembl")
tss_all_fname <- file.path(data_dir, "tss-all.rds")
if (file.exists(tss_all_fname)) {
  tss_all <- readRDS(tss_all_fname)
} else {
  tss_all <- getBM(attributes = c("ensembl_gene_id", "chromosome_name",
                                  "transcription_start_site", "strand"),
                   filters = "ensembl_gene_id",
                   values = gene_names,
                   mart = ensembl)
  saveRDS(tss_all, file = tss_all_fname)
}

# One tss per gene
tss <- tss_all %>%
  group_by(ensembl_gene_id) %>%
  summarize(chr = chromosome_name[1],
            strand = strand[1],
            tss = if (strand == 1) min(transcription_start_site) else max(transcription_start_site))

# Window is +/- 50 kb
window <- 50000
tss$start <- tss$tss - window
tss$end <- tss$tss + window
tss$strand <- ifelse(tss$strand == 1, "+", "-")
tss_gr <- makeGRangesFromDataFrame(tss, keep.extra.columns = TRUE)
seqlevels(tss_gr) <- paste0("ch", seqlevels(tss_gr))

# 3. Assign the SNP with the lowest p-value

# Ghana GWAS
# Obtain snp coordinates
snp_coords_fname_ghana <- file.path(data_dir, "/snp-coords-ghana.rds")
if (file.exists(snp_coords_fname_ghana)) {
  snp_coords_ghana <- readRDS(snp_coords_fname_ghana)
} else {
  snp_coords_ghana <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                               gwas_thye_ghana$SNP,
                               ifnotfound = "drop")
  saveRDS(snp_coords_ghana, file = snp_coords_fname_ghana)
}
stopifnot(mcols(snp_coords_ghana)$RefSNP_id %in% gwas_thye_ghana$SNP)

# Overlap snps and genes
overlaps_ghana <- findOverlaps(snp_coords_ghana, tss_gr, ignore.strand = TRUE)

# Add back names
results_ghana <- data.frame(as.matrix(overlaps_ghana))
colnames(results_ghana) <- c("rsID", "gene")
results_ghana$rsID <- mcols(snp_coords_ghana)$RefSNP_id[results_ghana$rsID]
results_ghana$gene <- mcols(tss_gr)$ensembl_gene_id[results_ghana$gene]

# Add gwas pval
rownames(gwas_thye_ghana) <- gwas_thye_ghana$SNP
results_ghana$gwas_p <- gwas_thye_ghana[results_ghana$rsID, "PVAL"]
stopifnot(!is.na(results_ghana$gwas_p))

# Gambia GWAS
# Obtain snp coordinates
snp_coords_fname_gambia <- file.path(data_dir, "/snp-coords-gambia.rds")
if (file.exists(snp_coords_fname_gambia)) {
  snp_coords_gambia <- readRDS(snp_coords_fname_gambia)
} else {
  snp_coords_gambia <- snpsById(SNPlocs.Hsapiens.dbSNP144.GRCh38,
                               gwas_thye_gambia$SNP,
                               ifnotfound = "drop")
  saveRDS(snp_coords_gambia, file = snp_coords_fname_gambia)
}
stopifnot(mcols(snp_coords_gambia)$RefSNP_id %in% gwas_thye_gambia$SNP)

# Overlap snps and genes
overlaps_gambia <- findOverlaps(snp_coords_gambia, tss_gr, ignore.strand = TRUE)

# Add back names
results_gambia <- data.frame(as.matrix(overlaps_gambia))
colnames(results_gambia) <- c("rsID", "gene")
results_gambia$rsID <- mcols(snp_coords_gambia)$RefSNP_id[results_gambia$rsID]
results_gambia$gene <- mcols(tss_gr)$ensembl_gene_id[results_gambia$gene]

# Add gwas pval
rownames(gwas_thye_gambia) <- gwas_thye_gambia$SNP
results_gambia$gwas_p <- gwas_thye_gambia[results_gambia$rsID, "PVAL"]
stopifnot(!is.na(results_gambia$gwas_p))

results <- merge(results_ghana, results_gambia, by = c("rsID", "gene"),
                 suffixes = c("_ghana", "_gambia"))

# Assign minimum
results <- results %>%
  group_by(gene) %>%
  summarize(gwas_p_ghana = min(gwas_p_ghana),
            gwas_p_gambia = min(gwas_p_gambia),
            n_snps = n())

# Test for relationship with linear regression ---------------------------------

# Add DE stats
limma_coef <- fit$coefficients
results <- merge(results, limma_coef, by.x = "gene", by.y = "row.names")

# Absolute values
for (test in colnames(limma_coef)) {
  results[, test] <- abs(results[, test])
}

# Add mean expression level
mean_expression_level <- fit$Amean[results$gene]
stopifnot(results$gene == names(mean_expression_level))
results$mean_expression_level <- mean_expression_level

# Linear regression
lm_out <- list()
for (test in c(colnames(limma_coef), "mean_expression_level")) {
  lm_out[[test]][["ghana"]] <- lm(results[, "gwas_p_ghana"] ~ results[, test])
  lm_out[[test]][["gambia"]] <- lm(results[, "gwas_p_gambia"] ~ results[, test])
}
# Extract results
extract_lm_stats <- function(lmo) {
  # lmo - lm object
  lmo_sum <- summary(lmo)
  slope <- lmo_sum$coefficients[2, 1]
  slope_se <- lmo_sum$coefficients[2, 2]
  t <- lmo_sum$coefficients[2, 3]
  p <- lmo_sum$coefficients[2, 4]
  return(c(t = t, p = p, slope = slope, slope_se = slope_se))
}

# Step 1: Convert to list of data frames
lm_out_tmp <- lapply(lm_out, function(x) {
  ldply(x, extract_lm_stats, .id = "population")
})
# Step 2: Convert to data frame
lm_out_df <- ldply(lm_out_tmp, .id = "test")

# Save results -----------------------------------------------------------------

write.table(results, file = file.path(data_dir, "results-gwas.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

write.table(lm_out_df, file = file.path(data_dir, "results-gwas-lm.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

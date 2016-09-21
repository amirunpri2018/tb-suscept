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
ensembl <- useMart(host = "dec2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
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

# Linear regression
lm_out <- list()
for (test in c(colnames(limma_coef), "n_snps")) {
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

# Enrichment analysis ----------------------------------------------------------

# Set seed for permutations
set.seed(12345)

enrich <- function(x, y, cutoff, xmin, xmax, breaks = 10,
                   x_direction = "greater", cutoff_direction = "greater") {
  intervals <- seq(xmin, xmax, length.out = breaks)
  enrichment <- numeric(length = breaks)
  sizes <- numeric(length = breaks)
  for (i in seq_along(intervals)) {
    if (x_direction == "greater") {
      y_sub <- y[x > intervals[i]]
    } else {
      y_sub <- y[x < intervals[i]]
    }
    sizes[i] <- length(y_sub)
    # browser()
    if (cutoff_direction == "greater") {
      enrichment[i] <- sum(y_sub > cutoff) / sizes[i]
    } else {
      enrichment[i] <- sum(y_sub < cutoff) / sizes[i]
    }
  }
  return(data.frame(enrichment = enrichment, intervals = intervals, sizes = sizes))
}

enrich_full <- function(x, y, cutoff, xmin, xmax, breaks = 10,
                        x_direction = "greater", cutoff_direction = "greater",
                        iterations = 100) {
  mat_enrichment <- matrix(nrow = breaks, ncol = iterations + 1)
  mat_intervals <- matrix(nrow = breaks, ncol = iterations + 1)
  mat_sizes <- matrix(nrow = breaks, ncol = iterations + 1)
  main <- enrich(x = x, y = y, cutoff = cutoff, xmin = xmin, xmax = xmax,
                 breaks = breaks, x_direction = x_direction,
                 cutoff_direction = cutoff_direction)
  mat_enrichment[, 1] <- main$enrichment
  mat_intervals[, 1] <- main$intervals
  mat_sizes[, 1] <- main$sizes
  # browser()
  for (iter in 1:iterations) {
    permuted <- enrich(x = sample(x), y = y, cutoff = cutoff, xmin = xmin, xmax = xmax,
                       breaks = breaks, x_direction = x_direction,
                       cutoff_direction = cutoff_direction)
    mat_enrichment[, iter + 1] <- permuted$enrichment
    mat_intervals[, iter + 1] <- permuted$intervals
    mat_sizes[, iter + 1] <- permuted$sizes
  }
  return(list(enrichment = mat_enrichment, intervals = mat_intervals,
              sizes = mat_sizes))
}

x <- enrich_full(x = results$gwas_p_gambia,
            y = results$status_ni,
            cutoff = 1,
            xmin = 0, xmax = 1,
            breaks = 25,
            x_direction = "lesser",
            cutoff_direction = "greater")

# http://stackoverflow.com/a/12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall=k)

for (gwas in c("gambia", "ghana")) {
  for (test in colnames(fit$coef)[1:4]) {
    fname_base <- file.path(data_dir, paste(gwas, test, sep = "-"))
    print(fname_base)
    if (gwas == "gambia") {
      enrich_result <- enrich_full(x = results$gwas_p_gambia,
                                   y = results[, test],
                                   cutoff = 1,
                                   xmin = 0, xmax = 1,
                                   breaks = 25,
                                   x_direction = "lesser",
                                   cutoff_direction = "greater")
    } else if (gwas == "ghana") {
      enrich_result <- enrich_full(x = results$gwas_p_ghana,
                                   y = results[, test],
                                   cutoff = 1,
                                   xmin = 0, xmax = 1,
                                   breaks = 25,
                                   x_direction = "lesser",
                                   cutoff_direction = "greater")
    }
    # enrichment values
    write.table(specify_decimal(enrich_result$enrichment, k = 5),
                file = paste0(fname_base, "-enrichment.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    # interval values
    write.table(specify_decimal(enrich_result$intervals, k = 5),
                file = paste0(fname_base, "-intervals.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
    # size values
    write.table(enrich_result$sizes,
                file = paste0(fname_base, "-sizes.txt"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

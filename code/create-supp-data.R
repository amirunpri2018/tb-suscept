#!/usr/bin/env Rscript

# Create Supplementary Data files.

# Need to be < 50 MB
# http://www.nature.com/srep/publish/guidelines#supplementary-info

suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("data.table"))

if(interactive()) {
  data_dir <- "../data"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
}
stopifnot(dir.exists(data_dir))

# S1 - counts ------------------------------------------------------------------

# No need to create a new file. Just symlink to counts.txt. Note: The path to
# the file is relative to where the symlink is created, not the current
# directory.
s1_fname <- file.path(data_dir, "Supplementary_Data_S1.tds")
if (!file.exists(s1_fname)) {
  file.symlink(from = "counts.txt", to = s1_fname)
}

# S2 - limma results -----------------------------------------------------------

results <- readRDS(file.path(data_dir, "results-limma-stats.rds"))
anno_gene <- read.delim(file.path(data_dir, "gene-annotation.txt"),
                        stringsAsFactors = FALSE)

s2 <- createWorkbook()
for (name in names(results)[1:4]) {
  d <- results[[name]]
  d <- d[order(rownames(d)), ]
  stopifnot(rownames(d) == anno_gene$ensembl_gene_id)
  d_final <- data.frame(id = anno_gene$ensembl_gene_id,
                        gene = anno_gene$external_gene_name,
                        d[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "qvalue")],
                        chr = anno_gene$chromosome_name,
                        anno_gene[, c("description", "phenotype", "go_id")],
                        go_description = anno_gene$go_descrip)
  addWorksheet(s2, sheetName = name)
  writeData(s2, sheet = name, d_final)
}
saveWorkbook(s2, file = file.path(data_dir, "Supplementary_Data_S2.xlsx"),
             overwrite = TRUE)

# Queries for paper
results[["treat_resist"]] %>% filter(qvalue < .1, abs(logFC) > 1) %>% nrow()
results[["treat_suscep"]] %>% filter(qvalue < .1, abs(logFC) > 1) %>% nrow()
results[["status_ni"]] %>% filter(qvalue < .1) %>% nrow()
results[["status_ii"]] %>% filter(qvalue < .1) %>% nrow()

# S3 - gwas results ------------------------------------------------------------

gwas <- read.delim(file.path(data_dir, "results-gwas.txt"),
                   stringsAsFactors = FALSE)
gwas <- gwas %>% select(-interact)

# Add gene annotation
gwas_anno <- merge(gwas, anno_gene,
                   by.x = "gene", by.y = "ensembl_gene_id")
stopifnot(nrow(gwas_anno) == nrow(gwas))
gwas_anno <- gwas_anno %>% rename(id = gene,
                                  chr = chromosome_name,
                                  gene = external_gene_name,
                                  go_description = go_descrip)

# For the paper:
gwas_anno %>%
  filter(gwas_p_gambia < .01, gwas_p_ghana < .01, status_ni > 2) %>%
  select(gene, chr, status_ni, gwas_p_gambia, gwas_p_ghana, description)

# http://stackoverflow.com/a/12135122
specify_decimal <- function(x, k) format(round(x, k), nsmall = k)

# Use different cutoffs
de_effect_size <- c(2, 1.5, 1)
pval <- c(.005, .01, .05)
sig_genes <- matrix(NA, nrow = length(de_effect_size) * length(pval), ncol = 4)
colnames(sig_genes) <- c("GWAS P cutoff", "Effect size cutoff",
                         "Number of genes", "Names")
row_counter = 1
for (effect in de_effect_size) {
  for (p in pval) {
    genes <- gwas_anno %>%
      filter(gwas_p_gambia < p,
             gwas_p_ghana < p,
             status_ni > effect) %>%
      select(gene) %>%
      unlist %>%
      sort %>%
      paste(., collapse = ", ")
    sig_genes[row_counter, "GWAS P cutoff"] <- specify_decimal(p, 3)
    sig_genes[row_counter, "Effect size cutoff"] <- specify_decimal(effect, 1)
    sig_genes[row_counter, "Number of genes"] <- length(strsplit(genes, ",")[[1]])
    sig_genes[row_counter, "Names"] <- genes
    row_counter = row_counter + 1
  }
}

s3 <- createWorkbook()
addWorksheet(s3, sheetName = "input-data")
writeData(s3, sheet = "input-data", gwas_anno)
addWorksheet(s3, sheetName = "top-genes")
writeData(s3, sheet = "top-genes", sig_genes)
saveWorkbook(s3, file = file.path(data_dir, "Supplementary_Data_S3.xlsx"),
             overwrite = TRUE)

# S4 - classifier results ------------------------------------------------------------

# Input data for modeling
training <- fread(file.path(data_dir, "training-input.txt"), data.table = FALSE)
testing <- fread(file.path(data_dir, "testing-input.txt"), data.table = FALSE)
combined_results <- readRDS(file.path(data_dir, "combined-limma.rds"))
genes_used_in_final_classifer <- rownames(combined_results[["status_ni"]])[combined_results[["status_ni"]]$qvalue < .05]
training <- training %>% rename(id = V1) %>%
  select(id, status, one_of(genes_used_in_final_classifer))
testing <- testing %>% rename(id = V1) %>%
  select(id, one_of(genes_used_in_final_classifer))

# Results
predictions <- readRDS(file.path(data_dir, "classifier-predictions.rds"))
predictions_final <- predictions[["svmLinear"]][["q0.05"]]$pred %>%
  select(id, obs, resist) %>%
  rename(status = obs,
         prob_tb_resist = resist)

predictions_lbb <- readRDS(file.path(data_dir, "classifier-predictions-lbb.rds"))
predictions_lbb_final <- predictions_lbb[["svmLinear"]][["q0.05"]] %>%
  mutate(id = rownames(.)) %>%
  select(id, resist) %>%
  rename(prob_tb_resist = resist)

# Gene information
class_anno <- anno_gene %>% filter(ensembl_gene_id %in% genes_used_in_final_classifer) %>%
  rename(id = ensembl_gene_id,
         gene = external_gene_name,
         chr = chromosome_name,
         go_description = go_descrip)

s4 <- createWorkbook()
addWorksheet(s4, sheetName = "gene-list")
writeData(s4, sheet = "gene-list", class_anno)
addWorksheet(s4, sheetName = "training-data")
writeData(s4, sheet = "training-data", training)
addWorksheet(s4, sheetName = "training-results")
writeData(s4, sheet = "training-results", predictions_final)
addWorksheet(s4, sheetName = "testing-data")
writeData(s4, sheet = "testing-data", testing)
addWorksheet(s4, sheetName = "testing-results")
writeData(s4, sheet = "testing-results", predictions_lbb_final)
saveWorkbook(s4, file = file.path(data_dir, "Supplementary_Data_S4.xlsx"),
             overwrite = TRUE)

#!/usr/bin/env Rscript

# Create Supplementary Data files.

# Need to be < 50 MB
# http://www.nature.com/srep/publish/guidelines#supplementary-info

suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("dplyr"))

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

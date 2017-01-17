#!/usr/bin/env Rscript

# Create Supplementary Data S4.

suppressPackageStartupMessages(library("openxlsx"))
suppressPackageStartupMessages(library("dplyr"))

if(interactive()) {
  data_dir <- "../data"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  data_dir <- args[1]
}
stopifnot(dir.exists(data_dir))

# Input data -------------------------------------------------------------------

# Differential expression effect sizes
fit <- readRDS(file.path(data_dir, "results-limma-fit.rds"))

# GWAS p-values
result_russia <- readRDS(file.path(data_dir, "gwas-russia.rds"))
result_gambia <- readRDS(file.path(data_dir, "gwas-gambia.rds"))
result_ghana <- readRDS(file.path(data_dir, "gwas-ghana.rds"))
result_sobota <- readRDS(file.path(data_dir, "gwas-sobota.rds"))
result_height <- readRDS(file.path(data_dir, "gwas-height.rds"))

# Gene annotation information
anno_gene <- read.delim(file.path(data_dir, "gene-annotation.txt"),
                        stringsAsFactors = FALSE)

# Extract and combine data sources ---------------------------------------------

# Differential expression effect sizes
df_de <- as.data.frame(fit$coefficients) %>% select(-interact)
n_genes <- nrow(df_de)
df_de$id <- row.names(df_de)

# GWAS p-values
# 
# Every GWAS assigned a SNP to a different number of genes. Merge them 
# one-by-one. For the GWAS lists, the y (referring to the y-axis in the 
# enrichment plots) is the GWAS SNP p-values and is  always the same across the
# differential expression contrasts.
df_de_p <- merge(df_de, as.data.frame(result_russia$status_ni$input$y),
                 by.x = "id", by.y = "row.names", all = TRUE)
df_de_p <- merge(df_de_p, as.data.frame(result_gambia$status_ni$input$y),
                 by.x = "id", by.y = "row.names", all = TRUE)
df_de_p <- merge(df_de_p, as.data.frame(result_ghana$status_ni$input$y),
                 by.x = "id", by.y = "row.names", all = TRUE)
df_de_p <- merge(df_de_p, as.data.frame(result_sobota$status_ni$input$y),
                 by.x = "id", by.y = "row.names", all = TRUE)
df_de_p <- merge(df_de_p, as.data.frame(result_height$status_ni$input$y),
                 by.x = "id", by.y = "row.names", all = TRUE)
stopifnot(nrow(df_de_p) == n_genes)
colnames(df_de_p)[(ncol(df_de_p) - 4):ncol(df_de_p)] <- c("gwas_p_russia",
                                                          "gwas_p_gambia",
                                                          "gwas_p_ghana",
                                                          "gwas_p_uganda",
                                                          "gwas_p_height")

# Gene annotation information
stopifnot(anno_gene$ensembl_gene_id == df_de_p$id)
df_full <- cbind(df_de_p, anno_gene)
df_full <- df_full %>% select(-ensembl_gene_id) %>%
  rename(chr = chromosome_name,
         gene = external_gene_name,
         go_description = go_descrip) %>%
  select(id, gwas_p_russia:gwas_p_height, status_ni:treat_suscep,
         gene:go_description)

df_top <- df_full %>% filter(abs(status_ni) >= 2) %>%
  arrange(gene)
stopifnot(nrow(df_top) == 17)

# Numbers ----------------------------------------------------------------------

# How many genes were assigned a GWAS SNP for each study?
df_full %>% select(starts_with("gwas_p")) %>%
  sapply(., function(x) sum(!is.na(x))) %>%
  sort
## gwas_p_gambia  gwas_p_ghana gwas_p_russia gwas_p_height gwas_p_uganda 
## 10265         10268         10827         10893         11060 

# Save data --------------------------------------------------------------------

s4 <- createWorkbook()
addWorksheet(s4, sheetName = "input-data")
writeData(s4, sheet = "input-data", df_full, keepNA = TRUE)
addWorksheet(s4, sheetName = "top-genes")
writeData(s4, sheet = "top-genes", df_top, keepNA = TRUE)
saveWorkbook(s4, file = file.path(data_dir, "Supplementary_Data_S4.xlsx"),
             overwrite = TRUE)

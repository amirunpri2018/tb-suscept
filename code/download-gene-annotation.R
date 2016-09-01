#!/usr/bin/env Rscript

# Add gene annotation for interpreting DE results.

suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))

# Input ------------------------------------------------------------------------

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

counts <- read.delim(file.path(data_dir, "counts.txt"),
                     check.names = FALSE, row.names = 1)
genes <- rownames(counts)

# Download from Ensembl  -------------------------------------------------------

# Download human exons (Ensembl 83, Dec 2015, grch38.p5, hg38)

ensembl <- useMart(host = "dec2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
# attributePages(ensembl)
# [1] "feature_page" "structure" "homologs" "sequences" "snp" "snp_somatic"
atts <- listAttributes(ensembl, page = "feature_page")
# atts[grep("escription", atts$description, ignore.case = TRUE), ]

# For speed purposes, perform separate queries. Otherwise the query multiplies
# with each new column that maps many to one.

# Obtain annotation results that are one per gene
basic <- getBM(attributes = c("ensembl_gene_id",
                            "external_gene_name",
                            "chromosome_name",
                            "description"),
             filters = "ensembl_gene_id",
             values = genes,
             mart = ensembl)
stopifnot(sort(genes) == sort(basic$ensembl_gene_id))

# Obtain phenotype descriptions
pheno <- getBM(attributes = c("ensembl_gene_id",
                              "phenotype_description"),
               filters = "ensembl_gene_id",
               values = genes,
               mart = ensembl)
stopifnot(genes %in% unique(pheno$ensembl_gene_id))

# Obtain GO categories
go <- getBM(attributes = c("ensembl_gene_id",
                           "go_id",
                           "name_1006"),
            filters = "ensembl_gene_id",
            values = genes,
            mart = ensembl)
# Not all genes have GO annotations
# stopifnot(genes %in% unique(go$ensembl_gene_id))

# Combine annotation per gene---------------------------------------------------

basic_gene <- basic %>% arrange(ensembl_gene_id)
pheno_gene <- pheno %>% group_by(ensembl_gene_id) %>%
  summarize(phenotype = paste(phenotype_description, collapse = ";")) %>%
  arrange(ensembl_gene_id)
stopifnot(basic_gene$ensembl_gene_id == pheno_gene$ensembl_gene_id)
basic_pheno <- cbind(basic_gene, pheno_gene[, "phenotype"])

go_gene <- go %>% group_by(ensembl_gene_id) %>%
  summarize(go_id = paste(go_id, collapse = ";"),
            go_descrip = paste(name_1006, collapse = ";")) %>%
  arrange(ensembl_gene_id)

combined <- merge(basic_pheno, go_gene, by = "ensembl_gene_id",
                  sort = TRUE, all.x = TRUE)
stopifnot(sort(genes) == combined$ensembl_gene_id)

# Save annotation --------------------------------------------------------------

write.table(combined, file.path(data_dir, "gene-annotation.txt"),
            quote = FALSE, sep = "\t", row.names = FALSE)

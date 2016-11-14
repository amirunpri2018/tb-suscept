#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("Biobase"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))

if(interactive()) {
  data_dir <- "../data"
} else {
  data_dir <- commandArgs(trailingOnly = TRUE)
}
stopifnot(dir.exists(data_dir))

# Input ------------------------------------------------------------------------

geo_fname <- file.path(data_dir, "../data/GSE11199.rds")
gse <- readRDS(geo_fname)
stopifnot(class(gse) == "ExpressionSet")
x <- exprs(gse)

# Prepare sample annotation ----------------------------------------------------

anno <- pData(gse)[, 1]
anno <- as.data.frame(anno)
anno <- separate(anno, col = anno, into = c("individual", "treatment"))
anno$treatment <- factor(anno$treatment, levels = c("unstim", "stim"),
                         labels = c("noninf", "infect"))
anno$status <- substr(anno$individual, 1, 1)
anno$tb <- factor(anno$status, levels = c("L", "P", "T"),
                  labels = c("latent", "pulmonary", "meningeal"))
anno$status <- ifelse(anno$status == "L", "resist", "suscep")
rownames(anno) <- paste(anno$individual, anno$status, anno$treatment, sep = "-")
colnames(x) <- rownames(anno)
# table(anno$treatment, anno$status)

# Map probes to genes ----------------------------------------------------------

# Download probe annotation.

ensembl <- useMart(host = "dec2015.archive.ensembl.org",
                   biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl")
probe2id <- getBM(attributes = c("affy_hg_u133_plus_2", "ensembl_gene_id",
                                 "chromosome_name", "status",
                                 "external_gene_name", "gene_biotype"),
                  filters = "affy_hg_u133_plus_2",
                  values = as.character(fData(gse)$ID),
                  mart = ensembl)
head(probe2id)

# Filtering by chromosome and biotype.

probes <- probe2id[probe2id$chromosome_name %in% c(1:22, "X", "Y", "MT") &
                     probe2id$gene_biotype == "protein_coding", ]

# Keeping only those probes that unambiguously tag one protein-coding gene.

genes_per_probe <- table(probes$affy_hg_u133_plus_2)
ambiguous_probes <- names(genes_per_probe)[genes_per_probe > 1]
good_probes <- names(genes_per_probe)[genes_per_probe == 1]
probes_1pergene <- probes[probes$affy_hg_u133_plus_2 %in% good_probes, ]
stopifnot(length(probes_1pergene$affy_hg_u133_plus_2) ==
            length(unique(probes_1pergene$affy_hg_u133_plus_2)))
probes_per_gene <- table(probes_1pergene$ensembl_gene_id)

# Many genes are tagged by more than one probe. I average them to get one
# value per gene.

x_filtered <- x[probes_1pergene$affy_hg_u133_plus_2, ]
x_filtered <- as.data.frame(x_filtered)
x_filtered$ensembl <- probes_1pergene$ensembl_gene_id
x_mean <- x_filtered %>%
  group_by(ensembl) %>%
  summarise_each(funs(mean)) %>%
  as.data.frame
stopifnot(length(x_mean$ensembl) == length(unique(x_mean$ensembl)))
rownames(x_mean) <- x_mean$ensembl
x_mean <- x_mean %>% select(-ensembl)

# Remove outlier ---------------------------------------------------------------

x_final <- x_mean[, colnames(x_mean) != "LTB2-resist-noninf"]
anno_final <- anno[rownames(anno) != "LTB2-resist-noninf", ]
stopifnot(colnames(x_final) == rownames(anno_final))

# Save data --------------------------------------------------------------------

write.table(x_final, file = file.path(data_dir, "thuong2008-exp.txt"),
            quote = FALSE, sep = "\t", col.names = NA)
write.table(anno_final, file = file.path(data_dir, "thuong2008-anno.txt"),
            quote = FALSE, sep = "\t", col.names = NA)

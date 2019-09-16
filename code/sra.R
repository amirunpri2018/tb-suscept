#!/usr/bin/env Rscript

# Connect samples to the SRA accesssions for downloading raw FASTQ data.
#
# Creates data/sra.txt, which contains accessions for GEO, SRX, and SRR
#
# Creates code/download-fastq-sra.sh, a bash script to download the 50 FASTQ files

library(glue)
library(dplyr)
library(GEOquery)
library(SRAdb)

# Obtain GEO accessions and links to SRA accessions (SRX) ----------------------
gse <- getGEO("GSE94116")
gse <- as.data.frame(gse[[1]], stringsAsFactors = FALSE)
sra_info <- gse %>%
  mutate(srx = substr(relation.1, start = 44, stop = 54)) %>%
  select(sample = title, geo = geo_accession, srx)
head(sra_info)

# Download SRA database --------------------------------------------------------
sqlfile <- "/tmp/SRAmetadb.sqlite"
if(!file.exists(sqlfile)) {
  sqlfile <- getSRAdbFile(destdir = "/tmp", method = "wget")
}
sra_con <- dbConnect(SQLite(), sqlfile)

# Obtain SRR -------------------------------------------------------------------
sra_info <- sra_info %>%
  mutate(srr = sraConvert(srx, "run", sra_con)$run)

# Save SRA metadata ------------------------------------------------------------
write.table(sra_info, "data/sra.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# Create download script -------------------------------------------------------
fastq_dump <- glue_data(sra_info,
                        "fastq-dump --stdout --gzip {srr} > {sample}.fastq.gz")
echo <- glue_data(sra_info,
                  "echo \"Downloading {srr}: {sample} ({seq_len(nrow(sra_info))} of {nrow(sra_info)})\"")
script <- paste(echo, fastq_dump, sep = "\n")
script <- c("#!/bin/bash",
            "set -eu",
            "\n# Download GSE94116 FASTQ files from SRA with sra-tools\n",
            script)
cat(script, file = "code/download-fastq-sra.sh", sep = "\n")

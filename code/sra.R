#!/usr/bin/env Rscript

# Connect samples to the SRA accesssions for downloading raw FASTQ data.

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

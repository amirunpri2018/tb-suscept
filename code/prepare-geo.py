#!/usr/bin/env python

# This script prepares the information about the raw and processed
# data files for submission to the Gene Expression Omnibus (GEO). It
# uses template version 2.1.

# Because this is a one time operation, the variables are hard-coded
# in the script. Full paths to the data are specified, so the script
# can be executed from anywhere in the filesystem.

import glob
import hashlib
import os
import re
import string
import subprocess as sp
import sys

sys.stdout.write("Running prepare-geo.py\n\n")

################################################################################
# Setup
################################################################################

# Specify directories (include the forward slash at the end)
# fastq files
raw_dir = "/scratch/midway/jdblischak/tb-suscept-geo/"
# processed data files
processed_dir = "/home/jdblischak/tb-suscept/data/"
# output files
outdir = processed_dir

sys.stdout.write("fastq directory:\t%s\nprocessed directory:\t%s\noutput directory:\t%s\ncurrent directory:\t%s\n\n"%(
                 raw_dir, processed_dir, outdir, os.getcwd()))

# Create output directory 
if not os.path.exists(outdir):
    os.mkdir(outdir)

# Read in names of fastq files
fastq_all_full = glob.glob(raw_dir + "*fastq.gz")
fastq_all = [f.split("/")[-1] for f in fastq_all_full]
fastq_all.sort()

################################################################################
# Samples
################################################################################

sys.stdout.write("Creating geo-samples.txt\n\n")

# For this section, each row is a sample. The columns include meta
# data about the sample and also all the raw and processed data files
# that correspond to that sample.

samples = open(outdir + "geo-samples.txt", "w")

# Static columns
source_name = "monocyte-derived dendritic cells"
organism = "Homo sapiens"
molecule = "polyA RNA"
description = ""
processed_data_1 = "Supplementary_Data_S1.tds"
processed_data_2 = "Supplementary_Data_S2.tds"

for fq in fastq_all:
    id = fq.replace(".fastq.gz", "")
    title = id
    name_parts = id.split("-")
    individual = name_parts[0]
    status = "susceptible" if name_parts[1] == "suscep" else "resistant"
    treatment = "infected" if name_parts[2] == "infect" else "non-infected"
    raw_file = fq
    samples.write("\t".join([id, title, source_name, organism,
                             individual, status, treatment,
                             molecule, description, processed_data_1,
                             processed_data_2, raw_file]) + "\n")

samples.close()

################################################################################
# Processed data files
################################################################################

sys.stdout.write("Creating geo-processed-data-files.txt\n\n")

# For this section, each row is a processed data file. The columns are
# the filename, the file type, and the md5sum.

processed = open(outdir + "geo-processed-data-files.txt", "w")

# static column
file_type = "tab-delimited text"

processed_files = [processed_dir + "Supplementary_Data_S1.tds",
                   processed_dir + "Supplementary_Data_S2.tds"]

for p in processed_files:
    p_md5 = hashlib.md5(open(p, 'r').read().encode()).hexdigest()
    fname = p.split("/")[-1]
    processed.write("\t".join([fname, file_type, p_md5]) + "\n")

processed.close()

################################################################################
# Raw files
################################################################################

sys.stdout.write("Creating geo-raw-files.txt\n\n")

# For this section, each row is a raw data file. The columns are the
# filename, information about the sequencing, and the
# md5sum.

raw = open(outdir + "geo-raw-files.txt", "w")

# static columns
file_type = "fastq"
instrument_model = "HiSeq 4000"
read_length = "50"
single_or_paired_end = "single"

for f in fastq_all:
    sys.stdout.write("Calculating md5sum for %s\n"%(f))
    md5result = sp.run("md5sum %s"%(raw_dir + f),
                       shell = True, stdout = sp.PIPE)
    md5sum = md5result.stdout.decode().split()[0]
    raw.write("\t".join([
        f,
        file_type,
        md5sum,
        instrument_model,
        read_length,
        single_or_paired_end]) + "\n")

raw.close()

#!/usr/bin/env python3

# Move and rename fastq files downloaded from FGF's FTP site.
#
# Usage:
#
# python3 arrange-fastq.py indir outdir
#
# indir - Highly nested directory structure from downloading data
# outdir - New output directory (created if does not exist)
#
# Ex:
#
# python3 arrange-fastq.py /project/gilad/jdblischak/tb-suscept/fgfftp.uchicago.edu/Genomics_Data/NGS/160511_K00242_0067_AHCMF5BBXX-SR-YG-MM1-3/FastQ /project/gilad/jdblischak/tb-suscept/fastq
#
# Some examples:
#
# orig: fgfftp.uchicago.edu/Genomics_Data/NGS/160525_K00242_0071_AHCN77BBXX-YG-MM1-3-(2)/FastQ/YG-MM1-02-contact-infected_S6_L002_R1_001.fastq.gz
# new:  /project/gilad/jdblischak/tb-suscept/fastq/02-contact-infected-HCN77BBXX-l2.fastq.gz
#
# orig: fgfftp.uchicago.edu/Genomics_Data/NGS/160525_K00242_0071_AHCN77BBXX-YG-MM1-3-(2)/FastQ/YG-MM2-19-contact-none_S39_L003_R1_001.fastq.gz
# new:  /project/gilad/jdblischak/tb-suscept/fastq/19-contact-none-HCN77BBXX-l3.fastq.gz
#

import glob
import os
import shutil
import sys

# Input arguments
args = sys.argv
assert len(args) == 3, "Incorrect number of arguments.\nUsage: python3 arrange-fastq.py indir outdir"
indir = args[1]
outdir = args[2]
assert os.path.exists(indir), "Input directory does not exist: %s"%(indir)
if not os.path.exists(outdir):
    os.mkdir(outdir)
# Add final forward slash if necessary
if indir[-1] != "/":
    indir = indir + "/"
if outdir[-1] != "/":
    outdir = outdir + "/"

# Obtain file names
files = glob.glob(indir + "/*fastq.gz")[:]

# Rename and move files
for f in files:
    path = f.rstrip('fastq.gz').split('/')
    flow_cell = path[-3].split("_")[-1].split("-")[0][1:]
    file_parts = path[-1].split('_')[:-1]
    lane = "l" + file_parts[2][-1]
    if file_parts[0] == "Undetermined":
        sample_name = file_parts[0].lower()
    else:
        sample_name = "-".join(file_parts[0].split("-")[-3:])
    new_name = outdir + sample_name + '-' + flow_cell + "-" + lane + '.fastq.gz'
    sys.stderr.write("Moving:\n%s\n%s\n\n"%(new_name, f))
    shutil.move(f, new_name)

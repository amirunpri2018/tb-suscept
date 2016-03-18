#!/usr/bin/env python

# Create labels for tubes.

# For Python2&3 compatibility
# https://docs.python.org/3/howto/pyporting.html
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

# Input
data_fname = "../data/experiment-info.txt"
assert os.path.exists(data_fname), "File %s does not exist"%(data_fname)
data_handle = open(data_fname, "r")
# skip first line
data_handle.readline()

# Output
infection_fname = "../data/infection-labels.txt"
infection_handle = open(infection_fname, "w")

extraction_fname = "../data/extraction-labels.txt"
extraction_handle = open(extraction_fname, "w")


extraction_dates = ["2016-03-21", "2016-03-22", "2016-03-23",
                    "2016-03-24", "2016-03-28"]


for line in data_handle:
    #print(line)
    cols = line.strip().split("\t")
    id = cols[0]
    status = cols[1]
    treatment = cols[2]
    infection = cols[3]
    arrival = cols[4]
    extraction = cols[5]
    extraction_int = int(extraction[-1]) - 1
    master_mix = cols[6]
    me = cols[7]
    individual = cols[8]
    birthdate = cols[9]

    infection_handle.write("%s\nQIAzol 1 mL\n%s\n%s; %s\n\n"%(
                           id, infection, extraction, individual))

    for i in range(2):
        # Top label
        extraction_handle.write("%s\n%s -> %s\n\n"%(
                                id, extraction, master_mix))
        # Side label
        extraction_handle.write("%s\nRNA ~20 uL #%d\n%s\n%s -> %s\n\n"%(
                                id, i + 1, extraction_dates[extraction_int], extraction, master_mix))

data_handle.close()
infection_handle.close()
extraction_handle.close()


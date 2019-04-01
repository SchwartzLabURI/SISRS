#!/usr/bin/env python3

# This script prepares data for a SISRS run by setting the data up and creating mapping scripts
# Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
# The composite genome is indexed by Bowtie2 and Samtools
# SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder

import os
from os import path
import sys
from glob import glob

#Set cwd to script location
script_dir = sys.path[0]

#Set TrimRead + SISRS directories based off of script folder location
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome')]

sisrs_output_template = """
SCRIPT_DIR/get_alignment.py TWOTAXA SISRS_DIR

SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment.nex MISSING
grep -oe "SISRS_[^/]*" SISRS_DIR/alignment_locs_mMISSING.txt | uniq -c | sort -k1 -nr | awk '{print $2}' > SISRS_DIR/alignment_locs_mMISSING_Clean.txt

SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment_bi.nex MISSING
grep -oe "SISRS_[^/]*" SISRS_DIR/alignment_bi_locs_mMISSING.txt | uniq -c | sort -k1 -nr | awk '{print $2}' > SISRS_DIR/alignment_bi_locs_mMISSING_Clean.txt

SCRIPT_DIR/filter_nexus_for_missing.py SISRS_DIR/alignment_pi.nex MISSING
grep -oe "SISRS_[^/]*" SISRS_DIR/alignment_pi_locs_mMISSING.txt | uniq -c | sort -k1 -nr | awk '{print $2}' > SISRS_DIR/alignment_pi_locs_mMISSING_Clean.txt
"""

keyList = ['TWOTAXA','SISRS_DIR','SCRIPT_DIR','MISSING']
keyDict = {'TWOTAXA':str(len(sisrs_tax_dirs) - 2),
           'SISRS_DIR':sisrs_dir,
           'SCRIPT_DIR':script_dir,
           'MISSING':str(int(sys.argv[1]))}
for key in keyList:
    new_sisrs = new_sisrs.replace(key,keyDict[key])
with open(sisrs_dir+"/Output_Alignment.sh", "w") as text_file:
    print(new_sisrs, file=text_file)
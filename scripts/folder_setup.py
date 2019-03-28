#!/usr/bin/env python3

# This script preps the folder architecture for a SISRS run.
# Taxon ID file (text file with Taxon IDs on new lines) must be in base SISRS folder
# Arguments: (1) Path to Taxon ID file (text file with Taxon IDs on new lines)
# Example: python scripts/folder_setup.py TaxonIDs
# Output: Script will create taxon folders in the RawReads, TrimReads, and SISRS_Run folders

import sys
import os
from os import path

taxon_ID_file = path.abspath(sys.argv[1])
sisrs_dir = path.dirname(taxon_ID_file)
taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

for x in taxa_list:
    os.mkdir(sisrs_dir+"/Reads/RawReads/"+x)
    os.mkdir(sisrs_dir+"/Reads/TrimReads/"+x)
    os.mkdir(sisrs_dir+"/SISRS_Run/"+x)

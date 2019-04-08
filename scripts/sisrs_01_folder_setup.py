#!/usr/bin/env python3

# This script preps the folder architecture for a SISRS run.
# Taxon ID file (text file with Taxon IDs on new lines) must be in base SISRS folder
# Arguments: (1) Path to Taxon ID file (text file with Taxon IDs on new lines)
# Example: python scripts/folder_setup.py TaxonIDs
# Output: Script will create lots of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders

import sys
import os
from os import path

taxon_ID_file = path.abspath(sys.argv[1])
sisrs_dir = path.dirname(taxon_ID_file)
taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

os.mkdir(sisrs_dir+"/Reads")
os.mkdir(sisrs_dir+"/Reads/RawReads")
os.mkdir(sisrs_dir+"/Reads/TrimReads")
os.mkdir(sisrs_dir+"/Reads/SubsetReads")
os.mkdir(sisrs_dir+"/SISRS_Run")
os.mkdir(sisrs_dir+"/Reference_Genome")
os.mkdir(sisrs_dir+"/Reference_Genome/Annotations")
os.mkdir(sisrs_dir+"/Post_SISRS_Processing")
os.mkdir(sisrs_dir+"/R_Analyses")
os.mkdir(sisrs_dir+"/Reference_Topology")

for x in taxa_list:
    os.mkdir(sisrs_dir+"/Reads/RawReads/"+x)
    os.mkdir(sisrs_dir+"/Reads/TrimReads/"+x)
    os.mkdir(sisrs_dir+"/SISRS_Run/"+x)

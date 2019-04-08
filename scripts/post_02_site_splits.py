import os
from os import path
import sys
from glob import glob
import subprocess
import stat
from itertools import islice

script_dir = sys.path[0]
processors = str(int(sys.argv[1]))
ref_species = sys.argv[2]
site_id = 'SISRS_Biallelic_NoMissing'

#Set directories based off of script folder location
post_processing_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing"
ref_topology_dir = path.dirname(path.abspath(script_dir))+"/Reference_Topology"
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"

if(not path.isdir(post_processing_dir+"/Site_Splits")):
    os.mkdir(post_processing_dir+"/Site_Splits")
site_split_dir = post_processing_dir+"/Site_Splits"

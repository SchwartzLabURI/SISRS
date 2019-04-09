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

if(not path.isdir(post_processing_dir+"/Annotations")):
    os.mkdir(post_processing_dir+"/Annotations")
annotation_dir = post_processing_dir+"/Annotations"

if(not path.isdir(annotation_dir+"/Composite")):
    os.mkdir(annotation_dir+"/Composite")
composite_annotation_dir = annotation_dir+"/Composite"

if(not path.isdir(annotation_dir+"/SISRS")):
    os.mkdir(annotation_dir+"/SISRS")
sisrs_annotation_dir = annotation_dir+"/SISRS"

site_output_dir = post_processing_dir+"/SISRS_Sites"
good_split_dir = site_output_dir + "/Good_Splits"
bad_split_dir = site_output_dir + "/Bad_Splits"

ref_topology_dir = path.dirname(path.abspath(script_dir))+"/Reference_Topology"
ref_annotation_dir = path.dirname(path.abspath(script_dir))+"/Reference_Genome/Annotations"
ref_annotation_files = glob(ref_annotation_dir+"/*.bed")

sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"

post_log_dir = post_processing_dir+"/logFiles"

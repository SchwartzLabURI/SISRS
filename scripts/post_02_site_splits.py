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
site_output_dir = post_processing_dir+"/SISRS_Sites"
ref_topology_dir = path.dirname(path.abspath(script_dir))+"/Reference_Topology"
sisrs_dir = path.dirname(path.abspath(script_dir))+"/SISRS_Run"
ref_tree = glob(ref_topology_dir+"/*.nwk")[0]
post_log_dir = post_processing_dir+"/logFiles"

if(not path.isdir(post_processing_dir+"/Site_Splits")):
    os.mkdir(post_processing_dir+"/Site_Splits")
site_split_dir = post_processing_dir+"/Site_Splits"

if(not path.isdir(site_split_dir+"/Good_Splits")):
    os.mkdir(site_split_dir+"/Good_Splits")
good_split_dir = site_split_dir+"/Good_Splits"

if(not path.isdir(site_split_dir+"/Bad_Splits")):
    os.mkdir(site_split_dir+"/Bad_Splits")
bad_split_dir = site_split_dir+"/Bad_Splits"

with open(post_log_dir+"/out_04_Site_Splits","w") as outfile:
    split_command = ['python','{}/site_splits.py'.format(script_dir),
            '{}'.format(ref_tree),
            '{siteoutputdir}/{taxa}_{siteid}_NoGaps.phylip-relaxed'.format(siteoutputdir=site_output_dir,taxa=ref_species,siteid=site_id),
            '{siteoutputdir}/{taxa}_{siteid}_NoGaps_LocList.txt'.format(siteoutputdir=site_output_dir,taxa=ref_species,siteid=site_id),
            ref_species+"_"+site_id,
            site_split_dir,
            bad_split_dir,
            good_split_dir]
    subprocess.call(split_command,stdout=outfile)

bad_split_alignment = ['python',
        '{}/alignment_slicer.py'.format(script_dir),
        '{}/alignment_bi_locs_m0.txt'.format(sisrs_dir),
        '{}/Bad_Split_LocList.txt'.format(bad_split_dir),
        '{}/alignment_bi_m0.phylip-relaxed'.format(sisrs_dir),
        ref_species+"_"+site_id+"_BadSplits",
        bad_split_dir]
subprocess.call(bad_split_alignment)

bad_split_define = ['python',
        '{}/bad_splits.py'.format(script_dir),
        '{}'.format(ref_tree),
        '{badsplitdir}/{refspecies}_{siteid}_BadSplits_NoGaps.phylip-relaxed'.format(badsplitdir=bad_split_dir,refspecies=ref_species,siteid=site_id),
        '{badsplitdir}/{refspecies}_{siteid}_BadSplits_NoGaps_LocList.txt'.format(badsplitdir=bad_split_dir,refspecies=ref_species,siteid=site_id),
        ref_species+"_"+site_id+"_BadSplits",
        bad_split_dir]

subprocess.call(bad_split_define)

#Number bad splits
# bad_split_number = [
#         'cut',
#         '-f1',
#         '{badsplitdir}/{refspecies}_{siteid}_BadSplits_Counts.tsv'.format(badsplitdir=bad_split_dir,refspecies=ref_species,siteid=site_id),
#         '|',
#         'nl',
#         '-v',
#         '1',
#         '-n',
#         'rn',
#         '-',
#         '|',
#         'sed',
#         "'s/\.\.\./\t/g'",
#         '-',
#         '>',
#         '{}/Numbered_BadSplits'.format(bad_split_dir)
# subprocess.call(bad_split_number)

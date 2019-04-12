#!/usr/bin/env python3

import os
from os import path
import sys
from glob import glob
import subprocess
import stat
from itertools import islice
import pandas as pd

script_dir = sys.path[0]
processors = str(int(sys.argv[1]))
ref_species = sys.argv[2]

site_id = 'SISRS_Biallelic_NoMissing'

post_processing_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing"
sisrs_site_dir = path.dirname(path.abspath(script_dir))+"/Post_SISRS_Processing/SISRS_Sites"

if(not path.isdir(post_processing_dir+"/Chronos")):
    os.mkdir(post_processing_dir+"/Chronos")
chronos_dir = post_processing_dir+"/Chronos"

all_sisrs=  pd.read_csv(sisrs_site_dir+'/'+ref_species+'_'+site_id+'_NoGaps_LocList.txt',header=None)[0].str.split("/", n = 1, expand = True)
all_sisrs.columns=['Contig','Site']

sisrs_10000 = pd.DataFrame(all_sisrs.groupby(['Contig']).size()).sort_values(0,ascending=False).head(10000).index.tolist()
sisrs_50000 = pd.DataFrame(all_sisrs.groupby(['Contig']).size()).sort_values(0,ascending=False).head(50000).index.tolist()

with open(chronos_dir+'/sisrs_10000_ContigList.txt','w') as outfile:
    for contig in sisrs_10000:
        outfile.write(contig+"\n")

with open(chronos_dir+'/sisrs_50000_ContigList.txt','w') as outfile:
    for contig in sisrs_50000:
        outfile.write(contig+"\n")

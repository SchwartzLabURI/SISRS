#!/usr/bin/env python3

"""

Input: -d/--directory, -m/--missing
"""

import argparse
import os
from os import path
import sys
import glob
import subprocess
from subprocess import Popen
from itertools import islice
from get_alignment import *
from filter_nexus_for_missing import *

def getData(outPath):
    '''
    This function is designed to get the appropriate directories needed to finish the
    sisrs run. 

    Arguments: path to the output directory.

    Returns: 
    string: directory containing Composite Genome
    list: paths to dir for each taxon
    string: sisrs run folder
    '''

    sisrs_dir = outPath+"/SISRS_Run"
    composite_dir = sisrs_dir + '/Composite_Genome'
    with open(outPath +'/TaxonList.txt') as f:
        sisrs_tax_dirs = f.readlines()
    sisrs_tax_dirs = [sisrs_dir+'/'+x.rstrip() for x in sisrs_tax_dirs]
    return composite_dir,sisrs_tax_dirs,sisrs_dir


def count_sites_by_contig(sisrs, ms):
    '''
    get counts of sites per contig (sorted most to least)

    arguments:
    sisrs (string): sisrs run folder
    ms (int): number of species missing

    returns: none
    '''

    files = ["alignment_pi_locs_m", "alignment_bi_locs_m", "alignment_locs_m"]
    gaps = ["", "_nogap"]
    for f in files:
        for g in gaps:
            alignment_count_command = f"grep -oe \"SISRS_[^/]*\" {sisrs}/{f}{ms}{g}.txt | uniq -c | sort -k1 -nr | awk \"{{print $2}}\" > {sisrs}/{f}{ms}{g}_Clean.txt"

            print(alignment_count_command)
            os.system(alignment_count_command)   


if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-m', '--missing',type=int, action='store',default=[0],nargs="*")
    args = my_parser.parse_args()

    sis = args.directory
    ms = args.missing #contains a list of missing

    composite_dir,sisrs_tax_dirs,sisrs = getData(sis)

    alignment=get_phy_sites(sisrs, composite_dir, len(sisrs_tax_dirs) - 2)

    numbi=alignment.numsnps() #prints numbers of snps, biallelic snps, and singletons

    alignment = write_alignment(sisrs + '/alignment.nex', alignment) #write nexus alignments


    for i in ms:
        i = int(i)
        for f in ["/alignment.nex", "/alignment_bi.nex", "/alignment_pi.nex"]:
            filter_nexus(sisrs + f, i, True)
            filter_nexus(sisrs + f, i, False)
        
        #get counts of sites per contig (sorted most to least)
        count_sites_by_contig(sisrs, i)

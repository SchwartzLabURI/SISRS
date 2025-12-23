#!/usr/bin/env python3

import argparse
import os
from pathlib import Path
from collections import Counter
import re
import subprocess
from get_alignment import *
from filter_nexus_for_missing import *
import gc

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
            input_file = Path(sisrs) / f'{f}{ms}{g}.txt'
            output_file = Path(sisrs) / f'{f}{ms}{g}_Clean.txt'
        
            # Process line by line (memory efficient)
            matches = []
            with open(input_file, 'r') as inf:
                for line in inf:
                    matches.extend(re.findall(r'SISRS_[^/]*', line))
        
            # Count and sort
            counts = Counter(matches)
            sorted_items = sorted(counts.items(), key=lambda x: x[1], reverse=True)
        
            # Write output
            with open(output_file, 'w') as outf:
                for item, count in sorted_items:
                    outf.write(f'{count:6d} {item}\n')  # Format like uniq -c

def makelinks(newdir, olddir):
    filestocopy = [f for f in os.listdir(olddir + "/SISRS_Run/") if f.startswith('alignment')]
    for f in filestocopy:
        # Creates the soft link to the files
        os.link(olddir + "/SISRS_Run/" + f,
                newdir + "/SISRS_Run/" + f)

def run7(sis, ms, link):

    if link:
        makelinks(sis, link)
    else:
        composite_dir, sisrs_tax_dirs, sisrs = getData(sis)
        site_dict_labels, species_data = get_phy_sites(sisrs, composite_dir, len(sisrs_tax_dirs) - 2)
        numsnps(site_dict_labels)  # prints numbers of snps, biallelic snps, and singletons
        write_alignment(sisrs + '/alignment.nex', site_dict_labels, species_data)

        del site_dict_labels
        del species_data
        gc.collect()

        for f in ["/alignment.nex", "/alignment_bi.nex", "/alignment_pi.nex"]:
            filter_nexus(sisrs + f, ms)
            filter_nexus(sisrs + f, ms)

        # get counts of sites per contig (sorted most to least)
        for i in ms:
            count_sites_by_contig(sisrs, i)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory', action='store', nargs="?")
    my_parser.add_argument('-m', '--missing', type=int, action='store', default=[0],nargs="*")
    my_parser.add_argument('--link', action='store', nargs="?", default=None)
    args = my_parser.parse_args()

    sis = args.directory
    ms = args.missing #contains a list of missing

    run7(sis, ms, args.link)

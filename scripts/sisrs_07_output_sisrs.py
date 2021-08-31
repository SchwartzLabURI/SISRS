#!/usr/bin/env python3

"""

Input: -m/--missing
"""

import argparse
import os
from os import path
import sys
import glob
import subprocess
from subprocess import Popen
from cmdCheck import *
from itertools import islice
from get_alignment import *
from filter_nexus_for_missing import *
from filter_nexus_for_missing_nogap import *

"""

This function is designed to get the appropriate directories needed to finish the
sisrs run. The parameters that is needed for it to properly run is just the
current sisrs working directory.
"""
def getData(outPath):
    """ This function gets the appropriate directories needed to finish the SISRS run. """

    sisrs_dir = outPath+"/SISRS_Run"
    composite_dir = sisrs_dir + '/Composite_Genome'
    with open(outPath +'/TaxonList.txt') as f:
        sisrs_tax_dirs = f.readlines()
    sisrs_tax_dirs = [sisrs_dir+'/'+x.rstrip() for x in sisrs_tax_dirs]
    return composite_dir,sisrs_tax_dirs,sisrs_dir


def run_alignment_locs_m(sisrs, ms):

    alignment_locs_m_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_locs_m" + str(ms) + ".txt",

    ' | uniq -c ',

    ' | sort -k1 -nr ',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_locs_m" + str(ms) + "_Clean.txt"

    ]

    print(alignment_locs_m_command)
    os.system("".join(alignment_locs_m_command))

def run_alignment_locs_m_nogap(sisrs, ms):

    alignment_locs_m_nogap_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_locs_m" + str(ms) + "_nogap.txt",

    ' | uniq -c',

    ' | sort -k1 -nr',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_locs_m" + str(ms) + "_nogap_Clean.txt"

    ]

    print(alignment_locs_m_nogap_command)
    os.system("".join(alignment_locs_m_nogap_command))

def run_alignment_bi_locs(sisrs, ms):

    alignment_bi_locs_m_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_bi_locs_m" + str(ms) + ".txt",

    ' | uniq -c',

    ' | sort -k1 -nr',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_bi_locs_m" + str(ms) + "_Clean.txt"

    ]

    print(alignment_bi_locs_m_command)
    os.system("".join(alignment_bi_locs_m_command))

def run_alignment_bi_locs_nogap(sisrs, ms):

    alignment_bi_locs_m_nogap_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_bi_locs_m" + str(ms) + "_nogap.txt",

    ' | uniq -c',

    ' | sort -k1 -nr',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_bi_locs_m" + str(ms) + "_nogap_Clean.txt"

    ]

    print(alignment_bi_locs_m_nogap_command)
    os.system("".join(alignment_bi_locs_m_nogap_command))

def run_alignment_pi_locs_m(sisrs, ms):

    alignment_pi_locs_m_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_pi_locs_m" + str(ms) + ".txt",

    ' | uniq -c',

    ' | sort -k1 -nr',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_pi_locs_m" + str(ms) + "_Clean.txt"

    ]
    print(alignment_pi_locs_m_command)
    os.system("".join(alignment_pi_locs_m_command))


def run_alignment_pi_locs_m_nogap(sisrs, ms):

    alignment_pi_locs_m_nogap_command = [

    'grep -oe "SISRS_[^/]*" ',

    sisrs + "/alignment_pi_locs_m" + str(ms) + "_nogap.txt",

    ' | uniq -c',

    ' | sort -k1 -nr',

    ' | awk "{print $2}" > ',

    sisrs + "/alignment_pi_locs_m" + str(ms) + "_nogap_Clean.txt"

    ]
    print(alignment_pi_locs_m_nogap_command)
    os.system("".join(alignment_pi_locs_m_nogap_command))




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

    alignment = write_alignment(sisrs + '/alignment.nex', alignment, numbi)


    for i in ms:
        i = int(i)

        filter_nexus(sisrs + "/alignment.nex", i)

        filter_nexus_no_gap(sisrs + "/alignment.nex", i)

        run_alignment_locs_m(sisrs, i)

        run_alignment_locs_m_nogap(sisrs, i)

        filter_nexus(sisrs + "/alignment_bi.nex", i)

        filter_nexus_no_gap(sisrs + "/alignment_bi.nex", i)

        run_alignment_bi_locs(sisrs, i)

        run_alignment_bi_locs_nogap(sisrs, i)

        filter_nexus(sisrs + "/alignment_pi.nex", i)

        filter_nexus_no_gap(sisrs + "/alignment_pi.nex", i)

        run_alignment_pi_locs_m(sisrs, i)

        run_alignment_pi_locs_m_nogap(sisrs, i)

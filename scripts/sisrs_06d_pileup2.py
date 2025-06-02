#!/usr/bin/env python3

'''
This script does a pileup for the specific species genome contigs then calls a base for each site
'''

import os
from get_pruned_dict import *
import argparse
from sisrs_06b_pileup import pileup

def sindex(outPath,sp):
    '''
    This function runs samtools index command.

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name 

    Returns: none.
    '''

    sin = f'samtools index {outPath}/SISRS_Run/{sp}/{sp}.bam'
    os.system(sin)


def prune(outPath, sp, minread, threshold):
    '''
    This function runs the scripts to output LocList file containing a called base for each site for a species relative the positions in the composite genome

    Arguments: path to the output directory, taxon name directory,
               threshold for the number of reads to call a site,
               threshold [<=1] for the proportion of sites required to be the same to call a site.

    Returns: none.
    '''
    #script_dir = os.path.dirname(os.path.abspath(__file__))
    #print(script_dir)
    f = "".join([outPath, '/SISRS_Run/', sp])
    print(f)
    getallbases_main(f, outPath+'/SISRS_Run/Composite_Genome', minread, threshold)

def run6d(sis, minread, folder, threshold):
    print(sis, folder)

    sindex(sis, folder)
    pileup(sis,folder)
    prune(sis, folder, minread, threshold)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-m','--minread',action='store',default=3,nargs="?")
    my_parser.add_argument('-s', '--species', action='store',nargs="?")
    my_parser.add_argument('-t', '--threshold', action='store',default=1,nargs="?")

    args = my_parser.parse_args()

    sis = args.directory
    minread = args.minread
    folder = args.species#.rstrip('/')
    threshold = args.threshold



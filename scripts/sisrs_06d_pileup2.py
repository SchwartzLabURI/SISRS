#!/usr/bin/env python3

'''

This script runs one sisrs alignment - specific taxon contigs
'''

import os
from os import path
import sys
from glob import glob
from cmdCheck import *
import argparse
import re

def sindex(outPath,sp):
    ''' This function runs samtools index command. '''

    outbam = "".join([outPath, '/SISRS_Run/', sp,
        '/',
        sp,
        '.bam'])
    sin = ['samtools index ', outbam]
    os.system("".join(sin))

'''
specific contigs
'''
def pileup(outPath,sp):
    ''' This function performs samtools mpileup on composite genome. '''

    outbam = "".join([outPath, '/SISRS_Run/', sp, #AotNan
        '/',
        sp,
        '.bam'])

    outpile = "".join([outPath, '/SISRS_Run/', sp, #AotNan
        '/',
        sp,
        '.pileups'])

    #samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups
    pileup = ['samtools mpileup -f ', outPath, '/SISRS_Run/Composite_Genome/contigs.fa ', outbam, ' > ', outpile ]

    os.system("".join(pileup))

def prune(outPath, sp, minread, threshold):
    ''' This function calls get_pruned_dict.py '''

    #put the script here instead of calling a new interpreter

    script_dir = os.path.dirname(os.path.abspath(__file__))
    f = "".join([outPath, '/SISRS_Run/', sp])
    specific = ['python3 ', script_dir, '/get_pruned_dict.py ', f, ' ', outPath,'/SISRS_Run/Composite_Genome ',
    minread, ' ', threshold ]# python SCRIPT_DIR/specific_genome.py SISRS_DIR/TAXA COMPOSITE_GENOME

    os.system("".join(specific))


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

    print(sis, folder)

    sindex(sis, folder)
    pileup(sis,folder)
    prune(sis, folder, minread, threshold)

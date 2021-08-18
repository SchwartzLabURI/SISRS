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

'''
specific contigs
'''
def pileup(outPath,sp):
    ''' This function runs samtools mpileup command. '''

    outbam = "".join([outPath, '/SISRS_Run/', sp, #AotNan
        '/',
        sp,
        '.bam'])

    outpile = "".join([outPath, '/SISRS_Run/', sp, #AotNan
        '/',
        sp,
        '.pileups'])

    pileup = ['samtools mpileup -f ', #samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups
        outPath,'/SISRS_Run/Composite_Genome/contigs.fa ',
        outbam, ' > ', outpile ]

    os.system("".join(pileup))

def specific_genome(outPath, sp):
    ''' This function calls specific_genome.py ''' #Ask Rachel about this function description

    #put the script here instead of calling a new interpreter

    script_dir = os.path.dirname(os.path.abspath(__file__))
    print(script_dir)
    f = "".join([outPath, '/SISRS_Run/', sp ])
    print(f)
    spec = ['python3 ', script_dir, '/specific_genome.py ', f, ' ', outPath,'/SISRS_Run/Composite_Genome/contigs.fa' ]
    # python SCRIPT_DIR/specific_genome.py SISRS_DIR/TAXA COMPOSITE_GENOME
    print(spec)
    os.system("".join(spec))

def faidx(outPath, sp):
    ''' This function indexes the files. '''
    
    f = "".join([outPath, '/SISRS_Run/', sp ])
    print(f)
    faid = ['samtools faidx ', f, '/contigs.fa' ] #samtools faidx SISRS_DIR/TAXA/contigs.fa
    os.system("".join(faid))

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-s', '--species', action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    sp = args.species
#    sp = sp.rstrip("/")
    print(sis, sp)

    pileup(sis,sp)
    specific_genome(sis, sp)
    faidx(sis, sp)

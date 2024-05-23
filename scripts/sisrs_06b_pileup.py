#!/usr/bin/env python3

'''
Pileup from the bowtie2 alignment
Make a genome specific to the species by replacing sites
Index new genome
'''

import os
from os import path
import sys
from glob import glob
from specific_genome import *
import argparse
import re

'''
specific contigs
'''
def pileup(outPath,sp):
    '''
    This function runs samtools mpileup command.

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name (not path) 

    Returns: none.
    '''

    outbam = "".join([outPath, '/SISRS_Run/', sp, #eg AotNan
        '/',
        sp,
        '.bam'])

    outpile = "".join([outPath, '/SISRS_Run/', sp, #eg AotNan
        '/',
        sp,
        '.pileups'])

    pileup = ['samtools mpileup -f ', #equivalent to: samtools mpileup -f COMPOSITE_GENOME SISRS_DIR/TAXA/TAXA.bam > SISRS_DIR/TAXA/TAXA.pileups
        outPath,'/SISRS_Run/Composite_Genome/contigs.fa ',
        outbam, ' > ', outpile ]

    os.system("".join(pileup))

def specific_genome(outPath, sp):
    '''
    Write a contigs.fa file specif to the species

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name (not path)

    Returns: none.
    '''
    #script_dir = os.path.dirname(os.path.abspath(__file__))
    #print(script_dir)
    f = "".join([outPath, '/SISRS_Run/', sp ])
    print(f)
    getbases_main(f, outPath+'/SISRS_Run/Composite_Genome/contigs.fa')

def faidx(outPath, sp):
    '''
    This function indexes the files with samtools faidx.

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name 

    Returns: none.
    '''

    f = "".join([outPath, '/SISRS_Run/', sp ])
    print(f)
    faid = ['samtools faidx ', f, '/contigs.fa' ] #samtools faidx SISRS_DIR/TAXA/contigs.fa
    os.system("".join(faid))

def run6b(sis, sp):
    #    sp = sp.rstrip("/")
    print(sis, sp)

    pileup(sis, sp)
    specific_genome(sis, sp)
    faidx(sis, sp)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-s', '--species', action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    sp = args.species

    run6b(sis, sp)
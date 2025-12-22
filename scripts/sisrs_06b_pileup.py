#!/usr/bin/env python3

'''
Pileup from the bowtie2 alignment
Make a genome specific to the species by replacing sites
Index new genome
'''

import subprocess
from os import path
from specific_genome import *
import argparse

def pileup(outPath,sp):
    '''
    This function runs samtools mpileup command.

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name (not path) 

    Returns: none.
    '''

    outbam = f'{outPath}/SISRS_Run/{sp}/{sp}.bam'
    outpile = f'{outPath}/SISRS_Run/{sp}/{sp}.pileups'
    contigs_fa = f'{outPath}/SISRS_Run/Composite_Genome/contigs.fa'

    pileup_command = [
        'samtools',
        'mpileup',
        '-f',
        contigs_fa,
        outbam]

    with open(outpile, 'w') as f:
        subprocess.run(pileup_command, stdout=f, check=True)

def specific_genome(outPath, sp):
    '''
    Write a contigs.fa file specif to the species

    Arguments: 
    outPath (string): path to the output directory
    sp (string): taxon name (not path)

    Returns: none.
    '''
    
    f = f'{outPath}/SISRS_Run/{sp}'
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

    faid_command = [
        'samtools',
        'faidx',
        f'{outPath}/SISRS_Run/{sp}/contigs.fa']

    subprocess.run(faid_command, check=True)

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

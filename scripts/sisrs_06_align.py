#!/usr/bin/env python3

'''

This script runs one sisrs alignment
'''

import os
from os import path
import sys
from glob import glob
from cmdCheck import *
import pandas as pd
from Bio import SeqIO
import argparse
import re

'''
This function runs bowtie2 on the reads in a folder treating all reads as unpaired
'''
def runBowtie(outPath,threads,readfolder):

    outbam = "".join([outPath, '/SISRS_run/', os.path.basename(readfolder), #AotNan
        '/',
        os.path.basename(readfolder),
        '_Temp.bam'])
    outbamb = "".join([outPath, '/SISRS_run/', os.path.basename(readfolder), #AotNan
        '/',
        os.path.basename(readfolder),
        '.bam'])

    bowtie_command = [
        'bowtie2 -p ',
        str(threads),
        ' -N 1 --local -x ',
        outPath,'/SISRS_Run/Composite_Genome/contigs -U ', #contigs base of filename
        ",".join(glob(readfolder+"*.fastq.gz")), #files in the readfolder
        ' | samtools view -Su -@ ', 
        str(threads),
        ' -F 4 - | samtools sort -@ ',
        str(threads),
        ' - -o ',
        outbam]
    os.system("".join(bowtie_command))

    samtools1 = [
        'samtools view -@ ',
        str(threads),
        ' -H ', outbam,
        ' > ', os.path.splitext(outbam)[0], '_Header.sam' ]
    samtools2 = [ 
        'samtools view -@ ', 
        str(threads),
        ' ', outbam, ' | grep -v "XS:" | cat ', os.path.splitext(outbam)[0], '_Header.sam - | samtools view -@ ',
        str(threads), ' -b - > ', outbamb]
    
    os.system("".join(samtools1))
    os.system("".join(samtools2)) #why is this command necessary?


if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-p','--processors',action='store',default=1,nargs="?")
    my_parser.add_argument('-f', '--folder', ,action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors
    folder = args.folder

    runBowtie(sis,proc,folder)

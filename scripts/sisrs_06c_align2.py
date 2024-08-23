#!/usr/bin/env python3

'''
This script outputs a sorted bam by aligning reads to the contigs.fa for the specific species
'''

import os
from os import path
from glob import glob
import argparse
import psutil

def bbuild(outPath, readfolder, proc):
    '''
    This function runs bowtie2-build on contigs.fa in the specified folder in the SISRS run

    Arguments: 
    outPath (string): path to the output directory
    readfolder (string): taxon name
    proc (int): number of processors to use.

    Returns: none.
    '''

    b = f"bowtie2-build {outPath}/SISRS_Run/{readfolder}/contigs.fa {outPath}/SISRS_Run/{readfolder}/contigs -p {proc}"
    # bowtie2-build SISRS_DIR/TAXA/contigs.fa SISRS_DIR/TAXA/contigs -p PROCESSORS

    os.system(b)

def runBowtie(outPath,threads,sp, ref):
    '''
    This function runs  bowtie2 on all fastq.gz in a folder and outputs a sorted bam.
    The function treats all reads as unpaired.

    Arguments: 
    outPath (string): path to the output directory
    threads (int): number of processors to use
    sp: taxon name (folder name not path)
    ref (string): reference to use for alignment

    Returns: none.
    '''

    outbam = f"{outPath}/SISRS_Run/{sp}/{sp}_Temp.bam"
    outbamb = f"{outPath}/SISRS_Run/{sp}/{sp}.bam"
    print(outbam)

    if path.isdir(outPath + "/Reads/TrimReads_nocl"):
        trim_read_dir = outPath + "/Reads/TrimReads_nocl/"
    else:
        trim_read_dir = outPath + "/Reads/TrimReads/"

    mem = psutil.virtual_memory()[1]/threads

    fastqs = ",".join(glob(os.path.expanduser("".join([trim_read_dir, sp, '/*.fastq.gz']))))
    bowtie_command = f"bowtie2 -p {threads} -N 1 --local -x {ref} -U {fastqs} | samtools view -Su -@ {threads} -F 4 - | samtools sort -@ {threads} -m {mem} - -o {outbam}"
    print(bowtie_command)
    os.system(bowtie_command)

    samtools1 = f"samtools view -@ {threads} -H {outbam} > {outbam}_Header.sam"
    samtools2 = f'samtools view -@ {threads} {outbam} | grep -v "XS:" | cat {outbam}_Header.sam - | samtools view -@ {threads} -b - > {outbamb}'

    print(samtools1)
    print(samtools2)

    os.system(samtools1)
    os.system(samtools2)

    os.remove(outbam)  #rm SISRS_DIR/TAXA/TAXA_Temp.bam
    os.remove(outbam+'_Header.sam') #rm SISRS_DIR/TAXA/TAXA_Header.sam

def run6c(sis, proc, f2):
    bbuild(sis, f2, proc)

    ref = f"{sis}/SISRS_Run/{f2}/contigs"
    runBowtie(sis, proc, f2, ref)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--directory', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', action='store', default=1, nargs="?", type=int)
    my_parser.add_argument('-f', '--folder', action='store', nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors
    f2 = args.folder

    run6c(sis, proc, f2)

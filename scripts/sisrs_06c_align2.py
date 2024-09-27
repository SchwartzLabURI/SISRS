#!/usr/bin/env python3

'''
This script outputs a sorted bam by aligning reads to the contigs.fa for the specific species
'''

import os
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

    outsam = f"{outPath}/SISRS_Run/{sp}/{sp}.sam"
    outbamt = f"{outPath}/SISRS_Run/{sp}/{sp}_temp.bam"
    outbam = f"{outPath}/SISRS_Run/{sp}/{sp}.bam"

    if os.path.isdir(outPath + "/Reads/TrimReads_nocl"):
        trim_read_dir = outPath + "/Reads/TrimReads_nocl/"
    else:
        trim_read_dir = outPath + "/Reads/TrimReads/"

    mem = 0.8*psutil.virtual_memory()[1]/threads

    fastqs = ",".join(glob(os.path.expanduser("".join([trim_read_dir, sp, '/*.fastq.gz']))))
    bowtie_command = f'bowtie2 -p {threads} -N 1 --local -x {ref} -U {fastqs} -S {outsam}'
    samview_command = f'samtools view -@ {threads} -F 4 -e "mapq >= 13" -b -o {outbamt} {outsam}'
    samsort_command = f'samtools sort -@ {threads} -m {mem} -O bam -o {outbam} {outbamt}'
    print(bowtie_command)
    os.system(bowtie_command)
    os.system(samview_command)
    os.system(samsort_command)

    os.remove(outsam)
    os.remove(outbamt)

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

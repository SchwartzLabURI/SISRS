#!/usr/bin/env python3

'''
This script outputs a sorted bam by aligning reads to the contigs.fa for the specific species
'''

import os
from pathlib import Path
from glob import glob
import subprocess
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

    bowtie_command = [
        'bowtie2-build',
        f'{outPath}/SISRS_Run/{readfolder}/contigs.fa',
        f'{outPath}/SISRS_Run/{readfolder}/contigs',
        '-p',
        f'{proc}']

    subprocess.run(bowtie_command, check=True)

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

    # Get fastq files using pathlib
    fastq_dir = Path(trim_read_dir).expanduser() / sp
    fastqs = ",".join([str(f) for f in fastq_dir.glob('*.fastq.gz')])

    # Convert commands to lists
    bowtie_command = [
        'bowtie2',
        '-p', f'{threads}',
        '-N', '1',
        '--local',
        '-x', ref,
        '-U', fastqs,
        '-S', outsam]

    samview_command = [
        'samtools',
        'view',
        '-@', f'{threads}',
        '-F', '4',
        '-e', 'mapq >= 13',
        '-b',
        '-o', outbamt,
        outsam]

    samsort_command = [
        'samtools',
        'sort',
        '-@', f'{threads}',
        '-m', mem,
        '-O', 'bam',
        '-o', outbam,
        outbamt]

    print(bowtie_command)
    subprocess.run(bowtie_command, check=True)
    subprocess.run(samview_command, check=True)
    subprocess.run(samsort_command, check=True)

    Path(outsam).unlink()
    Path(outbamt).unlink()

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

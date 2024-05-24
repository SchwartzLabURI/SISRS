#!/usr/bin/env python3

'''
This script outputs a sorted bam by aligning reads to the contigs.fa for the specific species
'''

import os
from os import path
from glob import glob
import argparse


def bbuild(outPath, readfolder, proc):
    '''
    This function runs bowtie2-build on contigs.fa in the specified folder in the SISRS run

    Arguments: 
    outPath (string): path to the output directory
    readfolder (string): taxon name
    proc (string? int?): number of processors to use.

    Returns: none.
    '''

    f = "".join([outPath, '/SISRS_Run/', readfolder ])
    b = ['bowtie2-build ', f, '/contigs.fa ', f, '/contigs -p ',proc ]
    # bowtie2-build SISRS_DIR/TAXA/contigs.fa SISRS_DIR/TAXA/contigs -p PROCESSORS

    os.system("".join(b))

def runBowtie(outPath,threads,sp):
    '''
    This function runs  bowtie2 on all fastq.gz in a folder and outputs a sorted bam.
    The function treats all reads as unpaired.

    Arguments: 
    outPath (string): path to the output directory
    threads (int? string?): number of processors to use
    sp: taxon name (folder name not path)

    Returns: none.
    '''

    outbam = "".join([outPath, '/SISRS_Run/', sp,
        '/',
        sp,
        '_Temp.bam'])
    outbamb = "".join([outPath, '/SISRS_Run/', sp,
        '/',
        sp,
        '.bam'])
    print(outbam)

    if path.isdir(outPath + "/Reads/TrimReads_nocl"):
        trim_read_dir = outPath + "/Reads/TrimReads_nocl/"
    else:
        trim_read_dir = outPath + "/Reads/TrimReads/"

    bowtie_command = [
        'bowtie2 -p ',
        str(threads),
        ' -N 1 --local -x ',
        outPath,'/SISRS_Run/', sp, '/contigs -U ', #contigs base of filename
        ",".join(glob(os.path.expanduser("".join([trim_read_dir, sp, '/*.fastq.gz'])))), #files in the readfolder
        ' | samtools view -Su -@ ',
        str(threads),
        ' -F 4 - | samtools sort -@ ',
        str(threads),
        ' - -o ',
        outbam]
    print(bowtie_command)
    os.system("".join(bowtie_command))

    samtools1 = [
        'samtools view -@ ',
        str(threads),
        ' -H ', outbam,
        ' > ', outbam, '_Header.sam' ]
    samtools2 = [
        'samtools view -@ ',
        str(threads),
        ' ', outbam, ' | grep -v "XS:" | cat ', outbam, '_Header.sam - | samtools view -@ ',
        str(threads), ' -b - > ', outbamb]

    print(samtools1)
    print(samtools2)

    os.system("".join(samtools1))
    os.system("".join(samtools2)) #why is this command necessary?

    os.remove(outbam)  #rm SISRS_DIR/TAXA/TAXA_Temp.bam
    os.remove(outbam+'_Header.sam') #rm SISRS_DIR/TAXA/TAXA_Header.sam

def run6c(sis, proc, f2):
    bbuild(sis, f2, proc)
    runBowtie(sis,proc,f2)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-p','--processors',action='store',default=1,nargs="?")
    my_parser.add_argument('-f', '--folder', action='store',nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors
    f2 = args.folder

    run6c(sis, proc, f2)

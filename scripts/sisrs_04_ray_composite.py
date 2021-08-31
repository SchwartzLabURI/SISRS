#!/usr/bin/env python3

'''

This script calls a Ray genome assembly on all reads in the SubsetReads directory
Ray requires working MPI via mpirun, even if running on  one node
This script also calls Bowtie2 and Samtools for genome indexing
Ensure Ray and Bowtie2 are compiled  with gzip support
Input:
(1) -p/--processors Number of processors
(2) -dir Path to output directory
Output: Ray assembly will be built in <basedir>/Ray_Composite_Genome (Contigs.fasta)
'''

import os
from os import path
import sys
from glob import glob
from cmdCheck import *
import subprocess
from subprocess import check_call
import argparse

'''
This function is designed to obtain the neccessary file paths and data that is
needed for the Ray command. It only requires the working directory for where all
of the sisrs information is stationed.
'''
def getDirs(sisrs_dir):
    ''' This function obtain the neccessary file paths and data that is needed for the Ray command. '''

    #Set SubsetRead and Genome directories based off of script folder location
    subset_read_dir = sisrs_dir+"/Reads/SubsetReads"
    ray_genome_dir = sisrs_dir+"/Ray_Composite_Genome"

    subset_reads = glob(subset_read_dir+"/*.gz")

    return ray_genome_dir, subset_reads

'''
This function is designed to run Ray on all files that are found in the subset
folders. Requieres the output folder developed, number of threads, and the files.
'''
def runRay(ray_genome_dir,subset_reads,threads):
    ''' This function runs Ray command on all files that are found in the subset folders. '''

    ray_command = [
        'mpirun',
        '-n','{}'.format(str(threads)),
        'Ray',
        '-k',
        '31',
        '{}'.format("-s " + " -s ".join(subset_reads)),
        '-o',
        '{}'.format(ray_genome_dir)]
    os.system(" ".join(ray_command))

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-dir','--directory',action='store',nargs="?")
    my_parser.add_argument('-p','--processors',action='store',default=1,nargs="?")
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors

    # obtain the directories and files needed for Ray
    ray_genome_dir, subset_reads = getDirs(sis)

    # Run the ray command
    runRay(ray_genome_dir,subset_reads,proc)

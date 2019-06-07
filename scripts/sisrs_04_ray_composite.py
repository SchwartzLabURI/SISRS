#!/usr/bin/env python3

# This script calls a Ray genome assembly on all reads in the SubsetReads directory
# Ray requires working MPI via mpirun, even if running on  one node
# This script also calls Bowtie2 and Samtools for genome indexing
# Ensure Ray and Bowtie2 are compiled  with gzip support
# Input: (1) Number of nodes for Ray. (2) Number of processors per node.
# Output: Ray assembly will be built in <basedir>/Ray_Composite_Genome (Contigs.fasta)

import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import check_call


'''
This function is designed to obtain the neccessary file paths and data that is
needed for the Ray command. It only requeires the working directory for were all
of the sisrs information is stationed.
'''
def getDirs(sisrs_dir):
    #Set SubsetRead and Genome directories based off of script folder location
    subset_read_dir = sisrs_dir+"/Reads/SubsetReads"
    ray_genome_dir = sisrs_dir+"/Ray_Composite_Genome"

    subset_reads = glob(subset_read_dir+"/*.gz")

    return ray_genome_dir, subset_reads

'''
This function is designed to run ray on all files that are found in the subset
folders. Requieres the output folder developed, number of threads, and the files.
'''
def runRay(ray_genome_dir,subset_reads,threads):
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

    cmd = sys.argv
    sis = path.dirname(path.abspath(path.dirname(cmd[0])))

    proc = 1
    if '-th' in cmd:
        try:
            proc = int(cmd[1])
        except:
            proc = 1
    else:
        print("SWITCHING TO DEFAULT THREAD OF 1")

    # obtain the directories and files needed for Ray
    ray_genome_dir, subset_reads = getDirs(sis)

    # Run the ray command
    runRay(ray_genome_dir,subset_reads,proc)

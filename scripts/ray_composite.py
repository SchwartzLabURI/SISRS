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


#Set cwd to script location
script_dir = sys.path[0]

node_count = int(sys.argv[1])
processors_per_node = int(sys.argv[2])

mpi_processor = node_count*processors_per_node

#Set SubsetRead and Genome directories based off of script folder location
subset_read_dir = path.dirname(path.abspath(script_dir))+"/Reads/SubsetReads"
ray_genome_dir = path.dirname(path.abspath(script_dir))+"/Ray_Composite_Genome"

subset_reads = glob(subset_read_dir+"/*.gz")

ray_command = [
    'mpirun',
    '-n','{}'.format(str(mpi_processor)),
    'Ray',
    '-k',
    '31',
    '{}'.format("-s " + " -s ".join(subset_reads)),
    '-o',
    '{}'.format(ray_genome_dir)]
os.system(" ".join(ray_command))

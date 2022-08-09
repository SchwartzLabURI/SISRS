#!/usr/bin/env python3

import os
import sys
import subprocess
from subprocess import call
from os import listdir
from os.path import isfile, join
import argparse

def align_contigs(output_path, thr):
    '''
    This function is designed to alighn contigs with MAFFT.

    Arguments: path to the output directory, numper of processors to use.

    Returns: none.
    '''

    path_in = output_path + "/SISRS_Run/contigs_outputs/"
    path_out = output_path + "/SISRS_Run/aligned_contigs/"

    #create if not created
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    # list of contigs
    contig_list = os.listdir(path_in)

    for con in contig_list:
        con_p_in = os.path.join(path_in, con)
        con_p_out = os.path.join(path_out, con)
        mafft_command = f"mafft --auto --thread {thr} {con_p_in} > {con_p_out}"
        os.system(mafft_command)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--dir', action='store',nargs="?")
    my_parser.add_argument('-p', '--proc', action='store',default=1,nargs="?", type = int)

    args = my_parser.parse_args()

    output_path = args.dir
    proc = args.proc

    align_contigs(output_path, proc)

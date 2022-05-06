#!/usr/bin/env python3

import os
import sys
import subprocess
from subprocess import call
from cmdCheck import *
from os import listdir
from os.path import isfile, join
import argparse

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--dir', action='store',nargs="?")
    my_parser.add_argument('-p', '--proc', action='store',default=1,nargs="?", type = int)

    args = my_parser.parse_args()

    output_path = args.dir
    proc = args.proc

    #first mask the ref sequence
    os.chmod('align_contigs.sh', 0o755)
    path_in = output_path + "SISRS_Run/contigs_outputs/ "

    path_out = output_path + "SISRS_Run/" + "aligned_contigs"
    #create if not created
    if not os.path.exists(path_out):
        os.mkdir(path_out)

    subprocess.call("./align_contigs.sh " + path_in + path_out + '/ ' + str(proc), shell=True)

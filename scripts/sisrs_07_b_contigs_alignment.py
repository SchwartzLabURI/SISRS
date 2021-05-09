#!/usr/bin/env python3


'''
Last edit: Yana Hrytsenko May 9th, 2021
'''
import os
import sys
import subprocess

if __name__ == '__main__':

    # Store the command line in a seperate argument
    cmd = sys.argv

    if len(cmd) < 4:
        print("THIS SCRIPT REQUIERS 4 ARGUMENTS (-dir <path to output data directory> -p <num processors to use>)")
        exit()

    out_dir = ""
    sisrs = ""

    if '-dir' in cmd or '--directory' in cmd:
        out_dir = isFound('-dir','--directory',cmd)
        sisrs = out_dir
        if (sisrs[-1] != '/'):
            output_path = sisrs + '/'
        else:
            output_path = sisrs

        contigs_from_consensus = output_path + "SISRS_Run/" + "contigs_outputs"

        if not os.path.exists(contigs_from_consensus):
            os.mkdir(contigs_from_consensus)
    else:
        print("SPECIFY THE OUTPUT PATH (-dir, --directory). PROGRAM EXITING.")
        exit()

    proc = 1
    if '-p' in cmd or '--processors' in cmd:
        try:
            proc = int(isFound('-p','--processors',cmd))
        except:
            proc = 1
    else:
        print("SWITCHING TO DEFAULT 1 PROCESSOR (-p,--processors)")


    #first mask the ref sequence
    os.chmod('align_contigs.sh', 0o755)
    path_in = output_path + 'SISRS_Run/contigs_outputs/ '
    path_out = output_path + 'SISRS_Run/contigs_alignments/ ' #TODO: create if not created
    subprocess.call("./align_contigs.sh " + path_in + path_out + proc, shell=True)

#!/usr/bin/env python3

'''

This script aligns data from one taxon to the composite genome
'''

import os
from glob import glob
import argparse
from sisrs_06c_align2 import runBowtie

def makelinks(output_folder, sp, old_folder):
    for f in os.listdir(old_folder + '/SISRS_Run/' + sp): #hard links
        os.link(old_folder + '/SISRS_Run/' + sp + '/' + f, output_folder + '/SISRS_Run/' + sp + '/' + f)

def run6(sis, proc, folder, link):
    print(sis, folder)

    if link:
        makelinks(sis, folder, link) #hard links
    else:
        ref = f"{sis}/SISRS_Run/Composite_Genome/contigs"
        runBowtie(sis, proc, folder, ref)

if __name__ == '__main__':

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--directory', action='store', nargs="?")
    my_parser.add_argument('-p', '--processors', type=int, action='store', default=1, nargs="?")
    my_parser.add_argument('-f', '--folder', action='store', nargs="?")
    my_parser.add_argument('--link', action='store', nargs="?", default=None)
    args = my_parser.parse_args()

    sis = args.directory
    proc = args.processors
    folder = args.folder

    run6(sis, proc, folder, args.link)

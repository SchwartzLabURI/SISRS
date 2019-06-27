#!/usr/bin/env python3

'''
This script prepares data for a SISRS run by setting the data up and creating mapping scripts
Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
The composite genome is indexed by Bowtie2 and Samtools
SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder
'''

import os
from os import path
import sys
from glob import glob
import subprocess

'''
Get the correct bash scripts based on the working sisrs directory. The only
parameter neeeded for this function is the current working sisrs directory.
'''
def sisrsSetup(sisrs_dir):

    #Set TrimRead + SISRS directories based off of script folder location
    sisrs_dir = sisrs_dir+"/SISRS_Run"
    sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))

    sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]

    return sisrs_tax_dirs

'''
This function is designed to run the bash scripts that are gathered by sisrsSetup
it only requiers the output that is given by sisrsSetup function as input.
'''
def runSisrs(sisrs_tax_dirs):

    #Create links to trimmed read files in SISRS_Run directory
    for tax_dir in sisrs_tax_dirs:
        taxa = path.basename(tax_dir[:-1])
        with open(tax_dir+'out_'+taxa+'_SISRS',"w") as file:
            with open(tax_dir+'err_'+taxa+'_SISRS',"w") as file2:
                cmd = tax_dir+taxa+'.sh'
                subprocess.call(['sh',cmd],stdout=file, stderr=file2)
        print("Completed SISRS filtering for "+taxa + "...\n")

if __name__ == '__main__':

    cmd = sys.argv

    sis = os.path.dirname(sys.path[0])

    sisrs_tax_dirs = sisrsSetup(sis)
    runSisrs(sisrs_tax_dirs)

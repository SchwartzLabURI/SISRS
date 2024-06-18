#!/usr/bin/env python3

'''
This script prepares the folder architecture for a SISRS run.

Arguments:
-d, --directory : Path to data folder
-dir, --outputdir : Path to the output directory

Output: Script will create necessary folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders.
'''

import os
from os import listdir
from os.path import isdir, isfile, join
import argparse

def readTaxaFromFolders(inpath,outpath):

    '''
    Makes a file containing the list of taxa (based on the names of taxon folders in the input folder). Makes a folder for each taxon in various output folders. 

    Arguments: 
    inpath (string): folder containing the folders that contain reads
    outpath (string): folder where output will go

    Returns: 
    list: contains taxon list
    '''

    # Error check to make sure it is a path that exists
    if not isdir(inpath):
        print("The path provided does not exist")
        exit()

    # Make list of taxa based on the folder names
    taxa_list = [folder for folder in listdir(inpath) if isdir(join(inpath,folder))]
    taxa_list = sorted(taxa_list)

    #write taxon list to file
    with open(outpath+'/'+'TaxonList.txt', 'w') as f:
        f.write("\n".join(taxa_list))

    #make folder for each taxon
    for x in taxa_list:
        os.mkdir(outpath+"/Reads/RawReads/"+x)
        os.mkdir(outpath+"/Reads/TrimReads/"+x)
        os.mkdir(outpath+"/SISRS_Run/"+x)

    return taxa_list

def fileStructure(sisrs_dir):
    '''
    Function designed to build the necessary file structure for SISRS to properly run.

    Arguments: 
    sisrs_dir (string): the folder for output

    Returns: none.
    '''

    os.mkdir(sisrs_dir+"/Reads")
    os.mkdir(sisrs_dir+"/Reads/RawReads")
    os.mkdir(sisrs_dir+"/Reads/TrimReads")
    os.mkdir(sisrs_dir+"/Reads/SubsetReads")
    os.mkdir(sisrs_dir+"/SISRS_Run")

def run1(d, outputdir):
    # check if dir exists and dir is a dir and dir has files - stop
    if os.path.exists(outputdir):
        print("This directory already exists. Please specify a name for a new directory.")
        exit()

    # Create the file structure
    fileStructure(outputdir)

    # Grab the taxon names
    tl = readTaxaFromFolders(d, outputdir)

if __name__ == "__main__":

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--directory', action='store', nargs="?")
    my_parser.add_argument('-dir', '--outputdir', action='store', nargs="?")
    args = my_parser.parse_args()

    run1(args.directory, args.outputdir)

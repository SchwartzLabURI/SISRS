#!/usr/bin/env python3

'''
Last edit: Devin McConnell May 23, 2019
'''

# This script preps the folder architecture for a SISRS run.
# Taxon ID file (text file with Taxon IDs on new lines) must be in base SISRS folder
# Arguments: (1) Path to Taxon ID file (text file with Taxon IDs on new lines)
# Example: python scripts/folder_setup.py TaxonIDs
# Output: Script will create lots of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders

import sys
import cProfile
import os
from os import listdir,path
from os.path import isdir, join


def devTaxaList(taxon_ID_file):
    # Error check to make sure the file can be opened
    try:
        taxonFile = open(taxon_ID_file)
    except:
        print("ERROR WHILE PROCESSING FILE!")
        exit()

    taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

    return taxa_list

'''
This function is relient on the fact that the folders inside the path are labeled
according to the correct taxon. Takes as the mainm argumemt the path it will be
search for the arguments to the files
'''
def readTaxaFromFolders(path):
    # Error check to make sure it is a path that exists
    if isdir(path) == False:
        print("ERROR WHILE PROCESSING TAXON IDS!")
        exit()

    # Make list of taxa based on the foler names
    taxa_list = [folder for folder in listdir(path) if isdir(join(path,folder))]

    return taxa_list

'''
This function was desinged to make soft links to all the files and make links to
them into the proper folders. This function requiers the path to the data, a list
of tha taxons/folder names, and a boolean value determining if the reads are
trimmer already or if it is rawdata.
'''
def makeLinks(path, taxa_list, trim):
    dest = ""
    if trim:
        dest = "TrimReads"
    else:
        dest = "RawReads"

    for x in taxa_list:
        # Only batins the files that are not starting with '.'
        l = [f for f in listdir(path + '/' + x) if not f.startswith('.')]
        for i in l:
            # Creates the soft link to the files
            os.symlink(path + '/' + x + '/' + i,
                sisrs_dir+"/Reads/%s/"%dest +x + '/' + i)


'''
Function that is designed to build the necessay file structure
for sisrs to properly run. It will take in as input the taxa list
and the working sisrs directory
'''
def fileStructure(sisrs_dir,taxa_list):
    os.mkdir(sisrs_dir+"/Reads")
    os.mkdir(sisrs_dir+"/Reads/RawReads")
    os.mkdir(sisrs_dir+"/Reads/TrimReads")
    os.mkdir(sisrs_dir+"/Reads/SubsetReads")
    os.mkdir(sisrs_dir+"/SISRS_Run")
    os.mkdir(sisrs_dir+"/Reference_Genome")
    os.mkdir(sisrs_dir+"/Reference_Genome/Annotations")
    os.mkdir(sisrs_dir+"/Post_SISRS_Processing")
    os.mkdir(sisrs_dir+"/R_Analyses")
    os.mkdir(sisrs_dir+"/Reference_Topology")

    for x in taxa_list:
        os.mkdir(sisrs_dir+"/Reads/RawReads/"+x)
        os.mkdir(sisrs_dir+"/Reads/TrimReads/"+x)
        os.mkdir(sisrs_dir+"/SISRS_Run/"+x)


if __name__ == '__main__':
    taxon_ID_file = path.abspath(sys.argv[1])
    sisrs_dir = path.dirname(taxon_ID_file)
    taxa_list = readTaxaFromFolders(taxon_ID_file)
    fileStructure(sisrs_dir,taxa_list)
    makeLinks(taxon_ID_file, taxa_list,False)

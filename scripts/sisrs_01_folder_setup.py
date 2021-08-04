#!/usr/bin/env python3

'''

This script preps the folder architecture for a SISRS run.
Arguments: -d: Path to data folder
TO DO: Optional Argument: -trm/--trimmed (lets system know the data has been pretrimmed)
Output: Script will create a number of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders
'''

import sys
import cProfile
import os
from os import listdir,path
from cmdCheck import *
from os.path import isdir, isfile, join
import argparse

def devTaxaList(taxon_ID_file):
    ''' Returned list containes taxon list. '''

    # Error check to make sure the file can be opened
    if isfile(taxon_ID_file) == False:
        print("ERROR WHILE PROCESSING FILE!")
        exit()

    taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]

    return taxa_list

'''
This function is relient on the fact that the folders inside the path are labeled
according to the correct taxon. Takes as the main argumemt the path it will be
searching for the arguments to the files.
'''
def readTaxaFromFolders(inpath,outpath):
    ''' Returned list containes taxon list. '''

    # Error check to make sure it is a path that exists
    if isdir(inpath) == False:
        print("ERROR WHILE PROCESSING PATH to TAXON FOLDERS!")
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

'''
This function was desinged to make soft links to all the files and make links to
them into the proper folders. This function requiers the path to the data, a list
of tha taxons/folder names, and a boolean value determining if the reads are
trimmer already or if it is raw data.
'''
def makeLinks(data_path, sisrs_dir, taxa_list, trim):
    ''' This function makes soft links to the data files. '''

    dest = ""
    if trim:
        dest = "TrimReads"
    else:
        dest = "RawReads"

    for x in taxa_list:
        # Only batins the files that are not starting with '.'
        l = [f for f in listdir(data_path + '/' + x) if not f.startswith('.')]
        for i in l:
            # Creates the soft link to the files
            os.link(data_path + '/' + x + '/' + i,
                sisrs_dir+"/Reads/%s/"%dest +x + '/' + i)

'''
Function that is designed to build the necessary file structure
for SISRS to properly run. It will take in as input the taxa list
and the working sisrs directory
'''
def fileStructure(sisrs_dir):
    ''' This function builds necessary file structure to run SISRS. '''

    if isdir(sisrs_dir) == False:
        os.mkdir(sisrs_dir)

#    try:
        os.mkdir(sisrs_dir+"/Reads")
        os.mkdir(sisrs_dir+"/Reads/RawReads")
        os.mkdir(sisrs_dir+"/Reads/TrimReads")
        os.mkdir(sisrs_dir+"/Reads/SubsetReads")
        os.mkdir(sisrs_dir+"/SISRS_Run")
        #os.mkdir(sisrs_dir+"/Reference_Genome") #uncomment to have this dir
        #os.mkdir(sisrs_dir+"/Reference_Genome/Annotations") #uncomment to have this dir

        #os.mkdir(sisrs_dir+"/Post_SISRS_Processing")
        #os.mkdir(sisrs_dir+"/R_Analyses")
        #os.mkdir(sisrs_dir+"/Reference_Topology")
#    except:
#        print("SISRS RUN ALREADY EXISTS HERE! IF YOU ARE ADDING NEW DATA USE FLAG -aD or -aT.\n")
#        print("OTHERWISE CHANGE DIRECTORIES AND TRY AGAIN. PROGRAM EXITING!")
#        exit()


if __name__ == "__main__":

    #run as python3 sisrs_01_folder_setup.py -d $D -dir $DIR
    
    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d','--directory',action='store',nargs="?")
    my_parser.add_argument('-dir','--outputdir',action='store',nargs="?")
    args = my_parser.parse_args()

    d = args.directory
    outputdir = args.outputdir

    # Create the file structure
    fileStructure(outputdir)

    # Grab the taxon names
    #tl = []
    #if id != "":
    #    tl = devTaxaList(id)
    #else:
    tl = readTaxaFromFolders(d, outputdir)

    # Create the file structure
    fileStructure(outputdir)

    #fix given new argparse
#    if rd != "":
#        if '-trm' in cmd or '--trimmed' in cmd:
#            makeLinks(rd, sisrs, tl, True)
#        else:
#    makeLinks(d, outputdir, tl, False)

#!/usr/bin/env python3

'''

This script preps the folder architecture for a SISRS run.
Arguments: (-rd/--rawData or -id) Path to data or Path to Taxon ID file
Optional Argument: -trm/--trimmed (lets system know the data has been pretrimmed)
Output: Script will create a number of folders, including taxon folders in the RawReads, TrimReads, and SISRS_Run folders
'''

import sys
import cProfile
import os
from os import listdir,path
from cmdCheck import *
from os.path import isdir, isfile, join

'''
This function will be developing the taxon list based on an input file of taxon
ids.
'''
def devTaxaList(taxon_ID_file):
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
def readTaxaFromFolders(path):
    # Error check to make sure it is a path that exists
    if isdir(path) == False:
        print("ERROR WHILE PROCESSING PATH to TAXON FOLDERS!")
        exit()

    # Make list of taxa based on the folder names
    taxa_list = [folder for folder in listdir(path) if isdir(join(path,folder))]

    return taxa_list

'''
This function was desinged to make soft links to all the files and make links to
them into the proper folders. This function requiers the path to the data, a list
of tha taxons/folder names, and a boolean value determining if the reads are
trimmer already or if it is raw data.
'''
def makeLinks(data_path, sisrs_dir, taxa_list, trim):
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
def fileStructure(sisrs_dir,taxa_list):

    if isdir(sisrs_dir) == False:
        os.mkdir(sisrs_dir)

    try:
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
    except:
        print("SISRS RUN ALREADY EXISTS HERE! IF YOU ARE ADDING NEW DATA USE FLAG -aD or -aT.\n")
        print("OTHERWISE CHANGE DIRECTORIES AND TRY AGAIN. PROGRAM EXITING!")
        exit()


    for x in taxa_list:
        os.mkdir(sisrs_dir+"/Reads/RawReads/"+x)
        os.mkdir(sisrs_dir+"/Reads/TrimReads/"+x)
        os.mkdir(sisrs_dir+"/SISRS_Run/"+x)

if __name__ == "__main__":

    # Store the command line in a seperate argument
    cmd = sys.argv
    #sisrs = os.path.dirname(sys.path[0]) #use for a default path up one dir

    if len(cmd) < 5:
        print("THIS SCRIPT REQUIERS 4 ARGUMENTS (-d <path to input data directory> -dir <path to output data directory>)")
        exit()

    id = ""
    rd = ""
    out_dir = ""
    sisrs = ""

    if '-dir' in cmd or '--directory' in cmd:
        out_dir = isFound('-dir','--directory',cmd)
        sisrs = out_dir
    else:
        print("SPECIFY THE OUTPUT PATH (-dir, --directory).PROGRAM EXITING.")
        exit()


    if '-id' in cmd:
        id = cmd[cmd.index('-id') + 1]

    elif '-d' in cmd or '--data' in cmd:
        rd = isFound('-d','--data',cmd)

    else:
        print("MISSING TAXA INFORMATION (-d,--data)")
        exit()

    # Grab the taxon names
    tl = []
    if id != "":
        tl = devTaxaList(id)
    else:
        tl = readTaxaFromFolders(rd)

    # Create the file structure
    fileStructure(sisrs, tl)

    if rd != "":
        if '-trm' in cmd or '--trimmed' in cmd:
            makeLinks(rd, sisrs, tl, True)
        else:
            makeLinks(rd, sisrs, tl, False)

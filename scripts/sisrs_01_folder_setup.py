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

#def devTaxaList(taxon_ID_file):
#    '''
#    Gets the list of taxa from the taxon ID file
#
#    Parameters: 
#    taxon_ID_file (string): taxon ID file name
#    
#    Returns: 
#    list: contains taxon list
#
#    '''
#
#    # Error check to make sure the file can be opened
#    if isfile(taxon_ID_file) == False:
#        print("ERROR WHILE PROCESSING FILE!")
#        exit()
#
#    taxa_list = [line.rstrip('\n') for line in open(taxon_ID_file)]
#
#    return taxa_list

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


def makeLinks(data_path, sisrs_dir, taxa_list, trim):

    '''
    This function was desinged to make soft links to all the files and makes links to
    them into the proper folders.

    Arguments: This function requiers the path to the data, a list of tha taxons/folder names,
    and a boolean value determining if the reads are trimmed already or if it is raw data.

    Returns: none.

    '''

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


def fileStructure(sisrs_dir):
    '''
    Function designed to build the necessary file structure for SISRS to properly run.

    Arguments: 
    sisrs_dir (string): the folder containing the folders containing the input data

    Returns: none.
    '''



    if not isdir(sisrs_dir):
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

def run1(d, outputdir):
    # check if dir exists and dir is a dir and dir has files - stop
    if os.path.exists(outputdir):
        print("This directory already exists. Please specify a name for a new directory.")
        exit()

    # Create the file structure
    fileStructure(outputdir)

    # Grab the taxon names
    # tl = []
    # if id != "":
    #    tl = devTaxaList(id)
    # else:
    tl = readTaxaFromFolders(d, outputdir)

    # fix given new argparse


#    if rd != "":
#        if '-trm' in cmd or '--trimmed' in cmd:
#            makeLinks(rd, sisrs, tl, True)
#        else:
#    makeLinks(d, outputdir, tl, False)

if __name__ == "__main__":

    # Get arguments
    my_parser = argparse.ArgumentParser()
    my_parser.add_argument('-d', '--directory', action='store', nargs="?")
    my_parser.add_argument('-dir', '--outputdir', action='store', nargs="?")
    args = my_parser.parse_args()

    run1(args.directory, args.outputdir)

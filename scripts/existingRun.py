'''
This script needs the -dir flag as input.
'''

import os
import sys
import subprocess
from os import path
from sisrs_01_folder_setup import *

# ****************DETERMINING IF SISRS RUN EXISTS******************************
'''
Fucntion to examine the folder structure of where the user specfies the SISRS_DIR
will beself.
'''
def folderStruct(sisrs_dir):

    folders = ["/Ray_Composite_Genome","/Reads","/SISRS_Run"]
    rtn = ""
    for f in folders:
        if not path.isdir(sisrs_dir+f):
            print("SISRS FOLDER STRUCTURE NOT COMPLETE.\n NO PREVIOUS SISRS RUN FOUND.\n EXITING")
            exit()

'''
Function to check for the .nex files and read the alignment.nex file.
This function will be used to extract the taxon list for further data exploration

Return the list of old taxon names
'''
def getOldTaxons(sisrs_dir):
    sis_dir = sisrs_dir + '/SISRS_Run/'

    files = ["alignment_bi.nex", "alignment.nex", "alignment_pi.nex"]

    for f in files:
        if not path.isfile(sis_dir+f):
            print(f + " : is missing.\n PREVIOUS SISRS RUN NOT COMPLETE! EXITING")
            exit()

    # Examine alignment.nex file
    taxons = []
    i = 0
    nexFile = open(sis_dir+"alignment.nex")

    for line in nexFile:
        if i >= 7:
            line = line.strip()
            if line == ';':
                break
            data = line.split()
            taxons += [data[0]]

            if len(data[1]) == 0:
                print("alignment.nex file is missing data.\n PREVIOUS SISRS RUN NOT COMPLETE! EXITING")
                exit()

        i += 1

    return taxons

'''
This function is deisgned to move the old data to a seperate location and zip it
in order to reduce space.

Return the list of nex taxons
'''
def moveFiles(taxonList, sisrs_dir, data_list, addTaxon, addData):

    run = sisrs_dir+'/SISRS_Run/'
    newTaxon = readTaxaFromFolders(data_list)

    taxon = []
    if addTaxon:
        taxon = [n for n in newTaxon if n not in taxonList]
    else:
        for n in newTaxon:
            if n in taxonList:
                taxon += [n]
            else:
                print("NEW TAXON DETECTED. PLEASE CHECK DATA! EXITING.")
                exit()

    #Create folder for BBDuk StdOut
    if not path.isdir(sisrs_dir+"/backup"):
        os.mkdir(sisrs_dir+"/backup")

    subprocess.call("mv {0}alignment* {1}/backup".format(run,sisrs_dir),shell=True)
    subprocess.call("mv {0}out_SISRS* {1}/backup".format(run,sisrs_dir),shell=True)
    subprocess.call("mv {0}Output_Alignment.sh {1}/backup".format(run,sisrs_dir),shell=True)
    if addData:
        for item in taxon:
            subprocess.call("mv {0}{1} {2}/backup".format(run,item,sisrs_dir),shell=True)
            command = [
                'zip -1',
                '{0}/backup/{1}.zip'.format(sisrs_dir,item),
                '{0}/backup/{1}'.format(sisrs_dir,item),
                ';',
                'rm -r',
                '{0}/backup/{1}'.format(sisrs_dir,item)]
            os.system(" ".join(command))

    return taxon

'''
Running the trimmer if needed
'''
def existingTrimmer(sisrs_dir, taxonList):
    print("NOT IMPLEMENTED YET")


def previousRun(cmd):
    folderStruct(cmd[0])
    taxons = getOldTaxons(cmd[0])
    newTaxons = moveFiles(taxons, cmd[0], cmd[1], cmd[8], cmd[9])

    existingTrimmer()

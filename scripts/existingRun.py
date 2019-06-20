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

Return True or false depending on if the structre is there
'''
def folderStruct(sisrs_dir):

    folders = ["/Post_SISRS_Processing", "/R_Analyses", "/Ray_Composite_Genome",
                "/Reads", "/Reference_Genome", "/Reference_Topology", "/SISRS_Run"]
    rtn = ""
    for f in folders:
        if path.isdir(sisrs_dir+f):
            rtn = True
        else:
            print("SISRS FOLDER STRUCTURE NOT COMPLETE.\n NO PREVIOUS SISRS RUN FOUND.\n EXITING")
            exit()

    return rtn

'''
Function to check for the .nex files and read the alignment.nex file.
This function will be used to extract the taxon list for further data exploration

Return true or false depdning on the presence of the files
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

# ***************** ADDING A TAXA MODE **************************************


# ***************** ADDING DATA MODE **************************************


if __name__ == '__main__':
    sis = sys.argv[1]
    data = sys.argv[2]
    addingTaxon = True if '-aT' in sys.argv else False
    addingData = True if '-aD' in sys.argv else False
    tList = []
    booVal = folderStruct(sis)

    if booVal:
        tList = getOldTaxons(sis)

    moveFiles(tList,sis,data, addingTaxon, addingData)

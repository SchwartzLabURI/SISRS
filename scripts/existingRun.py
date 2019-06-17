'''
This script needs the -dir flag as input.
'''

import os
import sys
from os import path

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
            return False
    return rtn

'''
Helper function to grab the previous files that were used in the SISRS_Run. This
function will be used by getOldTaxons???

Return a list of old files
'''

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
            return False

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
                return False

        i += 1

    return taxons

'''
Function to verify/ grab the RawReads/TrimReads data
This function will need to determine if there is duplicate data.
'''

'''
Function to check for subset reads

Return True or False if they are missing
'''

'''
Function to look at SH file from sisrs_05 and all of its output
'''

'''
FUnction to look at sh file from sisrs_07 and all of its output
'''

if __name__ == '__main__':
    sis = sys.argv[1]
    tList = []

    booVal = folderStruct(sis)

    if booVal == False:
        print("NO SISRS RUN FOUND/COMPLETE SISRS RUN FOUND!")

    if booVal:
        b = getOldTaxons(sis)
        if b == False:
            print("NO COMPLETE SISRS RUN FOUND!")
        else:
            tList = b

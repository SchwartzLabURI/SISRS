'''
This script needs the -dir flag as input.
'''

import os
import sys
import subprocess
from os import path
from sisrs_01_folder_setup import *
from sisrs_02_read_trimmer import *
from sisrs_05_setup_sisrs import *
from sisrs_06_run_sisrs import *
from sisrs_07_output_sisrs import *

# ****************DETERMINING IF SISRS RUN EXISTS******************************
'''
Fucntion to examine the folder structure of where the user specfies the SISRS_DIR
will beself.
'''
def folderStruct(sisrs_dir):

    folders = ["/Ray_Composite_Genome","/Reads", "/SISRS_Run"]
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

    # Add all of the empty output directiers back into the SISRS_Run folder
    for item in taxon:
        if isdir(sisrs_dir+"/SISRS_Run/"+item) == False:
            os.mkdir(sisrs_dir+"/SISRS_Run/"+item)

    return taxon

# *************************** RUNNING TRIMMER **********************************
'''
Move new data to the taxon folder
'''
def moveData(taxonList,sisrs_dir,data):
    for item in taxonList:
        if isdir(sisrs_dir+"/Reads/RawReads/"+item) == False:
            os.mkdir(sisrs_dir+"/Reads/RawReads/"+item)

    for x in taxonList:
        # Only batins the files that are not starting with '.'
        l = [f for f in listdir(data + '/' + x) if not f.startswith('.')]
        for i in l:
            # Creates the soft link to the files
            os.link(data + '/' + x + '/' + i,
                sisrs_dir+"/Reads/RawReads/%s/"%x + '/' + i)

'''
Running the trimmer if needed
'''
def existingTrimmer(sisrs_dir, taxonList,data_path,threads):
    moveData(taxonList,sisrs_dir,data_path)
    bb = findAdapter()
    rtn = setup(sisrs_dir)
    rtn[5] = [item for item in rtn[5] for item1 in taxonList if item1 in item]

    # raw_fastqc_command
    newdFastqc(threads,rtn[3],rtn[5])

    trim(rtn[5],rtn[1],bb,rtn[2])

    # trim_fastqc_command
    newdFastqc(threads,rtn[4],rtn[5])

# ********************** RUNNING SISRS SETUP ***********************************
'''
This function is designed to go through and and setup excatly what is needed for
a sisrs run. It will only go through and setup up the new data that has been
added.
'''
def runSetup(outPath,threads,minread,threshold,newTaxons):
    trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir = obtainDir(outPath)
    trim_read_tax_dirs = [item for item in trim_read_tax_dirs for item1 in newTaxons if item1 in item]
    fileChanges(ray_dir,composite_dir)
    indexCompGenome(composite_dir,threads)
    sisrs_template = beginSetUp(composite_dir,sys.path[0])
    copyShFile(trim_read_tax_dirs,sisrs_dir,sisrs_template,composite_dir,outPath,threads,minread,threshold,sys.path[0])

# ************************** Running SISRS *************************************
'''
This function is desinged to run all of the sisrs 6 script. It will only run the
new data tha has been added.
'''
def runSISRS(sisrs_dir,newTaxons):
    sisrs_tax_dirs = sisrsSetup(sisrs_dir)
    sisrs_tax_dirs = [item for item in sisrs_tax_dirs for item1 in newTaxons if item1 in item]
    runSisrs(sisrs_tax_dirs)

# ***************************** Output SISRS ***********************************
'''
This function is desinged to run all of the sisrs 7 script. It will only run the
new data tha has been added.
'''
def outputSISRS(outPath,missing,newTaxons):
    composite_dir,sisrs_tax_dirs,sisrs_dir = getData(outPath)
    sisrs_tax_dirs = [item for item in sisrs_tax_dirs for item1 in newTaxons if item1 in item]
    createBash(composite_dir,sisrs_tax_dirs,sisrs_dir,outPath,missing,sys.path[0])
    runBash(sisrs_dir,sisrs_tax_dirs)

# ********************** PREVIOUS RUN CALLER ***********************************
def previousRun(cmd):
    folderStruct(cmd[0])
    taxons = getOldTaxons(cmd[0])
    newTaxons = moveFiles(taxons, cmd[0], cmd[1], cmd[8], cmd[9])

    if not cmd[2]:
        existingTrimmer(cmd[0],newTaxons,cmd[1],cmd[3])

    runSetup(cmd[0],cmd[3],cmd[6],cmd[5],newTaxons)
    runSISRS(cmd[0],newTaxons)
    outputSISRS(cmd[0],cmd[7],newTaxons)

'''
This script is designed to determine if the a previous SISRS run has been done in
the directory that is going to be used for output. This script will determine if
there has been a previous run based on the folder structure in that directory and
if '.nex' files exist with appropriate corresponding data. After this script has
determined if there is a previous run it will backup all of the old alignment*,
out_SISRS*, Output_Alignment.sh, and if you were adding new data to a taxon
the taxon folder will be moved and zipped. After everything is backedup appropriatly
when adding data we will: trim reads if needed, generate <SPP>.sh scripts that map
to the previous runs refrence genome, and then rerun sisrs_07_output_sisrs.py
'''


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

def oldTaxonHelper(directory):
    files = ["alignment_bi.nex", "alignment.nex", "alignment_pi.nex"]

    for f in files:
        if not path.isfile(directory+f):
            print(f + " : is missing.\n PREVIOUS SISRS RUN NOT COMPLETE! EXITING")
            exit()

    # Examine alignment.nex file
    taxons = []
    i = 0
    nexFile = open(directory+"alignment.nex")

    # Pullout all of the old taxons that were used
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
Function to check for the .nex files and read the alignment.nex file.
This function will be used to extract the taxon list for further data exploration
Return the list of old taxon names
'''
def getOldTaxons(sisrs_dir):
    sis_dir = sisrs_dir + '/SISRS_Run/'
    data = listdir(sis_dir)

    miss = [x for x in data if 'missing_' in x]

    if miss == []:
        return oldTaxonHelper(sis_dir)
    else:
        return oldTaxonHelper(sis_dir + miss[0] +'/')

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

    #Create folder for backup data
    if not path.isdir(sisrs_dir+"/backup"):
        os.mkdir(sisrs_dir+"/backup")

    files = listdir(run)
    f = [x for x in files if 'missing_' in x]

    # If the files were moves into subdircotries of missing taxons move by folders
    if f != []:
        for item in f:
            subprocess.call("mv {0}{1} {2}/backup".format(run,item,sisrs_dir),shell=True)
    else:
        subprocess.call("mv {0}alignment* {1}/backup".format(run,sisrs_dir),shell=True)
        subprocess.call("mv {0}out_SISRS* {1}/backup".format(run,sisrs_dir),shell=True)
        subprocess.call("mv {0}Output_Alignment*.sh {1}/backup".format(run,sisrs_dir),shell=True)

    if addData:
        for item in taxon:
            os.mkdir(sisrs_dir+"/backup/"+item)
            subprocess.call("mv {0}{1}/contigs* {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)
            subprocess.call("mv {0}{1}/err* {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)
            subprocess.call("mv {0}{1}/out* {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)
            subprocess.call("mv {0}{1}/*.sh {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)
            subprocess.call("mv {0}{1}/*.bam* {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)
            subprocess.call("mv {0}{1}/*LocList {2}/backup/{1}".format(run,item,sisrs_dir),shell=True)

    # Add all of the empty output directiers back into the SISRS_Run folder
    for item in taxon:
        if isdir(sisrs_dir+"/SISRS_Run/"+item) == False:
            os.mkdir(sisrs_dir+"/SISRS_Run/"+item)

    return taxon

# *************************** RUNNING TRIMMER **********************************
'''
Move new data to the taxon folder
'''
def moveData(taxonList,sisrs_dir,data,trim):
    dest = ""
    if trim:
        dest = "TrimReads"
    else:
        dest = "RawReads"

    for item in taxonList:
        if isdir(sisrs_dir+"/Reads/%s/"%dest+item) == False:
            os.mkdir(sisrs_dir+"/Reads/%s/"%dest+item)

    newFiles = []
    for x in taxonList:
        # Only batins the files that are not starting with '.'
        l = [f for f in listdir(data + '/' + x) if not f.startswith('.')]
        newFiles += l
        for i in l:
            # Creates the soft link to the files
            os.link(data + '/' + x + '/' + i,
                sisrs_dir+"/Reads/%s/%s/"%(dest,x) + '/' + i)

    return newFiles

'''
Running the trimmer if needed
'''
def existingTrimmer(sisrs_dir,taxonList,data_path,threads,notTrimmed):
    newFiles = moveData(taxonList,sisrs_dir,data_path,notTrimmed)
    bb = findAdapter()
    rtn = setup(sisrs_dir)
    rtn[5] = [item for item in rtn[5] for item1 in taxonList if item1 in item]

    # raw_fastqc_command
    newdFastqc(threads,rtn[3],rtn[5],newFiles)

    trim(rtn[5],rtn[1],bb,rtn[2],newFiles)

    # trim_fastqc_command
    newdFastqc(threads,rtn[4],rtn[5],newFiles)

# ********************** RUNNING SISRS SETUP ***********************************
'''
This function is designed to go through and and setup excatly what is needed for
a sisrs run. It will only go through and setup up the new data that has been
added.
'''
def runSetup(outPath,threads,minread,threshold,newTaxons):

    # Grab the folders and files that are needed to setup sisrs
    trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir = obtainDir(outPath)

    # Only the new data
    trim_read_tax_dirs = [item for item in trim_read_tax_dirs for item1 in newTaxons if item1 in item]

    # Make the neccessary changes to the files
    fileChanges(ray_dir,composite_dir)

    # index the genomes
    indexCompGenome(composite_dir,threads)

    # Create the template bash script
    sisrs_template = beginSetUp(composite_dir,sys.path[0])

    # Make a bash script for all Taxons
    copyShFile(trim_read_tax_dirs,sisrs_dir,sisrs_template,composite_dir,outPath,threads,minread,threshold,sys.path[0])

# ************************** Running SISRS *************************************
'''
This function is desinged to run all of the sisrs 6 script. It will only run the
new data tha has been added.
'''
def runSISRS(sisrs_dir,newTaxons):

    # Get all folders and files needed to run sisrs
    sisrs_tax_dirs = sisrsSetup(sisrs_dir)

    # Only the new data
    sisrs_tax_dirs = [item for item in sisrs_tax_dirs for item1 in newTaxons if item1 in item]

    # Run all the SISRS bash scripts that were created
    runSisrs(sisrs_tax_dirs)

# ***************************** Output SISRS ***********************************
'''
This function is desinged to run all of the sisrs 7 script. It will only run the
new data tha has been added.
'''
def outputSISRS(outPath,missing,newTaxons):
    composite_dir,sisrs_tax_dirs,sisrs_dir = getData(outPath)
    sisrs_tax_dirs = [item for item in sisrs_tax_dirs for item1 in newTaxons if item1 in item]

    # If there is a range of allowed missing species do them all
    if '-' in str(missing):
        ms = missing.split('-')
        ms[0] = int(ms[0])
        ms[1] = int(ms[1])

        for i in range(ms[0],ms[1]+1):
            # Create the bash script to run sisrs as we need it to
            createBash(composite_dir,sisrs_tax_dirs,sisrs_dir,outPath,i,sys.path[0])

            #RunSisrs
            runBash(sisrs_dir,sisrs_tax_dirs,i)

            subprocess.call("mkdir {0}/missing_{1}".format(sisrs_dir,i),shell=True)
            subprocess.call("mv {0}/alignment* {0}/missing_{1}".format(sisrs_dir,i),shell=True)
            subprocess.call("mv {0}/out_SISRS* {0}/missing_{1}".format(sisrs_dir,i),shell=True)
            subprocess.call("mv {0}/Output_Alignment_m{1}.sh {0}/missing_{1}".format(sisrs_dir,i),shell=True)

    else:
        # Create the bash script to run sisrs as we need it to
        createBash(composite_dir,sisrs_tax_dirs,sisrs_dir,outPath,missing,sys.path[0])

        # Run sisrs
        runBash(sisrs_dir,sisrs_tax_dirs,missing)

# ********************** PREVIOUS RUN CALLER ***********************************
def previousRun(cmd):
    folderStruct(cmd[0])
    taxons = getOldTaxons(cmd[0])
    newTaxons = moveFiles(taxons, cmd[0], cmd[1], cmd[8], cmd[9])

    if not cmd[2]:
        existingTrimmer(cmd[0],newTaxons,cmd[1],cmd[3],cmd[2])
    else:
        moveData(newTaxons, cmd[0], cmd[1], cmd[2])

    runSetup(cmd[0],cmd[3],cmd[6],cmd[5],newTaxons)
    runSISRS(cmd[0],newTaxons)

    taxons += [x for x in newTaxons if x not in taxons]
    outputSISRS(cmd[0],cmd[7],taxons)

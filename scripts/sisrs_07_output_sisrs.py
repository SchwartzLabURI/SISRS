#!/usr/bin/env python3

'''
This script prepares data for a SISRS run by setting the data up and creating mapping scripts
Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
The composite genome is indexed by Bowtie2 and Samtools
SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder
Input: -ms/--missing
'''


import os
from os import path
import sys
from glob import glob
import subprocess
from subprocess import Popen
from cmdCheck import *
from itertools import islice

'''
This function is designed to get the appropriate directories needed to finish the
sisrs run. The parameters that is needed for it to properly run is just the
current sisrs working directory.
'''
def getData(outPath):

    #Set TrimRead + SISRS directories based off of script folder location
    sisrs_dir = outPath+"/SISRS_Run"
    composite_dir = sisrs_dir + '/Composite_Genome'
    sisrs_tax_dirs = sorted(glob(sisrs_dir+"/*/"))
    sisrs_tax_dirs = [x for x in sisrs_tax_dirs if not x.endswith('Composite_Genome/')]

    return composite_dir,sisrs_tax_dirs,sisrs_dir

'''
This fucntion is designed to create all the bash scripts that are still needed
to finish off the sisrs run. It requiers as input the current working directory
for sisrs, the location of the "SISRS_Run" folder, the missing data parameter,
and the list of composite genome paths.
'''
def createBash(composite_dir,sisrs_tax_dirs,sisrs_dir,outPath,missing,dir):

    sisrs_output_template = ""

    shFile = open(dir + '/sisrs_07_template.sh')
    for line in shFile:
        sisrs_output_template += line

    keyList = ['TWOTAXA','SISRS_DIR','SCRIPT_DIR','MISSING','COMPOSITE_DIR']
    keyDict = {'TWOTAXA':str(len(sisrs_tax_dirs) - 2),
               'SISRS_DIR':sisrs_dir,
               'SCRIPT_DIR':dir,
               'MISSING':str(missing),
               'COMPOSITE_DIR':composite_dir}
    for key in keyList:
        sisrs_output_template = sisrs_output_template.replace(key,keyDict[key])
    with open(sisrs_dir+"/Output_Alignment_m{}.sh".format(missing), "w") as text_file:
        print(sisrs_output_template, file=text_file)

'''
This function is designed to run all of the newly created bash scripts. It
requiers as input the path to the SISRS_Run directory and the path to all of
the folders that contain sh scripts.
'''
def runBash(sisrs_dir,sisrs_tax_dirs,missing):
    with open(sisrs_dir+"/out_SISRS_Alignment","w") as file:
        cmd = sisrs_dir+'/Output_Alignment_m{}.sh'.format(missing)
        p = Popen(['sh', cmd], stdout=file, stderr=subprocess.PIPE)
        output, err = p.communicate()
        rc = p.returncode
        print(output)
        print(err)
        print(rc)
        #subprocess.call(['sh',cmd],stdout=file, stderr=subprocess.PIPE)

    with open(sisrs_dir+"/out_SISRS_Log","w") as file:
        file.write("\nRead Mapping and SISRS Site Selection:\n")
        for tax_dir in sisrs_tax_dirs:
            taxa = path.basename(tax_dir[:-1])
            bowtie1 = ['grep','-A4',"'of these'",'{}'.format(tax_dir + "err_" + taxa + "_SISRS"),'|','sed','-n','2,6p']
            bowtie2 = ['grep','-A4',"'of these'",'{}'.format(tax_dir + "err_" + taxa + "_SISRS"),'|','sed','-n','9,13p']
            file.write("\n"+ taxa + " Composite Genome Mapping:\n\n")
            file.write((subprocess.check_output(' '.join(bowtie1),shell=True).decode("UTF8")))
            file.write("\n"+taxa+" Specific Genome Mapping:\n\n")
            file.write((subprocess.check_output(' '.join(bowtie2),shell=True).decode("UTF8")))

            with open(tax_dir + "out_" + taxa + "_SISRS") as f:
                file.write("\n"+taxa+" SISRS Site Selection:\n\n")
                for line in f:
                    if(str.startswith(line,'Of ')):
                        file.write(line)
        with open(sisrs_dir + "/out_SISRS_Alignment") as f2:
            file.write("\nSISRS Alignment Filtering:\n\n")
            for line in f2:
                file.write(line)

if __name__ == '__main__':
    cmd = sys.argv
    sis = os.path.dirname(sys.path[0])

    md = 0

    # Gets the allowed amount of missing data from the files
    if '-ms' in cmd or '--missing' in cmd:
        bool = isInt(isFound('-ms','--missing',cmd))
        md = int(isFound('-ms','--missing',cmd)) if bool else 0
        if bool == False:
            bool2 = isRange(isFound('-ms','--missing',cmd))
            md = isFound('-ms','--missing',cmd) if bool2 else 0
            if bool2 == False:
                print("THE ALLOWED AMOUNT OF MISSING DATA (-ms/--missing) MUST BE A VALID NUMBER --> SWITCHED TO 0")
                md = 0
    else:
        print("USING DEFAULT MISSING (-ms,--missing) DATA: 0")

    composite_dir,sisrs_tax_dirs,sisrs = getData(sis)

    if '-' in str(md):
        ms = md.split('-')
        ms[0] = int(ms[0])
        ms[1] = int(ms[1])

        for i in range(ms[0],ms[1]+1):
            # Create the bash script to run sisrs as we need it to
            createBash(composite_dir,sisrs_tax_dirs,sisrs,sis,i,sys.path[0])

            #RunSisrs
            runBash(sisrs,sisrs_tax_dirs,i)

            subprocess.call("mkdir {0}/missing_{1}".format(sisrs,i),shell=True)
            subprocess.call("mv {0}/alignment* {0}/missing_{1}".format(sisrs,i),shell=True)
            subprocess.call("mv {0}/out_SISRS* {0}/missing_{1}".format(sisrs,i),shell=True)
            subprocess.call("mv {0}/Output_Alignment_m{1}.sh {0}/missing_{1}".format(sisrs,i),shell=True)

    else:
        # Create the bash script to run sisrs as we need it to
        createBash(composite_dir,sisrs_tax_dirs,sisrs,sis,md,sys.path[0])

        # Run sisrs
        runBash(sisrs,sisrs_tax_dirs,md)

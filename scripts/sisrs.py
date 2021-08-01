#!/usr/bin/env python3
'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''
from existingRun import *
from sisrs_01_folder_setup import *
from sisrs_02_read_trimmer import *
from sisrs_03_read_subsetter import *
from sisrs_04_ray_composite import *
from sisrs_05_setup_sisrs import *
from sisrs_06_run_sisrs import *
from sisrs_07_output_sisrs import *
from cmdCheck import *

'''
Function to run all of the first script
'''
def sisrs01(data_path,sisrs_dir,trimmed):
    ''' This function runs step 1 of SISRS. '''

    taxa_list = []

    taxa_list = readTaxaFromFolders(data_path)

    fileStructure(sisrs_dir, taxa_list)


    if trimmed:
        makeLinks(data_path, sisrs_dir, taxa_list, True)
    else:
        makeLinks(data_path, sisrs_dir, taxa_list, False)

'''
Function to run all of the second script
'''
def sisrs2(trimmed,processors,sisrs_dir):
    ''' This function runs step 2 of SISRS. '''

    if not trimmed:
        bbduk_adapter = findAdapter()
        out = setup(sisrs_dir)

        # raw_fastqc_command
        fastqcCommand(processors,out[3],out[0])

        trim(out[5],out[1],bbduk_adapter,out[2],[])

        # trim_fastqc_command
        fastqcCommand(processors,out[4],out[1])
    else:
        print("FILES ALREADY TRIMMED")

'''
Function to run all of the thrid script
'''
def sisrs3(genomeSize, sisrs_dir):
    ''' This function runs step 3 of SISRS. '''

    setupInfo = setupDir(sisrs_dir,genomeSize)

    df, compiled_paired, compiled_single_end = countBasePair(setupInfo[3],setupInfo[6], setupInfo[7],setupInfo[5])

    print("Based on a genome size estimate of " + str(genomeSize) + " bp, and with " + str(len(setupInfo[3])) + " species, the requested subset depth is " + str(setupInfo[4]) + " bp per species")

    out = checkCoverage(df,setupInfo[4],setupInfo[1],compiled_paired, compiled_single_end)

    # Subset Single End
    subset(compiled_single_end,out[2],setupInfo[0],setupInfo[1],setupInfo[2],False)

    #Subset paired ends
    subset(compiled_paired,out[2],setupInfo[0],setupInfo[1],setupInfo[2],True)

'''
Function to run all of the fourth script
'''
def sisrs4(sisrs_dir,threads):
    ''' This function runs step 4 of SISRS. '''

    # obtain the directories and files needed for Ray
    ray_genome_dir, subset_reads = getDirs(sisrs_dir)

    # Run the ray command
    runRay(ray_genome_dir,subset_reads,threads)

'''
Function to run all of the fifth script
'''
def sisrs05(outPath,threads,minread,threshold):
    ''' This function runs step 5 of SISRS. '''

    # Grab the folders and files that are needed to setup sisrs
    trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir = obtainDir(outPath)

    # Make the neccessary changes to the files
    fileChanges(ray_dir,composite_dir)

    # index the genomes
    indexCompGenome(composite_dir,threads)

    # Create the template bash script
    sisrs_template = beginSetUp(composite_dir,sys.path[0])

    # Make a bash script for all Taxons
    copyShFile(trim_read_tax_dirs,sisrs_dir,sisrs_template,composite_dir,outPath,threads,minread,threshold,sys.path[0])

'''
Function to run all of the sixth script
'''
def sisrs06(sisrs_dir):
    ''' This function runs step 6 of SISRS. '''

    # Get all folders and files needed to run sisrs
    sisrs_tax_dirs = sisrsSetup(sisrs_dir)

    # Run all the SISRS bash scripts that were created
    runSisrs(sisrs_tax_dirs)

def sisrs07(outPath,missing):
    ''' This function runs step 7 of SISRS. '''

    # Grab the neccessary folders and files
    composite_dir,sisrs_tax_dirs,sisrs_dir = getData(outPath)

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

if __name__ == '__main__':
    cmdln = sys.argv
    rtn = commandLine(cmdln,os.path.dirname(sys.path[0]))

    if rtn[8] == True or rtn[9] == True:
        previousRun(rtn)
        exit()

    sisrs01(rtn[1],rtn[0],rtn[2])
    sisrs2(rtn[2],rtn[3],rtn[0])
    sisrs3(rtn[4],rtn[0])
    sisrs4(rtn[0],rtn[5])
    sisrs05(rtn[0],rtn[3],rtn[6],rtn[5])
    sisrs06(rtn[0])
    sisrs07(rtn[0],rtn[7])

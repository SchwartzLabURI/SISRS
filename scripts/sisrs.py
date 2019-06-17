#!/usr/bin/env python3
'''
Script Developed by Devin J. McConnell
Schwartz Lab May 2019
'''

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
def sisrs01(data_path,sisrs_dir,trimed):
    taxa_list = []

    taxa_list = readTaxaFromFolders(data_path)

    fileStructure(sisrs_dir, taxa_list)

    if trimed:
        makeLinks(data_path, sisrs_dir, taxa_list, True)
    else:
        makeLinks(data_path, sisrs_dir, taxa_list, False)

'''
Function to run all of the second script
'''
def sisrs2(trimed,processors,sisrs_dir):
    if not trimed:
        bbduk_adapter = findAdapter()
        out = setup(sisrs_dir)

        # raw_fastqc_command
        fastqcCommand(processors,out[3],out[0])

        trim(out[5],out[1],bbduk_adapter,out[2])

        # trim_fastqc_command
        fastqcCommand(processors,out[4],out[1])
    else:
        print("FILES ALREADY TRIMMED")

'''
Function to run all of the thrid script
'''
def sisrs3(genomeSize, sisrs_dir):
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

    # obtain the directories and files needed for Ray
    ray_genome_dir, subset_reads = getDirs(sisrs_dir)

    # Run the ray command
    runRay(ray_genome_dir,subset_reads,threads)

'''
Function to run all of the fifth script
'''
def sisrs05(outPath,threads,minread,threshold):
    trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir = obtainDir(outPath)
    fileChanges(ray_dir,composite_dir)
    indexCompGenome(composite_dir,threads)
    sisrs_template = beginSetUp(composite_dir,sys.path[0])
    copyShFile(trim_read_tax_dirs,sisrs_dir,sisrs_template,composite_dir,outPath,threads,minread,threshold,sys.path[0])

'''
Function to run all of the sixth script
'''
def sisrs06(sisrs_dir):
    sisrs_tax_dirs = sisrsSetup(sisrs_dir)
    runSisrs(sisrs_tax_dirs)

def sisrs07(outPath,missing):
    composite_dir,sisrs_tax_dirs,sisrs_dir = getData(outPath)
    createBash(composite_dir,sisrs_tax_dirs,sisrs_dir,outPath,missing,sys.path[0])
    runBash(sisrs_dir,sisrs_tax_dirs)

if __name__ == '__main__':
    cmdln = sys.argv
    rtn = commandLine(cmdln,os.path.dirname(sys.path[0]))
    sisrs01(rtn[1],rtn[0],rtn[3])
    sisrs2(rtn[2],rtn[3],rtn[0])
    sisrs3(rtn[4],rtn[0])
    sisrs4(rtn[0],rtn[3])
    sisrs05(rtn[0],rtn[3],rtn[6],rtn[5])
    sisrs06(rtn[0])
    sisrs07(rtn[0],rtn[7])

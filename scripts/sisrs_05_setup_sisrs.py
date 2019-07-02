#!/usr/bin/env python3

'''
This script prepares data for a SISRS run by setting the data up and creating
mapping scripts
Contigs are renamed and moved to the SISRS_Run/Composite_Genome directory
The composite genome is indexed by Bowtie2 and Samtools
SISRS scripts are generated from a template and saved to the SISRS_Run/TAXA folder
Input: -p/--processors Number of processors, -mr/--minRead, and -thresh/--threshold
'''

import os
from os import path
import sys
from glob import glob
from cmdCheck import *
import pandas as pd

'''
This function is designed to obtain all of the directories needed to finish the
setup of sisrs. This function only needs the working sisrs directory as a
parameter.
'''
def obtainDir(outPath):

    #Set TrimRead + SISRS directories based off of script folder location
    trim_read_dir =  outPath+"/Reads/TrimReads"
    trim_read_tax_dirs = sorted(glob(trim_read_dir+"/*/"))
    ray_dir = outPath+"/Ray_Composite_Genome"
    sisrs_dir = outPath+"/SISRS_Run"

    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('trimOutput/')]
    trim_read_tax_dirs = [x for x in trim_read_tax_dirs if not x.endswith('fastqcOutput/')]
    trim_read_tax_dirs = sorted([x for x in trim_read_tax_dirs if not x.endswith('subsetOutput/')])

    #Create composite genome folder
    if(not path.isdir(sisrs_dir+"/Composite_Genome")):
        os.mkdir(sisrs_dir+"/Composite_Genome")
    composite_dir =sisrs_dir+"/Composite_Genome"

    return trim_read_tax_dirs,ray_dir,sisrs_dir,composite_dir

'''
This function is designed to raname the contig files and make copies of the
length files that are a product of the output from Ray. It requieres as input
the directory were ray output all of its dated and the composite directory.
'''
def fileChanges(ray_dir,composite_dir):

    #Create renamed contig file in SISRS_Run/Composite_Genome
    rename_command = [
        'rename.sh',
        'in={}'.format(ray_dir+"/Contigs.fasta"),
        'out={}'.format(composite_dir+'/contigs.fa'),
        'prefix=SISRS',
        'addprefix=t',
        'ow=t']
    os.system(" ".join(rename_command))

    #Copy sequence length file from Ray
    contig_length = pd.read_csv(ray_dir+'/ContigLengths.txt',sep="\t",header=None)
    contig_length.columns = ['Contig','Length']
    contig_length.Contig = ("SISRS_" + contig_length.Contig)
    contig_length.to_csv(composite_dir+"/contigs_SeqLength.tsv",sep="\t",header=None,index=False)

'''
This function is designed to call the bowtie comman and the samtools command
in order to index the composite genomes correctly. This function requiers two
parameters the composite genome directory and the number of threads that it can
run with
'''
def indexCompGenome(composite_dir,threads):
    #Index composite genome using Bowtie2 and Samtools
    bowtie_command = [
        'bowtie2-build',
        '{}'.format(composite_dir+'/contigs.fa'),
        '{}'.format(composite_dir+'/contigs'),
        '-p',
        '{}'.format(str(threads))]
    os.system(" ".join(bowtie_command))

    samtools_command = [
        'samtools',
        'faidx',
        '{}'.format(composite_dir+'/contigs.fa')]
    os.system(" ".join(bowtie_command))

'''
This function is designed to start the file setups requiered for sisrs to run
properly. It requires as parameters the composite_dir and the location of were
all of the scripts are located
'''
def beginSetUp(composite_dir,script_dir):
    #Create file with an entry for each individual site in the alignment
    siteCount=0
    locListFile = open(composite_dir+'/contigs_LocList','a+')
    with open(composite_dir +"/contigs_SeqLength.tsv","r") as filein:
        for line in iter(filein):
            splitline=line.split()
            for x in range(1,(int(splitline[1])+1)):
                locListFile.write((splitline[0] +'/'+str(x)+'\n'))
                siteCount+=1
    locListFile.close()

    print("==== Site list created: " + str(siteCount) + " total sites ==== \n",flush=True)

    # Read in the bash script
    sisrs_template = ""

    shFile = open(script_dir + '/sisrs_05_template.sh')
    for line in shFile:
        sisrs_template += line

    return sisrs_template


'''
This function is designed to complete the sisrs setup and will be going through
to replace all of the non-specfic portions of the bash script with more specfic
portions of the scripts. This function requiers all of the directory information
alond with the number of threads, mindread, threshold, and the script directory.
'''
def copyShFile(trim_read_tax_dirs,sisrs_dir,sisrs_template,composite_dir,outPath,threads,minread,threshold,script):

    composite_dir = os.path.abspath(composite_dir)
    sisrs_dir = os.path.abspath(sisrs_dir)
    #Create links to trimmed read files in SISRS_Run directory
    for tax_dir in trim_read_tax_dirs:
        tax_dir = os.path.abspath(tax_dir)
        taxa = path.basename(tax_dir[:])
        read_link_command = [
            'cd',
            '{}'.format(sisrs_dir+"/"+taxa),
            ';',
            'cp',
            '-as',
            '{}/*.fastq.gz'.format(tax_dir),
            '.']
        os.system(" ".join(read_link_command))

        new_sisrs = sisrs_template

        keyList = ['PROCESSORS','BOWTIE2-INDEX','COMPOSITE_GENOME','SCRIPT_DIR','MINREAD','THRESHOLD','TAXA','SISRS_DIR','COMPOSITE_DIR','READS']
        keyDict = {'PROCESSORS':str(threads),
                   'BOWTIE2-INDEX':composite_dir+"/contigs",
                   'COMPOSITE_GENOME':composite_dir+"/contigs.fa",
                   'SCRIPT_DIR':script,
                   'MINREAD':str(minread),
                   'THRESHOLD':str(threshold),
                   'TAXA':taxa,
                   'SISRS_DIR':sisrs_dir,
                   'COMPOSITE_DIR':composite_dir,
                   'READS':",".join(glob(tax_dir+"/*.fastq.gz"))}
        for key in keyList:
            new_sisrs = new_sisrs.replace(key,keyDict[key])
        with open(sisrs_dir+"/"+taxa+"/"+taxa+".sh", "w") as text_file:
            print(new_sisrs, file=text_file)
        os.system('chmod +x '+sisrs_dir+"/"+taxa+"/"+taxa+".sh")

if __name__ == '__main__':
    cmd = sys.argv
    sis = os.path.dirname(sys.path[0])

    proc = 1
    mr = 3
    thres = 1

    if '-p' in cmd or '--processors' in cmd:
        try:
            proc = int(isFound('-p','--processors',cmd))
        except:
            print("INVALID NUMBER OF PROCESSORS (-p,--processors): DEFAULT PROCESSOR COUNT TO 1")
            proc = 1
    else:
        print("SETTING PROCESSORS COUNT (-p,--processors) TO 1")
        proc = 1

    if '-mr' in cmd or '--minread' in cmd:
        try:
            mr = int(isFound('-mr','--minread',cmd))
        except:
            print("INVALID NUMBER OF MINREAD (-mr,--minread): DEFAULT THREAD COUNT TO 3")
            mr = 3
    else:
        print("SETTING MINREAD (-mr,--minreads TO 3")
        mr = 3

    if '-thresh' in cmd or '--threshold' in cmd:
        try:
            thres = int(isFound('-thresh','--threshold',cmd))
        except:
            print("INVALID NUMBER OF THRESHOLD (-thresh,--threshold): DEFAULT THRESHOLD COUNT TO 1")
            thres = 1
    else:
        print("SETTING THRESHOLD COUNT (-thresh,--threshold) TO 1")
        thres = 1

    trim_read_tax_dirs,ray_dir,sisrs,composite_dir = obtainDir(sis)
    fileChanges(ray_dir,composite_dir)
    indexCompGenome(composite_dir,proc)
    sisrs_template = beginSetUp(composite_dir,sys.path[0])
    copyShFile(trim_read_tax_dirs,sisrs,sisrs_template,composite_dir,sis,proc,mr,thres,sys.path[0])
